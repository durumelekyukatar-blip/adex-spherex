"""Extract a screening spectrum for a sky position from SPHEREx Level-2 MEFs."""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS

SOURCE_MASK = 1 << 21  # SOURCE is an informational catalog-match flag, not a bad-pixel flag.
BAD_FLAG_MASK = 0xFFFFFFFF ^ SOURCE_MASK


def wavelength_at(hdus, x, y):
    spectral_wcs = WCS(hdus["IMAGE"].header, fobj=hdus, key="W")
    # The wavelength mapping is detector-fixed; do not apply the sky SIP distortion.
    spectral_wcs.sip = None
    wavelength, bandwidth = spectral_wcs.pixel_to_world(x, y)
    if not np.isfinite(wavelength.value) or not np.isfinite(bandwidth.value):
        return None, None
    return float(wavelength.to_value(u.um)), float(bandwidth.to_value(u.um))


def measure(path, ra, dec, radius):
    with fits.open(path, memmap=True) as hdus:
        image_hdu = hdus["IMAGE"]
        image = np.asarray(image_hdu.data, dtype=float)
        wcs = WCS(image_hdu.header, fobj=hdus).celestial
        px, py = (float(v) for v in wcs.world_to_pixel_values(ra, dec))
        if not np.isfinite(px + py):
            return None
        x, y = int(round(px)), int(round(py))
        if x < radius + 1 or y < radius + 1 or x >= image.shape[1] - radius - 1 or y >= image.shape[0] - radius - 1:
            return None
        flags_hdu = next((h for h in hdus if h.name.upper() == "FLAGS"), None)
        variance_hdu = next((h for h in hdus if h.name.upper() == "VARIANCE"), None)
        flags = np.asarray(flags_hdu.data) if flags_hdu is not None else np.zeros(image.shape, dtype=np.uint16)
        variance = np.asarray(variance_hdu.data, dtype=float) if variance_hdu is not None else None
        yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
        distance2 = (xx - px) ** 2 + (yy - py) ** 2
        aperture = distance2 <= radius ** 2
        annulus = (distance2 >= (radius + 2) ** 2) & (distance2 <= (radius + 6) ** 2)
        good = np.isfinite(image) & (np.bitwise_and(flags.astype(np.uint32), BAD_FLAG_MASK) == 0)
        if variance is not None:
            good &= np.isfinite(variance) & (variance >= 0)
        bg_pixels, ap_pixels = image[annulus & good], image[aperture & good]
        if bg_pixels.size < 8 or ap_pixels.size < 4:
            return None
        background = float(np.median(bg_pixels))
        brightness = float(np.mean(ap_pixels) - background)
        if variance is not None:
            varpix = variance[aperture & good]
            unc = float(np.sqrt(np.nansum(varpix)) / max(1, ap_pixels.size)) if varpix.size else None
        else:
            unc = None
        wavelength, bandwidth = wavelength_at(hdus, px, py)
        if wavelength is None:
            return None
        name = Path(path).name
        match = re.search(r"_(\d)D(\d)_", name)
        return {
            "source_file": name,
            "observation_id": name.split("_spx_")[0].removeprefix("level2_"),
            "date_obs": image_hdu.header.get("DATE-AVG", image_hdu.header.get("DATE-OBS")),
            "mjd_obs": image_hdu.header.get("MJD-AVG", image_hdu.header.get("MJD-OBS")),
            "exposure_seconds": image_hdu.header.get("XPOSURE"),
            "small_slew": int(match.group(1)) if match else None,
            "detector": int(image_hdu.header.get("DETECTOR", match.group(2) if match else 0)) or None,
            "wavelength_um": wavelength,
            "bandwidth_um": bandwidth,
            "brightness_mjy_sr": brightness,
            "uncertainty_mjy_sr": unc,
            "quality_ok": True,
            "x_pixel": float(px),
            "y_pixel": float(py),
        }


def auto_target(path):
    """Choose the highest unflagged local-contrast peak; review the result."""
    with fits.open(path, memmap=True) as hdus:
        hdu = hdus["IMAGE"]
        image = np.asarray(hdu.data, dtype=float)
        flags_hdu = next((h for h in hdus if h.name.upper() == "FLAGS"), None)
        var_hdu = next((h for h in hdus if h.name.upper() == "VARIANCE"), None)
        valid = np.isfinite(image)
        if flags_hdu is not None:
            valid &= np.bitwise_and(np.asarray(flags_hdu.data).astype(np.uint32), BAD_FLAG_MASK) == 0
        if not valid.any():
            raise ValueError("No unflagged pixels in first FITS image")
        replacement = float(np.median(image[valid]))
        residual = image - replacement
        if var_hdu is not None:
            variance = np.asarray(var_hdu.data, dtype=float)
            valid &= np.isfinite(variance) & (variance > 0)
            score = residual / np.sqrt(np.maximum(variance, np.finfo(float).tiny))
        else:
            scale = 1.4826 * np.median(np.abs(residual[valid] - np.median(residual[valid])))
            score = residual / max(float(scale), np.finfo(float).eps)
        valid[:16, :] = False
        valid[-16:, :] = False
        valid[:, :16] = False
        valid[:, -16:] = False
        score[~valid] = -np.inf
        score[~valid] = -np.inf
        y, x = np.unravel_index(np.argmax(score), score.shape)
        if not np.isfinite(score[y, x]) or score[y, x] < 5:
            raise ValueError("No unflagged local peak with S/N >= 5 found in the first image")
        ra, dec = WCS(hdu.header, fobj=hdus).celestial.pixel_to_world_values(float(x), float(y))
        return float(ra), float(dec), float(score[y, x])


def main():
    parser = argparse.ArgumentParser(description="Build a SPHEREx infrared screening spectrum at a sky coordinate.")
    parser.add_argument("fits", nargs="+", type=Path, help="SPHEREx Level-2 FITS files or folders containing them")
    parser.add_argument("--ra", type=float, help="Right ascension in ICRS degrees")
    parser.add_argument("--dec", type=float, help="Declination in ICRS degrees")
    parser.add_argument("--radius", type=float, default=2.0, help="Circular aperture radius in pixels (default: 2)")
    parser.add_argument("--out", type=Path, default=Path("data/spectrum.json"))
    parser.add_argument("--auto-target", action="store_true", help="Pick the highest local unflagged peak in first FITS (review it before use)")
    args = parser.parse_args()
    files = sorted({f for p in args.fits for f in (list(p.glob("*.fits")) + list(p.glob("*.fits.gz")) if p.is_dir() else [p])})
    if not files:
        parser.error("No FITS files found")
    if args.auto_target:
        if args.ra is not None or args.dec is not None:
            parser.error("Do not provide coordinates with --auto-target")
        ra, dec, snr = auto_target(files[0])
        print(f"Auto-selected review position: RA={ra:.7f} deg, Dec={dec:.7f} deg (local S/N {snr:.1f})")
    else:
        if args.ra is None or args.dec is None:
            parser.error("Provide --ra and --dec, or use --auto-target")
        ra, dec = args.ra, args.dec
    measurements = []
    for path in files:
        try:
            result = measure(path, ra, dec, args.radius)
            if result:
                measurements.append(result)
        except Exception as exc:
            print(f"Skipped {path.name}: {exc}")
    measurements.sort(key=lambda m: m["wavelength_um"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps({
        "product": "NASA SPHEREx DR1 Level-2 infrared screening spectrum",
        "mission": "NASA SPHEREx",
        "release": "DR1",
        "target": {"ra_deg": ra, "dec_deg": dec},
        "measurements": measurements,
        "units": {"wavelength": "um", "brightness": "MJy/sr"},
        "method": "circular-aperture mean minus median annulus; reject artifact flags but permit the informational SOURCE flag; variance-based uncertainty when available",
        "caveat": "Exploratory screening only. Validate with IRSA spectral-WCS calibration products, PSF/ePSF photometry, and pipeline documentation.",
    }, indent=2), encoding="utf-8")
    print(f"Wrote {len(measurements)} measurements to {args.out}")
    if len(measurements) < 2:
        print("Need at least two overlapping usable observations to draw a spectrum.")


if __name__ == "__main__":
    main()
