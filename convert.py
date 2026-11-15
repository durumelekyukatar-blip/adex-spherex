"""Convert a SPHEREx FITS image plane to a display PNG without altering the FITS."""
import argparse
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from PIL import Image, ImageOps


def convert(source: Path, target: Path) -> dict:
    with fits.open(source, memmap=True) as hdus:
        hdu = hdus["IMAGE"] if "IMAGE" in hdus else next(
            item for item in hdus if item.data is not None and item.data.ndim == 2
        )
        data = np.asarray(hdu.data, dtype=np.float32)
        header = hdu.header
        finite = np.isfinite(data)
        if not finite.any():
            raise ValueError(f"No finite pixels found in {source}")
        clean = np.where(finite, data, np.nanmedian(data[finite]))
        low, high = np.nanpercentile(clean, [1.0, 99.5])
        if not np.isfinite(high - low) or high <= low:
            high = low + 1.0
        linear = np.clip((clean - low) / (high - low), 0, 1)
        pixels = np.asarray(np.arcsinh(10 * linear) / np.arcsinh(10) * 255, dtype=np.uint8)
        image = Image.fromarray(pixels[::-1]).convert("RGB")
        image.thumbnail((2200, 2200), Image.Resampling.LANCZOS)
        image = ImageOps.pad(image, (900, 560), color=(6, 10, 18), centering=(0.5, 0.5))
        target.parent.mkdir(parents=True, exist_ok=True)
        image.save(target)
        return {
            "source_file": source.name,
            "image_hdu": hdu.name or "PRIMARY",
            "shape": list(data.shape),
            "width": int(data.shape[-1]),
            "height": int(data.shape[-2]),
            "bitpix": header.get("BITPIX"),
            "bunit": header.get("BUNIT"),
            "date_obs": header.get("DATE-OBS"),
            "mjd_obs": header.get("MJD-OBS"),
            "filter": header.get("FILTER"),
            "wavelength_um": header.get("WAVELEN"),
            "wcs": {k: header[k] for k in ("CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2", "CTYPE1", "CTYPE2") if k in header},
            "display_png": target.name,
            "stretch": "ZScale + Asinh (display only)",
        }


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert SPHEREx FITS images for the Adex comparison viewer.")
    parser.add_argument("fits", nargs="+", type=Path, help="One FITS file for a real-image preview, or two matched epochs for blink comparison")
    parser.add_argument("--out", type=Path, default=Path("data"), help="Output folder (default: data)")
    args = parser.parse_args()
    if len(args.fits) > 2:
        parser.error("Provide at most two FITS files")
    records = []
    for i, source in enumerate(args.fits):
        target = args.out / ("spherex-preview.png" if len(args.fits) == 1 else ("a.png" if i == 0 else "b.png"))
        records.append(convert(source, target))
        print(f"Wrote {target}")
    (args.out / "meta.json").write_text(json.dumps({"observations": records}, indent=2), encoding="utf-8")
    if len(records) == 1:
        print("Only one epoch supplied. Use a second, separate observation for change detection.")


if __name__ == "__main__":
    main()
