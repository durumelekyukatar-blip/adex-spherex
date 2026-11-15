# ADEX - SPHEREx infrared explorer

## Open the app

1. Extract the ZIP to a folder.
2. Open the extracted folder.
3. Double-click `index.html` (Chrome or Edge recommended).

The observation preview, measured spectrum, and reference images are included as local files. They load without a Python setup or a data download. If your browser restricts local files, open PowerShell in this folder and run `python -m http.server 8000`, then visit `http://localhost:8000`.

## Real infrared observation data

`data/raw/` contains four compressed NASA SPHEREx DR1 Level-2 FITS cutouts downloaded from the supplied IRSA URL list. Each is 0.1 degrees across and preserves image, flags, variance, PSF/ePSF, and spectral WCS products. The embedded graph uses four measurements across two small slews, at approximately 0.906, 0.928, 3.027, and 3.110 microns. Wavelength gaps are not connected.

The review coordinate is approximately RA 347.5607 degrees, Dec -3.6313 degrees. These nearby observations support infrared spectrum inspection; they do not establish source variability or an exoplanet detection.

## Exoplanet lab and visual references

The lab renders a hypothetical teaching exoplanet labelled `SIM-EXO-01 b`. This is not an official catalog name. The model parameters are user assumptions and are separate from the SPHEREx source.

- `data/spherex-cygnus-ice.jpg`: real SPHEREx infrared map of interstellar ice and PAH emission in Cygnus X. This is not an exoplanet image. Credit: NASA/JPL-Caltech/IPAC/Hora et al.
- `data/trappist-1-system-concept.jpg`: NASA artist concept of the real TRAPPIST-1 system, not a direct photograph.
- `data/wasp-39b-nirspec-spectrum.png`: NASA Webb NIRSpec transmission spectrum of WASP-39 b; the 4.3 micron feature indicates carbon dioxide absorption. The spectrum is measured; any planet artwork is illustrative.

The energy-balance climate model is educational. It omits clouds, ocean circulation, atmospheric chemistry, and weather. Its index is not a scientific measure of habitability.

## Add more URL-list observations

`data/spherex_urls.txt` contains the supplied URL list. Full-frame FITS files are about 63-70 MB. To download selected entries and re-extract a screening spectrum:

```powershell
python -m pip install -r requirements.txt
python download_spherex_list.py data/spherex_urls.txt --indices 1,2,3,4 --out data/raw
python extract_spherex_spectrum.py data/raw --auto-target --out data/spectrum.json
```

For a known target, use `--ra <RA-deg> --dec <Dec-deg>` instead of `--auto-target`.

## Scientific interpretation

- SPHEREx Level-2 images are calibrated surface brightness in MJy/sr. The extractor uses aperture and local-background estimates, per-pixel variance, and artifact-flag filtering.
- The in-file wavelength WCS is an approximation for visualization. This is exploratory screening, not precision photometry. Validate with official spectral-WCS calibration products, PSF/ePSF photometry, and pipeline documentation before scientific interpretation.

Official references: [IRSA SPHEREx data products](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/), [IRSA spectral image tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/).
