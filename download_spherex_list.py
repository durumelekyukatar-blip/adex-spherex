"""Download selected official IRSA FITS URLs from the supplied URL list."""
import argparse
from pathlib import Path
from urllib.request import urlopen


def main():
    parser = argparse.ArgumentParser(description="Download selected SPHEREx FITS files from a URL-list text file.")
    parser.add_argument("url_list", type=Path, help="Text file containing one IRSA URL per line")
    parser.add_argument("--indices", required=True, help="1-based URL indices, e.g. 1,3,5")
    parser.add_argument("--out", type=Path, default=Path("data/raw"))
    args = parser.parse_args()
    urls = [line.strip() for line in args.url_list.read_text(encoding="utf-8-sig").splitlines()
            if line.strip().startswith("https://")]
    indexes = [int(value) for value in args.indices.split(",")]
    args.out.mkdir(parents=True, exist_ok=True)
    for index in indexes:
        if index < 1 or index > len(urls):
            parser.error(f"URL index {index} is outside 1..{len(urls)}")
        url = urls[index - 1]
        target = args.out / url.rsplit("/", 1)[-1]
        if target.exists() and target.stat().st_size > 2880:
            print(f"Already present: {target}")
            continue
        print(f"Downloading {index}/{len(urls)}: {target.name}")
        with urlopen(url, timeout=90) as response, target.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
        print(f"Saved {target} ({target.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
