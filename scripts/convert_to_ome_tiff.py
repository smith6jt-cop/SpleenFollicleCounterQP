#!/usr/bin/env python3
"""Convert a directory of single-channel TIF files into a pyramidal OME-TIFF.

Produces a file matching the format used by existing QuPath project images:
- Pyramidal OME-TIFF with SubIFDs (5 resolution levels: 1x, 4x, 8x, 16x, 32x)
- 512x512 tiles, DEFLATE compression, uint16
- CYX axis order (one page per channel)

Usage:
    python convert_to_ome_tiff.py <input_dir> <output_path> [--pixel-size 0.508]
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import tifffile


def discover_channels(input_dir: Path) -> list[Path]:
    """Find channel TIF files, sort with DAPI first, rename DAPI-01 -> DAPI."""
    tifs = sorted(
        p for p in input_dir.glob("*.tif")
        if p.stem not in ("BLANK",) and not p.stem.startswith("BLANK")
    )
    if not tifs:
        sys.exit(f"No TIF files found in {input_dir}")

    # Sort: DAPI first, then alphabetical
    dapi = [p for p in tifs if p.stem in ("DAPI", "DAPI-01")]
    others = [p for p in tifs if p.stem not in ("DAPI", "DAPI-01")]
    return dapi + others


def channel_name(path: Path) -> str:
    """Extract channel name from filename, normalizing DAPI-01 to DAPI."""
    name = path.stem
    if name == "DAPI-01":
        return "DAPI"
    return name


def downsample_2x_cyx(img: np.ndarray) -> np.ndarray:
    """Block-mean 2x downsample of a 3D (C, Y, X) array, trimming odd dims."""
    c, h, w = img.shape
    img = img[:, : h - h % 2, : w - w % 2]
    return img.reshape(c, h // 2, 2, w // 2, 2).mean(axis=(2, 4)).astype(np.uint16)


def convert(input_dir: Path, output_path: Path, pixel_size: float) -> None:
    """Convert single-channel TIFs to pyramidal OME-TIFF."""
    channel_paths = discover_channels(input_dir)
    names = [channel_name(p) for p in channel_paths]
    n_channels = len(names)

    print(f"Found {n_channels} channels: {', '.join(names)}")

    # Load all channels into a single 3D (C, Y, X) array
    print("Loading channels ...")
    layers = []
    for i, (path, name) in enumerate(zip(channel_paths, names)):
        print(f"  [{i+1}/{n_channels}] {name}")
        layers.append(tifffile.imread(str(path)).astype(np.uint16))
    data = np.stack(layers)
    del layers
    print(f"Array shape: {data.shape} (C, Y, X)")

    # OME metadata — tifffile generates correct OME-XML from this
    metadata = {
        "axes": "CYX",
        "Channel": {"Name": names},
        "PhysicalSizeX": pixel_size,
        "PhysicalSizeXUnit": "um",
        "PhysicalSizeY": pixel_size,
        "PhysicalSizeYUnit": "um",
    }

    print(f"Writing {output_path} ...")
    with tifffile.TiffWriter(str(output_path), ome=True, bigtiff=True) as tw:
        # Full-resolution 3D array with 4 SubIFD levels
        tw.write(
            data,
            tile=(512, 512),
            compression="deflate",
            subifds=4,
            metadata=metadata,
        )

        # Generate and write 4 pyramid sub-levels as SubIFDs
        # Level 1: 4x (two 2x downsamples from full-res)
        sub = downsample_2x_cyx(downsample_2x_cyx(data))
        del data

        subifd_options = dict(
            tile=(512, 512),
            compression="deflate",
            subfiletype=1,
        )

        tw.write(sub, **subifd_options)

        # Levels 2-4: each 2x from previous
        for _ in range(3):
            sub = downsample_2x_cyx(sub)
            tw.write(sub, **subifd_options)

    print(f"Wrote {output_path} ({output_path.stat().st_size / 1e9:.2f} GB)")


def main():
    parser = argparse.ArgumentParser(
        description="Convert single-channel TIFs to pyramidal OME-TIFF"
    )
    parser.add_argument("input_dir", type=Path, help="Directory of channel TIF files")
    parser.add_argument("output_path", type=Path, help="Output OME-TIFF path")
    parser.add_argument(
        "--pixel-size",
        type=float,
        default=0.5077663810243286,
        help="Pixel size in µm (default: 0.5077663810243286)",
    )
    args = parser.parse_args()

    if not args.input_dir.is_dir():
        sys.exit(f"Input directory not found: {args.input_dir}")

    args.output_path.parent.mkdir(parents=True, exist_ok=True)

    convert(args.input_dir, args.output_path, args.pixel_size)


if __name__ == "__main__":
    main()
