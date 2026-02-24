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
from uuid import uuid4

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


def downsample_2x(img: np.ndarray) -> np.ndarray:
    """Block-mean 2x downsample, trimming odd dimensions."""
    h, w = img.shape
    img = img[: h - h % 2, : w - w % 2]
    return img.reshape(h // 2, 2, w // 2, 2).mean(axis=(1, 3)).astype(np.uint16)


def build_ome_xml(
    channel_names: list[str],
    size_y: int,
    size_x: int,
    pixel_size: float,
) -> str:
    """Build OME-XML metadata string."""
    uuid = f"urn:uuid:{uuid4()}"
    channels_xml = "\n".join(
        f'<Channel ID="Channel:0:{i}" Name="{name}" SamplesPerPixel="1"/>'
        for i, name in enumerate(channel_names)
    )
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06"
     xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
     UUID="{uuid}"
     xsi:schemaLocation="http://www.openmicroscopy.org/Schemas/OME/2016-06
     http://www.openmicroscopy.org/Schemas/OME/2016-06/ome.xsd">
  <Image ID="Image:0" Name="default">
    <Pixels BigEndian="false"
            DimensionOrder="XYCZT"
            ID="Pixels:0"
            Interleaved="false"
            PhysicalSizeX="{pixel_size}"
            PhysicalSizeXUnit="&#181;m"
            PhysicalSizeY="{pixel_size}"
            PhysicalSizeYUnit="&#181;m"
            SizeC="{len(channel_names)}"
            SizeT="1"
            SizeX="{size_x}"
            SizeY="{size_y}"
            SizeZ="1"
            Type="uint16">
      {channels_xml}
      <TiffData/>
    </Pixels>
  </Image>
</OME>"""
    return xml


def convert(input_dir: Path, output_path: Path, pixel_size: float) -> None:
    """Convert single-channel TIFs to pyramidal OME-TIFF."""
    channel_paths = discover_channels(input_dir)
    names = [channel_name(p) for p in channel_paths]
    n_channels = len(names)

    print(f"Found {n_channels} channels: {', '.join(names)}")

    # Read first channel to get dimensions
    first = tifffile.imread(str(channel_paths[0]))
    size_y, size_x = first.shape
    print(f"Image dimensions: {size_y} x {size_x} (Y x X)")
    del first

    # Number of sub-resolution levels (4x, 8x, 16x, 32x = 4 SubIFDs)
    n_subifds = 4

    ome_xml = build_ome_xml(names, size_y, size_x, pixel_size)

    print(f"Writing {output_path} ...")
    with tifffile.TiffWriter(str(output_path), bigtiff=True) as tw:
        for i, (path, name) in enumerate(zip(channel_paths, names)):
            print(f"  [{i+1}/{n_channels}] {name} ...", end=" ", flush=True)
            img = tifffile.imread(str(path)).astype(np.uint16)

            # Write full-resolution page
            options = dict(
                tile=(512, 512),
                compression="deflate",
                subifds=n_subifds,
            )
            # Attach OME-XML description only to the first page
            if i == 0:
                options["description"] = ome_xml
                options["metadata"] = None  # prevent tifffile auto-metadata

            tw.write(img, **options)

            # Generate and write pyramid levels as SubIFDs
            # Level 1: 4x downsample (2x twice)
            sub = downsample_2x(downsample_2x(img))
            del img  # free full-res memory

            subifd_options = dict(
                tile=(512, 512),
                compression="deflate",
                subfiletype=1,
            )

            tw.write(sub, **subifd_options)

            # Levels 2-4: each 2x from previous
            for _ in range(3):
                sub = downsample_2x(sub)
                tw.write(sub, **subifd_options)

            print("done")

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
