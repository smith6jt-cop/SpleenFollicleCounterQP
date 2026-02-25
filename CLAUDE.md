# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

QuPath (v0.6.0) project for quantitative analysis of multiplex Phenocycler (formerly CODEX) human spleen immunofluorescence images. The project segments tissue regions and small vessels via pixel classifiers, detects cells with InstanSeg deep learning, and exports measurements for downstream statistical analysis.

**Primary analysis goal:** Compare SmallVessel density across tissue regions (Parent column), grouped by genotype from `Groups.xlsx` column `rs3184504 (SH2B3)`. The primary region of interest is **Follicle**.

## Data Architecture

- **Images/** — 9 OME-TIFF multiplex images (19-33 channels, 16-bit, ~0.508 µm/pixel). Original samples on disk: HDL011, HDL043, HDL052, HDL055, HDL063, HDL086. HiperGator samples: 1901HBMP004 (29ch), HDL073_PC19 (19ch, pre-signal-isolation), HDL073_PC29 (29ch, post-signal-isolation)
- **FromHipergator/** — Raw single-channel TIF files from HiperGator processing (3 directories; `19-002_spleen_CC3-C` skipped)
- **data/[1-9]/** — QuPath per-image output (data.qpdata, summary.json, server.json, thumbnail.jpg)
- **Measurements/AllAnnotations.csv** — Primary analysis input (~422K rows). Columns: `Image, Object ID, Object type, Name, Classification, Parent, ROI, Centroid X µm, Centroid Y µm, Area µm², Perimeter µm, Num Detections, Length µm, Circularity, Solidity, Max diameter µm, Min diameter µm`
- **Groups.xlsx** — Sample metadata/genotype groupings; key column: `rs3184504 (SH2B3)`
- **classifiers/** — Pixel classifier JSONs (`SpleenRegions2.json`, `SmallVessels2.json`) and class definitions (`classes.json`)
- **scripts/Workflow.groovy** — QuPath analysis pipeline: region segmentation → vessel detection → shape measurements → InstanSeg cell detection → distance calculations
- **scripts/convert_to_ome_tiff.py** — Converts single-channel TIF directories into Bio-Formats-compatible pyramidal OME-TIFF (6 levels, 512x512 tiles, ADOBE_DEFLATE, big-endian, per-channel IFD+SubIFDs with manual OME-XML)
- **scripts/process_hdl73_channels.py** — Signal isolation for HDL73: autofluorescence subtraction using matched blank pairs via KINTSUGI `kintsugi.signal` module. Supports `--force` to re-process existing channels.
- **scripts/FastHierarchyAnnotationsDetections.groovy** — QuPath script for fast hierarchy export with annotations and detections
- **scripts/FastHierarchyAnnotationsOnly.groovy** — QuPath script for fast hierarchy export with annotations only
- **analysis/** — Downstream Python analysis outputs (vessel density notebook, figures, CSVs)
- **classifiers/pixel_classifiers/** — Pixel classifier model JSONs (e.g., `SmallVessel3.json`)
- **resources/display/** — Channel visualization configs

## Tissue Region Classes

Defined in `classifiers/classes.json`: **Follicle**, **PALS** (periarteriolar lymphoid sheath), **RedPulp**, **Trabeculae**, **LargeVessel**, **SmallVessel**, **Ignore***

## Key Data Relationships

- In `AllAnnotations.csv`, the `Parent` column encodes the annotation hierarchy. SmallVessel annotations are children of region annotations (e.g., `Annotation (Follicle)`, `Annotation (RedPulp)`).
- To compute vessel density per region: count SmallVessels per Parent region, normalize by parent region area.
- Image names in the CSV (e.g., `HDL011_PC33.ome.tiff`) must be mapped to sample IDs in `Groups.xlsx` by extracting the HDL### prefix.

## Sample ID Mapping (FromHipergator)

| HiperGator Directory | ALT ID | HANDEL ID | Output Filename | Channels | Status |
|---|---|---|---|---|---|
| `19-001_SP_CC2-A28` | 1901 | — | `1901HBMP004_PC29.ome.tiff` | 29 | Converted (1.74 GB) |
| `19-002_spleen_CC2-A_D200210` | 1902 | HDL073 | `HDL073_PC19.ome.tiff` | 19 | Converted (2.44 GB) |
| `HDL73_SPL_Processed/ImageJ` | 1902 | HDL073 | `HDL073_PC29.ome.tiff` | 29 | Signal-isolated (3.30 GB) |
| `19-002_spleen_CC3-C` | — | — | *skipped* | — | — |

## Analysis Stack

- **QuPath** (Groovy scripts) for image processing and segmentation
- **Python/Jupyter** for downstream data analysis and visualization
- **Key Python libraries:** pandas, openpyxl (for Groups.xlsx), scipy/statsmodels (statistics), matplotlib/seaborn/plotly (visualization), scikit-learn (clustering for single-cell data), tifffile/numpy (OME-TIFF conversion)

## OME-TIFF Format Specification

Target format for all images (Bio-Formats/QuPath compatible):
- **Pyramidal OME-TIFF** with SubIFDs (6 levels: 1x, 2x, 4x, 8x, 16x, 32x)
- **Tiles:** 512x512, **Compression:** ADOBE_DEFLATE (code 8), **Dtype:** uint16
- **Byte order:** big-endian (`byteorder='>'`)
- **Pixel size:** 0.5077663810243286 µm, **Axes:** CYX (one IFD page per channel)
- **Channel order:** DAPI first, then alphabetical
- **OME-XML:** Manual construction with `BigEndian="true"`, `Interleaved="false"`, `<LightPath/>`, per-channel `<TiffData>` blocks
- **µm encoding:** Raw UTF-8 µ bytes (`\xc2\xb5`), NOT XML entity `&#181;m`. Pass OME-XML as `bytes` to tifffile to bypass ASCII check.

## HDL73 Signal Isolation

**Blank position mapping** (from `FromHipergator/HDL73_SPL_meta/channelnames.txt`):

| Position | Blank Pair | Markers |
|----------|-----------|---------|
| a | Blank1a + Blank13a | CD20, CD31, CD34, CD35, Lyve1, PanCK, SMActin |
| b | Blank1b + Blank13b | CD8, CD15, CD21, CD44, CD45RO, CD5, CollagenIV, ECAD, FoxP3, Ki67, Podoplanin |
| c | Blank1c + Blank13c | CD3e, CD4, CD11c, CD107a, CD163, CD1c, CD45, CD68, HLADR, Vimentin |

- **Failed markers:** CD1c, CD5 (processed but may have poor signal)
- **Low signal preservation warnings:** PanCK (0.14), Podoplanin (0.19), SMActin (0.18)
- **Existing param files** in `FromHipergator/HDL73_SPL_Processed/Processing_parameters/` for: CD11c, CD15, CD1c, CD20, CD21, CD3e, CD4, CD5, CD8
- **DAPI** is copied directly (no subtraction)

## InstanSeg Cell Detection Channels (26)

DAPI, CD45RO, Ki67, FOXP3, CD38, CD20, CD4, CD44, CD31, CD11c, CD34, CD107a, PDL1, CD163, HLA-DR, CD68, CD8, CD21, CD66, CD141, CD57, CD3e, HLA-A, PD-1, CD45, Podoplanin

## Skills_Registry Submodule

A git submodule at `Skills_Registry/` providing cross-session agent memory and reusable skill templates.

- **Source:** https://github.com/smith6jt-cop/Skills_Registry
- **Relevant skill categories:**
  - `plugins/kintsugi/` — 34 microscopy and image analysis skills
  - `plugins/scientific/` — scientific computing and statistics skills
- **Commands:** `/advise` (get skill-informed guidance), `/retrospective` (post-session learning capture)
- **Update submodule:** `git submodule update --remote Skills_Registry`
