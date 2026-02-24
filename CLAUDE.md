# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

QuPath (v0.6.0) project for quantitative analysis of multiplex Phenocycler (formerly CODEX) human spleen immunofluorescence images. The project segments tissue regions and small vessels via pixel classifiers, detects cells with InstanSeg deep learning, and exports measurements for downstream statistical analysis.

**Primary analysis goal:** Compare SmallVessel density across tissue regions (Parent column), grouped by genotype from `Groups.xlsx` column `rs3184504 (SH2B3)`. The primary region of interest is **Follicle**.

## Data Architecture

- **Images/** — 9 OME-TIFF multiplex images (~270 GB total, 25-33 channels, 16-bit, ~0.508 µm/pixel). Sample IDs: HDL011, HDL018, HDL021, HDL043, HDL053, HDL055, HDL063, HDL086, HDL172
- **data/[1-9]/** — QuPath per-image output (data.qpdata, summary.json, server.json, thumbnail.jpg)
- **Measurements/AllAnnotations.csv** — Primary analysis input (~422K rows). Columns: `Image, Object ID, Object type, Name, Classification, Parent, ROI, Centroid X µm, Centroid Y µm, Area µm², Perimeter µm, Num Detections, Length µm, Circularity, Solidity, Max diameter µm, Min diameter µm`
- **Groups.xlsx** — Sample metadata/genotype groupings; key column: `rs3184504 (SH2B3)`
- **classifiers/** — Pixel classifier JSONs (`SpleenRegions2.json`, `SmallVessels2.json`) and class definitions (`classes.json`)
- **scripts/Workflow.groovy** — QuPath analysis pipeline: region segmentation → vessel detection → shape measurements → InstanSeg cell detection → distance calculations
- **resources/display/** — Channel visualization configs

## Tissue Region Classes

Defined in `classifiers/classes.json`: **Follicle**, **PALS** (periarteriolar lymphoid sheath), **RedPulp**, **Trabeculae**, **LargeVessel**, **SmallVessel**, **Ignore***

## Key Data Relationships

- In `AllAnnotations.csv`, the `Parent` column encodes the annotation hierarchy. SmallVessel annotations are children of region annotations (e.g., `Annotation (Follicle)`, `Annotation (RedPulp)`).
- To compute vessel density per region: count SmallVessels per Parent region, normalize by parent region area.
- Image names in the CSV (e.g., `HDL011_PC33.ome.tiff`) must be mapped to sample IDs in `Groups.xlsx` by extracting the HDL### prefix.

## Analysis Stack

- **QuPath** (Groovy scripts) for image processing and segmentation
- **Python/Jupyter** for downstream data analysis and visualization
- **Key Python libraries:** pandas, openpyxl (for Groups.xlsx), scipy/statsmodels (statistics), matplotlib/seaborn/plotly (visualization), scikit-learn (clustering for single-cell data)

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
