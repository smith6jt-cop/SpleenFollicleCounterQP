"""Shared data loading, genotype mapping, and analysis utilities.

Used by all H1–H5 hypothesis notebooks.
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial import cKDTree
from scipy.stats import kruskal, mannwhitneyu, spearmanr

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------
PROJECT = Path(__file__).resolve().parent.parent
DATA_CSV = PROJECT / "Measurements" / "AnnotationsFinal.csv"
GROUPS_XLSX = PROJECT / "Groups.xlsx"
FIGURES_DIR = PROJECT / "analysis" / "figures"
TABLES_DIR = PROJECT / "analysis" / "tables"
FIGURES_DIR.mkdir(exist_ok=True)
TABLES_DIR.mkdir(exist_ok=True)

GENO_ORDER = ["C/C", "C/T", "T/T"]
GENO_PALETTE = dict(zip(GENO_ORDER, sns.color_palette("Set2", 3)))
MAIN_REGIONS = ["Follicle", "PALS", "RedPulp", "Trabeculae"]
EXCLUDE_SAMPLES = {"HDL018", "HDL021", "HDL172"}

# ALT ID → HANDEL ID mapping (from Groups.xlsx)
_ALT_TO_HANDEL = {
    "1901": "1901HBMP004",
    "1902": "HDL073",
    "1903": "HDL075",
    "1904": "HDL070",
    "2006": None,
    "2007": None,
    "2008": "HDL098",
}


# ---------------------------------------------------------------------------
# Sample ID extraction
# ---------------------------------------------------------------------------
def extract_sample_id(image_name: str) -> str:
    """Extract a canonical sample ID from an image filename.

    Handles patterns like:
      - HDL011_PC33.ome.tiff → HDL011
      - HDL052SPLN_2025Aug6_Scan1.er.qptiff - resolution #1 → HDL052
      - 1901HBMP004_PC29.ome.tiff → 1901HBMP004
      - HDL073_PC29.ome.tiff → HDL073
    """
    # Try HDL### pattern first
    m = re.match(r"(HDL\d+)", image_name)
    if m:
        return m.group(1)
    # Try 1901HBMP### pattern
    m = re.match(r"(\d{4}HBMP\d+)", image_name)
    if m:
        return m.group(1)
    # Try bare 4-digit ALT ID
    m = re.match(r"(\d{4})", image_name)
    if m:
        return image_name.split("_")[0]
    return image_name


# ---------------------------------------------------------------------------
# Genotype mapping
# ---------------------------------------------------------------------------
def _build_genotype_map() -> dict:
    """Build sample_id → genotype dictionary from Groups.xlsx."""
    g = pd.read_excel(GROUPS_XLSX)
    geno_map = {}
    for _, row in g.iterrows():
        geno = row.get("rs3184504 (SH2B3)")
        if pd.isna(geno):
            continue
        # Map by HANDEL ID
        handel = row.get("HANDEL ID")
        if pd.notna(handel):
            geno_map[str(handel).strip()] = geno
        # Map by ALT ID patterns
        alt = row.get("ALT ID")
        if pd.notna(alt):
            alt_str = str(int(alt))
            mapped = _ALT_TO_HANDEL.get(alt_str)
            if mapped:
                geno_map[mapped] = geno
    return geno_map


GENOTYPE_MAP = _build_genotype_map()


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
def load_data() -> pd.DataFrame:
    """Load AnnotationsFinal.csv with Sample and Genotype columns.

    Drops excluded samples and rows without genotype.
    """
    df = pd.read_csv(DATA_CSV)
    df["Sample"] = df["Image"].apply(extract_sample_id)
    df["Genotype"] = df["Sample"].map(GENOTYPE_MAP)
    # Drop exclusions and unmapped
    df = df[~df["Sample"].isin(EXCLUDE_SAMPLES)]
    df = df.dropna(subset=["Genotype"])
    df["Genotype"] = pd.Categorical(df["Genotype"], categories=GENO_ORDER, ordered=True)
    return df


# ---------------------------------------------------------------------------
# Filtering helpers
# ---------------------------------------------------------------------------
def get_regions(df: pd.DataFrame) -> pd.DataFrame:
    """Return region-level annotations (Follicle, PALS, RedPulp, Trabeculae, LargeVessel)."""
    region_classes = {"Follicle", "PALS", "RedPulp", "Trabeculae", "LargeVessel"}
    return df[df["Classification"].isin(region_classes)].copy()


def get_vessels(df: pd.DataFrame) -> pd.DataFrame:
    """Return SmallVessel annotations with a Region column parsed from Parent."""
    v = df[df["Classification"] == "SmallVessel"].copy()
    v["Region"] = v["Parent"].str.extract(r"Annotation \((\w+)\)")
    return v


# ---------------------------------------------------------------------------
# Density computation
# ---------------------------------------------------------------------------
def compute_density(df: pd.DataFrame) -> pd.DataFrame:
    """Compute vessel density per image per region.

    Returns DataFrame with columns:
      Image, Sample, Genotype, Region, Vessel_Count, Region_Area_mm2,
      Density_per_mm2, RedPulp_Density, Density_Normalized
    """
    regions = get_regions(df)
    vessels = get_vessels(df)

    # Region areas
    area = (
        regions[regions["Classification"].isin(MAIN_REGIONS)]
        .groupby(["Image", "Sample", "Genotype", "Classification"], observed=True)["Area µm^2"]
        .sum()
        .reset_index()
        .rename(columns={"Classification": "Region", "Area µm^2": "Region_Area_um2"})
    )
    area["Region_Area_mm2"] = area["Region_Area_um2"] / 1e6

    # Vessel counts per region
    counts = (
        vessels[vessels["Region"].isin(MAIN_REGIONS)]
        .groupby(["Image", "Sample", "Genotype", "Region"], observed=True)
        .size()
        .reset_index(name="Vessel_Count")
    )

    # Merge — left join so regions with 0 vessels appear
    density = area.merge(counts, on=["Image", "Sample", "Genotype", "Region"], how="left")
    density["Vessel_Count"] = density["Vessel_Count"].fillna(0).astype(int)
    density["Density_per_mm2"] = density["Vessel_Count"] / density["Region_Area_mm2"]

    # RedPulp normalization
    rp = density[density["Region"] == "RedPulp"][["Image", "Density_per_mm2"]].rename(
        columns={"Density_per_mm2": "RedPulp_Density"}
    )
    density = density.merge(rp, on="Image", how="left")
    density["Density_Normalized"] = density["Density_per_mm2"] / density["RedPulp_Density"]

    return density


# ---------------------------------------------------------------------------
# Spatial assignment (for H5)
# ---------------------------------------------------------------------------
def assign_vessels_to_follicles(df: pd.DataFrame, image: str) -> pd.DataFrame:
    """Assign follicle-parented SmallVessels to nearest follicle centroid.

    Uses cKDTree for fast nearest-neighbor lookup.
    Returns vessel DataFrame with Follicle_ID and Follicle_Area columns.
    """
    img_df = df[df["Image"] == image]
    follicles = img_df[img_df["Classification"] == "Follicle"].copy()
    vessels = get_vessels(img_df)
    vessels = vessels[vessels["Region"] == "Follicle"].copy()

    if follicles.empty or vessels.empty:
        return pd.DataFrame()

    # Build tree from follicle centroids
    fol_coords = follicles[["Centroid X µm", "Centroid Y µm"]].values
    tree = cKDTree(fol_coords)

    # Query vessel centroids
    ves_coords = vessels[["Centroid X µm", "Centroid Y µm"]].values
    _, indices = tree.query(ves_coords)

    vessels["Follicle_ID"] = follicles["Object ID"].values[indices]
    vessels["Follicle_Area"] = follicles["Area µm^2"].values[indices]

    return vessels


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------
def rank_biserial(x, y):
    """Rank-biserial effect size: r = 1 - 2U/(n1*n2)."""
    u_stat, _ = mannwhitneyu(x, y, alternative="two-sided")
    n1, n2 = len(x), len(y)
    return 1 - 2 * u_stat / (n1 * n2)


def run_kruskal(data: pd.DataFrame, value_col: str, group_col: str = "Genotype"):
    """Run Kruskal-Wallis test across genotype groups.

    Returns (H_statistic, p_value) or (NaN, NaN) if <2 groups have data.
    """
    groups = [g[value_col].dropna().values for _, g in data.groupby(group_col, observed=True)]
    groups = [g for g in groups if len(g) > 0]
    if len(groups) < 2:
        return np.nan, np.nan
    h, p = kruskal(*groups)
    return h, p


def run_pairwise(data: pd.DataFrame, value_col: str, group_col: str = "Genotype"):
    """Run pairwise Mann-Whitney U tests with rank-biserial effect size.

    Returns list of dicts with keys: Comparison, U, p, r, n1, n2.
    """
    results = []
    groups = data.groupby(group_col, observed=True)[value_col]
    group_dict = {name: g.dropna().values for name, g in groups}
    pairs = [("C/C", "C/T"), ("C/C", "T/T"), ("C/T", "T/T")]
    for g1, g2 in pairs:
        x = group_dict.get(g1, np.array([]))
        y = group_dict.get(g2, np.array([]))
        if len(x) < 1 or len(y) < 1:
            results.append({"Comparison": f"{g1} vs {g2}", "U": np.nan, "p": np.nan,
                            "r": np.nan, "n1": len(x), "n2": len(y)})
            continue
        u, p = mannwhitneyu(x, y, alternative="two-sided")
        r = 1 - 2 * u / (len(x) * len(y))
        results.append({"Comparison": f"{g1} vs {g2}", "U": u, "p": p,
                        "r": r, "n1": len(x), "n2": len(y)})
    return results


def run_dosage_trend(data: pd.DataFrame, value_col: str, group_col: str = "Genotype"):
    """Spearman correlation with ordinal genotype (C/C=0, C/T=1, T/T=2).

    Returns (rho, p_value).
    """
    ordinal = data[group_col].map({"C/C": 0, "C/T": 1, "T/T": 2})
    valid = ordinal.notna() & data[value_col].notna()
    if valid.sum() < 3:
        return np.nan, np.nan
    return spearmanr(ordinal[valid], data.loc[valid, value_col])


def full_stats_table(data: pd.DataFrame, value_col: str, label: str = ""):
    """Run all three statistical tests and return a summary DataFrame."""
    h, kw_p = run_kruskal(data, value_col)
    pw = run_pairwise(data, value_col)
    rho, sp_p = run_dosage_trend(data, value_col)

    rows = [{"Test": "Kruskal-Wallis", "Metric": label, "Statistic": h, "p": kw_p, "Effect_Size": ""}]
    for r in pw:
        rows.append({"Test": f"Mann-Whitney ({r['Comparison']})", "Metric": label,
                      "Statistic": r["U"], "p": r["p"], "Effect_Size": f"r={r['r']:.3f}"})
    rows.append({"Test": "Spearman dosage", "Metric": label,
                  "Statistic": rho, "p": sp_p, "Effect_Size": f"rho={rho:.3f}" if not np.isnan(rho) else ""})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------
def setup_style():
    """Configure seaborn/matplotlib defaults."""
    sns.set_theme(style="whitegrid", font_scale=1.1)
    plt.rcParams["figure.dpi"] = 100
    plt.rcParams["savefig.dpi"] = 150
    plt.rcParams["figure.facecolor"] = "white"


def save_figure(fig, name: str, tight=True):
    """Save figure to analysis/figures/ as PNG."""
    path = FIGURES_DIR / f"{name}.png"
    if tight:
        fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    else:
        fig.savefig(path, dpi=150, facecolor="white")
    print(f"Saved: {path.relative_to(PROJECT)}")


def save_table(df: pd.DataFrame, name: str):
    """Save DataFrame to analysis/tables/ as CSV."""
    path = TABLES_DIR / f"{name}.csv"
    df.to_csv(path, index=False)
    print(f"Saved: {path.relative_to(PROJECT)}")
