"""Shared data loading, genotype mapping, and analysis utilities.

Used by all H1–H10 hypothesis notebooks.
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

CODEX_SAMPLES = {"1901HBMP004", "HDL073", "HBMP006", "HBMP007", "HDL098"}
PHENOCYCLER_SAMPLES = {"HDL011", "HDL043", "HDL052", "HDL053", "HDL055",
                       "HDL063", "HDL070", "HDL079", "HDL086", "HDL094"}
CODEX_CSV = PROJECT / "Measurements" / "ForSH2B3.csv"

# CODEX region/class harmonization
_CODEX_REGION_MAP = {"Red_Pulp": "RedPulp", "Sinusoid": "RedPulp",
                     "Trabecula": "Trabeculae", "Peripheral_White_Pulp": "PALS"}
_CODEX_CLASS_MAP = {**_CODEX_REGION_MAP, "Vessel": "SmallVessel",
                    "Vein": "LargeVessel", "Arteriole": "LargeVessel"}

# Clinical / SNP constants
CLINICAL_COLS = ["Age (yrs)", "Gender", "Ethnicity", "C-pep (ng/ml)", "HbA1c", "BMI"]
CONTINUOUS_CLINICAL = ["Age (yrs)", "C-pep (ng/ml)", "HbA1c", "BMI"]

# ALT ID → HANDEL ID mapping (from Groups.xlsx)
_ALT_TO_HANDEL = {
    "1901": "1901HBMP004",
    "1902": "HDL073",
    "1903": "HDL075",
    "1904": "HDL070",
    "2006": "HBMP006",
    "2007": "HBMP007",
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
    # Try bare 4-digit ALT ID — resolve through _ALT_TO_HANDEL
    m = re.match(r"(\d{4})", image_name)
    if m:
        prefix = image_name.split("_")[0]
        mapped = _ALT_TO_HANDEL.get(prefix)
        return mapped if mapped else prefix
    return image_name


# ---------------------------------------------------------------------------
# Genotype mapping
# ---------------------------------------------------------------------------
def _normalize_genotype(g):
    """Normalize a genotype string: T/C → C/T, G A → A/G."""
    if pd.isna(g):
        return np.nan
    s = str(g).strip()
    parts = re.split(r"[\s/]+", s)
    if len(parts) == 2:
        return "/".join(sorted(parts))
    return s


def _build_genotype_map() -> dict:
    """Build sample_id → genotype dictionary from Groups.xlsx."""
    g = pd.read_excel(GROUPS_XLSX)
    geno_col = "rs3184504" if "rs3184504" in g.columns else "rs3184504 (SH2B3)"
    geno_map = {}
    for _, row in g.iterrows():
        geno = row.get(geno_col)
        if pd.isna(geno):
            continue
        geno = _normalize_genotype(geno)
        # Map by HANDEL ID
        handel = row.get("HANDEL ID")
        if pd.notna(handel):
            geno_map[str(handel).strip()] = geno
        # Map by ALT ID patterns (Groups.xlsx uses 5-digit, e.g. 19001 → key 1901)
        alt = row.get("ALT ID")
        if pd.notna(alt):
            alt_str = str(int(alt))
            # Try direct lookup first, then 4-digit shorthand (first2 + last2)
            mapped = _ALT_TO_HANDEL.get(alt_str)
            if not mapped and len(alt_str) == 5:
                short = alt_str[:2] + alt_str[3:]
                mapped = _ALT_TO_HANDEL.get(short)
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


def load_all_data() -> pd.DataFrame:
    """Load AnnotationsFinal.csv + ForSH2B3.csv (CODEX), harmonized.

    Returns combined DataFrame with Sample, Genotype, Platform columns.
    """
    df = load_data()
    df["Platform"] = "Phenocycler"

    if CODEX_CSV.exists():
        codex = pd.read_csv(CODEX_CSV)
        codex["Classification"] = codex["Classification"].map(
            lambda c: _CODEX_CLASS_MAP.get(c, c))
        codex["Parent"] = codex["Parent"].str.replace(
            r"Annotation \((\w+)\)",
            lambda m: f"Annotation ({_CODEX_REGION_MAP.get(m.group(1), m.group(1))})",
            regex=True,
        )
        codex["Sample"] = codex["Image"].apply(extract_sample_id)
        codex["Genotype"] = codex["Sample"].map(GENOTYPE_MAP)
        codex["Platform"] = "CODEX"
        codex = codex[~codex["Sample"].isin(EXCLUDE_SAMPLES)]
        codex = codex.dropna(subset=["Genotype"])
        codex["Genotype"] = pd.Categorical(
            codex["Genotype"], categories=GENO_ORDER, ordered=True)
        df = pd.concat([df, codex], ignore_index=True)
        df["Genotype"] = pd.Categorical(
            df["Genotype"], categories=GENO_ORDER, ordered=True)

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


# ---------------------------------------------------------------------------
# Clinical data loading (H8, H10)
# ---------------------------------------------------------------------------
def _resolve_sample_from_row(row) -> str | None:
    """Resolve a Groups.xlsx row to a canonical sample ID."""
    handel = row.get("HANDEL ID")
    if pd.notna(handel):
        return str(handel).strip()
    alt = row.get("ALT ID")
    if pd.notna(alt):
        alt_str = str(int(alt))
        mapped = _ALT_TO_HANDEL.get(alt_str)
        if not mapped and len(alt_str) == 5:
            short = alt_str[:2] + alt_str[3:]
            mapped = _ALT_TO_HANDEL.get(short)
        if mapped:
            return mapped
    return None


def load_clinical() -> pd.DataFrame:
    """Load clinical metadata from Groups.xlsx.

    Returns DataFrame with columns: Sample, Genotype, Platform,
    Age (yrs), Gender, Ethnicity, C-pep (ng/ml), HbA1c, BMI.
    """
    g = pd.read_excel(GROUPS_XLSX)
    g["Sample"] = g.apply(_resolve_sample_from_row, axis=1)
    g = g.dropna(subset=["Sample"])
    g = g[~g["Sample"].isin(EXCLUDE_SAMPLES)]

    # Genotype
    geno_col = "rs3184504" if "rs3184504" in g.columns else "rs3184504 (SH2B3)"
    g["Genotype"] = g[geno_col].apply(_normalize_genotype)
    g = g.dropna(subset=["Genotype"])
    g["Genotype"] = pd.Categorical(g["Genotype"], categories=GENO_ORDER, ordered=True)

    # Platform
    g["Platform"] = g["Sample"].apply(
        lambda s: "CODEX" if s in CODEX_SAMPLES else "Phenocycler")

    # Clean clinical columns
    out = g[["Sample", "Genotype", "Platform"] + CLINICAL_COLS].copy()
    for col in CONTINUOUS_CLINICAL:
        out[col] = pd.to_numeric(
            out[col].where(out[col] != "No Data"), errors="coerce")
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# SNP panel loading (H9)
# ---------------------------------------------------------------------------
def load_snp_panel(project_samples=None):
    """Load SNP genotype and risk-dosage matrices from Groups.xlsx.

    Returns (geno_df, dosage_df):
      - geno_df: Sample × SNP matrix of normalized genotype strings
      - dosage_df: Sample × risk-annotated-SNP matrix of 0/1/2 dosage
    """
    g = pd.read_excel(GROUPS_XLSX)
    g["Sample"] = g.apply(_resolve_sample_from_row, axis=1)
    g = g.dropna(subset=["Sample"])
    g = g[~g["Sample"].isin(EXCLUDE_SAMPLES)]

    snp_cols = [c for c in g.columns
                if c.startswith("rs") and "_risk_het_prot" not in c]
    risk_cols = {c.replace("_risk_het_prot", ""): c
                 for c in g.columns if c.endswith("_risk_het_prot")}

    # Build genotype matrix
    geno_records = {}
    for _, row in g.iterrows():
        sid = row["Sample"]
        if project_samples and sid not in project_samples:
            continue
        rec = {}
        for snp in snp_cols:
            val = row[snp]
            rec[snp] = _normalize_genotype(val) if pd.notna(val) else np.nan
        geno_records[sid] = rec
    geno_df = pd.DataFrame.from_dict(geno_records, orient="index")
    geno_df.index.name = "Sample"

    # Build dosage matrix from risk annotations
    dosage_map = {"Protective": 0, "Het": 1, "Risk": 2}
    dosage_records = {}
    for _, row in g.iterrows():
        sid = row["Sample"]
        if project_samples and sid not in project_samples:
            continue
        rec = {}
        for snp, risk_col in risk_cols.items():
            val = row.get(risk_col)
            if pd.notna(val) and val in dosage_map:
                rec[snp] = dosage_map[val]
            else:
                rec[snp] = np.nan
        dosage_records[sid] = rec
    dosage_df = pd.DataFrame.from_dict(dosage_records, orient="index")
    dosage_df.index.name = "Sample"

    # Drop monomorphic SNPs from dosage
    polymorphic = dosage_df.columns[dosage_df.nunique() > 1]
    dosage_df = dosage_df[polymorphic]

    return geno_df, dosage_df


# ---------------------------------------------------------------------------
# Feature matrix construction (H8, H9, H10)
# ---------------------------------------------------------------------------
def build_feature_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Build per-sample morphological feature matrix from annotation data.

    Returns ~23 numeric columns indexed by Sample, plus Genotype and Platform.
    """
    density = compute_density(df)
    regions = get_regions(df)
    vessels = get_vessels(df)

    features = {}

    # --- Vessel density per region (raw + normalized) ---
    for region in MAIN_REGIONS:
        rd = density[density["Region"] == region]
        features[f"{region}_density"] = rd.set_index("Sample")["Density_per_mm2"]
        if region != "RedPulp":
            features[f"{region}_norm_density"] = rd.set_index("Sample")["Density_Normalized"]

    # --- Follicle metrics ---
    follicles = regions[regions["Classification"] == "Follicle"]
    fol_per_img = follicles.groupby(["Sample"], observed=True).agg(
        Follicle_count=("Object ID", "count"),
        Follicle_total_area=("Area µm^2", "sum"),
        Follicle_mean_area=("Area µm^2", "mean"),
    )
    for col in fol_per_img.columns:
        features[col] = fol_per_img[col]

    # Follicle fraction of total tissue
    total_area = regions[regions["Classification"].isin(MAIN_REGIONS)].groupby(
        "Sample", observed=True)["Area µm^2"].sum()
    fol_area = follicles.groupby("Sample", observed=True)["Area µm^2"].sum()
    features["Follicle_fraction"] = fol_area / total_area

    # --- Tissue proportions ---
    region_areas = (
        regions[regions["Classification"].isin(MAIN_REGIONS)]
        .groupby(["Sample", "Classification"], observed=True)["Area µm^2"]
        .sum()
        .unstack(fill_value=0)
    )
    region_total = region_areas.sum(axis=1)
    for r in MAIN_REGIONS:
        if r in region_areas.columns:
            features[f"{r}_fraction"] = region_areas[r] / region_total

    # --- White pulp metrics ---
    wp_area = pd.Series(0.0, index=region_total.index)
    for r in ["Follicle", "PALS"]:
        if r in region_areas.columns:
            wp_area = wp_area + region_areas[r]
    features["WP_fraction"] = wp_area / region_total
    if "Follicle" in region_areas.columns and "PALS" in region_areas.columns:
        features["Follicle_PALS_ratio"] = region_areas["Follicle"] / region_areas["PALS"].replace(0, np.nan)

    # --- Vessel morphology (median per sample) ---
    sv = df[df["Classification"] == "SmallVessel"].copy()
    if not sv.empty:
        sv["Elongation"] = sv["Max diameter µm"] / sv["Min diameter µm"].replace(0, np.nan)
        morph_agg = sv.groupby("Sample", observed=True).agg(
            Vessel_median_area=("Area µm^2", "median"),
            Vessel_median_circularity=("Circularity", "median"),
            Vessel_median_solidity=("Solidity", "median"),
            Vessel_median_elongation=("Elongation", "median"),
        )
        for col in morph_agg.columns:
            features[col] = morph_agg[col]

    # Combine all features
    feat_df = pd.DataFrame(features)
    feat_df.index.name = "Sample"

    # Add Genotype and Platform
    sample_meta = df.drop_duplicates("Sample").set_index("Sample")[["Genotype"]]
    feat_df = feat_df.join(sample_meta)
    feat_df["Platform"] = feat_df.index.map(
        lambda s: "CODEX" if s in CODEX_SAMPLES else "Phenocycler")

    return feat_df


# ---------------------------------------------------------------------------
# Polygenic risk score (H9)
# ---------------------------------------------------------------------------
def compute_prs(dosage_df: pd.DataFrame, min_snps: int = 10) -> pd.Series:
    """Compute unweighted polygenic risk score from risk dosage matrix.

    PRS = sum(dosage) / (2 × N_genotyped), normalized to [0,1].
    Returns NaN for samples with fewer than min_snps genotyped.
    """
    n_genotyped = dosage_df.notna().sum(axis=1)
    total_dosage = dosage_df.sum(axis=1)
    prs = total_dosage / (2 * n_genotyped)
    prs[n_genotyped < min_snps] = np.nan
    prs.name = "PRS"
    return prs


# ---------------------------------------------------------------------------
# Partial correlation (H8)
# ---------------------------------------------------------------------------
def partial_spearman(data: pd.DataFrame, x_col: str, y_col: str,
                     covariate_cols: list) -> tuple:
    """Partial Spearman correlation controlling for covariates.

    Rank-transforms all variables, regresses out covariates via OLS residuals,
    then computes Pearson correlation on residuals.
    Returns (rho, p_value).
    """
    from scipy.stats import pearsonr

    cols = [x_col, y_col] + list(covariate_cols)
    sub = data[cols].dropna()
    if len(sub) < 4:
        return np.nan, np.nan

    # Rank-transform
    ranked = sub.rank()

    # Regress out covariates
    if covariate_cols:
        C = ranked[covariate_cols].values
        C = np.column_stack([C, np.ones(len(C))])
        for col in [x_col, y_col]:
            y_vec = ranked[col].values
            beta, _, _, _ = np.linalg.lstsq(C, y_vec, rcond=None)
            ranked[col] = y_vec - C @ beta

    rho, p = pearsonr(ranked[x_col], ranked[y_col])
    return rho, p


# ---------------------------------------------------------------------------
# Platform diagnostic (H8, H10)
# ---------------------------------------------------------------------------
def platform_diagnostic(feature_df: pd.DataFrame,
                        feature_cols: list) -> pd.DataFrame:
    """Mann-Whitney test for CODEX vs Phenocycler on each feature.

    Returns DataFrame with Feature, U, p, rank_biserial, n_CODEX, n_PC.
    """
    results = []
    for feat in feature_cols:
        codex = feature_df.loc[feature_df["Platform"] == "CODEX", feat].dropna()
        pc = feature_df.loc[feature_df["Platform"] == "Phenocycler", feat].dropna()
        if len(codex) < 1 or len(pc) < 1:
            results.append({"Feature": feat, "U": np.nan, "p": np.nan,
                            "rank_biserial": np.nan, "n_CODEX": len(codex), "n_PC": len(pc)})
            continue
        u, p = mannwhitneyu(codex, pc, alternative="two-sided")
        r = 1 - 2 * u / (len(codex) * len(pc))
        results.append({"Feature": feat, "U": u, "p": p,
                        "rank_biserial": r, "n_CODEX": len(codex), "n_PC": len(pc)})
    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Platform residualization (H11)
# ---------------------------------------------------------------------------
def residualize_platform(feature_df: pd.DataFrame,
                         feature_cols: list) -> pd.DataFrame:
    """Remove platform effect from features via OLS residuals.

    For each feature: fit y ~ platform_indicator, return residuals + grand mean.
    Preserves genotype variance while removing systematic platform offset.
    """
    out = feature_df.copy()
    platform_ind = (out["Platform"] == "CODEX").astype(float).values
    X = np.column_stack([np.ones(len(platform_ind)), platform_ind])
    for col in feature_cols:
        y = out[col].values.astype(float)
        mask = np.isfinite(y)
        if mask.sum() < 3:
            continue
        beta, _, _, _ = np.linalg.lstsq(X[mask], y[mask], rcond=None)
        residuals = y.copy()
        residuals[mask] = y[mask] - X[mask] @ beta + beta[0]  # residual + intercept
        out[col] = residuals
    return out


# ---------------------------------------------------------------------------
# Cell region assignment by signed distance (H11)
# ---------------------------------------------------------------------------
def assign_region_by_distance(dist_follicle, dist_pals, dist_redpulp):
    """Assign cells to the region they are deepest inside.

    Takes three arrays of signed distances (negative = inside).
    Returns array of region labels: Follicle, PALS, RedPulp, or Other.
    """
    regions = np.array(["Follicle", "PALS", "RedPulp"])
    dists = np.column_stack([dist_follicle, dist_pals, dist_redpulp])
    min_idx = np.argmin(dists, axis=1)
    min_val = np.min(dists, axis=1)
    result = regions[min_idx]
    result[min_val >= 0] = "Other"
    return result


# ---------------------------------------------------------------------------
# KDE-based region boundary utilities (H14)
# ---------------------------------------------------------------------------
def compute_marker_kde(coords, x_range, y_range, grid_spacing=10, bandwidth=75):
    """Bin cell coordinates onto a grid and smooth with Gaussian filter.

    Parameters
    ----------
    coords : (N, 2) array of (x, y) positions in µm
    x_range : (x_min, x_max)
    y_range : (y_min, y_max)
    grid_spacing : bin size in µm (default 10)
    bandwidth : Gaussian sigma in µm (default 75)

    Returns
    -------
    density : 2-D array (ny, nx) of smoothed cell density
    x_edges : 1-D array of bin edges along x
    y_edges : 1-D array of bin edges along y
    """
    from scipy.ndimage import gaussian_filter

    x_edges = np.arange(x_range[0], x_range[1] + grid_spacing, grid_spacing)
    y_edges = np.arange(y_range[0], y_range[1] + grid_spacing, grid_spacing)
    counts, _, _ = np.histogram2d(
        coords[:, 0], coords[:, 1], bins=[x_edges, y_edges])
    # histogram2d returns (nx, ny); transpose to (ny, nx) for image convention
    counts = counts.T
    sigma = bandwidth / grid_spacing
    density = gaussian_filter(counts.astype(float), sigma=sigma)
    return density, x_edges, y_edges


def extract_region_polygons(density, threshold_frac, grid_spacing, x_min, y_min,
                            min_area=8400):
    """Extract contours from density grid and return Shapely polygons.

    Parameters
    ----------
    density : 2-D smoothed density array (ny, nx)
    threshold_frac : fraction of max density for contouring (e.g. 0.10)
    grid_spacing : µm per grid cell
    x_min, y_min : origin offset in µm
    min_area : discard polygons smaller than this (µm², default 8400)

    Returns
    -------
    list of shapely.geometry.Polygon in µm coordinates
    """
    from skimage.measure import find_contours
    from shapely.geometry import Polygon
    from shapely.validation import make_valid

    threshold = threshold_frac * density.max()
    if threshold <= 0:
        return []

    contours = find_contours(density, threshold)
    polygons = []
    for contour in contours:
        if len(contour) < 4:
            continue
        # contour is (row, col) → convert to (x_µm, y_µm)
        coords_um = np.column_stack([
            contour[:, 1] * grid_spacing + x_min,
            contour[:, 0] * grid_spacing + y_min,
        ])
        poly = Polygon(coords_um)
        if not poly.is_valid:
            poly = make_valid(poly)
        if poly.is_empty or poly.geom_type not in ('Polygon', 'MultiPolygon'):
            continue
        if poly.geom_type == 'MultiPolygon':
            for p in poly.geoms:
                if p.area >= min_area:
                    polygons.append(p)
        elif poly.area >= min_area:
            polygons.append(poly)
    return polygons


def assign_objects_to_polygons(coords, polygon_dict):
    """Assign point coordinates to regions via spatial index.

    Parameters
    ----------
    coords : (N, 2) array of (x, y) in µm
    polygon_dict : {region_name: list_of_polygons}

    Returns
    -------
    labels : array of region name strings (len N), "Unassigned" if outside all
    """
    from shapely.geometry import Point
    from shapely import STRtree

    n = len(coords)
    labels = np.full(n, "Unassigned", dtype=object)

    # Build flat list of (polygon, region_name) and STRtree
    all_polys = []
    poly_labels = []
    for region, polys in polygon_dict.items():
        for p in polys:
            all_polys.append(p)
            poly_labels.append(region)

    if not all_polys:
        return labels

    tree = STRtree(all_polys)

    # Vectorized query: predicate="within" means point.within(polygon)
    # (shapely 2.x evaluates query_geom.predicate(tree_geom))
    points = np.array([Point(x, y) for x, y in coords], dtype=object)
    pt_idx, tree_idx = tree.query(points, predicate="within")
    for pi, ti in zip(pt_idx, tree_idx):
        labels[pi] = poly_labels[ti]

    return labels
