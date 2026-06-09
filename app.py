
import os
import re
from pathlib import Path
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import streamlit.components.v1 as components

# ============================================================
# Page config
# ============================================================

st.set_page_config(
    page_title="HOF–Protein Compatibility Recommender",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============================================================
# Strong light-mode CSS
# ============================================================

st.markdown("""
<style>
:root { color-scheme: light !important; }

html, body, .stApp {
    background: #ffffff !important;
    color: #111111 !important;
}

html, body, p, div, span, label, h1, h2, h3, h4, h5, h6,
[class*="css"], [class*="st-"] {
    color: #111111 !important;
}

section[data-testid="stSidebar"] {
    background-color: #f7f7f7 !important;
    color: #111111 !important;
}

section[data-testid="stSidebar"] * {
    color: #111111 !important;
}

/* Selectbox/dropdown */
div[data-baseweb="select"],
div[data-baseweb="select"] *,
div[data-baseweb="select"] > div,
div[data-baseweb="popover"],
div[data-baseweb="popover"] *,
div[data-baseweb="menu"],
div[data-baseweb="menu"] *,
ul[role="listbox"],
div[role="listbox"],
li[role="option"],
div[role="option"] {
    background-color: #ffffff !important;
    color: #111111 !important;
    -webkit-text-fill-color: #111111 !important;
    fill: #111111 !important;
}

li[role="option"]:hover,
div[role="option"]:hover,
li[aria-selected="true"],
div[aria-selected="true"] {
    background-color: #e8f0fe !important;
    color: #111111 !important;
}

/* Inputs */
input, textarea {
    background-color: #ffffff !important;
    color: #111111 !important;
    -webkit-text-fill-color: #111111 !important;
    border: 1px solid #999999 !important;
}

button, button *, .stButton button, .stDownloadButton button {
    background-color: #ffffff !important;
    color: #111111 !important;
    -webkit-text-fill-color: #111111 !important;
    border-color: #999999 !important;
}

div[data-testid="stMetric"],
div[data-testid="stMetric"] *,
div[data-testid="stDataFrame"],
div[data-testid="stDataFrame"] *,
div[data-testid="stExpander"],
div[data-testid="stExpander"] * {
    background-color: #ffffff !important;
    color: #111111 !important;
}

button[data-baseweb="tab"] {
    background-color: #ffffff !important;
    color: #111111 !important;
}

button[data-baseweb="tab"][aria-selected="true"] {
    color: #d71920 !important;
    border-bottom: 2px solid #d71920 !important;
}

.card {
    background: #ffffff !important;
    border: 1px solid #dddddd !important;
    border-radius: 14px;
    padding: 18px;
    margin-bottom: 18px;
    box-shadow: 0 1px 5px rgba(0,0,0,0.06);
}

.small {
    color: #555555 !important;
    font-size: 0.9rem;
}

.tag {
    display: inline-block;
    background: #eef1f4 !important;
    color: #111111 !important;
    padding: 5px 10px;
    border-radius: 18px;
    margin: 4px 4px 4px 0;
    font-size: 0.85rem;
}

.good {
    background: #e8f4ec !important;
    color: #176d37 !important;
}

.warn {
    background: #fff3df !important;
    color: #8a5a00 !important;
}

.bad {
    background: #fdeaea !important;
    color: #982626 !important;
}

a { color: #0645ad !important; }
</style>
""", unsafe_allow_html=True)

# ============================================================
# File paths with fallbacks
# ============================================================

def first_existing(paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return None

WORKBOOK_PATH = first_existing([
    "HOF_MASTER_WORKBOOK_functional_group_typed_FIXED.xlsx",
    "HOF_MASTER_WORKBOOK_functional_group_typed_FIXED.xlsx",
    "/content/HOF_MASTER_WORKBOOK_functional_group_typed.xlsx",
    "/content/HOF_MASTER_WORKBOOK_enhanced.xlsx",
])

PROTEIN_PATH = first_existing([
    "proteins100_enhanced.csv",
    "proteins100_enhanced.csv",
])

RECLASSIFIED_WORKBOOK = first_existing([
    "HOF_DATABASE_framework_series_family_reclassified.xlsx",
    "HOF_DATABASE_framework_series_family_reclassified.xlsx",
])

HOF_DB_DIR = first_existing(["hof_db"])

if WORKBOOK_PATH is None:
    st.error("Missing HOF workbook. Upload HOF_MASTER_WORKBOOK_functional_group_typed_FIXED.xlsx.")
    st.stop()

if PROTEIN_PATH is None:
    st.error("Missing protein CSV. Upload proteins100_enhanced.csv.")
    st.stop()

# ============================================================
# Helpers
# ============================================================

def read_sheet(path, sheet):
    try:
        xl = pd.ExcelFile(path)
        if sheet in xl.sheet_names:
            return pd.read_excel(path, sheet_name=sheet)
    except Exception:
        pass
    return pd.DataFrame()

def num(s):
    return pd.to_numeric(s, errors="coerce")

def minmax(s):
    s = num(s)
    if s.notna().sum() == 0:
        return pd.Series(0.0, index=s.index)
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or pd.isna(hi) or hi == lo:
        return pd.Series(0.0, index=s.index)
    return (s - lo) / (hi - lo)

def percentile(s):
    s = num(s)
    if s.notna().sum() == 0:
        return pd.Series(0.0, index=s.index)
    return s.rank(pct=True).fillna(0)

def safe_float(x, default=np.nan):
    try:
        if pd.isna(x):
            return default
        return float(x)
    except Exception:
        return default

def fmt(x, n=3):
    try:
        if pd.isna(x):
            return "NA"
        return f"{float(x):.{n}f}"
    except Exception:
        return "NA"

def render_html_table(table_df, max_height=520):
    """Render a stable light-mode table inside an HTML component."""
    table_df = table_df.copy()
    for column in table_df.columns:
        if pd.api.types.is_float_dtype(table_df[column]):
            table_df[column] = table_df[column].round(4)

    html_table = table_df.to_html(index=False, escape=True, border=0)
    html_doc = f"""
    <!DOCTYPE html>
    <html>
    <head>
    <meta charset="utf-8">
    <style>
      html, body {{ margin: 0; padding: 0; background: #ffffff; color: #111111; }}
      .hof-table-wrap {{
        width: 100%; overflow-x: auto; overflow-y: auto;
        max-height: {int(max_height)}px; border: 1px solid #dddddd;
        border-radius: 10px; background: #ffffff;
        font-family: Arial, Helvetica, sans-serif;
      }}
      table {{
        border-collapse: collapse; width: max-content; min-width: 100%;
        background: #ffffff; color: #111111; font-size: 13px;
      }}
      th {{
        position: sticky; top: 0; background: #f1f3f5; color: #111111;
        padding: 9px; border: 1px solid #dddddd; text-align: left;
        white-space: nowrap; z-index: 2;
      }}
      td {{
        background: #ffffff; color: #111111; padding: 8px;
        border: 1px solid #e3e3e3; white-space: nowrap;
      }}
      tr:nth-child(even) td {{ background: #fafafa; }}
    </style>
    </head>
    <body>
      <div class="hof-table-wrap">{html_table}</div>
    </body>
    </html>
    """
    components.html(html_doc, height=min(int(max_height) + 20, 680), scrolling=True)

def canonical_fallback(x):
    s = str(x)
    s = Path(s).name.replace(".cif", "")
    s = re.sub(r"^\d+[_\-\s]+\d{6,8}[_\-\s]+", "", s)
    s = re.sub(r"^\d+[_\-\s]+", "", s)
    s = re.sub(r"[_\-\s]?\d{2,3}\s*K$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"[_\-\s]?(RT|roomtemp|roomtemperature)$", "", s, flags=re.IGNORECASE)
    s = re.sub(r"[^A-Za-z0-9]+", "", s).lower()
    return s if s else str(x)

def get_series_key(row):
    for c in ["framework_series_key", "framework_series"]:
        if c in row.index and pd.notna(row[c]) and str(row[c]).strip():
            return str(row[c])
    return canonical_fallback(row.get("cif_filename", ""))

def get_series_label(row):
    for c in ["framework_series", "framework_series_key"]:
        if c in row.index and pd.notna(row[c]) and str(row[c]).strip():
            return str(row[c])
    return get_series_key(row)

def get_family_label(row):
    for c in ["family_revised", "family_broad_revised", "family", "family_original"]:
        if c in row.index and pd.notna(row[c]) and str(row[c]).strip():
            return str(row[c])
    return "Unclassified/Other"

def has_literature(row):
    for c in [
        "final_doi", "discovery_publication_doi",
        "discovery_publication_citation", "doi", "DOI",
        "reference", "citation"
    ]:
        if c in row.index and pd.notna(row.get(c)) and str(row.get(c)).strip():
            return True
    return False

def literature_text(row):
    vals = []
    for c in [
        "final_doi", "discovery_publication_doi",
        "discovery_publication_citation", "doi", "DOI",
        "reference", "citation"
    ]:
        if c in row.index and pd.notna(row.get(c)) and str(row.get(c)).strip():
            vals.append(f"**{c}:** {row.get(c)}")
    return vals

def get_protein_id_candidates(protein_row):
    candidates = []
    for c in [
        "protein_name", "uniprot_id", "UniProtID",
        "alphafold_id", "pdb_id_or_alphafold_id",
        "protein_id", "name"
    ]:
        if c in protein_row.index and pd.notna(protein_row[c]) and str(protein_row[c]).strip():
            candidates.append(str(protein_row[c]).strip())
    return list(dict.fromkeys(candidates))

def find_html_views(cif_filename):
    if HOF_DB_DIR is None or not os.path.exists(HOF_DB_DIR):
        return []
    key = str(cif_filename).replace(".cif", "")
    hits = []
    for p in Path(HOF_DB_DIR).rglob("*.html"):
        ps = str(p)
        if key in ps:
            hits.append(ps)
    return sorted(hits)

# ============================================================
# Load data
# ============================================================

@st.cache_data(show_spinner=True)
def load_data():
    xl = pd.ExcelFile(WORKBOOK_PATH)

    core = read_sheet(WORKBOOK_PATH, "core_with_functional_groups")
    if core.empty:
        core = read_sheet(WORKBOOK_PATH, "core_master")
    if core.empty:
        core = pd.read_excel(WORKBOOK_PATH, sheet_name=xl.sheet_names[0])

    fg = read_sheet(WORKBOOK_PATH, "functional_group_scores")
    pairwise = read_sheet(WORKBOOK_PATH, "hof_protein_fg_match")
    proteins = pd.read_csv(PROTEIN_PATH)

    # Merge functional layer if needed
    if not fg.empty and "hof_functional_interaction_score" not in core.columns:
        fg_keep = [
            "cif_filename",
            "fg_typing_success", "fg_confidence", "fg_flags",
            "hof_functional_interaction_score",
            "hbond_donor_group_count", "hbond_acceptor_group_count",
            "carboxylate_or_carboxylic_acid_count",
            "hydroxyl_explicit_count", "amide_like_count",
            "amine_explicit_count", "pyridine_like_n_count",
            "aromatic_ring_count", "charged_or_strong_polar_group_count",
            "polar_functional_group_count",
            "functional_group_diversity_score",
            "raw_hbond_motif_score_scaled",
            "raw_polar_motif_score_scaled",
            "raw_charge_motif_score_scaled",
            "raw_aromatic_motif_score_scaled",
        ]
        fg_keep = [c for c in fg_keep if c in fg.columns]
        if "cif_filename" in fg_keep:
            core = core.merge(fg[fg_keep].drop_duplicates("cif_filename"), on="cif_filename", how="left")

    # Merge reclassified framework/family layer
    if RECLASSIFIED_WORKBOOK and os.path.exists(RECLASSIFIED_WORKBOOK):
        series_layer = read_sheet(RECLASSIFIED_WORKBOOK, "framework_series_layer")
        keep_cols = [
            "cif_filename",
            "family_original",
            "family_revised",
            "family_broad_revised",
            "framework_series",
            "framework_series_key",
            "related_parent_series",
            "variant_label",
            "temperature_K",
            "is_near_duplicate_variant",
            "series_size",
            "representative_cif",
            "is_representative",
        ]
        keep_cols = [c for c in keep_cols if c in series_layer.columns]
        if "cif_filename" in keep_cols:
            for c in keep_cols:
                if c != "cif_filename" and c in core.columns:
                    core = core.drop(columns=[c])
            core = core.merge(series_layer[keep_cols].drop_duplicates("cif_filename"), on="cif_filename", how="left")

    core["series_key_eval"] = core.apply(get_series_key, axis=1)
    core["series_label_eval"] = core.apply(get_series_label, axis=1)
    core["family_label_eval"] = core.apply(get_family_label, axis=1)
    core["series_size_eval"] = core.groupby("series_key_eval")["series_key_eval"].transform("size")
    core["has_literature_eval"] = core.apply(has_literature, axis=1)

    # Normalised HOF features
    for c in [
        "Df_A_best", "Di_A_best", "Dif_A_best",
        "ASA_m2_g_best", "AV_volume_fraction_best",
        "AV_cm3_g_best", "framework_openness_proxy",
        "hof_functional_interaction_score",
    ]:
        if c in core.columns:
            core[c + "_norm_eval"] = minmax(core[c])
        else:
            core[c + "_norm_eval"] = 0.0

    access_terms = []
    for c in [
        "AV_volume_fraction_best_norm_eval",
        "ASA_m2_g_best_norm_eval",
        "framework_openness_proxy_norm_eval",
    ]:
        access_terms.append(num(core[c]).fillna(0))

    core["accessibility_score_metric"] = np.nanmean(np.vstack(access_terms), axis=0)

    if "hof_functional_interaction_score" in core.columns:
        hfs = num(core["hof_functional_interaction_score"])
        if hfs.max(skipna=True) <= 1.05 and hfs.min(skipna=True) >= -0.05:
            core["hof_functional_score_metric"] = hfs.fillna(0).clip(0, 1)
        else:
            core["hof_functional_score_metric"] = minmax(hfs).fillna(0)
    else:
        core["hof_functional_score_metric"] = 0.0

    return core, pairwise, proteins

core, pairwise, proteins = load_data()

# ============================================================
# Pairwise matching and scoring
# ============================================================

def pairwise_for_protein(protein_row):
    if pairwise.empty or "protein_id" not in pairwise.columns:
        return pd.DataFrame()

    p = pairwise.copy()
    p["protein_id_str"] = p["protein_id"].astype(str)

    for cand in get_protein_id_candidates(protein_row):
        exact = p[p["protein_id_str"] == cand]
        if len(exact) > 0:
            return exact.drop(columns=["protein_id_str"])

        contains = p[p["protein_id_str"].str.contains(re.escape(cand), case=False, na=False)]
        if len(contains) > 0:
            return contains.drop(columns=["protein_id_str"])

    return pd.DataFrame()

def score_for_protein(protein_row):
    df = core.copy()

    p_eff = safe_float(protein_row.get("effective_diameter_A", np.nan))
    p_min = safe_float(protein_row.get("min_dimension_A", np.nan))

    if pd.isna(p_eff):
        p_eff = safe_float(protein_row.get("radius_gyration_A", np.nan)) * 2

    df["window_size_A"] = num(df.get("Df_A_best", np.nan))
    df["cavity_size_A"] = num(df.get("Di_A_best", np.nan))

    if not pd.isna(p_min) and p_min > 0:
        df["window_to_protein_ratio"] = df["window_size_A"] / p_min
    else:
        df["window_to_protein_ratio"] = np.nan

    if not pd.isna(p_eff) and p_eff > 0:
        df["cavity_to_protein_ratio"] = df["cavity_size_A"] / p_eff
    else:
        df["cavity_to_protein_ratio"] = np.nan

    df["surface_size_score_metric"] = np.clip(df["window_size_A"].fillna(0) / 10.0, 0, 1)

    df["encapsulation_size_score_metric"] = (
        np.clip(df["cavity_to_protein_ratio"].fillna(0), 0, 1)
        * np.clip(df["window_to_protein_ratio"].fillna(0), 0, 1)
    )

    df["pore_size_score_metric"] = (
        0.70 * df["surface_size_score_metric"]
        + 0.30 * np.sqrt(np.clip(df["encapsulation_size_score_metric"], 0, 1))
    ).clip(0, 1)

    pw = pairwise_for_protein(protein_row)

    if len(pw) > 0 and "cif_filename" in pw.columns:
        keep = [
            "cif_filename",
            "hbond_complementarity_score",
            "polar_charge_complementarity_score",
            "aromatic_hydrophobic_compatibility_score",
            "pairwise_functional_complementarity_score",
        ]
        keep = [c for c in keep if c in pw.columns]
        df = df.merge(pw[keep].drop_duplicates("cif_filename"), on="cif_filename", how="left")

    for c in [
        "hbond_complementarity_score",
        "polar_charge_complementarity_score",
        "aromatic_hydrophobic_compatibility_score",
    ]:
        if c not in df.columns:
            df[c] = np.nan

    df["hbond_score_metric"] = num(df["hbond_complementarity_score"]).fillna(0).clip(0, 1)
    df["electrostatic_score_metric"] = num(df["polar_charge_complementarity_score"]).fillna(0).clip(0, 1)
    df["hydrophobic_aromatic_score_metric"] = num(df["aromatic_hydrophobic_compatibility_score"]).fillna(0).clip(0, 1)

    df["general_hof_suitability"] = (
        0.35 * df["pore_size_score_metric"]
        + 0.30 * df["accessibility_score_metric"]
        + 0.35 * df["hof_functional_score_metric"]
    ).clip(0, 1)

    df["protein_specific_compatibility"] = (
        0.35 * df["hbond_score_metric"]
        + 0.35 * df["electrostatic_score_metric"]
        + 0.30 * df["hydrophobic_aromatic_score_metric"]
    ).clip(0, 1)

    df["general_percentile"] = percentile(df["general_hof_suitability"])
    df["protein_specific_percentile"] = percentile(df["protein_specific_compatibility"])

    df["score_rank_fusion_25_75"] = (
        0.25 * df["general_percentile"]
        + 0.75 * df["protein_specific_percentile"]
    ).clip(0, 1)

    df["score_rank_fusion_40_60"] = (
        0.40 * df["general_percentile"]
        + 0.60 * df["protein_specific_percentile"]
    ).clip(0, 1)

    df["score_additive_metric_only"] = (
        0.22 * df["pore_size_score_metric"]
        + 0.18 * df["accessibility_score_metric"]
        + 0.20 * df["hbond_score_metric"]
        + 0.16 * df["electrostatic_score_metric"]
        + 0.12 * df["hydrophobic_aromatic_score_metric"]
        + 0.12 * df["hof_functional_score_metric"]
    ).clip(0, 1)

    df["score_gated_protein_specific"] = (
        df["general_hof_suitability"]
        * (0.35 + 0.65 * df["protein_specific_compatibility"])
    ).clip(0, 1)

    df["score_protein_specific_only"] = df["protein_specific_compatibility"]
    df["score_general_only"] = df["general_hof_suitability"]

    # Functionality hypotheses
    df["pred_surface_engagement"] = (
        0.35 * df["surface_size_score_metric"]
        + 0.30 * df["accessibility_score_metric"]
        + 0.20 * df["protein_specific_compatibility"]
        + 0.15 * df["hof_functional_score_metric"]
    ).clip(0, 1)

    df["pred_encapsulation"] = (
        0.45 * df["encapsulation_size_score_metric"]
        + 0.25 * df["accessibility_score_metric"]
        + 0.15 * df["pore_size_score_metric"]
        + 0.15 * df["hbond_score_metric"]
    ).clip(0, 1)

    df["pred_interfacial_catalysis"] = (
        0.25 * df["accessibility_score_metric"]
        + 0.30 * df["hof_functional_score_metric"]
        + 0.25 * df["hbond_score_metric"]
        + 0.20 * df["electrostatic_score_metric"]
    ).clip(0, 1)

    df["pred_mesoporous_hosting"] = (
        0.30 * num(df.get("AV_volume_fraction_best_norm_eval", 0)).fillna(0)
        + 0.25 * num(df.get("ASA_m2_g_best_norm_eval", 0)).fillna(0)
        + 0.20 * df["pore_size_score_metric"]
        + 0.25 * df["hof_functional_score_metric"]
    ).clip(0, 1)

    # ------------------------------------------------------------
    # Physically distinct screening routes
    # ------------------------------------------------------------
    # Growth-mediated encapsulation does not require the final aperture
    # to admit the folded protein post-synthetically. It is intentionally
    # a soft, metric-only prioritisation route for in-situ incorporation.
    df["cavity_size_norm_metric"] = minmax(df["cavity_size_A"])
    df["accessibility_percentile"] = percentile(df["accessibility_score_metric"])
    df["hof_functional_percentile"] = percentile(df["hof_functional_score_metric"])

    df["growth_hosting_context"] = (
        0.35 * df["cavity_size_norm_metric"]
        + 0.30 * df["accessibility_score_metric"]
        + 0.20 * df["hof_functional_score_metric"]
        + 0.15 * df["hbond_score_metric"]
    ).clip(0, 1)
    df["growth_hosting_context_percentile"] = percentile(df["growth_hosting_context"])

    df["route_growth_mediated"] = (
        0.50 * df["protein_specific_percentile"]
        + 0.25 * df["general_percentile"]
        + 0.25 * df["growth_hosting_context_percentile"]
    ).clip(0, 1)

    # Surface-contact route: whole-protein entry is intentionally not required.
    df["route_surface_contact"] = (
        0.70 * df["protein_specific_percentile"]
        + 0.20 * df["accessibility_percentile"]
        + 0.10 * df["hof_functional_percentile"]
    ).clip(0, 1)

    # Advanced strict route: post-synthetic diffusion of the folded protein.
    strict_mask = pd.Series(False, index=df.index)
    if not pd.isna(p_min) and p_min > 0 and not pd.isna(p_eff) and p_eff > 0:
        strict_mask = (
            (df["window_size_A"] >= p_min)
            & (df["cavity_size_A"] >= p_eff)
        )
    df["strict_entry_feasible"] = strict_mask.fillna(False)
    df["route_strict_entry"] = (
        0.75 * df["protein_specific_percentile"]
        + 0.20 * df["general_percentile"]
        + 0.05 * df["accessibility_percentile"]
    ).clip(0, 1)

    # Informative feasibility tier used for display; not an outcome claim.
    df["entry_feasibility_ratio"] = np.minimum(
        num(df["window_to_protein_ratio"]).fillna(0),
        num(df["cavity_to_protein_ratio"]).fillna(0),
    )
    df["entry_feasibility_tier"] = np.select(
        [
            df["strict_entry_feasible"],
            df["entry_feasibility_ratio"] >= 0.75,
            df["entry_feasibility_ratio"] >= 0.25,
        ],
        [
            "A — strict post-synthetic entry feasible",
            "B — near-feasible / constrained entry",
            "C — growth-mediated incorporation only",
        ],
        default="D — surface-contact hypothesis only",
    )

    mode_cols = [
        "pred_surface_engagement",
        "pred_encapsulation",
        "pred_interfacial_catalysis",
        "pred_mesoporous_hosting",
    ]
    mode_names = {
        "pred_surface_engagement": "Surface engagement / immobilisation",
        "pred_encapsulation": "Protein encapsulation",
        "pred_interfacial_catalysis": "Interfacial catalysis / microenvironment effect",
        "pred_mesoporous_hosting": "Mesoporous hosting / confinement",
    }

    df["predicted_functionality"] = df[mode_cols].idxmax(axis=1).map(mode_names)
    df["predicted_functionality_score"] = df[mode_cols].max(axis=1)

    return df

# ============================================================
# Explanation and plotting
# ============================================================

def best_function(row):
    items = {
        "Surface engagement / immobilisation": safe_float(row.get("pred_surface_engagement", np.nan)),
        "Protein encapsulation": safe_float(row.get("pred_encapsulation", np.nan)),
        "Interfacial catalysis / microenvironment effect": safe_float(row.get("pred_interfacial_catalysis", np.nan)),
        "Mesoporous hosting / confinement": safe_float(row.get("pred_mesoporous_hosting", np.nan)),
    }
    items = {k: v for k, v in items.items() if not pd.isna(v)}
    if not items:
        return "Uncertain", np.nan, "weak"
    mode = max(items, key=items.get)
    score = items[mode]
    strength = "strong" if score >= 0.70 else "moderate" if score >= 0.45 else "weak"
    return mode, score, strength

def pros_cons(row):
    pros, cons = [], []

    family = row.get("family_label_eval", "Unclassified/Other")
    series = row.get("series_label_eval", row.get("series_key_eval", "unknown"))
    series_size = safe_float(row.get("series_size_eval", 1), 1)

    df_a = safe_float(row.get("Df_A_best", np.nan))
    di_a = safe_float(row.get("Di_A_best", np.nan))
    av = safe_float(row.get("AV_volume_fraction_best", np.nan))
    asa = safe_float(row.get("ASA_m2_g_best", np.nan))

    general_p = safe_float(row.get("general_percentile", np.nan))
    protein_p = safe_float(row.get("protein_specific_percentile", np.nan))
    protein_score = safe_float(row.get("protein_specific_compatibility", np.nan))
    general_score = safe_float(row.get("general_hof_suitability", np.nan))

    if series_size > 1:
        cons.append(f"This is part of a framework-series group with {int(series_size)} CIF variants; interpret ranking at series level.")

    pros.append(f"Assigned series: {series}; broad/revised family: {family}.")

    if not pd.isna(general_p) and not pd.isna(protein_p):
        if protein_p > general_p:
            pros.append("This candidate is promoted mainly by protein-specific compatibility rather than general HOF quality alone.")
        else:
            pros.append("This candidate combines good general HOF suitability with acceptable protein-specific matching.")
            if protein_p < 0.70:
                cons.append("Protein-specific percentile is not among the strongest candidates, so the match should be treated cautiously.")

    if not pd.isna(df_a):
        if df_a >= 12:
            pros.append(f"Large limiting aperture/window for this dataset (Df ≈ {df_a:.2f} Å).")
        elif df_a >= 6:
            pros.append(f"Moderate limiting aperture/window (Df ≈ {df_a:.2f} Å), supporting surface or local-region interactions.")
            cons.append("Window size is below whole-protein scale, so full encapsulation should not be overclaimed.")
        else:
            cons.append(f"Small limiting aperture/window (Df ≈ {df_a:.2f} Å), limiting whole-protein entry.")

    if not pd.isna(di_a):
        if di_a >= 20:
            pros.append(f"Relatively large cavity estimate (Di ≈ {di_a:.2f} Å), supporting confinement or hosting hypotheses.")
        elif di_a < 12:
            cons.append(f"Small cavity estimate (Di ≈ {di_a:.2f} Å), limiting protein-scale confinement.")

    if not pd.isna(av):
        if av >= 0.50:
            pros.append(f"High accessible void fraction (AV ≈ {av:.3f}), supporting open-framework screening.")
        elif av < 0.15:
            cons.append(f"Low accessible void fraction (AV ≈ {av:.3f}), suggesting limited internal accessibility.")

    if not pd.isna(asa):
        if asa >= 10000:
            pros.append(f"Very high accessible surface area in the descriptor table ({asa:.1f} m²/g).")
            cons.append("Unusually high ASA should be manually checked for descriptor or CIF artefacts.")
        elif asa >= 500:
            pros.append(f"Accessible surface area is favourable for surface contact ({asa:.1f} m²/g).")

    donor = safe_float(row.get("hbond_donor_group_count", np.nan))
    acceptor = safe_float(row.get("hbond_acceptor_group_count", np.nan))

    if not pd.isna(donor) and not pd.isna(acceptor):
        if donor + acceptor > 0:
            if acceptor > donor * 2:
                pros.append(f"Acceptor-rich H-bond profile: donors = {int(donor)}, acceptors = {int(acceptor)}.")
            elif donor > acceptor * 2:
                pros.append(f"Donor-rich H-bond profile: donors = {int(donor)}, acceptors = {int(acceptor)}.")
            else:
                pros.append(f"Mixed donor/acceptor profile: donors = {int(donor)}, acceptors = {int(acceptor)}.")
        else:
            cons.append("No typed H-bond donor/acceptor groups were detected.")

    motif_values = {
        "carboxylate/carboxylic-acid-like groups": safe_float(row.get("carboxylate_or_carboxylic_acid_count", 0), 0),
        "hydroxyl groups": safe_float(row.get("hydroxyl_explicit_count", 0), 0),
        "amide-like groups": safe_float(row.get("amide_like_count", 0), 0),
        "amine groups": safe_float(row.get("amine_explicit_count", 0), 0),
        "pyridine-like N motifs": safe_float(row.get("pyridine_like_n_count", 0), 0),
        "aromatic rings": safe_float(row.get("aromatic_ring_count", 0), 0),
    }
    present = {k: v for k, v in motif_values.items() if v > 0}

    if present:
        dom = max(present, key=present.get)
        pros.append(f"Dominant typed motif: {dom} ({int(present[dom])}).")
        if len(present) >= 4:
            pros.append("Chemically diverse motif profile, supporting multiple possible interaction routes.")
    else:
        cons.append("No major typed functional motif was detected in the available functional-group layer.")

    mode, mode_score, strength = best_function(row)
    pros.append(f"Most supported functionality hypothesis: {mode} ({strength}, score {fmt(mode_score, 3)}).")

    if protein_score < 0.12:
        cons.append("Absolute protein-specific compatibility score is low, even if percentile rank is favourable within this dataset.")

    if not has_literature(row):
        cons.append("No direct DOI/citation was available for this entry in the workbook.")
    else:
        pros.append("Literature/citation metadata is available for traceability.")

    if len(find_html_views(row.get("cif_filename", ""))) == 0:
        cons.append("No matching 3D HTML view was found.")
    else:
        pros.append("A matching 3D HTML view is available for inspection.")

    pros = list(dict.fromkeys(pros))
    cons = list(dict.fromkeys(cons))

    return pros[:8], cons[:8]

def radar_chart(row):
    labels = [
        "Pore",
        "Access",
        "HOF chem",
        "H-bond",
        "Electrostatic",
        "Hydrophobic",
        "Protein match",
    ]
    vals = [
        safe_float(row.get("pore_size_score_metric", 0), 0),
        safe_float(row.get("accessibility_score_metric", 0), 0),
        safe_float(row.get("hof_functional_score_metric", 0), 0),
        safe_float(row.get("hbond_score_metric", 0), 0),
        safe_float(row.get("electrostatic_score_metric", 0), 0),
        safe_float(row.get("hydrophobic_aromatic_score_metric", 0), 0),
        safe_float(row.get("protein_specific_compatibility", 0), 0),
    ]

    angles = np.linspace(0, 2*np.pi, len(labels), endpoint=False).tolist()
    vals += vals[:1]
    angles += angles[:1]

    fig = plt.figure(figsize=(5, 5))
    ax = plt.subplot(111, polar=True)
    ax.plot(angles, vals, linewidth=2)
    ax.fill(angles, vals, alpha=0.15)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylim(0, 1)
    ax.set_title("Metric profile", pad=18)
    return fig

# ============================================================
# Sidebar controls
# ============================================================

st.sidebar.title("HOF–Protein Recommender")

protein_labels = []
for _, r in proteins.iterrows():
    pname = str(r.get("protein_name", r.get("name", "Unknown protein")))
    uid = str(r.get("uniprot_id", r.get("UniProtID", "")))
    protein_labels.append(f"{pname} | {uid}")

selected_protein_label = st.sidebar.selectbox("Target protein", protein_labels, index=0)
protein_idx = protein_labels.index(selected_protein_label)
protein_row = proteins.iloc[protein_idx]

st.sidebar.markdown("---")
st.sidebar.subheader("Scoring model")

model_display = {
    "rank_fusion_25_75": "Rank fusion 25/75 — recommended",
    "rank_fusion_40_60": "Rank fusion 40/60 — conservative",
    "additive_metric_only": "Additive metric-only",
    "gated_protein_specific": "Gated protein-specific",
    "protein_specific_only": "Protein-specific only control",
    "general_only": "General HOF only control",
}

selected_model_label = st.sidebar.selectbox(
    "Ranking logic",
    list(model_display.values()),
    index=0
)

selected_model = [k for k, v in model_display.items() if v == selected_model_label][0]

score_col_map = {
    "rank_fusion_25_75": "score_rank_fusion_25_75",
    "rank_fusion_40_60": "score_rank_fusion_40_60",
    "additive_metric_only": "score_additive_metric_only",
    "gated_protein_specific": "score_gated_protein_specific",
    "protein_specific_only": "score_protein_specific_only",
    "general_only": "score_general_only",
}

score_col = score_col_map[selected_model]

st.sidebar.caption(
    "Default uses rank fusion: 25% general HOF suitability percentile + 75% protein-specific compatibility percentile."
)

st.sidebar.markdown("---")
st.sidebar.subheader("Screening objective")

screening_objective = st.sidebar.radio(
    "Choose the physical question",
    [
        "Growth-mediated encapsulation",
        "Interaction / surface-contact compatibility",
        "Strict infiltration / post-synthetic entry",
    ],
    index=0,
)

if screening_objective == "Growth-mediated encapsulation":
    st.sidebar.caption(
        "Default route. Protein incorporation may occur during HOF formation, so the final aperture is not required to admit the folded protein post-synthetically."
    )
elif screening_objective == "Interaction / surface-contact compatibility":
    st.sidebar.caption(
        "Ranks external and near-surface interaction hypotheses. Whole-protein pore entry is intentionally not required."
    )
else:
    st.sidebar.caption(
        "Advanced conservative route. Requires both Df ≥ protein minimum dimension and Di ≥ protein effective diameter."
    )

mode_filter = st.sidebar.selectbox(
    "Optional functionality emphasis",
    [
        "No additional emphasis",
        "Surface engagement / immobilisation",
        "Protein encapsulation",
        "Interfacial catalysis / microenvironment effect",
        "Mesoporous hosting / confinement",
    ],
    index=0,
)

top_n = st.sidebar.slider("Number of recommendations", 5, 50, 20, 5)

st.sidebar.markdown("---")
st.sidebar.subheader("Filters")

collapse_series = st.sidebar.checkbox("Collapse framework-series variants", value=True)
require_3d = st.sidebar.checkbox("Require 3D HTML view", value=False)
require_lit = st.sidebar.checkbox("Require literature/citation", value=False)
require_fg = st.sidebar.checkbox("Require functional-group typing success", value=False)
min_score = st.sidebar.slider("Minimum score", 0.0, 1.0, 0.0, 0.01)
search_text = st.sidebar.text_input("Search CIF / series / family", "")

# ============================================================
# Main app
# ============================================================

st.title("HOF–Protein Compatibility Recommender")

st.markdown(
    """
This interface ranks HOF candidates using **metric-only HOF–protein compatibility scoring**.
Literature, family and framework-series labels are used only for **filtering, grouping, display and validation**;
they are not used as scoring evidence.
"""
)

st.info(
    "Default physical route: growth-mediated encapsulation. The default does not require post-synthetic diffusion of the folded protein through the final aperture. Strict infiltration remains available as an advanced conservative screen."
)

if screening_objective == "Growth-mediated encapsulation":
    st.markdown(
        "**Growth-mediated encapsulation** ranks candidates for in-situ incorporation during HOF formation. "
        "The final aperture is not treated as a hard whole-protein entry requirement."
    )
elif screening_objective == "Interaction / surface-contact compatibility":
    st.markdown(
        "**Interaction / surface-contact compatibility** ranks surface and near-surface interaction hypotheses. "
        "Whole-protein entry is not required."
    )
else:
    st.markdown(
        "**Strict infiltration / post-synthetic entry** applies the conservative rule: "
        "Df ≥ protein minimum dimension and Di ≥ protein effective diameter."
    )

# Score
df = score_for_protein(protein_row)

if screening_objective == "Growth-mediated encapsulation":
    route_col = "route_growth_mediated"
    route_label = "Growth-mediated encapsulation"
elif screening_objective == "Interaction / surface-contact compatibility":
    route_col = "route_surface_contact"
    route_label = "Interaction / surface-contact compatibility"
else:
    route_col = "route_strict_entry"
    route_label = "Strict infiltration / post-synthetic entry"

# Optional function-conditioned reranking remains secondary and user-controlled.
# It does not use literature labels and therefore does not create circular validation.
function_col = None
if mode_filter == "Surface engagement / immobilisation":
    function_col = "pred_surface_engagement"
elif mode_filter == "Protein encapsulation":
    function_col = "pred_encapsulation"
elif mode_filter == "Interfacial catalysis / microenvironment effect":
    function_col = "pred_interfacial_catalysis"
elif mode_filter == "Mesoporous hosting / confinement":
    function_col = "pred_mesoporous_hosting"

if function_col is None:
    df["sort_score"] = df[route_col]
else:
    df["sort_score"] = 0.80 * df[route_col] + 0.20 * df[function_col]

strict_fallback = False
if screening_objective == "Strict infiltration / post-synthetic entry":
    strict_df = df[df["strict_entry_feasible"] == True].copy()
    if len(strict_df) == 0:
        strict_fallback = True
        st.warning(
            "No strict post-synthetic entry candidates were found for this protein. "
            "This does not exclude growth-mediated incorporation during HOF formation. "
            "The closest growth-mediated alternatives are shown below for comparison."
        )
        route_col = "route_growth_mediated"
        route_label = "Growth-mediated alternatives after strict-screen failure"
        if function_col is None:
            df["sort_score"] = df[route_col]
        else:
            df["sort_score"] = 0.80 * df[route_col] + 0.20 * df[function_col]
    else:
        df = strict_df

df = df[df["sort_score"] >= min_score].copy()

if require_lit:
    df = df[df["has_literature_eval"] == True].copy()

if require_fg and "fg_typing_success" in df.columns:
    df = df[df["fg_typing_success"] == True].copy()

if require_3d:
    df = df[df["cif_filename"].apply(lambda x: len(find_html_views(x)) > 0)].copy()

if search_text.strip():
    q = search_text.strip().lower()
    mask = (
        df["cif_filename"].astype(str).str.lower().str.contains(q, na=False)
        | df["series_label_eval"].astype(str).str.lower().str.contains(q, na=False)
        | df["family_label_eval"].astype(str).str.lower().str.contains(q, na=False)
    )
    df = df[mask].copy()

df = df.sort_values("sort_score", ascending=False)

if collapse_series:
    df = df.drop_duplicates("series_key_eval", keep="first").copy()

df = df.head(top_n).reset_index(drop=True)

# Header metrics
m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Protein", str(protein_row.get("protein_name", "Unknown"))[:18])
m2.metric("UniProt", str(protein_row.get("uniprot_id", protein_row.get("UniProtID", "")))[:12])
m3.metric("Ranked HOFs", len(df))
m4.metric("Objective", route_label[:18])
m5.metric("Total HOF rows", len(core))

with st.expander("Protein surface descriptor details"):
    pshow_cols = [
        "protein_name", "uniprot_id", "effective_diameter_A",
        "surface_net_charge_index", "surface_residue_count",
        "hydrophobic_surface_fraction", "aromatic_surface_fraction",
        "hbond_donor_density", "hbond_acceptor_density",
        "largest_positive_patch_area_A2", "largest_negative_patch_area_A2",
        "largest_hydrophobic_patch_area_A2",
    ]
    pshow_cols = [c for c in pshow_cols if c in proteins.columns]
    render_html_table(pd.DataFrame([protein_row[pshow_cols]]), max_height=260)

st.markdown("---")

# Table
st.subheader("Ranked recommendation table")

table_cols = [
    "cif_filename",
    "series_label_eval",
    "family_label_eval",
    "sort_score",
    route_col,
    score_col,
    "entry_feasibility_tier",
    "general_hof_suitability",
    "protein_specific_compatibility",
    "general_percentile",
    "protein_specific_percentile",
    "Df_A_best",
    "Di_A_best",
    "AV_volume_fraction_best",
    "ASA_m2_g_best",
    "hbond_score_metric",
    "electrostatic_score_metric",
    "hydrophobic_aromatic_score_metric",
    "hof_functional_score_metric",
    "predicted_functionality",
    "series_size_eval",
    "has_literature_eval",
]
table_cols = [c for c in table_cols if c in df.columns]
render_html_table(df[table_cols], max_height=620)

csv = df[table_cols].to_csv(index=False).encode("utf-8")
st.download_button(
    "Download ranked table as CSV",
    data=csv,
    file_name="hof_rank_fusion_recommendations.csv",
    mime="text/csv"
)

st.markdown("---")
st.subheader("Recommendation details")

for i, row in df.iterrows():
    st.markdown(f"<div class='card'>", unsafe_allow_html=True)

    c1, c2 = st.columns([1.15, 1])

    with c1:
        st.markdown(f"### #{i+1} — {row.get('cif_filename', 'Unknown CIF')}")
        st.markdown(f"<span class='tag'>Series: {row.get('series_label_eval', 'NA')}</span>", unsafe_allow_html=True)
        st.markdown(f"<span class='tag'>Family: {row.get('family_label_eval', 'NA')}</span>", unsafe_allow_html=True)

        if row.get("has_literature_eval", False):
            st.markdown("<span class='tag good'>Literature metadata: available</span>", unsafe_allow_html=True)
        else:
            st.markdown("<span class='tag warn'>Literature metadata: not found</span>", unsafe_allow_html=True)

        st.markdown(f"<span class='tag'>Screening route: {route_label}</span>", unsafe_allow_html=True)
        st.markdown(f"<span class='tag'>Feasibility tier: {row.get('entry_feasibility_tier', 'NA')}</span>", unsafe_allow_html=True)
        st.markdown(f"<span class='tag'>Predicted mode: {row.get('predicted_functionality', 'NA')}</span>", unsafe_allow_html=True)

    with c2:
        a, b, c, d = st.columns(4)
        a.metric("Route score", fmt(row.get("sort_score", np.nan), 3))
        b.metric("Model score", fmt(row.get(score_col, np.nan), 3))
        c.metric("Protein match", fmt(row.get("protein_specific_compatibility", np.nan), 3))
        d.metric("General HOF", fmt(row.get("general_hof_suitability", np.nan), 3))

    tab1, tab2, tab3, tab4, tab5 = st.tabs(
        ["Why recommended", "Metrics", "Polygon chart", "3D view", "Literature"]
    )

    with tab1:
        pros, cons = pros_cons(row)

        pc1, pc2 = st.columns(2)

        with pc1:
            st.markdown("#### Pros")
            for p in pros:
                st.markdown(f"- {p}")

        with pc2:
            st.markdown("#### Cons / cautions")
            for ctext in cons:
                st.markdown(f"- {ctext}")

        st.caption(
            "Functionality prediction is a metric-based screening hypothesis, not an experimentally validated claim."
        )

    with tab2:
        metric_cols = [
            "sort_score",
            route_col,
            score_col,
            "route_growth_mediated",
            "route_surface_contact",
            "route_strict_entry",
            "strict_entry_feasible",
            "entry_feasibility_ratio",
            "entry_feasibility_tier",
            "general_hof_suitability",
            "protein_specific_compatibility",
            "general_percentile",
            "protein_specific_percentile",
            "pore_size_score_metric",
            "accessibility_score_metric",
            "hof_functional_score_metric",
            "hbond_score_metric",
            "electrostatic_score_metric",
            "hydrophobic_aromatic_score_metric",
            "Df_A_best",
            "Di_A_best",
            "Dif_A_best",
            "AV_volume_fraction_best",
            "ASA_m2_g_best",
            "hbond_donor_group_count",
            "hbond_acceptor_group_count",
            "series_size_eval",
            "pred_surface_engagement",
            "pred_encapsulation",
            "pred_interfacial_catalysis",
            "pred_mesoporous_hosting",
        ]
        metric_cols = [c for c in metric_cols if c in row.index]
        render_html_table(pd.DataFrame({"metric": metric_cols, "value": [row.get(c) for c in metric_cols]}), max_height=520)

    with tab3:
        st.pyplot(radar_chart(row))

    with tab4:
        views = find_html_views(row.get("cif_filename", ""))
        if views:
            view_choice = st.selectbox(
                f"3D view for {row.get('cif_filename')}",
                views,
                key=f"view_{i}_{row.get('cif_filename')}"
            )
            try:
                html = Path(view_choice).read_text(encoding="utf-8", errors="ignore")
                components.html(html, height=620, scrolling=True)
            except Exception as e:
                st.warning(f"Could not load HTML view: {e}")
        else:
            st.info("No matching 3D HTML view found for this CIF.")

    with tab5:
        lits = literature_text(row)
        if lits:
            for item in lits:
                st.markdown(item)
        else:
            st.info("No DOI/citation metadata found for this entry.")

    st.markdown("</div>", unsafe_allow_html=True)

# Footer
st.markdown("---")
st.caption(
    "Scoring uses only HOF/protein metrics. Literature/family/series information is used for display, filtering, grouping, and validation only."
)
