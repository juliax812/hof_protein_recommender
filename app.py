import re
import unicodedata
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

st.set_page_config(
    page_title="Protein-oriented HOF prioritization",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
<style>
:root { color-scheme: light !important; }
html, body, .stApp { background:#ffffff !important; color:#111111 !important; }
section[data-testid="stSidebar"] { background:#f6f7f8 !important; }
section[data-testid="stSidebar"] * { color:#111111 !important; }
div[data-baseweb="select"] *, div[role="listbox"] * {
    color:#111111 !important;
    -webkit-text-fill-color:#111111 !important;
}
</style>
""",
    unsafe_allow_html=True,
)

PROTEIN_PATH = "protein_descriptors_FINAL_100_UI_MIN.csv"
WORKBOOK_PATH = "HOF_MASTER_WORKBOOK_functional_group_typed_FIXED.xlsx"
CURATED_PATH = "HOF_DATABASE_FINAL_CURATED_V1.xlsx"
HOF_DB_DIR = Path("hof_db")

for path in [PROTEIN_PATH, WORKBOOK_PATH, CURATED_PATH]:
    if not Path(path).exists():
        st.error(f"Missing required file: {path}")
        st.stop()


def num(values):
    return pd.to_numeric(values, errors="coerce")


def minmax(values):
    s = num(values)
    lo, hi = s.min(), s.max()
    if pd.isna(lo) or pd.isna(hi) or hi - lo <= 1e-15:
        return pd.Series(0.0, index=s.index, dtype=float)
    return ((s - lo) / (hi - lo)).fillna(0.0)


def percentile(values):
    return num(values).rank(method="average", pct=True).fillna(0.0)


def norm_key(value):
    s = unicodedata.normalize("NFKC", str(value)).casefold().replace(".cif", "")
    return re.sub(r"[^a-z0-9]+", "", s)


def fmt(value, digits=3):
    try:
        if pd.isna(value):
            return "NA"
        return f"{float(value):.{digits}f}"
    except Exception:
        return "NA"


def repair_label(value):
    """Repair legacy UTF-8-as-CP437 mojibake for display only.

    Identifiers, CIF keys, series keys, and 3D lookup keys remain untouched.
    """
    if pd.isna(value):
        return ""
    text = str(value)
    markers = ("┬", "╬", "ΓÇ", "╖", "╢", "╣", "╠", "╦", "╩", "╚", "╝")
    if not any(m in text for m in markers):
        return text
    try:
        repaired = text.encode("cp437").decode("utf-8")
    except Exception:
        return text
    old_score = sum(text.count(m) for m in markers)
    new_score = sum(repaired.count(m) for m in markers)
    return repaired if new_score < old_score else text


def framework_class(row):
    fn = str(row.get("cif_filename", ""))
    ccdc = int(row["ccdc_id"]) if pd.notna(row.get("ccdc_id")) else None
    metal = (
        float(row.get("metal_count_cif", 0) or 0)
        if pd.notna(row.get("metal_count_cif"))
        else 0.0
    )
    text = (fn + " " + str(row.get("discovery_publication_citation") or "")).lower()

    if ccdc == 1991981:
        return "retained_verified_hof"
    if ccdc in {1991980, 1991982, 1566036, 1566037, 1893340, 1893341}:
        return "excluded_confirmed_nonhof"
    if fn == "257-rub-hofs.cif":
        return "retained_verified_metal_hof"
    if metal <= 0:
        return "retained_source_hof_candidate"
    if any(x in text for x in ["hydrogen bonded", "hydrogen-bonded", "m-hof", "hof-", "hofs"]):
        return "retained_verified_metal_hof"
    if any(x in text for x in ["coordination polymer", "metal-organic framework", "porous 3d mof", " mofs "]):
        return "excluded_confirmed_nonhof"
    return "manual_review_excluded_main"


@st.cache_data(show_spinner="Building audited 651 × 100 matrix…")
def load_data():
    proteins = pd.read_csv(PROTEIN_PATH)
    core = pd.read_excel(WORKBOOK_PATH, sheet_name="core_with_functional_groups").drop_duplicates("cif_filename")
    series = pd.read_excel(CURATED_PATH, sheet_name="framework_series_layer")

    assert len(proteins) == 100
    assert proteins["uniprot_id"].nunique() == 100

    proteins["protein_id_eval"] = proteins["uniprot_id"].astype(str)
    for col in [
        "surface_hbond_donor_density",
        "surface_hbond_acceptor_density",
        "surface_pos_residue_frac",
        "surface_neg_residue_frac",
        "surface_hydrophobic_frac",
        "surface_aromatic_frac",
    ]:
        proteins[col + "_norm"] = minmax(proteins[col])

    core["framework_class_status"] = core.apply(framework_class, axis=1)
    core["main_analysis_included"] = core["framework_class_status"].str.startswith("retained")

    series_cols = [
        c
        for c in [
            "cif_filename",
            "framework_series_key",
            "framework_series",
            "family_revised",
            "broad_family_final",
            "display_group_final",
            "series_size",
            "representative_cif",
            "is_representative",
        ]
        if c in series.columns
    ]
    core = core.merge(series[series_cols].drop_duplicates("cif_filename"), on="cif_filename", how="left")

    missing_series = core["framework_series_key"].isna() | core["framework_series_key"].astype(str).str.strip().eq("")
    core["series_key_eval"] = core["framework_series_key"]
    core.loc[missing_series, "series_key_eval"] = core.loc[missing_series, "cif_filename"].map(norm_key)
    core["series_label_eval"] = core["framework_series"].where(
        core["framework_series"].notna(), core["series_key_eval"]
    )

    pentahof = core["ccdc_id"].eq(1991981)
    core.loc[pentahof, "series_key_eval"] = "pentahof1"
    core.loc[pentahof, "series_label_eval"] = "pentaHOF-1"

    for col in ["Di_A_best", "Df_A_best", "Dif_A_best", "ASA_m2_g_best", "AV_volume_fraction_best"]:
        core[col] = num(core[col]).fillna(0.0)

    di = core["Di_A_best"]
    core["relative_aperture_path_ratio"] = (
        0.5 * np.where(di > 0, (core["Df_A_best"] / di).clip(0, 1), 0.0)
        + 0.5 * np.where(di > 0, (core["Dif_A_best"] / di).clip(0, 1), 0.0)
    )

    core = core[core["main_analysis_included"] & core["fg_typing_success"].eq(True)].copy()
    assert len(core) == 651

    for col in [
        "hbond_donor_group_count_per_100_heavy_atoms_scaled",
        "hbond_acceptor_group_count_per_100_heavy_atoms_scaled",
        "raw_charge_motif_score_scaled",
        "raw_aromatic_motif_score_scaled",
    ]:
        core[col] = num(core[col]).fillna(0.0).clip(0, 1)

    core["Di_norm"] = minmax(core["Di_A_best"])
    core["AV_norm"] = minmax(core["AV_volume_fraction_best"])
    core["ASA_norm"] = minmax(core["ASA_m2_g_best"])
    core["O_norm"] = minmax(core["relative_aperture_path_ratio"])
    core["G_enc"] = 0.40 * core["Di_norm"] + 0.35 * core["AV_norm"] + 0.25 * core["O_norm"]
    core["G_surf"] = 0.70 * core["ASA_norm"] + 0.30 * core["O_norm"]
    core["P_G_enc"] = percentile(core["G_enc"])
    core["P_G_surf"] = percentile(core["G_surf"])

    hof_cols = [
        "cif_filename",
        "series_key_eval",
        "series_label_eval",
        "Di_A_best",
        "Df_A_best",
        "Dif_A_best",
        "ASA_m2_g_best",
        "AV_volume_fraction_best",
        "G_enc",
        "G_surf",
        "P_G_enc",
        "P_G_surf",
        "protein_screen_confidence_tier",
        "hbond_donor_group_count_per_100_heavy_atoms_scaled",
        "hbond_acceptor_group_count_per_100_heavy_atoms_scaled",
        "raw_charge_motif_score_scaled",
        "raw_aromatic_motif_score_scaled",
    ]
    protein_cols = [
        "protein_id_eval",
        "uniprot_id",
        "protein_name",
        "organism",
        "min_dimension_A",
        "effective_diameter_A",
        "mean_plddt",
        "surface_residue_fraction",
        "surface_hbond_donor_density_norm",
        "surface_hbond_acceptor_density_norm",
        "surface_pos_residue_frac_norm",
        "surface_neg_residue_frac_norm",
        "surface_hydrophobic_frac_norm",
        "surface_aromatic_frac_norm",
    ]

    h = core[hof_cols].copy()
    p = proteins[protein_cols].copy()
    h["_k"] = 1
    p["_k"] = 1
    pairs = h.merge(p, on="_k").drop(columns="_k")

    pairs["H"] = 0.5 * (
        pairs["hbond_donor_group_count_per_100_heavy_atoms_scaled"]
        * pairs["surface_hbond_acceptor_density_norm"]
        + pairs["hbond_acceptor_group_count_per_100_heavy_atoms_scaled"]
        * pairs["surface_hbond_donor_density_norm"]
    )
    pairs["Q"] = pairs["raw_charge_motif_score_scaled"] * 0.5 * (
        pairs["surface_pos_residue_frac_norm"] + pairs["surface_neg_residue_frac_norm"]
    )
    pairs["R"] = pairs["raw_aromatic_motif_score_scaled"] * 0.5 * (
        pairs["surface_hydrophobic_frac_norm"] + pairs["surface_aromatic_frac_norm"]
    )
    pairs["chemistry"] = 0.35 * pairs["H"] + 0.35 * pairs["Q"] + 0.30 * pairs["R"]
    pairs["Pp_C"] = pairs.groupby("protein_id_eval")["chemistry"].rank(method="average", pct=True)
    pairs["S_hp"] = pairs.groupby("cif_filename")["chemistry"].rank(method="average", pct=True)

    meta_cols = [
        c
        for c in [
            "cif_filename",
            "family",
            "family_original",
            "discovery_publication_year",
            "discovery_publication_doi",
            "discovery_publication_citation",
            "final_doi",
            "fg_confidence",
            "fg_flags",
            "family_revised",
            "broad_family_final",
            "display_group_final",
        ]
        if c in core.columns
    ]
    pairs = pairs.merge(core[meta_cols].drop_duplicates("cif_filename"), on="cif_filename", how="left")

    assert len(pairs) == 65100
    return pairs, proteins


PAIRS, PROTEINS = load_data()


@st.cache_data(show_spinner=False)
def build_3d_index():
    index = {}
    if HOF_DB_DIR.exists():
        for folder in HOF_DB_DIR.iterdir():
            if folder.is_dir():
                index.setdefault(norm_key(folder.name), []).append(folder)
    return index


VIEW_INDEX = build_3d_index()


def views_for_cif(cif):
    folders = VIEW_INDEX.get(norm_key(cif), [])
    if len(folders) != 1:
        return {}
    folder = folders[0]
    return {
        name: str(folder / f"{name}.html")
        for name in ["unit", "super", "wire"]
        if (folder / f"{name}.html").exists()
    }


def literature_items(row):
    out = []
    for col in ["final_doi", "discovery_publication_doi", "discovery_publication_citation"]:
        value = row.get(col)
        if pd.notna(value) and str(value).strip():
            out.append((col, str(value)))
    return out


def family_label(row):
    for col in ["broad_family_final", "family_revised", "family", "family_original", "display_group_final"]:
        value = row.get(col)
        if pd.notna(value) and str(value).strip():
            return repair_label(value)
    return "Unclassified"


def route_frame(protein_id, route, threshold):
    d = PAIRS[PAIRS["protein_id_eval"].eq(protein_id)].copy()

    if route == "Growth-mediated integration / encapsulation":
        d["PG"] = d["P_G_enc"]
        gate = (
            (d["AV_volume_fraction_best"] > 0)
            & (d["ASA_m2_g_best"] > 0)
            & (d["P_G_enc"] >= threshold)
        )
        note = (
            "Framework forms around/with the protein; final pore aperture is not treated "
            "as a whole-protein entry requirement."
        )
    elif route == "Accessible-interface contact":
        d["PG"] = d["P_G_surf"]
        gate = (d["ASA_m2_g_best"] > 0) & (d["P_G_surf"] >= threshold)
        note = (
            "Small-probe-accessible framework interfaces are compared; periodic ASA is not "
            "claimed to equal external particle area."
        )
    else:
        d["PG"] = d["P_G_enc"]
        gate = (
            (d["AV_volume_fraction_best"] > 0)
            & (d["ASA_m2_g_best"] > 0)
            & (d["P_G_enc"] >= threshold)
            & (d["Df_A_best"] >= d["min_dimension_A"])
            & (d["Di_A_best"] >= d["effective_diameter_A"])
        )
        note = (
            "Conservative static whole-protein entry filter; short or empty lists are expected."
        )

    d = d.loc[gate].copy()
    d["final_score"] = 0.65 * d["Pp_C"] + 0.30 * d["S_hp"] + 0.05 * d["PG"]
    return d, note


def profile_plot(row):
    labels = ["H", "Q", "R", "Pp_C", "S_hp", "Geometry"]
    values = [row["H"], row["Q"], row["R"], row["Pp_C"], row["S_hp"], row["PG"]]
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False)
    angles = np.r_[angles, angles[0]]
    values = np.r_[values, values[0]]
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111, polar=True)
    ax.plot(angles, values)
    ax.fill(angles, values, alpha=0.12)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1)
    return fig


st.sidebar.header("Screening controls")
protein_options = ["— Select protein —"] + [
    f"{row.protein_name} | {row.uniprot_id}"
    for _, row in PROTEINS.sort_values(["protein_name", "uniprot_id"]).iterrows()
]
selection = st.sidebar.selectbox("Target protein", protein_options)
route = st.sidebar.selectbox(
    "Screening route",
    [
        "Growth-mediated integration / encapsulation",
        "Accessible-interface contact",
        "Strict post-synthetic infiltration",
    ],
)
tau = st.sidebar.slider("Geometry feasibility threshold (τ)", 0.0, 0.75, 0.25, 0.05)
top_n = st.sidebar.slider("Recommendations", 5, 50, 20, 5)
collapse_series = st.sidebar.checkbox("Collapse framework-series variants", True)
require_3d = st.sidebar.checkbox("Require exact 3D view")
require_literature = st.sidebar.checkbox("Require literature metadata")
search_text = st.sidebar.text_input("Search CIF / series / family", "")

st.sidebar.markdown("---")
st.sidebar.caption("Frozen manuscript model")
st.sidebar.code("H/Q/R = 0.35/0.35/0.30\nFinal = 0.65 Pp_C + 0.30 S_hp + 0.05 P_G")
st.sidebar.caption("Literature, family labels and 3D availability never enter the score.")

st.title("Protein-oriented HOF prioritization")
st.markdown(
    "**Auditable decision support for the frozen 100-protein / 651-HOF analysis.** "
    "Scores are prioritization hypotheses, not experimental compatibility probabilities."
)

c1, c2, c3, c4 = st.columns(4)
c1.metric("Curated proteins", 100)
c2.metric("Typed HOFs", 651)
c3.metric("Pairwise chemistry rows", "65,100")
c4.metric("Typed framework series", 550)

if selection == "— Select protein —":
    st.info("Select a protein to generate a route-specific shortlist.")
    st.stop()

protein = PROTEINS[(PROTEINS["protein_name"] + " | " + PROTEINS["uniprot_id"]).eq(selection)].iloc[0]
protein_id = str(protein["uniprot_id"])
shortlist, route_note = route_frame(protein_id, route, tau)

if require_3d:
    shortlist = shortlist[shortlist["cif_filename"].map(lambda x: bool(views_for_cif(x)))]
if require_literature:
    shortlist = shortlist[shortlist.apply(lambda row: bool(literature_items(row)), axis=1)]
if search_text.strip():
    query = search_text.casefold()
    shortlist = shortlist[
        shortlist["cif_filename"].astype(str).map(repair_label).str.casefold().str.contains(query, na=False)
        | shortlist["series_label_eval"].astype(str).map(repair_label).str.casefold().str.contains(query, na=False)
        | shortlist.apply(lambda row: query in family_label(row).casefold(), axis=1)
    ]

shortlist = shortlist.sort_values(["final_score", "PG", "cif_filename"], ascending=[False, False, True])
if collapse_series:
    shortlist = shortlist.drop_duplicates("series_key_eval")
shortlist = shortlist.head(top_n).reset_index(drop=True)
shortlist.insert(0, "rank", np.arange(1, len(shortlist) + 1))

st.subheader(f"{protein.protein_name} — {protein_id}")
st.caption(route_note)
p1, p2, p3, p4 = st.columns(4)
p1.metric("Organism", str(protein.get("organism", "NA"))[:28])
p2.metric("Mean pLDDT", fmt(protein.get("mean_plddt"), 1))
p3.metric("Surface residues", int(protein.get("surface_residue_count", 0)))
p4.metric("Candidates shown", len(shortlist))

if (
    pd.notna(protein.get("mean_plddt"))
    and float(protein["mean_plddt"]) < 70
) or (
    pd.notna(protein.get("surface_residue_fraction"))
    and float(protein["surface_residue_fraction"]) >= 0.95
):
    st.warning(
        "This protein belongs to the lower-confidence structural sensitivity subset; "
        "interpret surface-specific ranking cautiously."
    )

display_table = shortlist.copy()
display_table["series"] = display_table["series_label_eval"].map(repair_label)
display_table["CIF"] = display_table["cif_filename"].map(repair_label)
show_cols = [
    "rank", "series", "CIF", "final_score", "Pp_C", "S_hp", "PG", "H", "Q", "R",
    "Df_A_best", "Di_A_best", "AV_volume_fraction_best", "ASA_m2_g_best",
]
st.dataframe(display_table[show_cols], use_container_width=True, hide_index=True)

download_cols = [
    "rank", "series_key_eval", "series_label_eval", "cif_filename", "final_score", "Pp_C",
    "S_hp", "PG", "H", "Q", "R", "Df_A_best", "Di_A_best",
    "AV_volume_fraction_best", "ASA_m2_g_best",
]
st.download_button(
    "Download shortlist",
    shortlist[download_cols].to_csv(index=False),
    f"HOF_{protein_id}_top{top_n}.csv",
    "text/csv",
)

for _, row in shortlist.iterrows():
    st.markdown("---")
    st.markdown(f"### #{int(row['rank'])} — {repair_label(row['series_label_eval'])}")
    st.caption(
        f"{repair_label(row['cif_filename'])} · {family_label(row)} · "
        f"HOF confidence: {row.get('protein_screen_confidence_tier', 'NA')}"
    )

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Final", fmt(row["final_score"]))
    m2.metric("Chemistry %ile", fmt(row["Pp_C"]))
    m3.metric("Differentiation", fmt(row["S_hp"]))
    m4.metric("Geometry", fmt(row["PG"]))

    tab1, tab2, tab3, tab4 = st.tabs(["Why shortlisted", "Score profile", "3D view", "Literature"])

    with tab1:
        supporting = []
        cautions = []
        supporting.append(
            "Protein-specific chemistry is in the top quartile."
            if row["Pp_C"] >= 0.75
            else "Candidate passes the selected route gate."
        )
        supporting.append(
            "Strong target-relative differentiation."
            if row["S_hp"] >= 0.75
            else "More generalist than target-specific."
            if row["S_hp"] < 0.5
            else "Moderate target-relative differentiation."
        )
        strongest = max(
            {"H-bond": row["H"], "polar/charged motif": row["Q"], "aromatic/hydrophobic": row["R"]},
            key={"H-bond": row["H"], "polar/charged motif": row["Q"], "aromatic/hydrophobic": row["R"]}.get,
        )
        supporting.append("Largest raw chemistry contribution: " + strongest)
        if row["PG"] < 0.35:
            cautions.append("Route geometry is relatively weak despite passing the threshold.")
        if route == "Strict post-synthetic infiltration":
            cautions.append("Static dimensions do not model flexibility, solvent or kinetics.")

        st.markdown("**Supporting evidence**")
        for item in supporting:
            st.markdown("- " + item)
        st.markdown("**Cautions**")
        if cautions:
            for item in cautions:
                st.markdown("- " + item)
        else:
            st.markdown("- No additional descriptor-level caution triggered.")
        st.caption("Descriptor explanation only; not an experimental outcome prediction.")

    with tab2:
        st.pyplot(profile_plot(row), use_container_width=False)

    with tab3:
        view_map = views_for_cif(row["cif_filename"])
        if not view_map:
            st.info("No unique exact-normalized 3D HTML view for this CIF.")
        else:
            view_name = st.radio(
                "View",
                list(view_map),
                horizontal=True,
                key=f"view_{protein_id}_{row['rank']}",
            )
            components.html(
                Path(view_map[view_name]).read_text(encoding="utf-8", errors="ignore"),
                height=620,
                scrolling=True,
            )

    with tab4:
        items = literature_items(row)
        if not items:
            st.info("No DOI/citation metadata attached to this structural entry.")
        else:
            for key, value in items:
                st.markdown(f"**{key}:** {value}")
        st.caption("Literature is traceability/context only and is not a score input.")

st.markdown("---")
st.caption(
    "Frozen analysis: 100 unique proteins, 651 functionally typed HOFs, "
    "65,100 pairwise records; default τ = 0.25."
)
