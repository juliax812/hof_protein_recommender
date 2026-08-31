# Protein-oriented HOF prioritization

Streamlit decision-support interface synchronized to the frozen manuscript analysis.

## Final analysis universe

- 100 unique UniProt protein accessions
- 651 functionally typed retained HOFs
- 550 typed framework-series ranking units
- 65,100 HOF–protein chemistry records reconstructed at app start

## Frozen scoring model

Pairwise chemistry uses three interpretable components:

- H: hydrogen-bond complementarity
- Q: polar/charged motif abundance co-occurrence proxy
- R: aromatic/hydrophobic compatibility

Chemistry score:

`0.35 H + 0.35 Q + 0.30 R`

Final route-aware prioritization score:

`0.65 Pp_C + 0.30 S_hp + 0.05 P_G`

where `Pp_C` is the protein-specific chemistry percentile, `S_hp` is the HOF-relative differentiation percentile across the 100-protein panel, and `P_G` is the selected route-geometry percentile.

## Routes

1. Growth-mediated integration / encapsulation
2. Accessible-interface contact
3. Optional strict post-synthetic infiltration

Default geometry threshold: `tau = 0.25`.

Literature metadata, family labels, and 3D availability are displayed for traceability and optional filtering; they do not enter the score.

## Data sources used by the live app

- `HOF_MASTER_WORKBOOK_functional_group_typed_FIXED.xlsx`
- `HOF_DATABASE_FINAL_CURATED_V1.xlsx`
- `protein_descriptors_FINAL_100_UI_MIN.csv`
- `hof_db/` exact-normalized HTML structure views

The app rebuilds the complete 651 x 100 pairwise layer from these audited source tables at startup so the live UI cannot silently drift from the manuscript equations.

## Interpretation

The interface is an auditable prioritization and hypothesis-generation tool. It is not an experimentally validated compatibility-probability model.
