# HOF–Protein Compatibility Recommender

Permanent Streamlit deployment package with 3D HOF HTML views.

## Entry point

`app.py`

## Default ranking model

`rank_fusion_25_75`

- 25% general HOF suitability percentile
- 75% protein-specific compatibility percentile

## Scoring scope

Literature, family and framework-series labels are used only for:

- filtering
- grouping
- duplicate handling
- interpretation
- post-hoc validation

They are not used as score inputs.

## 3D data

The `hof_db/` folder contains HTML visualisation files.
