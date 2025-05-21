# 🧬 SynTact: Synthetic Targeting AI for Cancer Therapy

## Overview
SynTact is an AI-powered research platform designed to accelerate the discovery of cancer-specific CAR-T cell therapy targets. Our mission is to **democratize access to precision immunotherapy tools** by providing an open-source, modular system that enables researchers to identify tumor-specific antigens and eventually generate CAR DNA constructs to target them.

Our initial focus is on **ovarian cancer**, with the broader goal of extending support to multiple solid and hematologic malignancies.

---

## 🚀 Project Goals

1. **Identify tumor-specific surface antigens** for a given cancer type using ML-powered analysis of public omics data
2. **Score and rank targets** based on expression, safety, localization, and immunogenicity
3. **Create a user-friendly web interface** to help researchers upload data and receive actionable CAR-T target suggestions
4. **(Future)** Generate or evaluate **DNA sequences** (CAR constructs) to optimize binding efficacy and persistence
5. Build a **collaborative tool** for research labs, biotech companies, and academic consortia

---

## 📍 Current Focus: Track A — Target Discovery

We are currently in **Track A** of development: identifying ideal antigen targets for CAR-T therapy in ovarian cancer.

### ✅ Stage 1: Research Foundation
- Conducted a **deep literature review** on CAR-T targets for ovarian cancer (MUC16, FOLR1, Mesothelin, B7-H3, etc.)
- Compiled **expression and localization databases**:
  - Tumor data: TCGA-OV (Ovarian Serous Cystadenocarcinoma)
  - Normal data: GTEx
  - Surface annotations: Human Protein Atlas (HPA), Cell Surface Protein Atlas (CSPA)
- Identified ML modeling strategies:
  - Classical classifiers (e.g., XGBoost) on expression-based features
  - Protein embedding models (e.g., ProtBERT, ESM-2)
  - Ranking models combining safety and efficacy features
- Selected **evaluation metrics**: ROC-AUC, PR-AUC, Mean Average Precision, expert validation

### 🧪 Stage 2: Development Pipeline

We are now beginning to build out our **modular ML pipeline and web interface**, including:

- Data ingestion: Parsing TCGA + GTEx matrices
- Feature extraction: Surface flag, log fold-change, immunogenicity signals
- Initial model: Ranking or binary classifier with explainability (e.g., SHAP)
- MVP interface (Streamlit or Gradio):
  - Upload or select cancer dataset
  - View ranked CAR-T targets
  - Inspect gene details (expression, safety, known trials)

---

## 📁 Repository Structure

```
SynTact/
├── research/               # Research notes, figures, and presentations
├── data/                   # Raw and processed datasets (TCGA, GTEx, HPA, etc.)
├── notebooks/              # Jupyter notebooks for EDA, modeling, prototyping
├── src/                    # Core pipeline code (data loading, scoring, modeling)
├── app/                    # Web interface (Streamlit, Gradio, etc.)
├── reports/                # Research briefs, experimental results
├── README.md               # Project overview and documentation
├── requirements.txt        # Dependencies
├── .gitattributes
└── .gitignore              # Git ignore patterns

```

---

## 🤝 Join the Mission
We are looking for collaborators — researchers, students, and engineers — who want to help:
- Expand SynTact to other cancers (glioblastoma, breast, etc.)
- Build the DNA sequence generation module (Track B)
- Contribute to modeling, UI, or experimental validation

If you’re interested, feel free to open an issue or email us directly!

---

## 📌 Citation & Acknowledgements
SynTact is built upon the work of the open cancer research community. We acknowledge datasets from TCGA, GTEx, HPA, and prior publications on CAR-T therapy.

> This project is a student-led research initiative with aspirations for clinical translation and academic collaboration.