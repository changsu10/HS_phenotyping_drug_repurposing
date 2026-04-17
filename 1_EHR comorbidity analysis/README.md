# EHR Comorbidity Analysis

Code for the EHR-based phenotypic association analysis of hidradenitis suppurativa (HS), as described in the manuscript. Identifies HS-associated phecodes and medications using multivariable logistic regression, with validation via propensity score matching (PSM).

## Requirements

- Python 3.12.4 
- R 4.1

## Data Access

This analysis was run on the **All of Us Researcher Workbench** (https://workbench.researchallofus.org/). Individual-level data cannot be shared due to data use agreements. Approved researchers can reproduce the analysis using these scripts on the Workbench.

Reference mapping files required:
- SNOMED CT to ICD-10-CM Map (UMLS: https://www.nlm.nih.gov/research/umls/)
- Phecode Map v1.2 (https://phewascatalog.org/phecodes)

## Workflow

Run notebooks in order:

| Notebook | Description |
|---|---|
| `01_code_mapping.ipynb` | Map SNOMED/ICD codes to phecodes |
| `02_create_cohort.ipynb` | Define HS cases and controls |
| `03_cohort_preprocessing.ipynb` | Apply eligibility criteria; add demographics |
| `04_prepare_input_matrix.ipynb` | Build binary feature matrix for regression |
| `05_logistic_regression.ipynb` | Run multivariable logistic regression per phecode |
| `06_PSM.ipynb` | 1:10 propensity score matching for validation |

