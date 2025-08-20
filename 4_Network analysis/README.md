# HS Network Analysis Workflow

Network analysis workflow for Hidradenitis Suppurativa research using knowledge graphs.

## Files

- `01_data_preparation_relationships_entities.ipynb` - Prepare relationship and entity data
- `02_outlier_detection_analysis.ipynb` - Detect outliers in HS data
- `03_entity_weight_calculation.ipynb` - Calculate entity weights from beta values
- `04_network_analysis_propagation.ipynb` - Network propagation analysis

## Setup

```bash
!pip install -r requirements.txt
```


### Required Files
- **Relationship Data**: `dis_dis.csv`, `dis_gene.csv`, `dis_sys.csv`, `gene_gene.csv`
- **Entity Vocabularies**: `disease_vocab.csv`, `gene_vocab.csv`, `symptom_vocab.csv`
- **UMLS Data**: `MRCONSO.RRF`, `umls_icdall_phe_CUI.csv`
- **ICD Mappings**: `ICD10_update_zero_leading.csv`, `ICD_9_concat.csv`
- **HS Data**: `logistic_regression_results.csv` or `plot_data_df.csv`
- **GPSnet Results**: Various threshold files (e.g., `run1004_0.001.txt`)
