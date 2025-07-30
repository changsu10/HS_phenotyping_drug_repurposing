# Integration of Real-World Data, Biomedical Knowledge Graph, and Single-cell Transcriptomics Reveals Comorbidities and Therapeutic Targets in Hidradenitis Suppurativa 
Hidradenitis Suppurativa is a chronic inflammatory skin disease with significant comorbidities, yet its underlying mechanisms remain poorly understood. There is a pressing need for effective therapeutic strategies. To identify comorbidities, molecular pathways, genetic factors, and potential repurposed drug candidates associated with HS using the All of Us Workbench database and knowledge graph.  We identified significant associations between HS and metabolic, dermatologic, mental, circulatory, and genitourinary conditions. Females were more severely affected than males, with strong associations to lifestyle factors such as smoking and alcohol consumption. Network analysis emphasized the role of ILB1 and TNF in HS genetics and highlighted the relationship between HS and neoplasms. Differential expression analysis revealed fibroblasts as the most strongly associated cell type. Additionally, we identified GLP-1 receptor agonists as promising repurposed drug candidates for HS. This study provides valuable insights into the genetic and clinical underpinnings of HS, highlighting key comorbidities, pathways, and genes. The findings also suggest that GLP-1 receptor agonists may be a promising therapeutic target for HS, warranting further investigation. Future research should validate these findings in diverse cohorts and explore targeted interventions to improve patient outcomes. 

![Pipeline](./figure1.png?raw=true "Title")

## EHR comorbidity analysis

## Single Cell RNA-seq Analysis
**1.integration_cellTypeAnnotation**:We integrate the NYU data and two public datasets (GSE173706 and GSE195452) using Harmony. We then annotate cell types based on marker genes.

**2.DEG.r**:We performed differential expressed gene analysis using pesudobulk approach in Seurat.

**3.GPSnet_module_generation**: We generated gene module using the GPSnet algorithm. The original code is also available at https://github.com/ChengF-Lab/GPSnet/tree/master

## Drug repurposing
**drug_repurposing**: Drug repurposing using network proximity approach. The code and data is also avaliable at https://github.com/ChengF-Lab/GPSnet/tree/master/code
