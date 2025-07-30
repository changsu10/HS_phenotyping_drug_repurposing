## Single Cell RNA-seq Analysis
1. **1.integration_cellTypeAnnotation**:We integrate the NYU data and two public datasets (GSE173706 and GSE195452) using Harmony. We then annotate cell types based on marker genes.
2. **2.DEG.r**:We performed differential expressed gene analysis using pesudobulk approach in Seurat.
3. **3.GPSnet_module_generation**: We generated gene module using the GPSnet algorithm. The original code is also available at https://github.com/ChengF-Lab/GPSnet/tree/master
4. **4.drug_repurposing**: Drug repurposing using network proximity approach. The code and data is also avaliable at https://github.com/ChengF-Lab/GPSnet/tree/master/code