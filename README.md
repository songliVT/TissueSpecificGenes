# TissueSpecificGenes

Identification of tissue specific genes using weighted entropy method, which is originally described in ROKU paper PMID:16764735.

**Input data**: matrix of gene expression data, each row is one gene, each column represent one tissue type. Gene expression average FPKM from normalized RNAseq read counts.

To run the script, use the following command:

> Rscript GetTissueSpecificGenes.R FPKM.example.csv [...]

Parameters:

**inputfn**: input file name, FPKM.example.csv is a sample tissue specific expression data. 

**outputfn**: output file name, this is a matrix of the same deminsion of the original expression matrix. entries in this matrix take values of (1,0,-1). 1 means this gene is signifiantly highly expressed in this tissue. -1 means significantly repressed in this tissue. 0 mean not significant. default value is 'TissueSpecificData.csv'
  
**lowexp**: This is the threshold used to filter lowly expressed gene.  If a gene expressed below this threshold in all conditions, this gene will not be included in the analysis.  Default value is 0.05. 

**bgfold**: Fold change threshold.  This is used to define background gene set. Default value is 2. 

**bgmedian**: Another threshold to filter relatively lowly expressed genes. Genes with median expression values smaller than this threshold will not be included in the results. Default value is 0.5. 

**pvalue**: This is the p value threshold. This is based on emperical distribution of the background gene set. Default value is 0.001

