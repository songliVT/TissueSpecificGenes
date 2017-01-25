# this is the working version 2 of the calling tissue specific genes.
# use expression FPKM

###########################################################
# step 1. load expression from gene FPKM in root data
###########################################################
#

args = commandArgs(trailingOnly=TRUE)
# usage: rscript GetTissueSpecificGenes_012517.R FPKM.example.csv

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  outputfn='TissueSpecificData.csv'
  lowexp<-0.05 # threshold for low expression values
  bgfold=2 # fold change used to define background gene set
  bgmedian=0.5 # threshold for median expression level
  pvalue=0.001 # p value
}

inputfn<-args[1]

# inputfn<-'FPKM.example.csv' # input data matrix for gene expression

geneFPKM<-read.table(inputfn,sep=',',as.is=T,header=T,row.names = 1)
ncelltype<-ncol(geneFPKM) # get number of tissue/cell types

# a cell type specific genes can only be specific to less than half of the cell types
nmax<-floor(ncelltype/2)

# remove gene with FPKM < 0.05 in all cell types. 
tmp<-geneFPKM;tmp[tmp>=lowexp]=1;tmp[tmp<lowexp]=0;tmpsum<-apply(tmp,1,sum)
geneFPKMkeep<-geneFPKM[tmpsum>0,]

###########################################################
# step 2. make the entropy calculation
###########################################################

source('ROKU_modified.R')

###########################################################
# get back ground genes for estimating the distribution of entropy.
# get fold change for genes
###########################################################
getfold<-function(x)
{
  x<-x+0.5 # add pseudo count to avoid 0s
  return(max(x)/min(x))
}
geneFold<-apply(geneFPKM,1,getfold)
geneMedian<-apply(geneFPKM,1,median)

###########################################################
# get background genes as those withing 2 fold change and median expressoin > 0.5 
# this filter is to not use really lowly expressed genes as inputs.
###########################################################


bggene<-names(geneFold[geneFold<bgfold&geneMedian>bgmedian]) # get background genes as those withing 2 fold change
bggeneFPKM<-geneFPKM[bggene,]
print(paste('number of genes in the background set:', 
      nrow(bggeneFPKM)))

# get the WH for for shaved exprssion
bgWHmat<-weightcenterentropy_log_shave(bggeneFPKM, ntoshave=nmax)

# get threshold at a pvalue, default pvalue = 0.001

TR_0001<-apply(bgWHmat,2,quantile,prob=pvalue)

# run the calculation for all genes
allWH_shave<-weightcenterentropy_log_shave(bggeneFPKM)

###########################################################
# Get tissue specific genes.
###########################################################

geneCellSpecific_P0001<-weightcenterentropy_log_shave_identify(geneFPKM,TR_0001)
#save(geneCellSpecific_P0001,file=outputfn)

write.table(geneCellSpecific_P0001$index,file=outputfn,sep=',')

