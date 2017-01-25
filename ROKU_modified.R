# weighted entropy for a matrix.
#

weightcenterentropy_log_shave<-function(x,ntoshave=7)
{#x is a matrix. 
  # use centered entropy + entropy without center.
  # sort the expression values from high to low and 
  # and calculate the entropy for the rest of the values
  # 
  # ntoshave is the number of samples to remove
  # total output diminsion is ntoshare+1 columns 
  # this is used to get the empirical distribution of the reads
  
  delta=0.05
  ngene<-nrow(x)
  genenames<-rownames(x)
  ncond<-ncol(x)
  x<-x+delta # penalize low expression data
  x1<-log2(x)
  tmp.mad<-sqrt(apply(x1,1,var)) # get standard deviation 
  reg.mad<-quantile(tmp.mad,0.05) # 5 percentile
  tmp.median<-apply(x1,1,median)
  cenexp<-abs((x1-tmp.median)/(tmp.mad+reg.mad))
  
  #sort the expression
  sort_cenexp<-matrix(0,ncol=ncond,nrow=ngene)
  sort_x<-matrix(0,ncol=ncond,nrow=ngene)
  allorder<-NULL
  
  for(i in 1:ngene){
    tmporder<-order(cenexp[i,])
    sort_x[i,]<-unlist(x[i,tmporder])
    sort_cenexp[i,]<-unlist(cenexp[i,tmporder])
  }
  
  ent<-function(y){
    y<-y+1e-06;p<-y/sum(y);e<--sum(p*log2(p))
  }

  ent_ori<-function(y){ #entroyp in original scale
    y <- y+delta # this is to penalize genes with low expression levels
    p <- y / sum(y)
    e <- - sum(p * log2(p))
  }
  
  outdim<-ntoshave+1 #output diminsion is the number of cell types minus one
  outmat<-matrix(0,ncol=outdim,nrow=ngene) 
  rownames(outmat)<-rownames(x)
  for(i in 1:outdim)
  { # caculate the WH for different conditions.
    nskip=i-1
    tmpcenexp<-sort_cenexp[,1:(ncond-nskip)]
    tmpx<-sort_x[,1:(ncond-nskip)]
    outent<-apply(tmpcenexp,1,ent)
    outent1<-apply(tmpx,1,ent_ori)
  
    maxH<-max(outent1); minH<-min(outent1)
    pH<-(outent1-minH)/(maxH-minH)
    outdat<-outent*(1-pH)+outent1*pH
    outmat[,i]<-outdat
  }
  colnames(outmat)<-c(1:outdim)-1
  
  return(outmat)
}


weightcenterentropy_log_shave_identify<-function(x,thre)
{#x is a matrix. 
 # calculate the WH by iterartively remove the most highly expressed tissue type.
 # thre is a threshold vector, 
 # names of the vector are the thresholds.
 # 
  ntoshave <- length(thre)-1
  # preprocessing the expression data
  delta=0.05
  ngene<-nrow(x)
  ncond<-ncol(x)
  x<-x+delta # penalize low expression data
  x1<-log2(x)
  tmp.mad<-sqrt(apply(x1,1,var))
  reg.mad<-quantile(tmp.mad[tmp.mad>0],0.05) # 5 percentile
  tmp.median<-apply(x1,1,median)
  oricenexp<-(x1-tmp.median)/(tmp.mad+reg.mad) # centered expression at original scale
  cenexp<-abs(oricenexp)
  
  #sort the expression
  sort_cenexp<-matrix(0,ncol=ncond,nrow=ngene)
  sort_x<-matrix(0,ncol=ncond,nrow=ngene)
  allorder<-matrix(0,ncol=ncond,nrow=ngene)
  celltypes<-colnames(cenexp)
  for(i in 1:ngene){
    tmporder<-order(cenexp[i,])
    celltypeorder<-celltypes[tmporder]
    allorder[i,]<-tmporder
    sort_x[i,]<-unlist(x[i,tmporder])
    sort_cenexp[i,]<-unlist(cenexp[i,tmporder])
  }
  
  ent<-function(y){
    y<-y+1e-06;p<-y/sum(y);e<--sum(p*log2(p))
  }
  
  ent_ori<-function(y){ #entroyp in original scale
    y <- y+delta # this is to penalize genes with low expression levels
    p <- y / sum(y)
    e <- - sum(p * log2(p))
  }
  
  out_index<-matrix(0,ncol=ncond,nrow=ngene)
  colnames(out_index)<-celltypes
  rownames(out_index)<-rownames(x)
  outdim<-ntoshave+1 #output diminsion is the number of cell types minus one
  
  outmat<-matrix(0,ncol=outdim,nrow=ngene)  # this is the matrix with shaved WH values
  rownames(outmat)<-rownames(x)

  for(i in 1:outdim)
  { # caculate the WH for different conditions.
    nskip=i-1
    tmpcenexp<-sort_cenexp[,1:(ncond-nskip)]
    tmpx<-sort_x[,1:(ncond-nskip)]
    outent<-apply(tmpcenexp,1,ent)
    outent1<-apply(tmpx,1,ent_ori)
    
    maxH<-max(outent1); minH<-min(outent1)
    pH<-(outent1-minH)/(maxH-minH)
    outdat<-outent*(1-pH)+outent1*pH
    outmat[,i]<-outdat
  }
  colnames(outmat)<-c(1:outdim)-1
  
  # make a matrix of threshold.
  thremat<-t(matrix(rep(thre,nrow(outmat)),nrow=length(thre))) 
  
  signmat<-thremat-outmat # outmat is smaller than thremat then the gene is significant
  signmat[signmat>=0]=1;signmat[signmat<=0]=0
  
  # go throuhg each gene and find the significant genes.
  outnsig<-NULL # number of significant genes
  for(i in 1:ngene)
  {
    #print(i)
    thisexporicen<-oricenexp[i,] # centered expression at original scale
    thisorder<-allorder[i,] # order of the samples.
    tmporiexp<-rev(thisexporicen[thisorder])
    thissig<-signmat[i,]
    if(sum(thissig)==ncol(signmat)){nsig = ncol(signmat)}# all significant.
    else{nsig<-min(which(thissig!=1))-1} #some are not significant
    outnsig<-c(outnsig,nsig)
    if (nsig==0){next} # do nothing is no significant gene
    else{
      tmpsigkeep<-tmporiexp[1:nsig]
      upnames<-names(tmpsigkeep)[tmpsigkeep>0]
      downnames<-names(tmpsigkeep)[tmpsigkeep<0]
      out_index[i,upnames]=1 # assign up regulated genes
      out_index[i,downnames]=-1 # assign down regulated genes
    }

  }
  out<-NULL
  out$index<-out_index
  out$oricenexp<-oricenexp
  out$signmat<-signmat
  out$WH<-outmat
  out$thre<-thre
  return(out)
}








