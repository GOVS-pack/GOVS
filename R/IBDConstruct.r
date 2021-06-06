################################ IBD analysis and map construct #############
################################ function ###################################
#' @importFrom grDevices hcl
#' @importFrom stats dist hclust
#' @importFrom utils setTxtProgressBar txtProgressBar
pmatrix=function(k,snp_parent,numparent=24,p1=0.97){
  pres=pblapply(1:k,function(i){
    # p1=0.97               
    marker<-vector(length=0)
    
    marker[1]<-"marker"
    marker[2]<-snp_parent[i,1]
    marker[3]<-3
    marker[4]<-snp_parent[i,3]
    marker[5:(2+numparent)]<-""
    
    allele1<-vector(length=0)
    allele1[1]<-"allele"
    allele1[2]<-"NN"
    allele1[3:(2+numparent)]<-round(1/numparent, 3)	
    
    allele2<-vector(length=0)
    allele3<-vector(length=0)
    
    a2<-unlist(strsplit(snp_parent[i,2],"/"))[1]
    a3<-unlist(strsplit(snp_parent[i,2],"/"))[2]
    aa2<-paste(a2,a2,sep = "")
    aa3<-paste(a3,a3,sep = "")
    allele2[1]<-"allele"
    allele2[2]<-a2
    
    allele3[1]<-"allele"
    allele3[2]<-a3
    
    n1<-(snp_parent[i,5:(4+numparent)]==aa2)
    n2<-(snp_parent[i,5:(4+numparent)]==aa3)
    
    
    if (sum(n1+n2)==numparent) 
    {allele2[2+which(snp_parent[i,5:(4+numparent)]==aa2)]<-p1/sum(n1)
    allele2[2+which(snp_parent[i,5:(4+numparent)]!=aa2)]<-(1-p1)/sum(n2)	
    allele3[2+which(snp_parent[i,5:(4+numparent)]==aa3)]<-p1/sum(n2)
    allele3[2+which(snp_parent[i,5:(4+numparent)]!=aa3)]<-(1-p1)/sum(n1)}
    if (sum(n1+n2)!=numparent)
    { p<-which(n1+n2==0)
    n1[1,p]<-1/2
    n2[1,p]<-1/2
    allele2_m<-n1*p1/sum(n1)
    allele3_m<-n2*p1/sum(n2)
    allele2_m[(allele2_m==0)]<-(1-p1)/sum(allele2_m==0)
    allele3_m[(allele3_m==0)]<-(1-p1)/sum(allele3_m==0)
    allele2[3:(2+numparent)]<-allele2_m
    allele3[3:(2+numparent)]<-allele3_m
    }
    
    pres<-data.frame(marker,allele2,allele3,allele1)
    
    return(t(pres))
  })
  return(do.call(rbind,pres))
}

omatrix=function(p,snp,nummarker,G,numparent=24,rou){
  # Z<-array(0,c(dim(snp)[1],numparent,2))
  P<-matrix(0,nummarker,ncol=numparent)        
  Q<-matrix(0,nummarker,ncol=numparent) 
  
  if(snp[1]==paste(p[2,2],p[2,2],sep="")){
    P[1,]<-as.numeric(p[2,3:(numparent+2)])
  } else if(snp[1]==paste(p[3,2],p[3,2],sep="")){
    P[1,]<-as.numeric(p[3,3:(numparent+2)])} else
    {P[1,]<-as.numeric(p[4,3:(numparent+2)])} 
  
  
  if(snp[nummarker]==paste(p[2+4*(nummarker-1),2],p[2+4*(nummarker-1),2],sep="")){
    Q[nummarker,]<-as.numeric(p[2+4*(nummarker-1),3:(numparent+2)])
  } else if(snp[nummarker]==paste(p[3+4*(nummarker-1),2],p[3+4*(nummarker-1),2],sep="")){
    Q[nummarker,]<-as.numeric(p[3+4*(nummarker-1),3:(numparent+2)])} else
    {Q[nummarker,]<-as.numeric(p[4+4*(nummarker-1),3:(numparent+2)])} 
  
  pb = txtProgressBar(min = 0, max = nummarker, initial = 0,style=3) 
  for (i in 1:(nummarker-1)){
    pai<-matrix(0,nrow=1,ncol=numparent)
    Qpai<-matrix(0,nrow=1,ncol=numparent)
    
    # determined by generations
    if (rou[i]>1E-6)     
    {lamda<-G*rou[i]/50} else  
    {lamda<-G*(2E-8)}
    #  lamda<-G*rou[i]/1E4   
    
    r<-matrix((1-exp(-lamda))/numparent,nrow=numparent,ncol=numparent)
    diag(r)<-exp(-lamda)+(1-exp(-lamda))/numparent
    
    if (rou[nummarker-i]>1E-6)
    {Qlamda<-G*(rou[nummarker-i]/50)} else
    {Qlamda<-G*(2E-8)}
    #  Qlamda<-G*rou[nummarker-i]/1E4   
    
    Qr<-matrix((1-exp(-Qlamda))/numparent,nrow=numparent,ncol=numparent)
    diag(Qr)<-exp(-Qlamda)+(1-exp(-Qlamda))/numparent
    
    
    if(snp[i+1]==paste(p[2+4*i,2],p[2+4*i,2],sep=""))
    {pai<-as.numeric(p[2+4*i,3:(numparent+2)])} else if
    (snp[i+1]==paste(p[3+4*i,2],p[3+4*i,2],sep=""))
    {pai<-as.numeric(p[3+4*i,3:(numparent+2)])} else
    {pai<-as.numeric(p[4+4*i,3:(numparent+2)])}
    
    if(snp[nummarker-i]==paste(p[2+4*(nummarker-1-i),2],p[2+4*(nummarker-1-i),2],sep=""))
    {Qpai<-as.numeric(p[2+4*(nummarker-1-i),3:(numparent+2)])} else if
    (snp[nummarker-i]==paste(p[3+4*(nummarker-1-i),2],p[3+4*(nummarker-1-i),2],sep=""))
    {Qpai<-as.numeric(p[3+4*(nummarker-1-i),3:(numparent+2)])} else
    {Qpai<-as.numeric(p[4+4*(nummarker-1-i),3:(numparent+2)])}
    
    
    f=apply(r,2,function(m){
      m*pai/as.vector(m%*%pai)})
    
    P[i+1,]<-P[i,]%*%t(f)
    
    Qf=apply(Qr,2,function(m){
      m*Qpai/as.vector(m%*%Qpai)})
    Q[nummarker-i,]<-Q[nummarker-i+1,]%*%t(Qf)
    
    setTxtProgressBar(pb,i)
  }
  # Z[,,1]<-P
  # Z[,,2]<-Q
  Z=matrix((P+Q)/2,nrow=nrow(P))
  return(Z)
}

binConstruct <- function(Z,markerInfo,threshold=NULL,omit = T){
  binRes <- c()
  for (i in 1:nrow(Z)) {
    x <- Z[i,]
    a <- which(x == max(x))
    if (!is.null(threshold)) {
      if (length(a) > 1 | max(x) < threshold) {
        binRes <- c(binRes,NA)
      }else{
        binRes <- c(binRes,a)
      }
      binRes[which(is.na(binRes))] <- 0
    }else{
      if (length(a) > 1) {
        binRes <- c(binRes,NA)
      }else{
        binRes <- c(binRes,a)
      }
      binRes[which(is.na(binRes))] <- 0
    }
  }
  
  if(omit & (length(which(binRes == 0) > 0))){
    idx.omit <- which(binRes == 0)
    binRes <- binRes[-idx.omit]
    markerInfo <- markerInfo[-idx.omit,]
  }
  
  dffRes <- diff(binRes)
  idx <- which(dffRes != 0)
  
  sta <- numeric()
  for (i in 1:length(idx)) {
    sta <- c(sta,mean(as.numeric(c(markerInfo[idx[i]+1,4],markerInfo[idx[i],4]))))
  }
  
  end <- c(sta,as.numeric(markerInfo[nrow(markerInfo),4]))
  sta <- c(as.numeric(markerInfo[1,4]),sta)
  
  binsInfo <- data.frame(bin = 1:length(sta),chr = as.numeric(markerInfo[length(sta),3]),sta = sta, end = end, len = end - sta)
  bin <- binRes[c(idx,idx[length(idx)]+1)]
  bin <- matrix(bin,ncol = 1)
  
  binsData <- list(bin = bin,binsInfo = binsInfo)
  binsData
}

######################################
#' @export IBDConstruct
IBDConstruct <- function(snpParents,snpProgeny,markerInfo,q,rou,G,threshold = NULL,omit = T){
  if (nrow(snpParents) != length(snpProgeny)) {
    stop("the markers number of parents  are not euqal to progeny!")
  }
  
  # compute priority probability matrix
  k <- nrow(snpParents)
  numparent <- ncol(snpParents)
  snp_parent <- cbind(markerInfo,snpParents)
  message("compute priority probability matrix...")
  pmatrixP <- pmatrix(k,snp_parent,numparent=numparent,p1=q)
  # compute posterior probalility matrix of progeny 
  message("compute posterior probalility matrix of progeny...")
  
  # if(is.null(rou)){
  #   rou <- rep(0,k-1)
  # }
  Zmatrix <-omatrix(p = pmatrixP,snp = snpProgeny,nummarker = k,numparent = numparent,G=G,rou = rou)
  ## construct bins map
  message("\nconstruct bins map...")
  binC <- binConstruct(Z = Zmatrix,markerInfo = markerInfo,threshold=threshold,omit = omit)
  
  return(binC)
}
