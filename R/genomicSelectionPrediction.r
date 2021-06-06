## GBLUP
#' @import rrBLUP
#' @export GBLUP
GBLUP <- function(amat,y,idx1,idx2,fix = NULL,model = FALSE){
  amat <- amat
  ytrain <- as.numeric(y[idx1])
  #yreal <- as.numeric(y[idx2])
  
  #misIdx <- which(is.na(yreal) == T)
  y[idx2] <- NA
  y <- as.numeric(y)
  if (is.null(fix)){
    data <- data.frame(y=y,gid=rownames(amat))
    ans <- kin.blup(data=data,geno="gid",pheno="y",K=amat)
  }else{
    class(fix) <- "numeric"
    data <- data.frame(y=y,gid=rownames(amat),fix = fix)
    ans <- kin.blup(data=data,geno="gid",pheno="y",K=amat,fixed = names(data)[3:ncol(data)])
  }
  ypred <- ans$pred[idx2]
  if (model){
    finalres <- list(model = ans,predRes = ypred)
  }else{
    finalres <- ypred
  }
  #yreal <- yreal[-misIdx]
  #yreal <- as.numeric(yreal)
  #ypred <- ypred[-misIdx]
  #names(yreal) <- names(ypred) <- rownames(amat)[idx2][-misIdx]
  
  #evalres <- evalfunc(yreal,ypred,length(yreal))
  #finalres <- list(resMat = cbind(yreal,ypred),evalres = evalres)
  finalres
  #finalres
}

## rrBLUP
#' @import rrBLUP
#' @export SNPrrBLUP
.trainModel_RRBLUP <- function(markerMat, phenVec,X = NULL){
  phen_answer<-mixed.solve(phenVec, Z=markerMat, K=NULL, SE = FALSE, return.Hinv=FALSE,X = X)
  beta <- phen_answer$beta
  phD <- phen_answer$u
  e <- as.matrix(phD)
  return( list(beta = beta, e = e, phD = phD) )
}
SNPrrBLUP <- function(x,y,idx1,idx2,fix = NULL,model = FALSE){
  
  trainG <- x[idx1,]
  testG <- x[idx2,]
  
  if(is.null(fix)){
    ytrain <- as.numeric(y[idx1])
    res <- .trainModel_RRBLUP(markerMat = trainG,phenVec = ytrain,X = NULL)
    ypred <- testG %*% res$e
    ypred <- ypred[,1] + as.numeric(res$beta)
  }else{
    class(fix) <- "numeric"
    trainfix <- fix[idx1,]
    testfix <- fix[idx2,]
    
    ytrain <- as.numeric(y[idx1])
    yreal <- y[idx2]
    #misIdx <- which(is.na(yreal) == T)
    
    res <- .trainModel_RRBLUP(markerMat = trainG,phenVec = ytrain,X = trainfix)
    ypred <- testG %*% res$e
    beta <- matrix(res$beta,nrow = ncol(fix))
    beta <- testfix %*% beta
    ypred <- ypred[,1] + beta
  }
  #yreal <- yreal[-misIdx]
  #yreal <- as.numeric(yreal)
  #ypred <- ypred[-misIdx]
  #names(yreal) <- names(ypred) <- rownames(x)[idx2][-misIdx]
  
  #evalres <- evalfunc(yreal,ypred,length(yreal))
  #finalres <- list(resMat = cbind(yreal,ypred),evalres = evalres)
  
  if (model){
    finalres <- list(model = res,predRes = ypred)
  }else{
    finalres <- ypred
  }
  
  finalres
}
