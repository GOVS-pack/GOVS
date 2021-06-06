\name{GBLUP}
\alias{GBLUP}

\title{
Genomic prediction by GBLUP (kinship matrix input)
}
\description{
Genotpye-to-phenotype prediciton via genomic best linear unbiased prediction (GBLUP) model. The inputs is genotypes.
}
\usage{
GBLUP(amat,y,idx1,idx2,fix = NULL,model = FALSE)
}
\arguments{
  \item{amat}{
  Additive relationship matrix, which compute from genetic matrix (See \code{\link{A.mat}}).
}
  \item{y}{
  An numeric array of phenotype.
}
  \item{idx1}{
  An array of index for training set.
}
  \item{idx2}{
  An array of index for testing (predicted) set.
}
  \item{fix}{
  A matrix containing other variables as fixed effects in mixed model.
}
  \item{model}{
  Boolean, if output the model, default FALSE.
}
}

\value{
A array for predicted value.
}
\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255. doi: 10.3835/plantgenome2011.08.0024
}
\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\seealso{
\code{\link{rrBLUP}}:\cr
\code{\link{A.mat}}:\cr
\code{\link{kin.blup}}
}
\examples{
## Not run! 
## load hapmap data (genomic data) of MZ hybrids
data(MZ)

## load phenotypic data of MZ hybrids
data(phe)

## pre-process for G2P prediction 
rownames(MZ) <- MZ[,1]
MZ <- MZ[,-c(1:11)]
MZ.t <- t(MZ)

## conversion
MZ.n <- transHapmap2numeric(MZ.t)
dim(MZ.t)
## Additive relationship matrix compute
library(rrBLUP)
amat <- A.mat(MZ.n)

## prediction
idx1 <- sample(1:1404,1000)
idx2 <- setdiff(1:1404,idx1)
predRes <- GBLUP(amat,phe$EW,idx1,idx2,fix = NULL,model = FALSE)
}