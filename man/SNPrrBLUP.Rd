\name{SNPrrBLUP}
\alias{SNPrrBLUP}

\title{
Genomic prediction by rrBLUP
}
\description{
Genotpye-to-phenotype prediciton via ridge regression best linear unbiased prediction (rrBLUP) model. The inputs is genotypes.
}
\usage{
SNPrrBLUP(x,y,idx1,idx2,fix = NULL,model = FALSE)
}

\arguments{
  \item{x}{
  Genotypic matrix in numberic format (See \code{\link{transHapmap2numeric}}), row represents sample well column represents feature (SNP).
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
  \item{model}{
  Boolean, if output the model, default FALSE.
}
  \item{fix}{
  A matrix containing other variables as fixed effects in mixed model.
}
}

\value{
A array for predicted value
}
\references{
Endelman, J.B. 2011. Ridge regression and other kernels for genomic selection with R package rrBLUP. Plant Genome 4:250-255. doi: 10.3835/plantgenome2011.08.0024
}
\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\seealso{
\code{\link{transHapmap2numeric}}\cr\cr

\code{\link{rrBLUP}}:\cr
\code{\link{mixed.solve}}
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

## prediction
idx1 <- sample(1:1404,1000)
idx2 <- setdiff(1:1404,idx1)
predRes <- SNPrrBLUP(MZ.n,phe$EW,idx1,idx2,fix = NULL,model = FALSE)
}
