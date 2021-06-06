\name{transHapmap2numeric}
\alias{transHapmap2numeric}

\title{
Convert genotypic data in character to number for training models
}
\description{
This function help users to transform genetic matrix from character format to numeric format. AA-0, Aa-1, aa-2, A is major allele and a is minor allele.
}
\usage{
transHapmap2numeric(G)
}

\arguments{
  \item{G}{
Genetic matrix of character, row represents sample and column represents SNP.
}
}
\details{
0: AA \cr
1: Aa \cr
2: aa \cr\cr
A is major allele and a is minor allele
}
\value{
A matrix for genotypic data in numeric format.
}

\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\examples{
## Not run! 
## load hapmap data (genomic data) of MZ hybrids
data(MZ)

## pre-process for input of transHapmap2numeric 
rownames(MZ) <- MZ[,1]
MZ <- MZ[,-c(1:11)]
MZ.t <- t(MZ)

## conversion
MZ.n <- transHapmap2numeric(MZ.t)
dim(MZ.t)
}
