\name{extractGenome}
\alias{extractGenome}

\title{
Extraction & assembly module
}
\description{
Extracting genome fragment from candidates based the results of genome optimization and then assembling all fragments so that produce optimized genome(virtual genome).
}
\usage{
extractGenome(hmp,binsInfo,ID = NULL,designInfo,output = NULL,
              bins,extractContent ="Genotype")
}
\arguments{
  \item{hmp}{
  The genetic data in hapmap format.
}
  \item{binsInfo}{
  Data frame, including bins index, start, end, length of bins locus.
}
  \item{ID}{
  A character array regarding sample IDs for `hmp`, if NULL, the `hmp` data must involve header.
}
  \item{designInfo}{
  Outputs of \code{\link{genomeOptimization}}, a matrix consists of sample IDs regarding the fragment source among condidates at each bin locus.
}
  \item{output}{
  The prefix of output files regarding the assembled genome.
}
  \item{bins}{
  Results of IBD analysis (bins matrix), each row represents a bin fragment as well as each column represents each sample.
}
  \item{extractContent}{
  Character, the content of virtual genome, "Bin" for bin source well "Genotype" for genetic data, default "Genotype".
}
}

\value{
A matrix involves of three optimal genomes (optimal: max, moderate: med and poor: min)
If output is defined, the assembled genome will be written to files with suffix ".G" or ".bin".
}

\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\seealso{
\code{\link{GOVS}}
}
\examples{
## Not run !
## load test data
# Phenotypic data:
data(phe)
# bins data:
data(bins)
# genomic data:
data(MZ)
# bins info 
data(binsInfo)

# 1.run example for EW trait (write result files to local)
genomeOptimization(pheno = phe,trait = "EW",bins = bins,output = "MZ_test_1404")

# 2.run example for EW trait (get the results of genomeOptimization in R)
GO_Res <- genomeOptimization(pheno = phe,trait = "EW",bins = bins,output = NULL)

# 3.extract and assemble virtual genome (genomic genome)
Vgenome <- extractGenome(hmp = MZ,binsInfo = binsInfo,designInfo = GO_Res$overall)

# 4.extract and assemble virtual genome (bin source genome)
Vbin <- extractGenome(hmp = MZ,binsInfo = binsInfo,designInfo = GO_Res$overall,
                      bins = bins,extractContent = "BINsource")
}