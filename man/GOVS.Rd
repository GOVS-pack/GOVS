\name{GOVS}
\alias{GOVS}

\title{
Genome optimization via virtual simulation
}
\description{
One-stop function for a complete progress of genome optimization.
}
\usage{
 GOVS(hmp,ID = NULL,pheno,trait,bins,binsInfo,which = "max",
      output = NULL,module = "DES",designInfo,extractContent = "Genotype")
}

\arguments{
  \item{hmp}{The genetic data in hapmap format.}
  \item{ID}{A character array regarding sample IDs for "hmp", if NULL, the "hmp" data must involve header.}
  \item{pheno}{Phenotypic data frame, the first column describes sample IDs.}
  \item{trait}{The names of interest trait (The trait must be included in "Pheno" data frame).}
  \item{bins}{Results of IBD analysis (bins matrix), each row represents a bin fragment as well as each column represents each sample.}
  \item{binsInfo}{Data frame, including bins index, start, end, length of bins locus.}
  \item{which}{A character, defining which virtual genome used for statistics, default "max".}
  \item{output}{The prefix of output files.}
  \item{module}{Character represents the module combination for analysis, default "DES", "D" for genome optimization module, "E" for extraction & assembly module. "S" for statistic module. "D","E","S","DE","DES" and "ES" are alternative for different module combinations. Note that different combination need different essential inputs.}
  \item{designInfo}{Data frame, the results of genome optimization module, it's necessary for "ES" and "S" module.}
  \item{extractContent}{Character, the content of virtual genome, "Bin" for bin source well "Genotype" for genetic data, default "Genotype".}
}

\value{
A list regarding defined module.\cr
If output is defined, all files will be written to file with prefix of "output".\cr\cr
Details see:\cr
\code{\link{genomeOptimization}}\cr
\code{\link{extractGenome}}\cr
\code{\link{statDesign}}
}

\author{
  Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\seealso{
\code{\link{genomeOptimization}}\cr
\code{\link{extractGenome}}\cr
\code{\link{statDesign}}
}
\examples{
## Not run !
## load test data
# Phenotypic data:
data(phe)
# genomic data:
data(MZ)
# bins data:
data(bins)
# bins information:
data(binsInfo)
# example for one-stop solution for GOVS
GOVS_res <- GOVS(MZ,pheno = phe,trait = "EW",which = "max",bins = bins,
                 binsInfo = binsInfo,module = "DES")
}