\name{statDesign}
\alias{statDesign}

\title{
Statistic of optimal genome
}
\description{
Statistic summary for analysis of the contribution of all candidates to optimal genome. The results directly guide the lines selection and population improvement route.
}
\usage{
statDesign(designInfo,which = "max",binsInfo,pheno,trait,output = NULL)
}
\arguments{
  \item{designInfo}{
  Outputs of \code{\link{genomeOptimization}}, a matrix consists of sample IDs regarding the fragment source among condidates at each bin locus.
}
  \item{binsInfo}{
  Data frame, including bins index, start, end, length of bins locus.
}
  \item{which}{
  A character, defining which virtual genome used for statistics, default "max".
}
  \item{pheno}{
  Phenotypic data frame, the first column describes sample names.
}
  \item{trait}{
  The names of interest trait (The trait must be included in `Pheno` data frame). 
}
  \item{output}{
  The prefix of output files regarding the summary and statistic infomation via the process of genome optimization.
}
}

\value{
A data frame regarding statistics results of virtual simulation.
  \item{Lines}{Lines}
  \item{Bins(#)}{The number (#) of bins that a line contributed to the simulated genome.}
  \item{Bins(\%)}{The number of bins that a line contributed accounting for the proportion (\%) of simulated genome.}
  \item{Fragments(\%)}{The total length of genomic fragments that a line contributed accounting for the proportion (\%) of simulated genome.}
  \item{phenotype}{The phenotypic value of the corresponding lines or their offspring}
  \item{phenotypeRank}{The phenotype rank}
  \item{Cumulative(\%)}{The cumulative percentage of fragments contributing to the simulated genome.}
If output is defined, the statistics results will be written to file with suffix "statRes.csv".
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

# 1.run example for EW trait (write result files to local)
genomeOptimization(pheno = phe,trait = "EW",bins = bins,output = "MZ_test_1404")

# 2.run example for EW trait (get the results of genomeOptimization in R)
GO_res <- genomeOptimization(pheno = phe,trait = "EW",bins = bins,output = NULL)

# 3.Statistics for genome optimization
sta_res <- statDesign(designInfo = GO_Res$overall,binInfo = binsInfo,pheno = phe,
                      trait = "EW")
}
