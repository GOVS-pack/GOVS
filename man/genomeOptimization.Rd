\name{genomeOptimization}
\alias{genomeOptimization}

\title{
Genome optimization function
}
\description{
Virtual genome optimization based IBD(bins) data
}
\usage{
genomeOptimization(pheno,bins,trait,output = NULL)
}

\arguments{
\item{pheno}{Phenotypic data frame, the first column describes sample names.}
\item{bins}{ Results of IBD analysis (bins matrix), each row represents a bin fragment as well as each column represents each sample.}
\item{trait}{The names of interest trait (The trait must be included in 'Pheno' data frame}
\item{output}{The prefix of output files regarding the scheme of virtual genome.}
}

\value{
Genome optimization results
  \item{max}{Details for "optimal" simulation}
  \item{med}{Details for "moderate" simulation}
  \item{min}{Details for "poor" simulation}
  \item{overall}{Lines' names for genome optimization results of three simulated genome. Used for \code{\link{extractGenome}}}
  \item{pvalue}{P-value for ANOVA}
  
If output is defined, all above five elements will be written to five files with suffix ".max", ".med",".min", ".merge",".pvalue".
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
GO_Res <- genomeOptimization(pheno = phe,trait = "EW",bins = bins,output = NULL)
}
