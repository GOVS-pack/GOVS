\name{IBDConstruct}
\alias{IBDConstruct}

\title{
Construct IBD map
}
\description{
A IBD map was constructed of contributions from the parents onto the progeny lines using a hidden Markov model (HMM). 
}
\usage{
IBDConstruct(snpParents,snpProgeny,markerInfo,q,
             rou,G,threshold = NULL,omit = T)
}

\arguments{
  \item{snpParents}{
  A matrix for the parents' genotype, lines in column and marker in row.
}
  \item{snpProgeny}{
  An array for the progeny' genotype, marker number must equal to snpParents.
}
  \item{markerInfo}{
  A matrix or dataframe with four cols(marker ID, allele, chromsome and physical position) regarding genotypic information.
}
  \item{q}{
  The quality of sequencing, range 0 to 1 to define the quality of marker.
}
  \item{rou}{
  Correlations between any pairs of flanking markers, that estimated with the offspring-LD level after corrected by parent-LD level, it can be obtained by genetic location.
}
  \item{G}{
  Generations that the offsprings decented from the parents.
}
  \item{threshold}{
  The threshold of posterior.
}
  \item{omit}{
  Whether to omit untraceable segments, default True.
}
}
\details{
The details see: \cr
https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1930-x#MOESM14 \cr
CUBIC: an atlas of genetic architecture promises directed maize improvement \cr
https://github.com/heroalone/HMM-IBD
}
\value{
A list regarding Constructed bin map.\cr
\item{bin}{Results of IBD analysis (bins matrix), each row represents a bin fragment.}
\item{binsInfo}{Data frame, including bins index, start, end, length of bins locus.}
}
\references{
Liu H J, Wang X, Xiao Y, et al. CUBIC: an atlas of genetic architecture promises directed maize improvement[J]. Genome biology, 2020, 21(1): 1-17. \cr

https://github.com/heroalone/HMM-IBD
}

\examples{
  ## load example data
  data(IBDTestData)
  ## compute rou from genetic position
  rou = IBDTestData$posGenetic
  rou = diff(rou)
  rou = ifelse(rou<0,0,rou)
  ## constract IBD map of chr10 for one progeny
  IBDRes <- IBDConstruct(snpParents = IBDTestData$snpParents,
  markerInfo = IBDTestData$markerInfo,
  snpProgeny = IBDTestData$snpProgeny,q = 0.97,G = 9,rou = rou)
}
