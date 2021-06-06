\name{binsPlot}
\alias{binsPlot}

\title{
Visualization of IBD map results
}
\description{
Visualization of IBD map results
}
\usage{
binsPlot(IBDRes,color,parentInfo,parentNum)
}

\arguments{
  \item{IBDRes}{
  The results of IBDConstruct, see \code{\link{IBDConstruct}}.
}
  \item{color}{
  A named vector for defining color of parents.
}
  \item{parentInfo}{
  A named vector for defining label of parents.
}
  \item{parentNum}{
  The number of parents.
}
}

\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\seealso{
\code{\link{IBDConstruct}}
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
  
  ## plot
  # color
  color <- c("#DA053F","#FC0393","#C50F84","#D870D4","#DCA0DC","#4A0380",
             "#9271D9","#0414FB","#2792FC","#4883B2","#2CFFFE","#138B8A",
             "#42B373","#9BFB9C","#84FF2F","#566B32","#FED62D","#FD8A21",
             "#F87E75","#B01D26","#7E0006","#A9A9A9","#FFFE34","#FEBFCB")
  names(color) <- 1:24

  # parent label
  parentInfo <- c("5237","E28","Q1261","CHANG7-2","DAN340","HUANGC","HYS",
                  "HZS","TY4","ZI330","ZONG3","LX9801","XI502","81515",
                  "F349","H21","JI853","JI53","LV28","YUANFH","SHUANG741",
                  "K12","NX110","ZONG31")
  names(parentInfo) <- 1:24
  
  # plot
  binsPlot(IBDRes,color,parentInfo,24)
}
