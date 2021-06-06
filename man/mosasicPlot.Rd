\name{mosaicPlot}
\alias{mosaicPlot}

\title{
Visualization of overall bins data
}
\description{
Mosasic plot for overall bins data
}
\usage{
mosaicPlot(bins,binsInfo,chr,resolution = 500,list,parentNum = 24,color,
           clust = T,methods = "ward.D2",
           dist_method = "euclidean")
}

\arguments{
  \item{bins}{
  Results of IBD analysis (bins matrix), each row represents a bin fragment as well as each column represents each sample.
}
  \item{binsInfo}{
  Data frame, including bins index, start, end, length of bins locus.
}
  \item{chr}{
  Which chromesome will be used to plot mosaic.
}
  \item{resolution}{
  To set the resolution of mosaic plot, default 500.
}
  \item{list}{
  The names of lines to visualize mosaic plot.
}
  \item{parentNum}{
  The number of parent, if color not defined, this parameter is used to auto generate color palette. 
}
  \item{color}{
  A array to define color palette.
}
  \item{clust}{ 
  Boolean values determining if lines should be hclust object.
}
  \item{methods}{
  Clustering method used. Accepts the same values as \code{\link{hclust}}.
}
  \item{dist_method}{
  The distance measure to be used for clustering. See \code{\link{dist}}
}
}

\author{
Qian Cheng, Shuqin Jiang, Xiangfeng Wang
}

\seealso{
\code{\link{bins}}
}
\examples{
## load data
data(bins)
data(binsInfo)
## color
color <- c("#DA053F","#FC0393","#C50F84","#D870D4","#DCA0DC","#4A0380",
            "#9271D9","#0414FB","#2792FC","#4883B2","#2CFFFE","#138B8A",
            "#42B373","#9BFB9C","#84FF2F","#566B32","#FED62D","#FD8A21",
            "#F87E75","#B01D26","#7E0006","#A9A9A9","#FFFE34","#FEBFCB")

mosaicPlot(bins = bins,binsInfo = binsInfo,chr = 1,resolution = 500,
                  color = color,
                  list = colnames(bins)[1:200])
}
