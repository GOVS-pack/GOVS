## plot bin 
#' @export
#' @importFrom scales pretty_breaks
#' @import ggplot2
#' @export binsPlot
binsPlot <- function(IBDRes,color,parentInfo,parentNum){
  if (missing(color)) {
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    color <- gg_color_hue(parentNum)
    names(color) <- 1:parentNum
  }
  
  if (missing(parentInfo)) {
    parentInfo <- 1:parentNum
  }
  
  df <- data.frame(sta=IBDRes$binsInfo$sta,end =IBDRes$binsInfo$end,
                   ysta=IBDRes$bin-0.5,yend=IBDRes$bin+0.5,
                   size=IBDRes$binsInfo$len,bin=IBDRes$bin
  )
  df[,c("sta","end","size")] <- df[,c("sta","end","size")]/1000000
  
  p <- ggplot(df) +
    geom_rect(aes(xmin = sta, xmax = end, ymin = ysta, ymax = yend,size = size,fill = as.character(bin)))+
    scale_fill_manual(values = color)+
    scale_y_continuous(breaks=seq(1,parentNum,1),labels = parentInfo)+
    scale_x_continuous(expand=c(0,1),breaks = scales::pretty_breaks(n = 5))+
    labs(
      x = "Physical position (Mb)",
      y = "Parents"
    )+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",color = "black"),
          axis.text = element_text(size = 12),axis.title=element_text(size=12,face = "bold"),
          legend.position = "none",legend.text = element_text(size = 10),
          legend.title = element_text(size = 12,face = "bold"),legend.background = element_blank(),
          plot.title = element_text(hjust = 0.5,size = 12,face = "bold"),
          # axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          strip.background = element_blank(),strip.text = element_text(size =12),
          panel.spacing = unit(0.5, "lines"))
  p
}



# mosaic plot 
#' @import ggplot2
#' @import pheatmap
#' @export mosaicPlot
mosaicPlot <- function(bins,binsInfo,chr,resolution = 500,list,parentNum = 24,color,
                       clust = T,methods = "ward.D2",
                       dist_method = "euclidean"){
  bins.tmp.chr <- bins[which(binsInfo$chr == chr),]
  bins.tmp <- bins.tmp.chr[,list]
  
  if(missing(color)){
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    color <- gg_color_hue(parentNum)
  }
  
  if (clust == T) {
    d <- dist(t(bins.tmp), method = dist_method) # distance matrix
    fit <- hclust(d, method=methods)
    orde_list <- fit$order
    bins.tmp <-bins.tmp[orde_list]
  }
  binsInfo.tmp <- binsInfo[which(binsInfo$chr == chr),]
  
  idx <- round((binsInfo.tmp$len1/sum(binsInfo.tmp$len1))*resolution)
  len <- length(which(round((binsInfo.tmp$len1/sum(binsInfo.tmp$len1))*resolution) < 1))
  # message(len,"|",length(idx),"|",len/length(idx))
  bins.tmp.mat <- matrix(NA,1,length(list))
  for (i in 1:length(idx)) {
    if (idx[i] == 0) {
      next
    }else{
      temp.mat <- matrix(NA,idx[i],length(list))
      for (j in 1:length(list)) {
        temp.mat[,j] <- bins.tmp[i,j] 
      }
      bins.tmp.mat <- rbind(bins.tmp.mat,temp.mat)
    }
  }
  
  bins.tmp.mat <- bins.tmp.mat[-1,]
  bins.tmp.mat <- t(bins.tmp.mat)

  pheatmap(bins.tmp.mat,cluster_cols = F,
           cluster_rows = F,show_rownames = F,show_colnames = F,
           treeheight_row = 0,treeheight_col =0,color = color,legend = F)
}
