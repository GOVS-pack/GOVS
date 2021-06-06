## function 3 : statistic for design res and pheno rank
#============================== statistic genome parameters
# designInfo: data frame, function1 design results
# binsInfo: data frame, bins index, sta, end, length info 
# pheno: data frame, pheno file including sample seqname and aim trait
# trait: character, aim trait name
# output: output file prefix
#' @export statDesign
statDesign <- function(designInfo,which = "max",binsInfo,pheno,trait,output = NULL){
  
  bin.info <- binsInfo
  trim.info <- designInfo[,which]
  pheno <- pheno
  file.name <- output
  
  sta.res <- matrix(NA,nrow(pheno),7)
  colnames(sta.res)[1:7] <- c("Lines","Bins(#)","Bins(%)","Fragments(%)","phenotype","phenotypeRank","Cumulative(%)")
  sta.res[,1] <- pheno[,1]
  sta.res[,5] <- pheno[,trait]
  sta.res[,6] <- (nrow(pheno) + 1 - rank(pheno[,trait],na.last = F,ties.method = "first"))

  pb = txtProgressBar(min = 0, max = nrow(pheno), initial = 0,style=3) 
  for(i in 1:nrow(pheno)){
    if(sta.res[i,1] %in% trim.info){
      idx <- which(trim.info == sta.res[i,1])
      bin_prob <- length(idx)/nrow(bin.info) * 100
      prob <- sum(bin.info[idx,5])/(sum(bin.info[,ncol(bin.info)])) * 100
      sta.res[i,2:4] <- c(length(idx),bin_prob,prob)
    }else{
      sta.res[i,2:4] <- 0 
    }
    setTxtProgressBar(pb,i)
  }
  sta.res[,1] <- pheno[,1]
  
  sta.res <- sta.res[order(as.numeric(sta.res[,4]),decreasing = T),]
  
  sum_prob <- c()
  for (i in 1:nrow(sta.res)) {
    sum_prob <- c(sum_prob,sum(as.numeric(sta.res[1:i,4])))
  }
  
  sta.res[,7] <- sum_prob
  sta.res <- sta.res
  
  sta.res <- as.data.frame(sta.res)
  sta.res[,2:7] <- sapply(sta.res[,2:7],as.numeric)
  
  message("\nStatistics is done!",'\n')
  message(Sys.time(),'\n')
  
  if(!is.null(output)){
    write.table(sta.res,file = paste0(file.name,"_","sta_res.csv"),col.names = T,row.names = F,quote = F,sep = ",")
  }else{
    sta.res
  }
}
