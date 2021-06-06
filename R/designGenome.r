## function 1 : design genome 
#============================== design genome function parameters
# pheno: data frame, first column is smaple seqname
# bins: matrix, bins matrix original
# trait: character, trait name
# output: output file prefix
#' @import lsmeans
#' @export genomeOptimization
genomeOptimization <- function(pheno,bins,trait,output = NULL){
  bins <- bins
  pheno <- pheno
  trait <- trait
  output <- output
  
  rownames(pheno) <- as.character(pheno[,1])
  bins <- t(bins[,(colnames(bins) %in% rownames(pheno))])
  ParsF <- pheno[match(rownames(bins),rownames(pheno)),]
  message("Data has been loaded\nStart running ...\n")
  message(Sys.time(),'\n')

  final_bin_max <- final_bin_med <- final_bin_min <- data.frame()

  pb = txtProgressBar(min = 0, max = dim(bins)[2], initial = 0,style=3) 
  for(i in 1:dim(bins)[2]){
    
    pheno_bin <- data.frame(line=rownames(ParsF), bin=bins[,i], pheno=ParsF[,trait])
    pheno_bin$bin <- as.factor(pheno_bin$bin)
    pheno_lm <- lm(pheno~bin, data=pheno_bin)
    pheno_lsm <- summary(lsmeans(pheno_lm, "bin"))
    bin_max <- pheno_lsm[which(pheno_lsm$lsmean==max(pheno_lsm$lsmean, na.rm=T)),]
    bin_min <- pheno_lsm[which(pheno_lsm$lsmean==min(pheno_lsm$lsmean, na.rm=T)),]
    if(median(pheno_lsm$lsmean, na.rm=T) %in% pheno_lsm$lsmean){
      bin_med <- pheno_lsm[which(pheno_lsm$lsmean==median(pheno_lsm$lsmean, na.rm=T)),]
    }else{
      pheno_lsm_sort <- pheno_lsm[order(pheno_lsm$lsmean, decreasing=TRUE), ]
      bin_med <- pheno_lsm_sort[dim(pheno_lsm)[1]/2,]    
    }
    
    #### ?1 bin with only one line; ?2 max/min/median value with greater than two lines 
    pheno_bin_max <- pheno_bin[pheno_bin$bin %in% bin_max$bin,]
    line_bin_max <- pheno_bin_max[which(pheno_bin_max$pheno==max(pheno_bin_max$pheno, na.rm=T)),]
    line_bin_max$lsmean <- pheno_lsm[match(line_bin_max$bin, pheno_lsm$bin), "lsmean"]
    line_bin_max$pvalue <- anova(pheno_lm)[[5]][1]
    no_bin_max <- as.data.frame(table(pheno_bin_max$bin))
    no_bin_max <- no_bin_max[no_bin_max$Var1 %in% line_bin_max$bin, ]
    line_bin_max$Freq <- no_bin_max[match(line_bin_max$bin, no_bin_max$Var1), "Freq"]
    line_bin_max$binID <- rep(colnames(bins)[i], dim(line_bin_max)[1])
    line_bin_max$bin <- as.character(line_bin_max$bin)
    final_bin_max <- rbind(final_bin_max,line_bin_max[1,])
    
    pheno_bin_min <- pheno_bin[pheno_bin$bin %in% bin_min$bin,]
    line_bin_min <- pheno_bin_min[which(pheno_bin_min$pheno==min(pheno_bin_min$pheno, na.rm=T)),]
    line_bin_min$lsmean <- pheno_lsm[match(line_bin_min$bin, pheno_lsm$bin), "lsmean"]
    line_bin_min$pvalue <- anova(pheno_lm)[[5]][1]
    no_bin_min <- as.data.frame(table(pheno_bin_min$bin))
    no_bin_min <- no_bin_min[no_bin_min$Var1 %in% line_bin_min$bin, ]
    line_bin_min$Freq <- no_bin_min[match(line_bin_min$bin, no_bin_min$Var1), "Freq"]
    line_bin_min$binID <- rep(colnames(bins)[i], dim(line_bin_min)[1])
    line_bin_min$bin <- as.character(line_bin_min$bin)
    final_bin_min <- rbind(final_bin_min,line_bin_min[1,])
    
    pheno_bin_med <- pheno_bin[pheno_bin$bin %in% bin_med$bin,]
    if(median(pheno_bin_med$pheno, na.rm=T) %in% pheno_bin_med$pheno){
      line_bin_med <- pheno_bin_med[which(pheno_bin_med$pheno==median(pheno_bin_med$pheno, na.rm=T)),]
    }else{
      pheno_bin_med_sort <- pheno_bin_med[order(pheno_bin_med$pheno, decreasing=TRUE), ]
      line_bin_med <- pheno_bin_med_sort[dim(pheno_bin_med)[1]/2,]    
    }  
    line_bin_med$lsmean <- pheno_lsm[match(line_bin_med$bin, pheno_lsm$bin), "lsmean"]
    line_bin_med$pvalue <- anova(pheno_lm)[[5]][1]
    no_bin_med <- as.data.frame(table(pheno_bin_med$bin))
    no_bin_med <- no_bin_med[no_bin_med$Var1 %in% line_bin_med$bin, ]
    line_bin_med$Freq <- no_bin_med[match(line_bin_med$bin, no_bin_med$Var1), "Freq"]
    line_bin_med$binID <- rep(colnames(bins)[i], dim(line_bin_med)[1])
    line_bin_med$bin <- as.character(line_bin_med$bin)
    final_bin_med <- rbind(final_bin_med,line_bin_med[1,])
    
    setTxtProgressBar(pb,i)
    # print(i)
  }

  
  message('\nComputation is done! Merge res ...\n')
  message(Sys.time(),'\n')
  
  final_merge <- cbind(final_bin_max[,1],final_bin_med[,1],final_bin_min[,1])
  final_pvalue <- cbind(final_bin_max[,5],final_bin_med[,5],final_bin_min[,5])
  colnames(final_merge) <- colnames(final_pvalue) <- c("max","med","min")
  
  if (!is.null(output)) {
    write.table(final_bin_max, file = paste0(output,".max"), row.names=F, quote=F, sep="\t", col.names=F)
    write.table(final_bin_med, file = paste0(output,".med"), row.names=F, quote=F, sep="\t", col.names=F)
    write.table(final_bin_min, file = paste0(output,".min"), row.names=F, quote=F, sep="\t", col.names=F)
    
    write.table(final_merge, file = paste0(output,".merge"), row.names=F, quote=F, sep="\t", col.names=T)
    write.table(final_pvalue, file = paste0(output,".pvalue"), row.names=F, quote=F, sep="\t", col.names=T)
    message('All optimal source information have been written to',output,'\n')
  } else {
    final_list <- list(max = final_bin_max,med = final_bin_med,min = final_bin_min, overall = final_merge, pvalue = final_pvalue)
    final_list
  }
}
