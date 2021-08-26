## function 2 : extract optimization genome according to function 1(design results) 
#============================== extract optimization genome parameters
# hmp: hapmap data for candidate population, SNP rs must coded with pattern "chr[0-9].s_[0-9]*" (eg: chr1.s_4831, "1" for chromosome, "4831" for locus). 
# ID: character array, if NULL, the hapmap including header, if not, provided ID
# designInfo: data frame,function1 design results
# output: output file prefix
#' @export extractGenome
extractGenome <- function(hmp,binsInfo,ID = NULL,designInfo,output = NULL,bins,extractContent = "Genotype"){
  
  if (extractContent == "BINsource") {
    message(paste0('Parameter "extractContent" is ', extractContent,". Bin source virtual genome will been extracted and assembled \n"))
    bins <- as.matrix(bins)
    design.info <- as.matrix(designInfo)
    trim.mat <- design.info
    res.G <- matrix(NA,1,ncol(trim.mat))
    names <- c()
    for (i in 1:nrow(trim.mat)) {
      oneBin <- design.info[i,]
      ## if all design info is NA,next
      if (sum(is.na(oneBin)) != 0) {
        next
      }else{
        res.G <- rbind(res.G,bins[i,oneBin])
        names <- c(names,rownames(bins)[i])
      }
    }
    
    res.G <- res.G[-1,]
    rownames(res.G) <- names
    colnames(res.G) <- colnames(trim.mat)
    if (!is.null(output)) {
      write.table(res.G,row.names = T,col.names = T,quote = F,file = paste0(output,".bin"),sep = "\t")
      message(paste0('Extract bin source file is written to PATH: "',output,'.bin"'))
    }else{
      message(paste0(dim(res.G)[1]," bins are used\nExtraction completed!"))
      res.G
    }
  }else{
    message(paste0('Parameter "extractContent" is ', extractContent,", virtual genome will been extracted and assembled \n"))
    G <- hmp 
    bin.info <- as.matrix(binsInfo)
    if (is.null(ID) == F) {
      G.header = ID
      message("Line IDs are redefined.")
      colnames(G)[12:ncol(G)] <- G.header
    }else{
      message("This message indicates that the genotypic matrix contains line IDs, if not, please check data!")
    }
    
    design.info <- designInfo
    G <- trimws(G)
    rownames(G) <- G[,1]
    trim.mat <- as.matrix(design.info)
    # message(paste0("Design info is: ",dim(trim.mat),"\n"))
    
    bin.info <- cbind(bin.info,trim.mat)
    bin.info <- na.omit(bin.info)
    dim(bin.info)
    
    check1 <- all(bin.info[,6:ncol(bin.info)] %in% colnames(G))
    if(check1 == T){
      message(paste0("Data input is: ",check1,"\n"))
    }else{
      stop(paste0("Data input is: ",check1,"\n","Please check your data!"))
    }
    
    ## SNP.info
    SNP.info <- G[,1]
    chr <- do.call(rbind,strsplit(SNP.info,"chr"))
    chr <- do.call(rbind,strsplit(chr,".s_"))
    SNP.info <- cbind(SNP.info,chr)
    colnames(SNP.info) <- c("ID","chr","pos")
    
    # bin.info <- cbind(bin.info,trim.mat)
    # bin.info <- na.omit(bin.info)
    # dim(bin.info)
    # 
    names.all <- c()
    res.G <- matrix(NA,1,ncol(trim.mat))
    
    pb = txtProgressBar(min = 0, max = nrow(binsInfo), initial = 0,style=3) 
    
    k <- 0
    for (i in as.numeric(names(table(binsInfo[,2])))) {
      snp.sub <- subset(SNP.info,SNP.info[,2] == i)
      bin.sub <- subset(bin.info,bin.info[,2] == i)
      G.sub <- subset(G,as.numeric(G[,3]) == i)
      
      idx <- which((as.numeric(snp.sub[,3]) >= as.numeric(bin.sub[1,3])) & (as.numeric(snp.sub[,3]) <= as.numeric(bin.sub[1,4])))
      res.G <- rbind(res.G,G.sub[idx,bin.sub[1,6:ncol(bin.info)]])
      names.all <- c(names.all,G.sub[idx,1])
      for(j in 2:nrow(bin.sub)){
        idx <- which((as.numeric(snp.sub[,3]) > as.numeric(bin.sub[j,3])) & (as.numeric(snp.sub[,3]) <= as.numeric(bin.sub[j,4])))
        res.tmp <- G.sub[idx,bin.sub[j,6:ncol(bin.info)]]
        res.G <- rbind(res.G,res.tmp)
        names.all <- c(names.all,G.sub[idx,1])
        
        # print(paste0(i,":",j))
        
        k <- k + 1
        setTxtProgressBar(pb,k)
      }
    }
    res.G <- res.G[-1,]
    # res.G <- t(res.G)  
    rownames(res.G) <- names.all
    colnames(res.G) <- colnames(trim.mat)
    if (!is.null(output)) {
      write.table(res.G,row.names = T,col.names = T,quote = F,file = paste0(output,".G"),sep = "\t")
      message(paste0('\nExtract genome file is written to PATH: "',output,'.G"'))
    }else{
      message("\n",paste0(dim(bin.info)[1],"bins are used\nExtraction completed!"))
      res.G
    }
  }
}
