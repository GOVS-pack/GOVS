#' @import grid
#' @import readr
#' @export transHapmap2numeric
#' @export gcaCompute
#' @export scaCompute
#' @export reviseFunc
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
## trans 
#' @import pbapply
transHapmap2numeric <- function(G){
  rownamesG <- rownames(G)
  homo <- c("AA","TT","CC","GG")
  hete <- c("AC","AG","AT","CG","TC","TG","CA","GA","TA","GC","CT","GT")
  res <- pbapply(G,2,function(x){
    a <- table(x)
    patern <- names(a)
    lg <- length(a) 
    if(lg == 1){
      if(patern %in% homo){
        SNP <- rep(0,nrow(G))
      }else{
        SNP <- rep(1,nrow(G))
      }
    }else{
      SNP <- rep(2,nrow(G))
      SNP[which((x %in% hete) == TRUE)] <- 1
      b <- table(x[which((x %in% homo) == TRUE)])
      SNP[which(x == names(which(b == max(b)))[1])] <- 0
    }
    SNP
  })
  res
  rownames(res) <- rownamesG
  res
}
## cor func 
cor_func <- function(x,y,dight = 3,method = "pearson"){
  mat <- cbind(x,y)
  mat <- na.omit(mat)
  cor <- round(cor(mat[,1],mat[,2],method = method),digits = 3)
  cor.t <- formatC(cor.test(mat[,1],mat[,2],method = method)$p.value,format = "e",digits = 3)
  
  text <- switch (method,
                  pearson =  paste0("Pearson's product-moment correlation:\n",
                                    "R=",cor,"\n",
                                    "pvalue=",cor.t),
                  kendall =  paste0("Kendall's rank correlation tau:\n",
                                    "R=",cor,"\n",
                                    "pvalue=",cor.t),
                  spearman =  paste0("Spearman's rank correlation rho:\n",
                                     "R=",cor,"\n",
                                     "pvalue=",cor.t)
  )
  
  res <- list(R=cor,pvalue=cor.t,text = text)
  res
}

## gca compute function 
gcaCompute <- function(phe_df,which,trait){
  
  dimension <- length(unique(phe_df[,which]))
  gca <- matrix(NA,dimension,length(trait))
  rownames(gca) <- unique(phe_df[,which])
  colnames(gca) <- paste0("gca_",trait)
  for(i in 1:dimension){
    sub <- subset(phe_df,phe_df[,which] == unique(phe_df[,which])[i])
    for (j in trait) {
      eval(parse(text = paste0("mean_",j,"<- mean(as.numeric(sub[,'",j,"']),na.rm = T)")))
      eval(parse(text = paste0("allmean_",j,"<- mean(as.numeric(phe_df[,'",j,"']),na.rm = T)")))
      eval(parse(text = paste0("gca[",i,",'",paste0('gca_',j),"'] <- mean_",j,"-allmean_",j))) 
    }
  }
  gca[which(gca == "NaN")] <- NA
  class(gca) <- "numeric"
  gca
}

## sca compute func
scaCompute <- function(phe_df,which_male,which_female,trait,seqname){
  gca_male <- gcaCompute(phe_df = phe_df,which = which_male,trait = trait)
  gca_female <- gcaCompute(phe_df = phe_df,which = which_female,trait = trait)

  dimension <- dim(phe_df)[1]
  sca <- matrix(NA,dimension,length(trait))
  rownames(sca) <- seqname
  colnames(sca) <- paste0("sca_",trait)
  
  for (i in trait) {
    eval(parse(text = paste0("allmean_",i,"<- mean(as.numeric(phe_df[,'",i,"']),na.rm = T)")))
  }
 
  for (i in 1:dim(phe_df)[1]) {
    for (j in trait) {
      eval(parse(text = paste0("sca[",i,",'",paste0('sca_',j),"'] <- phe_df[",i,",'",j,"'] -allmean_",j,"- gca_female[phe_df[",i,",'",which_female,"'],'gca_",j,"'] - gca_male[phe_df[",i,",'",which_male,"'],'gca_",j,"']"))) 
    }
  }
  sca
}

## revese pred res function
reviseFunc <- function(ori,aim,cut = 10,sample_names){
  a <- cut(ori,breaks = cut)
  b <- cut(aim,breaks = cut)
  
  ori.mat <- data.frame(phe = ori,interval = a,rank = as.numeric(a))
  aim.mat <- data.frame(phe = aim,interval = b,rank = as.numeric(b),seqname = sample_names)
  ori.mean <- mean(ori.mat$phe,na.rm = T)
  
  res.all <- c()
  res.all.name <- c()
  for (i in 1:cut) {
    tmp.ori <- subset(ori.mat,ori.mat$rank == i)
    tmp.aim <- subset(aim.mat,aim.mat$rank == i)
    
    a.dim <- dim(tmp.ori)[1]
    b.dim <- dim(tmp.aim)[1]
    if(a.dim == 0){
      mean.tmp <- mean(parse_number(strsplit(names(table(ori.mat[,2]))[i],",")[[1]]))
      sd.tmp <- var(parse_number(strsplit(names(table(ori.mat[,2]))[i],",")[[1]]))
      res.tmp <- mean.tmp + tmp.aim$phe * sd.tmp
      res.all <- c(res.all,res.tmp)
      res.all.name <- c(res.all.name,as.character(tmp.aim$seqname))
    }else if(a.dim == 1){
      mean.tmp <- mean(c(tmp.ori$phe,parse_number(strsplit(names(table(ori.mat[,2]))[i],",")[[1]])))
      sd.tmp <- var(c(tmp.ori$phe,parse_number(strsplit(names(table(ori.mat[,2]))[i],",")[[1]])))
      # res.tmp <- if(tmp.ori$phe >= mean.tmp){mean.tmp - tmp.aim$phe * sd.tmp}else{mean.tmp + tmp.aim$phe * sd.tmp} 
      res.tmp <- mean.tmp + tmp.aim$phe * sd.tmp
      res.all <- c(res.all,res.tmp)
      res.all.name <- c(res.all.name,as.character(tmp.aim$seqname))
    }else{
      sd.tmp <- var(tmp.ori$phe)
      mean.tmp <- mean(tmp.ori$phe)
      res.tmp <- mean.tmp + tmp.aim$phe * sd.tmp 
      res.all <- c(res.all,res.tmp)
      res.all.name <- c(res.all.name,as.character(tmp.aim$seqname))
    }
  }
  res.all <- cbind(res.all.name,res.all)
  res.all <- as.numeric(res.all[match(sample_names,res.all[,1]),2])
  # print(cor(res.all,aim))
  res.all
}
