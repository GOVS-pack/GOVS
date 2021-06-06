## GOVS ########################################
#' @import lsmeans
#' @importFrom stats anova cor cor.test lm median na.omit var
#' @importFrom utils read.table write.table
#' @importFrom readr parse_number
#' @import grid
#' @export GOVS
GOVS <- function(hmp,ID = NULL,pheno,trait,bins,binsInfo,which = "max",output = NULL,module = "DES",designInfo = NULL,extractContent = "Genotype"){
  output <- output
  
  hmp <- hmp
  pheno <- pheno
  trait <- trait
  binsInfo <- binsInfo
  bins <- bins
  ID = ID
  which = which
  module = module
  extractContent = extractContent
  designInfo <- designInfo
  
  if (!is.null(output)) {
    switch(module,
           D = {genomeOptimization(pheno,bins,trait,output)},
           E = {extractGenome(hmp,binsInfo,ID = ID,designInfo = designInfo,output,bins = bins,extractContent)},
           S = {statDesign(designInfo,which,binsInfo,pheno,trait,output)},
           DE = {genomeOptimization(pheno,bins,trait,output)
             DI <- as.matrix(read.table(paste0(output,".merge"),header = T,sep = "\t", stringsAsFactors = F))
             extractGenome(hmp,binsInfo,ID = ID,bins = bins,designInfo = DI,output=output,extractContent = extractContent)},
           DS = {genomeOptimization(pheno,bins,trait,output)
             DI <- as.matrix(read.table(paste0(output,".merge"),header = T,sep = "\t", stringsAsFactors = F))
             statDesign(designInfo = DI,which,binsInfo,pheno,trait,output)},
           ES = {extractGenome(hmp,binsInfo,ID = ID,bins,designInfo,output,extractContent)
             statDesign(designInfo = designInfo,which,binsInfo,pheno,trait,output)},
           DES = {genomeOptimization(pheno,bins,trait,output)
             DI <- as.matrix(read.table(paste0(output,".merge"),header = T,sep = "\t", stringsAsFactors = F))
             extractGenome(hmp,binsInfo,ID = ID,bins,designInfo,output,extractContent)
             statDesign(designInfo = DI,which,binsInfo,pheno,trait,output)}
    )
  }else{
    switch(module,
           D = {D.res <- genomeOptimization(pheno,bins,trait,output)
                final_Res <- list(GORes = D.res)
                final_Res},
           E = {E.res <- extractGenome(hmp,binsInfo,ID = ID,designInfo = designInfo,output,bins = bins,extractContent)
                final_Res <- list(virtualGenome = E.res)
                final_Res},
           S = {S.res <- statDesign(designInfo = designInfo,which = which,binsInfo = binsInfo,pheno = phe,trait = trait,output)
                final_Res <- list(statRes = S.res)
                final_Res},
           DE = {D.res <- genomeOptimization(pheno,bins,trait,output)
                 E.res <- extractGenome(hmp,binsInfo,ID = ID,designInfo = D.res$overall,output,bins = bins,extractContent)
                 final_Res <- list(GORes = D.res,virtualGenome = E.res)
                 final_Res},
           DS = {D.res <- genomeOptimization(pheno,bins,trait,output)
                 S.res <- statDesign(designInfo = D.res$overall,which,binsInfo,pheno,trait,output)
                 final_Res <- list(GORes = D.res,statRes = S.res)
                 final_Res},
           ES = {E.res <- extractGenome(hmp,binsInfo,ID = ID,bins,designInfo = D.res,output,extractContent)
                 S.res <- statDesign(designInfo = designInfo,which,binsInfo,pheno,trait,output)
                 final_Res <- list(GORes = D.res,statRes = S.res)
                 final_Res},
           DES = {D.res <- genomeOptimization(pheno,bins,trait,output)
                  E.res <- extractGenome(hmp,binsInfo,ID = ID,designInfo = D.res$overall,output,bins = bins,extractContent)
                  S.res <- statDesign(designInfo = D.res$overall,which,binsInfo,pheno,trait,output)
                  final_Res <- list(GORes = D.res,virtualGenome = E.res,statRes = S.res)
                  final_Res}
    )
  }
}
