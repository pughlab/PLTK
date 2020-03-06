
#' List clone counts for diversity measures
#' @description This funtion lists read counts (column "cloneCount" on mixcr output) 
#' for dowstream diversity calculations. Note that this function does not include files with only one clone.
#' This is just to remove the highly shallow samples which happens if using RNAseq data. This should not cause any issues 
#' when using the function with captured data.
#' 
#' @param datapath path to mixcr files
#' @param chain any of: TRA, TRB, TRD, TRG, IGH, IGL, IGK
#'
#' 
#' 
immunelistfx <- function(datapath, chain){
  file_list<- list.files(datapath, pattern = paste("CLONES", chain, 
                                                   sep = "_"))
  readlist = list()
  i <- 1
  for(f in file_list){
    mixcrfle <- read.table(paste(datapath, f, 
                                 sep = ""), 
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE,
                           na.strings = c("", "NA"))
    f <- substr(f, 11, nchar(f)-4)
    if(nrow(mixcrfle) <= 1){next()}
    readlist[[i]] <- mixcrfle$cloneCount
    names(readlist)[i] <- f
    i <- i + 1
  }
  return(readlist)
}
