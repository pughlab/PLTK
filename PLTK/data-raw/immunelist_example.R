
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


celllines_list <- immunelistfx("~/git/PLTK/PLTK/data-raw/", "TRB")
usethis::use_data(celllines_list)
