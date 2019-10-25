require(ggalluvial)
library(randomcoloR)


#' Tracking recurring immune cell clonotypes
#' @description Hello pugh lab!
#' @param datapath path to mixcr clones files 
#' @param plotpath path to plot directory
#' @param chain could be any of: TRA, TRB, TRD, TRG, IGH, IGL, IGK
#' @param filelist a list of files to track clonotypes. Output from list.files()
#' @param countfrac plots either clonal fraction of absolut counts. 
#' cloneCount or cloneFraction
#' @param clnefrc specify a cut-off from 0 to 1 to track and plot only a subset of clonotypes.
#' Useful when you have too many clonotypes to plot.
#' Clonal fraction of 0.001 is usually a good starting point.



clontrack.fx <- function(datapath, plotpath, chain, filelist, countfrac, clnefrc){
  
  if (!(countfrac %in% c("cloneFraction", "cloneCount"))) {
    stop("Error: unknown argument ", countfrac, ". Please provide either cloneFraction or cloneCount.")
  }  
  if (!(chain %in% c("TRA", "TRB", "TRD", "TRG"))) {
    stop("Error: unknown argument ", chain, ". Please provide one of the following: TRA, TRB, TRD, TRG.")
  }   
  
  message("list of files to track clones: ")
  print(filelist)
  
  #Compile a big file with patient's mixcr files loaded in
  i <- 1
  for (f in filelist){
    mixcrfle <- read.table(paste(datapath, f, sep = ""), 
                           header = TRUE, sep = "\t",
                           stringsAsFactors = FALSE,
                           na.strings = c("", "NA"))
    if(i == 1){
      compldfle <- mixcrfle[!duplicated(mixcrfle$aaSeqCDR3),]
      compldfle <- cbind(cloneno = row.names(compldfle), 
                         filename = f, 
                         compldfle)
      i <- i + 1   
    }
    else{
      compldfle1 <- mixcrfle[!duplicated(mixcrfle$aaSeqCDR3),]
      compldfle1 <- cbind(cloneno = row.names(compldfle1), 
                          filename = f, 
                          compldfle1)
      compldfle <- rbind(compldfle, compldfle1)
      rm(compldfle1)
    }
  }
  
  #Clean the sample name. It should be in this format: CHP_XXX-0X
  compldfle$samplename <- gsub(paste(".*",chain, sep = ""), "", compldfle$filename)
  compldfle$samplename <- gsub("-PBMC-DNA_2000000.txt", "", compldfle$samplename)  
  # Subset
  CDR3_fraction <- compldfle[, c("samplename","aaSeqCDR3","cloneFraction", "cloneCount")]
  
  # Subset df
  CDR3_fraction <- compldfle[, c("samplename","aaSeqCDR3","cloneFraction", "cloneCount")]
  # Subset to include only clonotypes with more than specified clonal fraction    
  CDR3_fraction <- CDR3_fraction[CDR3_fraction$cloneFraction > clnefrc,] 
  ## append the empty clonotypes after here.   
  
  
  # Number of samples
  mysamples <- unique(CDR3_fraction$samplename)
  
  #Assign colors to recurring clonotypes
  recurring <- unique(CDR3_fraction$aaSeqCDR3[duplicated(CDR3_fraction$aaSeqCDR3)])
  notrecurring <- CDR3_fraction$aaSeqCDR3[!CDR3_fraction$aaSeqCDR3 %in% recurring]
  
  message("Total number of recurring clonotypes: ")     
  print(length(recurring))
  
  if(length(recurring) > 50){
    recurring_df <- CDR3_fraction[CDR3_fraction$aaSeqCDR3 %in% recurring,]
    recurringcdr3_ordered <- unique(recurring_df$aaSeqCDR3[order(recurring_df$cloneCount, decreasing = TRUE)])
    message("Total number of recurring clonotypes > 50 ")   
    message("Tracking top 10 recurring clonotypes ")  
    myColors <- distinctColorPalette(10)
    
    myColors <- c(myColors, rep("white",length(recurring)-10),
                  rep("white",length(notrecurring)))
    names(myColors) <- c(recurringcdr3_ordered, notrecurring)
    
    message("these are what we color: ")  
    print(myColors[myColors != "white"])         
  }
  else{
    myColors <- distinctColorPalette(length(recurring))
    myColors <- c(myColors, rep("white",length(notrecurring)))
    names(myColors) <- c(recurring, notrecurring)
    
    message("these are what we color: ")  
    print(myColors[myColors != "white"]) 
  }
  
  
  # Generate a row for each sample that doesnot have recurring clonotype
  ## This ensures alluvia are colored
  
  for(c in recurring){
    tmp <- CDR3_fraction[CDR3_fraction$aaSeqCDR3 == c,]
    nonexsiting <- mysamples[!mysamples %in% tmp$samplename]
    if(length(nonexsiting) > 0){
      newentries <- data.frame("samplename" = nonexsiting, "aaSeqCDR3" = c, 
                               "cloneFraction" = 0, "cloneCount" = 0)
      CDR3_fraction <- rbind(CDR3_fraction, newentries)
    }
  }
  
  p <-  ggplot(CDR3_fraction, aes(x = samplename, 
                                  y = eval(as.name(countfrac)),
                                  fill = aaSeqCDR3,
                                  stratum = aaSeqCDR3,
                                  alluvium = aaSeqCDR3,
                                  label = aaSeqCDR3))
  
  myp <- p + geom_alluvium(decreasing = FALSE) + 
    geom_stratum(decreasing = FALSE, stat = "alluvium") + 
    scale_fill_manual(breaks = names(myColors[myColors != "white"]),
                      values = myColors) +
    theme(axis.title.y = element_text(size = 50),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.position = "bottom",
          plot.margin = unit(c(0.2,0,0,0),"cm")) + 
    labs(y = countfrac) 
  
  pdf(paste(plotpath, "clonetracking_", 
            chain, countfrac, ".pdf", sep = ""),
      width = 15, 
      height = 20,
      useDingbats = FALSE,
      onefile = FALSE)       
  print(myp)  
  dev.off()      
  
}