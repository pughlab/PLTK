# This is the original clonetracking function

require(ggalluvial)
library(randomcoloR)
library(immunarch)


source("/Users/anabbi/OneDrive - UHN/R_src/ggplot2_theme.R")

datapath <- "/Users/anabbi/OneDrive - UHN/Documents/INTERCEPT/Data/"
plotpath <- "/Users/anabbi/OneDrive - UHN/Documents/INTERCEPT/Plots/"

clontrack.fx <- function(datapath, plotpath, chain, patient_id, countfrac){
    
  if (!(countfrac %in% c("cloneFraction", "cloneCount"))) {
    stop("Error: unknown argument ", countfrac, ". Please provide either cloneFraction or cloneCount.")

  }  
  if (!(chain %in% c("TRA", "TRB", "TRD", "TRG"))) {
    stop("Error: unknown argument ", chain, ". Please provide one of the following: TRA, TRB, TRD, TRG.")

  }   

  flelst <- list.files(datapath, recursive = TRUE,
                       pattern = paste("CLONES", chain, sep = "_"))

# subset to include only downsampled files
  ds_flelst <- flelst[grep("200000", flelst)]
# subset to patient_id: in CHP_XXX format
  ds_flelst_pt <- ds_flelst[grepl(patient_id, ds_flelst)]
    
    message("list of available files for patient: ", patient_id)
    print(ds_flelst_pt)

#Compile a big file with patient's mixcr files loaded in
  i <- 1
  for (f in ds_flelst_pt){
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
# Number of samples
    mysamples <- unique(CDR3_fraction$samplename)
    
#Assign colors to recurring clonotypes
    colorthem <- unique(CDR3_fraction$aaSeqCDR3[duplicated(CDR3_fraction$aaSeqCDR3)])
    norcolorthem <- CDR3_fraction$aaSeqCDR3[!CDR3_fraction$aaSeqCDR3 %in% colorthem]

    myColors <- distinctColorPalette(length(colorthem))
    myColors <- c(myColors, rep("white",length(norcolorthem)))
    names(myColors) <- c(colorthem, norcolorthem)

    message("these are what we color: ")  
    print(myColors[myColors != "white"])
    
# Generate a row for each sample that doesnot have recurring clonotype
## This ensures alluvia are colored
    
    for(c in colorthem){
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
    scale_fill_manual(values = myColors) +
    theme(axis.title.y = element_text(size = 50),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.position = "none",
        plot.margin = unit(c(0.2,0,0,0),"cm")) + 
    labs(y = countfrac) 

    
# return(myp)
    pdf(paste(plotpath, "clonetracking_", patient_id, 
              chain, countfrac, ".pdf", sep = ""),
        width = 15, 
        height = 20,
        useDingbats = FALSE,
        onefile = FALSE)       
    print(myp)  
    dev.off()      

}



library(immunarch)
Clonetracking.fx <- function(datapath, immunarchobj, mypatient, chain , reporder, plotpath){
    
#Subset to the patient of interest. Should be in this format: CHP_XXX
    mypt <- immunarchobj[grep(mypatient,names(immunarchobj))]

#Subset to the chain of interest. 
    mypt_chain <- mypt

# Reorder to reporder
    mynames <- paste(mypatient, reporder, sep = "-0")
    mypt_chain <- mypt_chain[mynames]
    target = mypt_chain[[1]] %>% select(CDR3.aa) #%>% head(10)
    tc <- trackClonotypes(mypt_chain, target,.norm = T)
    
#This is from immunarch vis function    
    melted = reshape2::melt(tc) %>% rename(Count = value, Sample = variable)
    last_obj_column_i = match("Sample", names(melted)) - 1
    melted[, Clonotype := do.call(paste, .SD), .SDcols = 1:last_obj_column_i]
#Which col other than the first (CDR3) and second (first ref rep) hass the minimum number of 0s? To order and catch the most recurrent clonotypes
    whichcol <- colnames(tc[,3:ncol(tc)])[colSums(tc[,3:ncol(tc)] == 0) == min(colSums(tc[,3:ncol(tc)] == 0))]
    tcorder <- tc[order(tc[[whichcol]], decreasing = TRUE),]
    
    myColors <- brewer.pal(length(tcorder$CDR3.aa[rowSums(tcorder == 0) < 3]),"Paired")
    #myColors <- myColors[-3]
    myColors <- c(myColors, rep("white",length(tcorder$CDR3.aa[rowSums(tcorder == 0) >= 3])))
    names(myColors) <- tcorder$CDR3.aa
    
    p <-  ggplot(melted, aes(x = Sample, y = Count,
                       fill = Clonotype,
                       stratum = Clonotype,
                       alluvium = Clonotype,
                       label = Clonotype))
    myp <- p + geom_alluvium(decreasing = FALSE) + geom_stratum(decreasing = FALSE, stat = "alluvium") + 
    scale_fill_manual(values = myColors) +
    theme(legend.position = "none") + mytheme
    
    return(myp)
#    pdf(paste("clonetracking_", mytitle, sep = ""),
#        width = 15, 
#        height = 20,
#        useDingbats = FALSE,
#        onefile = FALSE)       
#    myplot  
#    dev.off()    
}



#' This clonetrack function is for INSPIRE
#' Hello people!
#' @param datapath path to mixcr TRB clones
#' @param plotpath path to plot directory
#' @param patient_id in this format: INS-X-XXX
#' @param cycleorder order the plot. In this format. c("SB", "C3B", "C6B", etc)
#' @param countfrac either cloneCount or cloneFraction
#'
#' @return
#' @export
#'
#' @examples
clontrack.fx <- function(datapath, plotpath, patient_id, 
                         cycleorder , countfrac){
  
  if (!(countfrac %in% c("cloneFraction", "cloneCount"))) {
    stop("Error: unknown argument ", countfrac, ". Please provide either cloneFraction or cloneCount.")}
  
  
  flelst <- list.files(datapath, recursive = TRUE,
                       pattern = "CLONES_TRB")
  
  # subset to patient_id
  flelst_pt <- flelst[grepl(patient_id, flelst)]
  
  message("list of available files for patient: ", patient_id)
  print(flelst_pt)     
  
  #Compile a big file with patient's mixcr files loaded in
  i <- 1
  for (f in flelst_pt){
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
  
  #Clean up the sample name. 
  compldfle$samplename <- gsub(paste(".*","INSPIRE_", sep = ""), "", compldfle$filename)
  compldfle$samplename <- gsub("_TCR.*", "", compldfle$samplename)  
  
  # Subset df
  CDR3_fraction <- compldfle[, c("samplename","aaSeqCDR3","cloneFraction", "cloneCount")]
  # Subset to include only clonotypes with more than 0.001 clonal fraction    
  CDR3_fraction <- CDR3_fraction[CDR3_fraction$cloneFraction > 0.001,]
  
  # Number of samples
  mysamples <- unique(CDR3_fraction$samplename)
  # Reorder sample name with cycleorder
  mysamples <-  mysamples[sapply(cycleorder, function(x) { grep(x, mysamples) })]
  CDR3_fraction$samplename <- factor(CDR3_fraction$samplename, 
                                     levels = mysamples)
  
  #Find recurring clonotypes and order them based on clone count
  recurring <- unique(CDR3_fraction$aaSeqCDR3[duplicated(CDR3_fraction$aaSeqCDR3)])
  notrecurring <- CDR3_fraction$aaSeqCDR3[!CDR3_fraction$aaSeqCDR3 %in% recurring]
  
  recurring_df <- CDR3_fraction[CDR3_fraction$aaSeqCDR3 %in% recurring,]
  recurringcdr3_ordered <- unique(recurring_df$aaSeqCDR3[order(recurring_df$cloneCount, decreasing = TRUE)])
  
  message("Total number of recurring clonotypes: ")     
  print(length(recurring))
  
  #Assign colors to top 10 recurring clonotypes  
  myColors <- distinctColorPalette(10)
  myColors <- c(myColors, rep("white",length(recurring)-10),
                rep("white",length(notrecurring)))
  names(myColors) <- c(recurringcdr3_ordered, notrecurring)
  
  message("these are what we color: ")  
  print(myColors[myColors != "white"])
  
  p <-  ggplot(CDR3_fraction, aes(x = samplename, 
                                  y = eval(as.name(countfrac)),
                                  fill = aaSeqCDR3,
                                  stratum = aaSeqCDR3,
                                  alluvium = aaSeqCDR3,
                                  label = aaSeqCDR3))
  
  myp <- p + geom_alluvium(decreasing = FALSE) + 
    geom_stratum(decreasing = FALSE, stat = "alluvium") + 
    scale_fill_manual(values = myColors) +
    theme(axis.title.y = element_text(size = 50),
          axis.title.x = element_blank(),
          axis.line = element_line(color = "black"),
          axis.text = element_text(size = 50),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          legend.key = element_rect(fill = "white", colour = "white"),
          legend.position = "none",
          plot.margin = unit(c(0.2,0,0,0),"cm")) + 
    labs(y = countfrac) 
  
  pdf(paste(plotpath, "clonetracking_", patient_id, 
            countfrac, ".pdf", sep = ""),
      width = 15, 
      height = 20,
      useDingbats = FALSE,
      onefile = FALSE)       
  print(myp)  
  dev.off()       
}


