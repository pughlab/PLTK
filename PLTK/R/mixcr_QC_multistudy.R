

#' Plot Cap-TCRseq QC metrics across multiple studies
#' @description This function generates series of plots as QC metrics for Cap-TCRseq data from multiple studies.
#' Samples are ordered by total sequencing reads within each study. Numbers on the right are just plot numbers.
#' barplots with y-axis labeled as "% values" are values from the plot above normalized by aligned TCR reads. For example,
#' barplot below plot number 12 is: ("Clonotypes eliminated by PCR error correction"/"Successfully aligned reads")*100
#' See example_qc_ms.csv in raw-data
#' for headers in align_assemble_stats_df.
#'
#'
#' @param align_assemble_stats_df a merged file generated from align_stats and assemble_stats file.
#' @param plotname name for output plot
#' @param plotpath path to plot directory
#' @param sampletype type of sample eg: cfDNA, Blood or Tumor
#'
#' @examples

mixcrQC.multistudy.fx <- function(align_assemble_stats_df, plotname , plotpath, sampletype){

    # Function to align plots (from stackoverflow)
    align_plots1 <- function (...) {
      pl <- list(...)
      stopifnot(do.call(all, lapply(pl, inherits, "gg")))
      gl <- lapply(pl, ggplotGrob)
      bind2 <- function(x, y) gtable:::rbind_gtable(x, y, "first")
      combined <- Reduce(bind2, gl[-1], gl[[1]])
      wl <- lapply(gl, "[[", "widths")
      combined$widths <- do.call(grid::unit.pmax, wl)
      grid::grid.newpage()
      grid::grid.draw(combined)
    }

    # Theme for plots
    mytheme <- theme(axis.title.y = element_text(size = 25),
                     axis.title.x = element_blank(),
                     axis.line = element_line(color = "black"),
                     axis.text.y = element_text(size = 20),
                     axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "transparent",colour = NA),
            legend.key = element_rect(fill = "white", colour = "white"),
            plot.margin = unit(c(3,3,0,2),"cm"))

    #subset to sampletype
    align_assemble_stats_type <- align_assemble_stats_df[which(align_assemble_stats_df$sampletype == sampletype),]

    message("Total number of samples in the analysis:", nrow(align_assemble_stats_type))

    # order by study and total reads
    align_assemble_stats_type$samplename <- factor(align_assemble_stats_type$samplename,
                                                   levels = align_assemble_stats_type$samplename[order(align_assemble_stats_type$study,
                                                                                                       align_assemble_stats_type$Total.sequencing.reads)])

    studyfreq <- as.data.frame(table(align_assemble_stats_type$study), stringsAsFactors = F)
    message("Studies included in the analysis:")
    print(studyfreq)

    cumsum_numbers <- cumsum(studyfreq$Freq)

    # plots
    myplot_totalseq <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Total.sequencing.reads), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank()) +
      # This is to avoid printing min, max and median on y-axis
      geom_point(aes(y = Successfully.aligned.reads), color = "red",size = 4) +
      theme(axis.title.y.right = element_text(color = "red")) +
      geom_vline(xintercept = cumsum_numbers, color = "blue")

    myplot_totalseq <- myplot_totalseq + annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10",
                         sec.axis = dup_axis(name = "Aligned reads"))

    # barplot percentage aligned reads
    myplot_percaligned <- ggplot(aes(x = samplename,
                                     y = 100*(Successfully.aligned.reads/Total.sequencing.reads)), data = align_assemble_stats_type) +
      geom_bar(stat = "identity") +
      mytheme + theme(axis.text.x = element_blank())   +
      geom_vline(xintercept = cumsum_numbers, color = "blue")

    myplot_failedcdr3 <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Alignment.failed.because.of.absence.of.CDR3.parts), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "2", x = Inf,
                    y = max(Alignment.failed.because.of.absence.of.CDR3.parts)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    myplot_failedscore <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Alignment.failed.because.of.low.total.score), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank()) +
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "3", x = Inf,
                    y = max(Alignment.failed.because.of.low.total.score)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    myplot_failedhits <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Alignment.failed..no.hits..not.TCR.IG..), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "4", x = Inf,
                    y = max(Alignment.failed..no.hits..not.TCR.IG..)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    myplot_count <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Final.clonotype.count), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "5", x = Inf,
                    y = max(Final.clonotype.count)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    myplot_averagereads <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Average.number.of.reads.per.clonotype), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "6", x = Inf,
                    y = max(Average.number.of.reads.per.clonotype)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    myplot_readsinclons <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Reads.used.in.clonotypes..percent.of.total), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue")  +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "7", x = Inf,
                    y = max(Reads.used.in.clonotypes..percent.of.total)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    barplot_readsinclons <- ggplot(aes(x = samplename,
                                       y = 100*(Reads.used.in.clonotypes..percent.of.total/Successfully.aligned.reads)), data = align_assemble_stats_type) +
      geom_bar(stat = "identity") +
      mytheme + theme(axis.text.x = element_blank())   +
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10")

    myplot_lowreads <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Mapped.low.quality.reads..percent.of.used), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "8", x = Inf,
                    y = max(Mapped.low.quality.reads..percent.of.used)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    barplot_lowreads <- ggplot(aes(x = samplename,
                                   y = 100*(Mapped.low.quality.reads..percent.of.used/Successfully.aligned.reads)), data = align_assemble_stats_type) +
      geom_bar(stat = "identity") +
      mytheme + theme(axis.text.x = element_blank())   +
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10")

    myplot_pcrerror <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Reads.clustered.in.PCR.error.correction..percent.of.used), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "9", x = Inf,
                    y = max(Reads.clustered.in.PCR.error.correction..percent.of.used)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    barplot_pcrerror <- ggplot(aes(x = samplename,
                                   y = 100*(Reads.clustered.in.PCR.error.correction..percent.of.used/Successfully.aligned.reads)), data = align_assemble_stats_type) +
      geom_bar(stat = "identity") +
      mytheme + theme(axis.text.x = element_blank())   +
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10")

    myplot_readsdropped <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y =  Reads.dropped.due.to.the.lack.of.a.clone.sequence), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "10", x = Inf,
                    y = max(Reads.dropped.due.to.the.lack.of.a.clone.sequence)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    barplot_readsdropped <- ggplot(aes(x = samplename,
                                       y = 100*(Reads.dropped.due.to.the.lack.of.a.clone.sequence/Successfully.aligned.reads)), data = align_assemble_stats_type) +
      geom_bar(stat = "identity") +
      mytheme + theme(axis.text.x = element_blank())   +
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10")

    myplot_failedmapping <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Reads.dropped.due.to.failed.mapping), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank())+
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "11", x = Inf,
                    y = max(Reads.dropped.due.to.failed.mapping)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    barplot_failedmapping <- ggplot(aes(x = samplename,
                                        y = 100*(Reads.dropped.due.to.failed.mapping/Successfully.aligned.reads)), data = align_assemble_stats_type) +
      geom_bar(stat = "identity") +
      mytheme + theme(axis.text.x = element_blank())   +
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10")

    myplot_eliminatedPCRerror <- ggplot(data = align_assemble_stats_type, aes(x = samplename)) +
      geom_point(aes(y = Clonotypes.eliminated.by.PCR.error.correction), color = "black",size = 4) +
      mytheme + theme(axis.text.x = element_blank()) +
      geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10") +
      geom_text(aes(label = "12", x = Inf,
                    y = max(Clonotypes.eliminated.by.PCR.error.correction)),
                hjust = 0, size = 30) + coord_cartesian(clip = 'off')

    barplot_eliminatedPCRerror <- ggplot(aes(x = samplename,
                                             y = 100*(Clonotypes.eliminated.by.PCR.error.correction/Successfully.aligned.reads)), data = align_assemble_stats_type) +
      geom_bar(stat = "identity") +
      mytheme + geom_vline(xintercept = cumsum_numbers, color = "blue") +
      annotation_logticks(sides = "l") +
      scale_y_continuous(trans = "log10")

    pdf(file = paste0(plotpath, plotname, ".pdf"),
        width = 0.1 * nrow(align_assemble_stats_type), #play with this to get the right size of your plot
        height = 3 * 19,#total number of plots
        useDingbats = FALSE,
        onefile = FALSE)

    align_plots1(myplot_totalseq + ylab("Total\nsequencing reads"),
                 myplot_percaligned + ylab("Percentage\naligned reads"),
                 myplot_failedcdr3 + ylab("Alignment failed:\nno CDR3 parts"),
                 myplot_failedscore + ylab("Alignment failed:\nlow total score"),
                 myplot_failedhits + ylab("Alignment failed:\nno TCR hits"),
                 myplot_count + ylab("Final \nclonotype count"),
                 myplot_averagereads + ylab("Average number of\nreads per clonotype"),
                 myplot_readsinclons + ylab("Reads used\nin clonotypes"),
                 barplot_readsinclons + ylab("% Aligned reads"),
                 myplot_lowreads + ylab("Mapped\nlow quality reads"),
                 barplot_lowreads + ylab("% Aligned reads"),
                 myplot_pcrerror + ylab("Reads clustered in\nPCR error\n correction"),
                 barplot_pcrerror + ylab("% Aligned reads"),
                 myplot_readsdropped + ylab("Reads dropped:\nlack of \nclone sequence"),
                 barplot_readsdropped + ylab("% Aligned reads"),
                 myplot_failedmapping + ylab("Reads dropped:\nfailed mapping"),
                 barplot_failedmapping + ylab("% Aligned reads"),
                 myplot_eliminatedPCRerror + ylab("Clonotypes eliminated by\nPCR error \ncorrection"),
                 barplot_eliminatedPCRerror + ylab("% Aligned reads"))
    dev.off()
  }

