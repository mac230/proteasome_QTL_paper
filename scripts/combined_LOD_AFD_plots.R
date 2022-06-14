data_dir <- "~/data/illumina/2021.10.26_FPFA005_UBI_replicate_analyses/peaks/merged_peak_tables"

dat <- list()

for(i in 1:length(dir(data_dir))) {
    dat[[i]] <- read.table(dir(data_dir, full.names = T)[i],
                           sep = ",", header = T,
                           stringsAsFactors = F)
}

dat_f <- do.call("rbind", dat)

dat_avg <- data.frame(reporter = dat_f$reporter,
                      chr = dat_f$chr,
                      type = dat_f$type,
                      LOD = (dat_f$rep_1_LOD + dat_f$rep_2_LOD) / 2,
                      delta_AF = (dat_f$rep_1_delta_AF + dat_f$rep_2_delta_AF) / 2,
                      left_index = (dat_f$rep_1_left_Index + dat_f$rep_2_left_Index) / 2,
                      max_index = (dat_f$rep_1_max_Index + dat_f$rep_2_max_Index) / 2,
                      right_index = (dat_f$rep_1_right_Index + dat_f$rep_2_right_Index) / 2)

## common annotations, functions, etc ----------------
## check for Bioconductor and install if not available
ifelse(!requireNamespace("BiocManager", quietly = TRUE),
       install.packages("BiocManager",
                        dependencies = TRUE,
                        repos = "http://cran.wustl.edu/",
                        quiet = TRUE),
       paste0("Bioconductor available"))
require("BiocManager")

bioc_package_installer <- function(x) {
    if (!requireNamespace(x))
        BiocManager::install(x, INSTALL_opts = "--no-lock")
}

bioc_package_installer("VariantAnnotation")

library("VariantAnnotation")
source("~/QTL_scripts/gTest.R")
source("~/QTL_scripts/x_qtl_seq_functions_170831.R")
source("~/QTL_scripts/mp_JB_170901.R")
source("~/QTL_scripts/peaksFromVector.R")


## -----
## getGcoords function
getGcoords <- function (chr,     ## which chromosome
                        pos,     ## which position
                        spacing, ## how much spacing betw. chromosomes (for plotting)
                        sgd.table = "~/QTL_scripts/sacCer3ChromLengths.txt") {

    ## read in the lengths of the s. cer chromosomes,
    ## produce a cumulative sum of each chromosome's
    ## length plus that of the preceding chromosomes
    ## plus the spacing arg (typically 1e5)
    ## e.g., start pos. of chrII = 230218 (length chrI) + 1e5 (spacing)
    ## note that 'spacing' arg gets added to each chromosome,
    ## so that, e.g., end of chrIII is (230218 + 1243402) + 2e5, not 1e5
    offind <- as.vector(cumsum(read.table(sgd.table,
                                          header=FALSE,
                                          sep = "\t")[,2] + spacing))

    ## drop the last length from 'offind',which
    ## comes from adding the length of the
    ## mitochondrial chromosome.  uses negative
    ## subsetting via the 'length' arg
    offind <- offind[-length(offind)]

    ## add 0, since this is where we start.
    ## importantly, the cumulative sum for
    ## chr I becomes 0 following this, so
    ## that, e.g., pos 500 in chrI = 500
    offind <- c(0, offind)

    ## name the chromosomes, which is really,
    ## where they end plus the spacing
    names(offind) <- as.character(read.table(sgd.table,
                                             header = FALSE,
                                             sep = "\t")[,1])

    ## return the cumulative sum of positions for
    ## the chromosome supplied via the 'chr' arg
    ## this is wrapped in an 'sapply' so that you
    ## can supply mulitple positions simultaneously
    chr_off <- as.numeric(sapply(chr, function(x) {
                                     offind[x]
                                 }))

    ## now, add the position to the cumulative sum
    ## for the chromosome your position falls in,
    ## so:
    ## chr I position  500 = 500
    ## chr II position 500 = 330718:
    ## (230218 [chrI length] + 500 [pos.] + 1e5 [spacing]) => 330718
    return(chr_off + pos)
}

## padding for the plots
spacing       <- 1e5
trimFromEnd   <- 15e3
obsMin        <- 10
LoessSpan     <- 0.1

## AF and multipool line threshold, determined by null sorts
af_thresh       <- 0.07606
multi_thresh    <- 4.5

## read in chr lengths
## need this for building the heatmap
chr_lengths <- read.table("~/QTL_scripts/sacCer3ChromLengths.txt",
                          header = F)

## don't need mitochondrial chromosome
chr_lengths <- chr_lengths[1:16, ]

## vertical lines that separate chromosomes
## have to manually add the length of chrXVI
## through the second step below:
chr_dividers <- getGcoords(chr = 1:16,
                           pos = rep(1, 16),
                           spacing = spacing)

## 2nd step
chr_dividers <- c(chr_dividers,
                  getGcoords(chr = 16,
                             pos = chr_lengths$V2[16],
                             spacing = spacing))

## where to place the chromosome labels
chr_labels <- sapply(1:(length(chr_dividers) - 1),
                     function(i) {
                         (chr_dividers[i] + chr_dividers[i + 1]) / 2
                         })

## text labels for plot
chr_text <- as.roman(1:16)

## convert peaks from base to gcoords
dat_avg$coords <- getGcoords(chr = dat_avg$chr,
                             pos = dat_avg$max_index,
                             spacing = spacing)

## assign colors based on reporter
cols <- c("#D07900", "#B01013")
dat_avg$col <- ifelse(dat_avg$reporter == "ODC_TFT",
                      cols[1],
                      cols[2])

## delta_AF plot
delta_AF_plot <-  xyplot(dat_avg$delta_AF[1] ~ dat_avg$coords[1],
                         xlim = c(min(chr_dividers) - spacing,
                                  max(chr_dividers) + spacing),
                         xlab = "Chromosome",
                         ylab = "(High - Low UPS Activity Pools)",
                         ylim = c(-0.26, 0.26),
                         scales = list(x = list(at = chr_labels,
                                                labels = chr_text,
                                                cex = 1.75),
                                       y = list(at = c(-0.25, -0.12, 0, 0.12, 0.25),
                                                labels = c("-0.25", "-0.12", "0.0",
                                                           "0.12", "0.25"),
                                                cex = 1.75),
                                       tck = c(1, 0)),
                         par.settings = list(par.xlab.text = list(cex = 1.75),
                                             par.ylab.text = list(cex = 1.75),
                                             layout.widths = list(left.padding = 5),
                                             clip = list(panel = F)),
                         key = list(corner = c(0.01, 0.98),
                                    lines = list(col = cols,
                                                 lineheight = 2,
                                                 lwd = 4,
                                                 size = 3),
                                    padding.text = 3,
                                    text = list(labels = c("ODC TFT",
                                                           "Rpn4 TFT"),
                                                cex = 1.75),
                                    background = gray(1, alpha = 0.5),
                                    border = gray(0)),
                         type = "1",
                         panel = function(...) {
                             panel.xyplot(...)
                             panel.abline(v = chr_dividers,
                                          lty = 1,
                                          col = gray(0.9),
                                          lwd = 1.5)
                             for(r in 1:nrow(dat_avg)) {
                                 panel.lines(x = rep(dat_avg[r, ]$coords, 2),
                                             y = c(0, dat_avg[r, ]$delta_AF),
                                             col = dat_avg[r, ]$col,
                                             lty = 1,
                                             lwd = 4)
                                 panel.abline(h = 0, col = gray(0.4), lty = 1, lwd = 4)
                                 ## empirical 99th percentile
                                 ## panel.abline(h = c(af_thresh, -af_thresh),
                                 ##              lty = 2,
                                 ##              col = "red")
                             }
                         }
                         )

pdf(file = paste0("~/emacs/ubi_QTL_paper/results/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_delta_AF_line_plot.pdf"),
                  width = 14, height = 7)
print(delta_AF_plot)
grid.text(label = expression(paste(Delta, "RM Allele Frequency")),
          x = 0.0125,
          y = 0.52,
          rot = 90,
          default.units = "npc",
          gp = gpar(cex = 1.75))
dev.off()

LOD_plot <-  xyplot(dat_avg$LOD[1] ~ dat_avg$coords[1],
                    xlim = c(min(chr_dividers) - spacing,
                             max(chr_dividers) + spacing),
                    ylim = c(-1, 41),
                    xlab = "Chromosome",
                    ylab = "LOD",
                    scales = list(x = list(at = chr_labels,
                                           labels = chr_text,
                                           cex = 1.75),
                                  y = list(at = seq(from = 0, to = 40, by = 10),
                                           cex = 1.75),
                                  tck = c(1, 0)),
                    par.settings = list(par.xlab.text = list(cex = 1.75),
                                        par.ylab.text = list(cex = 1.75),
                                        layout.widths = list(left.padding = 5),
                                        clip = list(panel = F)),
                    key = list(corner = c(0.01, 0.98),
                               lines = list(col = cols,
                                            lineheight = 2,
                                            lwd = 4,
                                            size = 3),
                               padding.text = 3,
                               text = list(labels = c("ODC TFT",
                                                      "Rpn4 TFT"),
                                           cex = 1.75),
                               background = gray(1, alpha = 0.5),
                               border = gray(0)),
                    type = "1",
                    panel = function(...) {
                        panel.xyplot(...)
                        panel.abline(v = chr_dividers,
                                     lty = 1,
                                     col = gray(0.9),
                                     lwd = 1.5)
                        for(r in 1:nrow(dat_avg)) {
                            panel.lines(x = rep(dat_avg[r, ]$coords, 2),
                                        y = c(0, dat_avg[r, ]$LOD),
                                        col = dat_avg[r, ]$col,
                                        lty = 1,
                                        lwd = 4)
                            panel.abline(h = 0, col = gray(0.4), lty = 1, lwd = 4)
                            ## empirical 99th percentile
                            ## panel.abline(h = c(af_thresh, -af_thresh),
                            ##              lty = 2,
                            ##              col = "red")
                        }
                    }
                    )

pdf(file = paste0("~/emacs/ubi_QTL_paper/results/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_LOD_line_plot.pdf"),
                  width = 14, height = 7)
print(LOD_plot)
dev.off()
