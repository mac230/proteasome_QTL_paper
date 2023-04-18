## check for Bioconductor and install if not available
ifelse(!requireNamespace("BiocManager", quietly = TRUE),
       install.packages("BiocManager",
                        dependencies = TRUE,
                        repos = "http://cran.wustl.edu/",
                        quiet = TRUE),
       paste0("Bioconductor available"))
require("BiocManager")

## -----
## load packages or install if not available
## have to split these out by bioconductor vs. non-bioconductor
## non-bioconductor
package_installer <- function(x){
    if(!requireNamespace(x, quietly = TRUE))
        install.packages(x, dependencies = TRUE,
                         repos = "http://cran.wustl.edu/",
                         quiet = TRUE, INSTALL_opts = '--no-lock')}
packages <- c("colorspace", "lattice",
              "ggvis", "dygraphs", "grid",
              "DescTools", "viridis",
              "latticeExtra", "session")
sapply(X = packages, FUN = package_installer)
sapply(X = packages, FUN = require, character.only = TRUE)


## -----
## bioconductor
bioc_package_installer <- function(x){if(!requireNamespace(x))
                                          BiocManager::install(x, INSTALL_opts = '--no-lock')}
bioc_packages <-  c("Gviz", "GenomicRanges", "rtracklayer", "org.Sc.sgd.db")
sapply(X = bioc_packages, FUN = bioc_package_installer)
sapply(X = bioc_packages, FUN = require, character.only = TRUE)


## -----
## load the QTL intervals
base_dir <- "~/emacs/proteasome_QTL_paper/scripts_data/chr_VII_QTL_range_plot/"
ranges <- read.csv(file = paste0(base_dir,
                                "chr_VII_QTL_peaks.csv"),
                   sep = ",",
                   header = T,
                   stringsAsFactors = F)
QTL_peak <- median(ranges$peak)

## calculate the confidence interval for the peak
err_peak <- qt(0.95, df = length(ranges$peak) - 1) *
    sd(ranges$peak) / sqrt(length(ranges$peak))

QTL_from    <- round((QTL_peak - err_peak - 3e3) / 1e5, digits = 2) * 1e5
QTL_to      <- round((QTL_peak + err_peak + 3e3) / 1e5, digits = 2) * 1e5
QTL_ci_down <- QTL_peak - err_peak
QTL_ci_up   <- QTL_peak + err_peak


## -----
## load variants and make a BY/RM region
snps <- read.table(paste0(base_dir,
                          "/gdata_42k_VEP_wExtraMarkersAttached.txt"),
                   sep = "\t",
                   header = FALSE,
                   quote = "",
                   stringsAsFactors = FALSE)

snps$chr <- unlist(strsplit(x = snps$V1,
                            split = ":.+"))

colnames(snps) <- c("snp", "pos", "alt", "class", "impact", "gene_name", "sys_name",
                    "type", "V9", "V101", "V11", "V12", "AA", "change", "V15", "V16",
                    "V17", "chr")

snps$variant <- unlist(lapply(X = strsplit(x = unlist(
                                               lapply(X = strsplit(x = snps[, 2],
                                                                   split = ":"),
                                                      FUN = function(x) {x[[2]]})),
                                           split = "-"),
                              FUN = function(x) {x[[1]]}))

snps$var <- as.numeric(snps$variant)

plot_snps <- snps[snps$chr == "chrVII", ]

## final frame of all SNPs in the region (111 SNPs)
plot_snps <- plot_snps[plot_snps$var > QTL_from - 500 & plot_snps$var < QTL_to + 500, ]
plot_snps$val <- rep(1, nrow(plot_snps))

snp_cols <- rep(gray(0.5), nrow(plot_snps))

snp_cols[plot_snps$impact == "LOW"]      <- gray(0.5)
snp_cols[plot_snps$impact == "MODERATE"] <- "orange"
snp_cols[plot_snps$impact == "HIGH"]     <- "red"

## -----
## now plot the interval
## availableDisplayPars("IdeogramTrack")
i_track <- IdeogramTrack(genome = "sacCer3",
                         chromosome = "chrVII")

## availableDisplayPars("GenomeAxisTrack")
a_track <- GenomeAxisTrack()

v_track <- AnnotationTrack(start = plot_snps$var,
                           end   = plot_snps$var + 1,
                           genome = "sacCer3",
                           chromosome = "chrVII",
                           strand = rep("*", nrow(plot_snps)),
                           from = QTL_from - 500,
                           to = QTL_to + 500)

f1_track <- AnnotationTrack(start = plot_snps$var[1],
                            end   = plot_snps$var[2],
                            genome = "sacCer3",
                            chromosome = "chrVII",
                            from = QTL_from,
                            to = QTL_to)

f2_track <- AnnotationTrack(start = plot_snps$var[1],
                            end   = plot_snps$var[2],
                            genome = "sacCer3",
                            chromosome = "chrVII",
                            from = QTL_from,
                            to = QTL_to)

g_track <- UcscTrack(track      = "SGD Gene",
                       genome     = "sacCer3",
                       chromosome = "chrVII",
                       from       = QTL_from,
                       to         = QTL_to,
                       trackType  = "GeneRegionTrack",
                       name       = "SGD Genes",
                       rstarts    = "exonStarts",
                       rends      = "exonEnds",
                       gene       = "name",
                       strand     = "strand",
                       symbol     = "name",
                       id         = "name",
                       transcript = "name")

## availableDisplayPars("GeneRegionTrack")
displayPars(a_track) <- list(cex = 2,
                             col = gray(0.4),
                             fontcolor = "black",
                             ## size of the tick marks
                             distFromAxis = 2.5,
                             labelPos = "above",
                             ## 3'/5' indicators
                             add35 = F,
                             lwd = 2.5,
                             scale = NULL,
                             frame = F,
                             size = 0.3,
                             stackHeight = 1,
                             ticksAt = seq(from = 4.04e5, to = 4.18e5, by = 2e3))

displayPars(g_track) <- list(stacking = "squish",
                             showTitle = F,
                             shape = "arrow",
                             transcriptAnnotation = "gene",
                             fill = "grey",
                             just.group = "above",
                             fontsize = 20,
                             fontcolor.group = "black",
                             ## arrowHeadMaxWidth = 50,
                             ## arrowHeadWidth = 25,
                             ## border color for all track items
                             col = "green",
                             stackHeight = 0,
                             col.border.title = "black",
                             col.title = "black",
                             ## size of the gene name font
                             cex.group = 1.75,
                             ## size = 0.2
                             lwd = 2,
                             frame = F,
                             fontcolor = "black",
                             stackHeight = 1
                             )

displayPars(i_track) <- list(cex = 2,
                             showTitle = F,
                             col.title = "black",
                             outline = F,
                             showAxis = F,
                             reverseStrand = F,
                             fontcolor = "black",
                             ## for the highlighted region
                             lty = 1,
                             lwd = 2,
                             v = 4,
                             stackHeight = 1)

displayPars(v_track) <- list(showTitle = F,
                             col.line = "red",
                             stacking = "dense",
                             ## make the plot area small
                             ## for the SNPs we're plotting
                             size = 0.1,
                             lwd = 1,
                             min.height = 20,
                             ## cut the whitespace above and below
                             stackHeight = 1)

displayPars(g_track) <- list(stacking = "squish",
                             shape = "arrow",
                             fontcolor.group = "black",
                             transcriptAnnotation = "gene",
                             just.group = "above",
                             arrowHeadMaxWidth = 50,
                             col = "black",
                             fill = gray(0.8),
                             showTitle = F,
                             cex.group = 2,
                             stackHeight = 0.5,
                             size = 0.2)

displayPars(g_track) <- list(cex.group = 1.5,
                             arrowHeadMaxWidth = 35,
                             stackHeight = 0.8,
                             lwd = 2)
## have to combine these or it double plots the genes
h_track <- HighlightTrack(trackList = list(g_track),
                           start = c(QTL_ci_down, QTL_peak - 20),
                           end = c(QTL_ci_up, QTL_peak + 20),
                           chromosome = "chrVII",
                           genome = "sacCer3")
## can combine here as well
displayPars(h_track) <- list(col = c("red", "red"),
                             fill = c("#FFECEE", "red"),
                             size = 0.1)
displayPars(i_track) <- list(size = 0.04)
displayPars(a_track) <- list(size = 0.12,
                             labelPos = "above")
displayPars(v_track) <- list(size = 0.07,
                             fill = gray(0.65),
                             col =  gray(0.65),
                             lwd = 0.5)
displayPars(h_track) <- list(size = 0.5)
displayPars(f1_track) <- list(size = 0.25,
                              showTitle = F)
displayPars(f2_track) <- list(size = 0.25,
                              showTitle = F)
pdf(file = paste0("~/emacs/ubi_QTL_paper/results/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_gviz_plot.pdf"),
    height = 7, width = 12)
plotTracks(list(i_track, a_track, v_track, h_track, f1_track, f2_track))
dev.off()
