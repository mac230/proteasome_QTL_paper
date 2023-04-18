#############
## USER INPUT
#############
## set the specific directory you'll work
## in and name the comparison table
## TRAILING SLASH AT END OF DIR
## below, your project, e.g.,
## "2020.08.17_FPFA002_TDH3pr_Arg_N-end_TFT_sorts/"
base_dir       <- "~/emacs"
proj           <- "/proteasome_QTL_paper/scripts_data/QTL_plots/"
proj_dir       <- paste0(base_dir, proj)
data_dir       <- paste0(proj_dir, "rdata/")
comp_table     <- paste0(proj_dir,
                         "FPFA005_ubi_reporters_comparison_table.txt")
mpr            <- "_multipoolResults"
rd             <- ".RData"
pop            <- "_pop_0"
#################
## END USER INPUT
#################

needed.dirs <- c("results/", "rdata/", "/peaks")

dir.maker <- function(x){if(!dir.exists(paths = paste0(proj_dir, x)))
                             dir.create(path = paste0(proj_dir, x))}

sapply(X = needed.dirs, FUN = dir.maker)

results_dir     <- paste0(proj_dir, "results/")
output_dir      <- paste0(proj_dir, "r_output/")
peaks_dir       <- paste0(proj_dir, "peaks/")
experiment_file <- as.data.frame(read.table(comp_table,
                                            stringsAsFactors=FALSE,
                                            head=TRUE))
## stopping point
## SNPs is a giant table w/ SNP positions
SNPs <-  read.table(file = paste0(proj_dir,
                           "SNPs_170809_BY_positions.txt"),
                    stringsAsFactors = FALSE,
                    head = FALSE)

## as of 8/31/17, the SNPs seem not to be fully
## filtered and are out of sorting order
## I think the next section duplicates this code
## w/o the loop.  Frank's comment suggests he
## thinks the 'for' loop isn't working, but I
## think it is.  Does Maggie's code fix it?
## I guess not, given the dates.
for (thisChr in unique(SNPs[,1])){
    SNPs[SNPs[,1] == thisChr, 2] <- sort(SNPs[SNPs[,1] == thisChr, 2])
}

withMultipool <- TRUE

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

library("lattice")
library("VariantAnnotation")
source(paste0(proj_dir, "gTest.R"))
source(paste0(proj_dir, "x_qtl_seq_functions_170831.R"))
source(paste0(proj_dir, "mp_JB_170901.R"))
source(paste0(proj_dir, "peaksFromVector.R"))

## data frame with all yeast genes, plus
## chr, pos., strand, and names
geneInfo <- read.table(paste0(proj_dir,
                              "ensemblGenes_ensembl83_160307_MOD.txt"),
                       stringsAsFactors = FALSE,
                       sep = "\t",
                       header = TRUE)

## rownames become systemtatic names
rownames(geneInfo) <- geneInfo[, "geneID"]

## "geneName" is the common name, e.g., 'HOG1'
## for some (many?) rows of 'geneInfo', there
## is no 'geneName', so it's just an empty string
## e.g., head(allNames)
allNames <- geneInfo[, "geneName"]

names(allNames) <- geneInfo[, 1]
allNames[which(allNames == "")] <- names(allNames)[which(allNames == "")]
allNamesInv <- names(allNames)
names(allNamesInv) <- allNames

## padding for the plots
spacing       <- 1e5
trimFromEnd   <- 15e3
obsMin        <- 10
LoessSpan     <- 0.1

## AF line threshold, determined by null sorts
AFThres       <- 0.07606
## multipool LOD threshold, determined by null sorts
multiThres    <- 4.5
mpr           <- "_multipoolResults"
rd            <- ".RData"
p             <- 1

## read in chr lengths
## need this for building the heatmap
chr_lengths <- read.table(paste0(proj_dir,
                                 "sacCer3ChromLengths.txt"),
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


## -----
## ODC results 
## load 'multipoolResults.RData' for ODC replicate 1
## adds 'multiPeaks' and 'multipoolOutput',
## which we re-assign
experiment_file <- experiment_file[c(1, 7), ]

load(paste0(data_dir,
            experiment_file[1, 1], "_",
            experiment_file[1, 3], "_",
            experiment_file[1, 1], "_",
            experiment_file[1, 4],
            mpr, rd))

## load ".RData" for ODC replicate 1
## adds 'theseCounts', which we re-assign
load(paste0(data_dir,
            experiment_file[1, 1], "_",
            experiment_file[1, 3], "_",
            experiment_file[1, 1], "_",
            experiment_file[1, 4],
            rd))

odc_one_peaks    <- multiPeaks
odc_one_output   <- multipoolOutput
odc_one_counts   <- theseCounts

## load 'multipoolResults.RData' for ODC replicate 2
load(paste0(data_dir,
            experiment_file [1, 2], "_",
            experiment_file [1, 3], "_",
            experiment_file [1, 2], "_",
            experiment_file [1, 4],
            mpr, rd))

## load ".RData" for ODC replicate 2
load(paste0(data_dir,
            experiment_file [1, 2], "_",
            experiment_file [1, 3], "_",
            experiment_file [1, 2], "_",
            experiment_file [1, 4],
            rd))

odc_two_peaks    <- multiPeaks
odc_two_output   <- multipoolOutput
odc_two_counts   <- theseCounts


## -----
## rpn4 results
## load 'multipoolResults.RData' for rpn4 replicate 1
load(paste0(data_dir,
            experiment_file[2, 1], "_",
            experiment_file[2, 3], "_",
            experiment_file[2, 1], "_",
            experiment_file[2, 4],
            mpr, rd))

## load ".RData" for ODC replicate 1
## adds 'theseCounts', which we re-assign
load(paste0(data_dir,
            experiment_file[2, 1], "_",
            experiment_file[2, 3], "_",
            experiment_file[2, 1], "_",
            experiment_file[2, 4],
            rd))

rpn4_one_peaks    <- multiPeaks
rpn4_one_output   <- multipoolOutput
rpn4_one_counts   <- theseCounts

## load 'multipoolResults.RData' for ODC replicate 2
load(paste0(data_dir,
            experiment_file [2, 2], "_",
            experiment_file [2, 3], "_",
            experiment_file [2, 2], "_",
            experiment_file [2, 4],
            mpr, rd))

## load ".RData" for ODC replicate 2
load(paste0(data_dir,
            experiment_file [2, 2], "_",
            experiment_file [2, 3], "_",
            experiment_file [2, 2], "_",
            experiment_file [2, 4],
            rd))

rpn4_two_peaks    <- multiPeaks
rpn4_two_output   <- multipoolOutput
rpn4_two_counts   <- theseCounts    


## compute BY allele frequencies
## and sequencing coverage
allele_f_calc <- function(sample, pool, ...) {
    n  <- sample[, paste0(pool, "_ref")]
    d1 <- sample[, paste0(pool, "_ref")]
    d2 <- sample[, paste0(pool, "_alt")]
    n / (d1 + d2)
}

cover_calc <- function(sample, pool, ...) {
    a <- sample[, paste0(pool, "_ref")]
    r <- sample[, paste0(pool, "_alt")]
    a + r
}


## -----
## allele freq. calculations
## ODC rep. 1
odc_one_counts$h_BY_AF <- allele_f_calc(odc_one_counts, "high")
odc_one_counts$l_BY_AF <- allele_f_calc(odc_one_counts, "low")

## ODC rep. 2 
odc_two_counts$h_BY_AF <- allele_f_calc(odc_two_counts, "high")
odc_two_counts$l_BY_AF <- allele_f_calc(odc_two_counts, "low")

## rpn4 rep. 1
rpn4_one_counts$h_BY_AF <- allele_f_calc(rpn4_one_counts, "high")
rpn4_one_counts$l_BY_AF <- allele_f_calc(rpn4_one_counts, "low")

## rpn4 rep. 2
rpn4_two_counts$h_BY_AF <- allele_f_calc(rpn4_two_counts, "high")
rpn4_two_counts$l_BY_AF <- allele_f_calc(rpn4_two_counts, "low")


## -----
## sequencing coverage calculations
## ODC rep. 1
odc_one_counts$h_cover <- cover_calc(odc_one_counts, "high")
odc_one_counts$l_cover <- cover_calc(odc_one_counts, "low")

## ODC rep. 2
odc_two_counts$h_cover <- cover_calc(odc_two_counts, "high")
odc_two_counts$l_cover <- cover_calc(odc_two_counts, "low")

## rpn4 rep. 1
rpn4_one_counts$h_cover <- cover_calc(rpn4_one_counts, "high")
rpn4_one_counts$l_cover <- cover_calc(rpn4_one_counts, "low")

## rpn4 rep. 2
rpn4_two_counts$h_cover <- cover_calc(rpn4_two_counts, "high")
rpn4_two_counts$l_cover <- cover_calc(rpn4_two_counts, "low")


## -----
## coordinates for the plot 
gcoords <- getGcoords(odc_one_counts$chr,
                      odc_one_counts$pos,
                      spacing)



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


names(chr_labels)[16] = "chrXVI"




## make the ylimit of the plot the max of either:
## 1. the significance threshold
## 2. the highest peak of the first replicate
## 3. the highest peak of the second replicate
ylimMax <- max(c(multiThres,
                 sapply(c(odc_one_output,
                          odc_two_output,
                          rpn4_one_output,
                          rpn4_two_output),
                        function(x) {
                            max(x[[2]][, 2])}))) + 1


roll_high_odc_one <- rollLoessByChrWithWeights(data.frame(odc_one_counts[, "chr"],
                                                          odc_one_counts[, "h_BY_AF"],
                                                          gcoords,
                                                          median(odc_one_counts[, "h_cover"]),
                                                          stringsAsFactors = FALSE),
                                               LoessSpan)

roll_low_odc_one <- rollLoessByChrWithWeights(data.frame(odc_one_counts[, "chr"],
                                                         odc_one_counts[, "l_BY_AF"],
                                                         gcoords,
                                                         median(odc_one_counts[, "l_cover"]),
                                                         stringsAsFactors = FALSE),
                                              LoessSpan)

roll_high_odc_two <- rollLoessByChrWithWeights(data.frame(odc_two_counts[, "chr"],
                                                          odc_two_counts[, "h_BY_AF"],
                                                          gcoords,
                                                          median(odc_two_counts[, "h_cover"]),
                                                          stringsAsFactors = FALSE),
                                               LoessSpan)

roll_low_odc_two <- rollLoessByChrWithWeights(data.frame(odc_two_counts[, "chr"],
                                                         odc_two_counts[, "l_BY_AF"],
                                                         gcoords,
                                                         median(odc_two_counts[, "l_cover"]),
                                                         stringsAsFactors = FALSE),
                                              LoessSpan)

roll_high_rpn4_one <- rollLoessByChrWithWeights(data.frame(rpn4_one_counts[, "chr"],
                                                           rpn4_one_counts[, "h_BY_AF"],
                                                           gcoords,
                                                           median(rpn4_one_counts[, "h_cover"]),
                                                           stringsAsFactors = FALSE),
                                                LoessSpan)

roll_low_rpn4_one <- rollLoessByChrWithWeights(data.frame(rpn4_one_counts[, "chr"],
                                                          rpn4_one_counts[, "l_BY_AF"],
                                                          gcoords,
                                                          median(rpn4_one_counts[, "l_cover"]),
                                                          stringsAsFactors = FALSE),
                                               LoessSpan)

roll_high_rpn4_two <- rollLoessByChrWithWeights(data.frame(rpn4_two_counts[, "chr"],
                                                           rpn4_two_counts[, "h_BY_AF"],
                                                           gcoords,
                                                           median(rpn4_two_counts[, "h_cover"]),
                                                           stringsAsFactors = FALSE),
                                                LoessSpan)

roll_low_rpn4_two <- rollLoessByChrWithWeights(data.frame(rpn4_two_counts[, "chr"],
                                                          rpn4_two_counts[, "l_BY_AF"],
                                                          gcoords,
                                                          median(rpn4_two_counts[, "l_cover"]),
                                                          stringsAsFactors = FALSE),
                                               LoessSpan)


## code below returns the height on the y axis
## for the '*''s on the allele frequency difference
## plot.  it's called for ea. of the 2 replicates
## first, initialize 2 empty lists to put
## (roll_high - roll_low) values into
## y1_afd = replicate 1
## y2_afd = replicate 2
odc_one_afd  <- lapply(1:16, function(x) {
                           NULL
                       })
odc_two_afd  <- odc_one_afd
rpn4_one_afd <- odc_one_afd
rpn4_two_afd <- odc_one_afd

## chromosome names for pattern matching
chr_indices <- sapply(1:16, function(x) {
                          paste0("chr",
                                 as.character(as.roman(x)))
                      })

## name the elements according to chromosomes
names(odc_one_afd) <- chr_indices
names(odc_two_afd) <- chr_indices
names(rpn4_one_afd) <- chr_indices
names(rpn4_two_afd) <- chr_indices


## -----
## ODC rep 1
for (y in seq_along(odc_one_peaks)) {
    if (!is.null(odc_one_peaks[[y]]))
        for (q in seq_along(odc_one_peaks[[y]][, 1])) {
            ## has to return an index you can feed to roll_high_one
            h_roll_chr <- roll_high_odc_one[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                 x = names(roll_high_odc_one))]
            l_roll_chr <- roll_low_odc_one[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                x = names(roll_low_odc_one))]
            t_snp <- SNPs[SNPs$V1 == chr_indices[[y]], 2]
            t_min <- which.min(abs(odc_one_peaks[[y]][q, "maxIndex"] - t_snp))
            odc_one_afd[[y]][q] <- (h_roll_chr[t_min] - l_roll_chr[t_min])
            odc_one_afd[[y]][q] <- ifelse(odc_one_afd[[y]][q] > 0,
                                   ifelse(odc_one_afd[[y]][q] > AFThres,
                                          odc_one_afd[[y]][q] + 0.05,
                                          AFThres + 0.05),
                                   ifelse(odc_one_afd[[y]][q] < -AFThres,
                                          odc_one_afd[[y]][q] - 0.05,
                                          -AFThres - 0.05))
        }
}


## -----
## ODC rep 2
for (y in seq_along(odc_two_peaks)) {
    if (!is.null(odc_two_peaks[[y]]))
        for (q in seq_along(odc_two_peaks[[y]][, 1])) {
            ## has to return an index you can feed to roll_high_one
            h_roll_chr <- roll_high_odc_two[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                 x = names(roll_high_odc_two))]
            l_roll_chr <- roll_low_odc_two[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                x = names(roll_low_odc_two))]
            t_snp <- SNPs[SNPs$V1 == chr_indices[[y]], 2]
            t_min <- which.min(abs(odc_two_peaks[[y]][q, "maxIndex"] - t_snp))
            odc_two_afd[[y]][q] <- (h_roll_chr[t_min] - l_roll_chr[t_min])
            odc_two_afd[[y]][q] <- ifelse(odc_two_afd[[y]][q] > 0,
                                   ifelse(odc_two_afd[[y]][q] > AFThres,
                                          odc_two_afd[[y]][q] + 0.05,
                                          AFThres + 0.05),
                                   ifelse(odc_two_afd[[y]][q] < -AFThres,
                                          odc_two_afd[[y]][q] - 0.05,
                                          -AFThres - 0.05))
        }
}

## -----
## Rpn4 rep 1
for (y in seq_along(rpn4_one_peaks)) {
    if (!is.null(rpn4_one_peaks[[y]]))
        for (q in seq_along(rpn4_one_peaks[[y]][, 1])) {
            ## has to return an index you can feed to roll_high_one
            h_roll_chr <- roll_high_rpn4_one[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                  x = names(roll_high_rpn4_one))]
            l_roll_chr <- roll_low_rpn4_one[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                 x = names(roll_low_rpn4_one))]
            t_snp <- SNPs[SNPs$V1 == chr_indices[[y]], 2]
            t_min <- which.min(abs(rpn4_one_peaks[[y]][q, "maxIndex"] - t_snp))
            rpn4_one_afd[[y]][q] <- (h_roll_chr[t_min] - l_roll_chr[t_min])
            rpn4_one_afd[[y]][q] <- ifelse(rpn4_one_afd[[y]][q] > 0,
                                    ifelse(rpn4_one_afd[[y]][q] > AFThres,
                                           rpn4_one_afd[[y]][q] + 0.05,
                                           AFThres + 0.05),
                                    ifelse(rpn4_one_afd[[y]][q] < -AFThres,
                                           rpn4_one_afd[[y]][q] - 0.05,
                                           -AFThres - 0.05))
        }
}


## -----
## Rpn4 rep 2
for (y in seq_along(rpn4_two_peaks)) {
    if (!is.null(rpn4_two_peaks[[y]]))
        for (q in seq_along(rpn4_two_peaks[[y]][, 1])) {
            ## has to return an index you can feed to roll_high_one
            h_roll_chr <- roll_high_rpn4_two[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                  x = names(roll_high_rpn4_two))]
            l_roll_chr <- roll_low_rpn4_two[grep(pattern = paste0(chr_indices[[y]], "\\."),
                                                 x = names(roll_low_rpn4_two))]
            t_snp <- SNPs[SNPs$V1 == chr_indices[[y]], 2]
            t_min <- which.min(abs(rpn4_two_peaks[[y]][q, "maxIndex"] - t_snp))
            rpn4_two_afd[[y]][q] <- (h_roll_chr[t_min] - l_roll_chr[t_min])
            rpn4_two_afd[[y]][q] <- ifelse(rpn4_two_afd[[y]][q] > 0,
                                    ifelse(rpn4_two_afd[[y]][q] > AFThres,
                                           rpn4_two_afd[[y]][q] + 0.05,
                                           AFThres + 0.05),
                                    ifelse(rpn4_two_afd[[y]][q] < -AFThres,
                                           rpn4_two_afd[[y]][q] - 0.05,
                                           -AFThres - 0.05))
        }
}


## make pdfs
pdf(file = paste0(results_dir,
                  experiment_file[p, 5], "_",
                  experiment_file[p, 3], "_",
                  experiment_file[p, 4],
                  "_AFD_LOD_combined.pdf"),
    width=11,
    height=11)
par(mfrow = c(2, 1))

cols <- c("#EE8866", "#77AADD")
sig_line_col <- "#981F1FDD"   ## red


## -----
## ODC QTL data
odc_peak_pos <- c(69800, 418100,  85150, 291350,
                  20000, 409000, 666850, 768150,
                  47800, 410900, 441750)

odc_afd_peak_height <- c(0.099, -0.121, -0.101, -0.149,
                         -0.150,  0.225,  0.182,  0.111,
                         0.192,  0.127, -0.110)

odc_afd_peak_height <- ifelse(test = odc_afd_peak_height < 0,
                              yes = odc_afd_peak_height - 0.05,
                              no = odc_afd_peak_height + 0.025)

odc_lod_peak_height <- c(9.76, 7.13, 5.64, 12.83, 8.14,
                         28.74, 16.36, 8.13, 18.96, 7.96, 8.81)

odc_peak_chr <- c(2, 2, 4, 5, 7, 7, 10, 12, 13, 13, 14)

odc_peaks <- data.frame(pos = odc_peak_pos,
                        afd = odc_afd_peak_height,
                        lod = odc_lod_peak_height,
                        chr = odc_peak_chr)

## -----
## Rpn4 QTL data
rpn4_peak_pos <- c(240600, 259650,  88550, 882500,
                   672850, 544150, 167400)

rpn4_afd_peak_height <- c(-0.133, -0.127, -0.148, -0.111,
                          0.233, 0.145, -0.224)

rpn4_afd_peak_height <- ifelse(test = rpn4_afd_peak_height < 0,
                              yes = rpn4_afd_peak_height - 0.05,
                              no = rpn4_afd_peak_height + 0.025)

rpn4_lod_peak_height <- c(12.64, 10.09, 10.21, 6.80,
                          40.11, 16.58, 30.00)

rpn4_peak_chr <- c(4, 5, 7, 7, 12, 14, 15)

rpn4_peaks <- data.frame(pos = rpn4_peak_pos,
                         afd = rpn4_afd_peak_height,
                         lod = rpn4_lod_peak_height,
                         chr = rpn4_peak_chr)


## -----
## allele frequency difference plot
afd_plot <- xyplot((roll_high_odc_one - roll_low_odc_one) ~ gcoords,
                   type = c("l"),
                   xlab = "Chromosome",
                   ylab = "(High - Low Proteasome Activity Pool)",
                   xlim = c(min(chr_dividers) - spacing,
                            max(chr_dividers) + spacing),
                   ylim = c(-0.35, 0.35),
                   col = gray(1, alpha = 1),
                   lwd = 1.25,
                   key = list(corner = c(0.02, 0.99),
                              lines = list(col = c(cols, sig_line_col),
                                           lty = 1,
                                           lwd = 3,
                                           size = 2.5),
                              text = list(labels = c("ODC TFT", "Rpn4 TFT", "99.9% Quantile"),
                                          cex = 1.25),
                              between = 1.5,
                              padding.text = 2,
                              background = gray(1, alpha = 0.7)),
                   scales = list(tck = c(0.75, 0),
                                 x = list(at = chr_labels,
                                          labels = as.roman(1:16),
                                          cex = 1.25),
                                 y = list(cex = 1.25)),
                   par.settings = list(par.xlab.text = list(cex = 1.25),
                                       par.ylab.text = list(cex = 1.25),
                                       layout.widths = list(left.padding = 5),
                                       clip = list(panel = F)),
                   panel = function(...) {
                       ## 0 AFD line
                       panel.abline(h = 0,
                                    lty = 1,
                                    col = gray(0.7),
                                    lwd = 1.5)
                       ## allele frequency difference thresholds
                       panel.abline(h = c(AFThres, -AFThres),
                                    lty = 1,
                                    lwd = 0.8,
                                    col = sig_line_col)
                       ## chromosome cutoffs
                       panel.abline(v = c(chr_dividers, 13566086),
                                    lty = 1,
                                    lwd = 1,
                                    col = gray(0.875))
                       ## ODC rep one
                       panel.xyplot(...)
                       for (j in unique(SNPs[,1])){
                           panel.points(gcoords[odc_one_counts[, "chr"] == j],
                                        (roll_high_odc_one - roll_low_odc_one)[odc_one_counts[, "chr"] == j],
                                        type = "l",
                                        lwd = 1.5,
                                        col = cols[1]) ## dark grey
                       }
                       for (j in unique(SNPs[,1])){
                           panel.points(gcoords[odc_two_counts[, "chr"] == j],
                                        (roll_high_odc_two - roll_low_odc_two)[odc_two_counts[, "chr"] == j],
                                        type = "l",
                                        lwd = 1.5,
                                        col = cols[1]) ## dark grey
                       }
                       for (j in unique(SNPs[,1])){
                           panel.points(gcoords[rpn4_one_counts[, "chr"] == j],
                                        (roll_high_rpn4_one - roll_low_rpn4_one)[rpn4_one_counts[, "chr"] == j],
                                        type = "l",
                                        lwd = 1.5,
                                        col = cols[2]) ## dark grey
                       }
                       for (j in unique(SNPs[,1])){
                           panel.points(gcoords[rpn4_two_counts[, "chr"] == j],
                                        (roll_high_rpn4_two - roll_low_rpn4_two)[rpn4_two_counts[, "chr"] == j],
                                        type = "l",
                                        lwd = 1.5,
                                        col = cols[2]) ## dark grey
                       }
                       panel.text(labels = "RM Allele Frequency Difference",
                                  x = -1.6e6,
                                  y = 0,
                                  srt = 90,
                                  cex = 1.25)
                       panel.text(x = getGcoords(chr = odc_peaks$chr,
                                                 pos = odc_peaks$pos,
                                                 spacing = 1e5),
                                  y = odc_peaks$afd,
                                  labels = rep("*", times = nrow(odc_peaks)),
                                  col = cols[1],
                                  cex = 2.5)
                       panel.text(x = getGcoords(chr = rpn4_peaks$chr,
                                                 pos = rpn4_peaks$pos,
                                                 spacing = 1e5),
                                  y = rpn4_peaks$afd,
                                  labels = rep("*", times = nrow(rpn4_peaks)),
                                  col = cols[2],
                                  cex = 2.5)
                   })


## -----
## LOD plot
lod_plot <- xyplot(gcoords ~ gcoords,
                   type = c("l"),
                   xlab = "Chromosome",
                   ylab = "      LOD",
                   xlim = c(min(chr_dividers) - spacing,
                            max(chr_dividers) + spacing),
                   ylim = c(-2, 47.5),
                   col = gray(1, alpha = 0.9),
                   lwd = 1.25,
                   key = list(corner = c(0.02, 0.99),
                              lines = list(col = c(cols, sig_line_col),
                                           lty = 1,
                                           lwd = 3,
                                           size = 2.5),
                              text = list(labels = c("ODC TFT", "Rpn4 TFT", "Significance Threshold"),
                                          cex = 1.25),
                              between = 1.5,
                              padding.text = 2,
                              background = gray(1, alpha = 0.7)),
                   scales = list(tck = c(0.75, 0),
                                 x = list(at = chr_labels,
                                          labels = as.roman(1:16),
                                          cex = 1.25),
                                 y = list(cex = 1.25)),
                   par.settings = list(par.xlab.text = list(cex = 1.25),
                                       par.ylab.text = list(cex = 1.25),
                                       layout.widths = list(left.padding = 7),
                                       clip = list(panel = F)),
                   panel = function(...) {
                       ## 0 AFD line
                       panel.abline(h = 0,
                                    lty = 1,
                                    col = gray(0.7),
                                    lwd = 1.5)
                       ## allele frequency difference thresholds
                       panel.abline(h = multiThres,
                                    lty = 1,
                                    lwd = 0.8,
                                    col = sig_line_col)
                       ## chromosome cutoffs
                       panel.abline(v = c(chr_dividers, 13566086),
                                    lty = 1,
                                    lwd = 1,
                                    col = gray(0.875))
                       ## ODC rep one
                       panel.xyplot(...)
                       for (j in 1:16) {
                           panel.points(getGcoords(paste0("chr", as.roman(j)),
                                                   rpn4_one_output[[j]][[2]][, 1],
                                                   spacing),
                                        rpn4_one_output[[j]][[2]][, 2],
                                        type = "l",
                                        lwd = 1.5,
                                        col = cols[2])
                       }
                       for (j in 1:16) {
                           panel.points(getGcoords(paste0("chr", as.roman(j)),
                                                   rpn4_two_output[[j]][[2]][, 1],
                                                   spacing),
                                        rpn4_two_output[[j]][[2]][, 2],
                                        type = "l",
                                        lwd = 1.5,
                                        col = cols[2])
                       }
                       for (j in 1:16) {
                           panel.points(getGcoords(paste0("chr", as.roman(j)),
                                                   odc_one_output[[j]][[2]][, 1],
                                                   spacing),
                                        odc_one_output[[j]][[2]][, 2],
                                        type = "l",
                                        lwd = 1.5,
                                        col = cols[1])
                       }
                       for (j in 1:16) {
                           panel.points(getGcoords(paste0("chr", as.roman(j)),
                                                   odc_two_output[[j]][[2]][, 1],
                                                   spacing),
                                        odc_two_output[[j]][[2]][, 2],
                                        type = "l",
                                        lwd = 1.5,
                                        col = cols[1])
                       }
                       panel.text(x = getGcoords(chr = odc_peaks$chr,
                                                 pos = odc_peaks$pos,
                                                 spacing = 1e5),
                                  y = odc_peaks$lod + 2,
                                  labels = rep("*", times = nrow(odc_peaks)),
                                  col = cols[1],
                                  cex = 2.5)
                       panel.text(x = getGcoords(chr = rpn4_peaks$chr,
                                                 pos = rpn4_peaks$pos,
                                                 spacing = 1e5),
                                  y = rpn4_peaks$lod + 2,
                                  labels = rep("*", times = nrow(rpn4_peaks)),
                                  col = cols[2],
                                  cex = 2.5)
                   })


## -----
## make the combined pdf
pdf(file = "~/emacs/ubi_QTL_paper/figures/all_reporters_AFD_LOD_combined.pdf",
    height = 10, width = 12)
print(afd_plot,
      split = c(1, 1, 1, 2),
      newpage = T)
print(lod_plot,
      split = c(1, 2, 1, 2),
      newpage = F)
grid.text(label = c("A", "B"),
          x = c(0.017, 0.017),
          y = c(0.98, 0.475),
          default.units = "npc",
          gp = gpar(cex = 2.5,
                    fontface = "bold"))
dev.off()
