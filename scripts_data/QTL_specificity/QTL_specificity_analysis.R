## -----
## start w/ all peaks, then subset per N-end rule branch

## -----
## setup and read in data
library("lattice")
base_dir    <- "~/emacs/proteasome_QTL_paper/scripts_data/QTL_specificity/"
peaks_table <- "merged_replicating_N-end_Rule_QTLs.csv"
peaks <- read.csv(file = paste0(base_dir, peaks_table),
                  header = T, stringsAsFactors = F)



## -----
## background
## consider peaks shared when:
## [1] they're on the same chromosome
## [2] they're w/in 100 kb of ea. other
## [3] they have the same direction of effect


## -----
## 'overlap_test'
## function to test whether two peaks overlap
overlap_test <- function(row_one, row_two, ...) {
    ifelse(test =
               ## QTL on same chromosome
               row_one$chr == row_two$chr &
               ## less than 100 kb apart
               abs(row_one$max_Index - row_two$max_Index) < 1e5 &
               ## same direction of effect 
               (row_one$delta_AF * row_two$delta_AF) > 0,
           yes = 1, 
           no = 0)
}

## make sure the function works as we intend:
## true positive
overlap_test(peaks[1, ], peaks[38, ])

## true negative
overlap_test(peaks[1, ], peaks[2, ])

## -----
## take only reporters w/ at least 7 QTLs:
out_list <- list()

i <- 1
for(i in 1:length(unique(peaks$reporter))) {

    curr_reporter <- unique(peaks$reporter)[i]
    curr_peaks    <- peaks[peaks$reporter == curr_reporter, ]

    ## skip to the next loop iteration if
    ## we don't have enough QTLs for a reporter
    if(nrow(curr_peaks) < 7) {
        next
    }

    out_list[[i]] <- curr_peaks

    }

out_peaks <- do.call("rbind", out_list)

length(unique(out_peaks$reporter))
## 14 reporters w/ at least 7 peaks


## -----
## compute overlap fractions:
out_sub_peaks <- out_peaks
overlap_results <- list()

i <- 1
j <- 1
for(i in 1:(length(unique(out_sub_peaks$reporter)))) {

    ref_reporter  <- unique(out_sub_peaks$reporter)[1]
    ref_peaks     <- out_sub_peaks[out_sub_peaks$reporter == ref_reporter, ]
    ref_n_peaks   <- nrow(ref_peaks)

    out <- data.frame(ref_reporter  = vector(),
                      ref_n_peaks   = vector(),
                      test_reporter = vector(),
                      test_n_peaks  = vector(),
                      total_peaks   = vector(),
                      overlaps      = vector(),
                      overlap_frac  = vector())

    for(j in 1:length(unique(out_sub_peaks$reporter))) {

        test_reporter   <- unique(out_sub_peaks$reporter)[j]
        test_peaks      <- out_sub_peaks[out_sub_peaks$reporter == test_reporter, ]
        test_n_peaks    <- nrow(test_peaks)
        total_peaks     <- ref_n_peaks + test_n_peaks
        overlap         <- sum(sapply(X = 1:nrow(ref_peaks),
                                 FUN = function(x) {
                                     sum(sapply(X = 1:nrow(test_peaks),
                                                FUN = function(y) {
                                                    overlap_test(ref_peaks[x, ], test_peaks[y, ])}))
                                 }))
        overlap_frac    <- overlap / (total_peaks - overlap)
        out[j, ] <- c(ref_reporter,
                      ref_n_peaks,
                      test_reporter,
                      test_n_peaks,
                      total_peaks,
                      overlap,
                      overlap_frac)
    }

    overlap_results[[i]] <- out
    out_sub_peaks <- out_sub_peaks[out_sub_peaks$reporter != ref_reporter, ]
 
}

overlap_results

overlap_frame  <- do.call("rbind", overlap_results)
overlap_subset <- overlap_frame$ref_reporter != overlap_frame$test_reporter
overlap_frame  <- overlap_frame[overlap_subset, ]

## now add the proteasome activity reporters:
overlap_frame[1 + nrow(overlap_frame), ] <- c("Rpn4_TFT", 7, "ODC_TFT", 11, 18, 3, (3 / (18 - 3)))

ord1 <- order(overlap_frame$overlap_frac,
              decreasing = F)

overlap_frame[ord1, ]

xyplot(as.numeric(overlap_frame$overlap_frac[ord1]) ~ 1:nrow(overlap_frame),
       type = c("p", "g"),
       xlim = c(-1, 93),
       ylim = c(-0.025, 0.825),
       ylab = "Fraction of Overlapping QTLs",
       xlab = "Reporter Pairs Ordered by Overlap Fraction",
       scales = list(x = list(at = c(-10, 101)),
                     y = list(at = seq(from = 0, to = 1, by = 0.1),
                              cex = 1.25),
                     tck = c(1, 0)),
       par.settings = list(par.ylab.text = list(cex = 1.25),
                           par.xlab.text = list(cex = 1.25)),
       panel = function(...) {
           panel.xyplot(...)
           panel.points(x = 51,
                        y = 0.2,
                        col = "Firebrick1",
                        pch = 19,
                        cex = 1.25)
})
## code end

## -----
## Arg/N-end rule only
library("lattice")
base_dir    <- "~/emacs/proteasome_QTL_paper/scripts/QTL_specificity/"
peaks_table <- "merged_replicating_N-end_Rule_QTLs.csv"
peaks <- read.csv(file = paste0(base_dir, peaks_table),
                  header = T, stringsAsFactors = F)

path_test <- function(x) {
    ifelse(x == "Ala" | x == "Cys" | x == "Gly" |
           x == "Met" | x == "Pro" | x == "Ser" |
           x == "Thr" | x == "Val",
           "Ac/N-end Pathway", "Arg/N-end Pathway")
}

path_subset <- path_test(gsub(pattern = "_TFT",
                              replacement = "",
                              x = peaks$reporter))

arg_subset <- path_subset == "Arg/N-end Pathway"

peaks <- peaks[arg_subset, ]

## -----
## take only reporters w/ at least 7 QTLs:
out_list <- list()

i <- 1
for(i in 1:length(unique(peaks$reporter))) {

    curr_reporter <- unique(peaks$reporter)[i]
    curr_peaks    <- peaks[peaks$reporter == curr_reporter, ]

    ## skip to the next loop iteration if
    ## we don't have enough QTLs for a reporter
    if(nrow(curr_peaks) < 7) {
        next
    }

    out_list[[i]] <- curr_peaks

    }


out_peaks <- do.call("rbind", out_list)

length(unique(peaks$reporter))
length(unique(out_peaks$reporter))
## 12 Arg/N-degrons; 7 w/ >= 7 QTLs 


## -----
## compute overlap fractions:
out_sub_peaks <- out_peaks
overlap_results <- list()

i <- 1
j <- 1
for(i in 1:(length(unique(out_sub_peaks$reporter)))) {

    ref_reporter  <- unique(out_sub_peaks$reporter)[1]
    ref_peaks     <- out_sub_peaks[out_sub_peaks$reporter == ref_reporter, ]
    ref_n_peaks   <- nrow(ref_peaks)

    out <- data.frame(ref_reporter  = vector(),
                      ref_n_peaks   = vector(),
                      test_reporter = vector(),
                      test_n_peaks  = vector(),
                      total_peaks   = vector(),
                      overlaps      = vector(),
                      overlap_frac  = vector())

    for(j in 1:length(unique(out_sub_peaks$reporter))) {

        test_reporter   <- unique(out_sub_peaks$reporter)[j]
        test_peaks      <- out_sub_peaks[out_sub_peaks$reporter == test_reporter, ]
        test_n_peaks    <- nrow(test_peaks)
        total_peaks     <- ref_n_peaks + test_n_peaks
        overlap         <- sum(sapply(X = 1:nrow(ref_peaks),
                                 FUN = function(x) {
                                     sum(sapply(X = 1:nrow(test_peaks),
                                                FUN = function(y) {
                                                    overlap_test(ref_peaks[x, ], test_peaks[y, ])}))
                                 }))
        overlap_frac    <- overlap / (total_peaks - overlap)
        out[j, ] <- c(ref_reporter,
                      ref_n_peaks,
                      test_reporter,
                      test_n_peaks,
                      total_peaks,
                      overlap,
                      overlap_frac)
    }

    overlap_results[[i]] <- out
    out_sub_peaks <- out_sub_peaks[out_sub_peaks$reporter != ref_reporter, ]
 
}

overlap_results

overlap_frame     <- do.call("rbind", overlap_results)
overlap_subset    <- overlap_frame$ref_reporter != overlap_frame$test_reporter
arg_overlap_frame <- overlap_frame[overlap_subset, ]
arg_overlap_frame$path <- "Arg/N-end"


## -----
## Ac/N-end rule only
library("lattice")
base_dir    <- "~/emacs/proteasome_QTL_paper/scripts/QTL_specificity/"
peaks_table <- "merged_replicating_N-end_Rule_QTLs.csv"
peaks <- read.csv(file = paste0(base_dir, peaks_table),
                  header = T, stringsAsFactors = F)

path_test <- function(x) {
    ifelse(x == "Ala" | x == "Cys" | x == "Gly" |
           x == "Met" | x == "Pro" | x == "Ser" |
           x == "Thr" | x == "Val",
           "Ac/N-end Pathway", "Arg/N-end Pathway")
}

path_subset <- path_test(gsub(pattern = "_TFT",
                              replacement = "",
                              x = peaks$reporter))

ac_subset <- path_subset == "Ac/N-end Pathway"

peaks <- peaks[ac_subset, ]

## -----
## take only reporters w/ at least 7 QTLs:
out_list <- list()

i <- 1
for(i in 1:length(unique(peaks$reporter))) {

    curr_reporter <- unique(peaks$reporter)[i]
    curr_peaks    <- peaks[peaks$reporter == curr_reporter, ]

    ## skip to the next loop iteration if
    ## we don't have enough QTLs for a reporter
    if(nrow(curr_peaks) < 7) {
        next
    }

    out_list[[i]] <- curr_peaks

    }


out_peaks <- do.call("rbind", out_list)

length(unique(peaks$reporter))
length(unique(out_peaks$reporter))
## 8 Ac/N-degrons; 7 w/ >= 7 QTLs 


## -----
## compute overlap fractions:
out_sub_peaks <- out_peaks
overlap_results <- list()

i <- 1
j <- 1
for(i in 1:(length(unique(out_sub_peaks$reporter)))) {

    ref_reporter  <- unique(out_sub_peaks$reporter)[1]
    ref_peaks     <- out_sub_peaks[out_sub_peaks$reporter == ref_reporter, ]
    ref_n_peaks   <- nrow(ref_peaks)

    out <- data.frame(ref_reporter  = vector(),
                      ref_n_peaks   = vector(),
                      test_reporter = vector(),
                      test_n_peaks  = vector(),
                      total_peaks   = vector(),
                      overlaps      = vector(),
                      overlap_frac  = vector())

    for(j in 1:length(unique(out_sub_peaks$reporter))) {

        test_reporter   <- unique(out_sub_peaks$reporter)[j]
        test_peaks      <- out_sub_peaks[out_sub_peaks$reporter == test_reporter, ]
        test_n_peaks    <- nrow(test_peaks)
        total_peaks     <- ref_n_peaks + test_n_peaks
        overlap         <- sum(sapply(X = 1:nrow(ref_peaks),
                                 FUN = function(x) {
                                     sum(sapply(X = 1:nrow(test_peaks),
                                                FUN = function(y) {
                                                    overlap_test(ref_peaks[x, ], test_peaks[y, ])}))
                                 }))
        overlap_frac    <- overlap / (total_peaks - overlap)
        out[j, ] <- c(ref_reporter,
                      ref_n_peaks,
                      test_reporter,
                      test_n_peaks,
                      total_peaks,
                      overlap,
                      overlap_frac)
    }

    overlap_results[[i]] <- out
    out_sub_peaks <- out_sub_peaks[out_sub_peaks$reporter != ref_reporter, ]
 
}

overlap_results

overlap_frame     <- do.call("rbind", overlap_results)
overlap_subset    <- overlap_frame$ref_reporter != overlap_frame$test_reporter
ac_overlap_frame  <- overlap_frame[overlap_subset, ]
ac_overlap_frame$path <- "Ac/N-end"

## -----
## now bind into a single dataframe
path_lists <- list()
path_lists[[1]] <- arg_overlap_frame
path_lists[[2]] <- ac_overlap_frame
path_frame <- do.call("rbind", path_lists)

path_frame[1 + nrow(path_frame), ] <- c("Rpn4_TFT", 7,
                                         "ODC_TFT", 11,
                                         18, 3, (3 / (18 - 3)),
                                         "Proteasome")


## -----
## where is 0.2 in the distribution?
mean(path_frame$overlap_frac < 0.2)
## [1] 0.2325581



## -----
## color setup
path_frame$col <- ifelse(test = path_frame$path == "Arg/N-end",
                         yes = "slategray1", 
                         no = ifelse(test = path_frame$path == "Ac/N-end",
                                     yes = "plum2", 
                                     no = "darkviolet")) 

ord1 <- order(path_frame$overlap_frac)
path_ord_frame <- path_frame[ord1, ]
path_ord_frame$n <- 1:nrow(path_ord_frame)

pdf(file = "~/emacs/ubi_QTL_paper/results/overlap_plot_by_branch.pdf")
xyplot(overlap_frac ~ n,
       groups = path,
       data = path_ord_frame,
       col = c("deeppink2", "navy", "firebrick1"),
       type = c("p"),
       xlim = c(-1, 44),
       ylim = c(-0.025, 0.825),
       cex = 1.5,
       key = list(corner = c(0.01, 0.99),
                  points = list(col = c("firebrick1", "navy", "deeppink2"),
                                pch = c(19, 1, 1),
                                cex = 2),
                  text = list(labels = c("Proteasome",
                                         "Arg/N-end Rule",
                                         "Ac/N-end Rule")),
                  between = 1,
                  background = gray(1, alpha = 0.5),
                  padding.text = 2.5),
       ylab = "Fraction of Overlapping QTLs",
       xlab = "",
       scales = list(x = list(at = c(-10, 101)),
                     y = list(at = seq(from = 0, to = 1, by = 0.1),
                              cex = 1.25),
                     tck = c(1, 0)),
       par.settings = list(par.ylab.text = list(cex = 1.25),
                           par.xlab.text = list(cex = 1.25),
                           clip = list(panel = F)),
       panel = function(...) {
           panel.abline(h = seq(from = 0, to = 1, by = 0.1),
                        v = seq(from = 0, to = 50, by = 5),
                        col = gray(0.9))
           panel.xyplot(...)
           panel.points(x = 11,
                        y = 0.2,
                        col = "Firebrick1",
                        pch = 19,
                        cex = 1.25)
           panel.text(x = 21.5,
                      y = -0.06,
                      fontface = "plain",
                      srt = 0,
                      cex = 1.25,
                      labels = "Reporter Pairs Ordered by Overlap Fraction")
       })
## code end 
dev.off()
