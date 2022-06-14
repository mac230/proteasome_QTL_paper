## -----
##############
## USER INPUT:
##############

## the only things that should need to be changed for this script to work
## no trailing '/' at the end!
base_dir <- "~/emacs/ubi_QTL_paper"
setwd(base_dir)
needed_dirs <- c("/fcs", "/results", "/tables", "/scripts")
dir_maker <- function(x){
    ifelse(!dir.exists(paths = paste0("./", x)),
           dir.create(path = paste0("./", x)),
           paste0("dir ", paste0(getwd(), x), " exists_"))
}
sapply(X = needed_dirs, FUN = dir_maker)
work_dir       <- paste0(base_dir, "/fcs")
results_dir    <- paste0(base_dir, "/results")
tables_dir     <- paste0(base_dir, "/tables")


## -----
## <<Required_Packages>>
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
packages <- c("colorspace", "lattice", "ggvis", "dygraphs", "DescTools", "viridis")
sapply(X = packages, FUN = package_installer)
sapply(X = packages, FUN = require, character.only = TRUE)


## -----
## bioconductor
source("~/emacs/R/functions/load_flow_packages.R")


## -----
## required for merging flowsets into a single flowframe
source(file = "https://raw.githubusercontent.com/mac230/flow_scripts/master/set2frame.R")


## -----
## read in the data
dat <- read.flowSet(path = work_dir,
                    min.limit = 1,
                    alter.names = T)


##-----
## <<TFT_Transformation>>
## use the transform function to get the TFT/PSV parameters we want
## start by converting 0's in fluors to 1's via truncate transform
trunc.trans   <- truncateTransform("Convert 0's to 1's.", a = 1)
trunc.fluors  <- function(x){
    transform(x,
              `GFP.A` = trunc.trans(`GFP.A`),
              `mCherry.A` = trunc.trans(`mCherry.A`))}
dat <- fsApply(x = dat, FUN = trunc.fluors)

dat_list <- list()
for(i in 1:length(dat)) {
    dat_in <- as.data.frame(exprs(dat[[i]]))
    dat_in$log_GFP   <- log(dat_in$GFP.A, 10)
    dat_in$log_RFP   <- log(dat_in$mCherry.A, 10)
    dat_in$TFT_ratio <- -1 * log(dat_in$mCherry.A / dat_in$GFP.A, 2)
    dat_list[[i]] <- dat_in[dat_in$log_GFP > 2.5, ]
}


## -----
## now gate / subset and plot 
## a function to gate the cells to include only haploids.
## we identify these as a sharp peak in the lower end of
## the fsc density plot.  I take 10% above and below the
## max density value
dat_sub       <- list()
high_list     <- list()
low_list      <- list()
high_val      <- list()
low_val       <- list()

for (i in 1:length(dat)) {

    ## gate on FSC and get fluorescence-positive cells
    high_FSC <- quantile(x = dat_list[[i]]$FSC.A, probs = 0.70)
    low_FSC  <- quantile(x = dat_list[[i]]$FSC.A, probs = 0.30)
    dat_sub[[i]] <- dat_list[[i]][dat_list[[i]]$FSC.A > low_FSC &
                                  dat_list[[i]]$FSC.A < high_FSC &
                                  dat_list[[i]]$log_GFP > 3, ]

    ## values for overplotting extreme phenotype pools
    high_val[[i]]   <- quantile(x = dat_sub[[i]]$TFT_ratio, probs = 0.975)
    low_val[[i]]    <- quantile(x = dat_sub[[i]]$TFT_ratio, probs = 0.025)
    high_list[[i]] <- dat_sub[[i]][dat_sub[[i]]$TFT_ratio >= high_val[[i]], ]
    low_list[[i]]  <- dat_sub[[i]][dat_sub[[i]]$TFT_ratio <= low_val[[i]], ]
}

## lapply(X = dat_list,
##        FUN = nrow)
## lapply(X = dat_sub,
##        FUN = nrow)

## plot all 4 SFA samples and see which looks best
## 'high_TFT' = high UPS activity, so low TFT ratio
low_col <- "#ec6ab7ff"
high_col <-  "#5AC168ff"

low_t_col <- "#ec6ab722"
high_t_col <-  "#5AC16811"

## 61 / 62 = ODC
## 65 / 66 = rpn4
dat[[1]]@description$GUID.original
## 1 / 2 = ODC
## 3 / 4 = rpn4


## -----
## RFP ~ GFP plots
pdf(file = paste0(results_dir, "/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_ODC_sort_xy.pdf"))


## -----
## ODC TFT
print(
xyplot(log_RFP ~ log_GFP,
       data = dat_sub[[1]],
       type = c("p", "g"),
       xlim = c(2.4, 4.6),
       ylim = c(2.4, 4.6),       
       ## xlim = c(2.75, 4.7),
       ## ylim = c(2.8, 4.75),
       ## ylim = c(2.9, 4.85),
       pch = 19,
       col = gray(0.3, alpha = 0.2),
       cex = 0.25,
       xlab = expression("log"["10"]*" GFP"),
       ylab = expression("log"["10"]*" RFP"),
       key = list(lines = list(col = c(low_col,
                                       high_col),
                               lwd = 10,
                               size = 2.5),
                  text = list(labels = c("2% Low UPS Activity Gate",
                                         "2% High UPS Activity Gate"),
                              cex = 1.75)),
       scales = list(tck = c(1, 0),
                     alternating = F,
                     x = list(cex = 1.75,
                              at = seq(from = 2, to = 7, by = 0.5)),
                     y = list(cex = 1.75,
                              at = seq(from = 2, to = 7, by = 0.5))),
       par.settings = list(par.ylab.text = list(cex = 1.75),
                           par.xlab.text = list(cex = 1.75)),
       panel = function(...) {
           panel.xyplot(...)
           panel.points(x = high_list[[1]]$log_GFP,
                        y = high_list[[1]]$log_RFP,
                        pch = 19, col = high_col, cex = 0.25)
           panel.points(x = low_list[[1]]$log_GFP,
                        y = low_list[[1]]$log_RFP,
                        pch = 19, col = low_col, cex = 0.25)
}
))
## code end
dev.off()


## -----
## Rpn4 TFT
pdf(file = paste0(results_dir, "/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_Rpn4_sort_xy.pdf"))

print(
xyplot(log_RFP ~ log_GFP,
       data = dat_sub[[3]],
       type = c("p", "g"),
       xlim = c(2.75, 4.7),
       ylim = c(2.8, 4.75),
       pch = 19,
       col = gray(0.3, alpha = 0.2),
       cex = 0.25,
       xlab = expression("log"["10"]*" GFP"),
       ylab = expression("log"["10"]*" RFP"),
       key = list(lines = list(col = c(low_col,
                                       high_col),
                               lwd = 10,
                               size = 2.5),
                  text = list(labels = c("2% Low UPS Activity Gate",
                                         "2% High UPS Activity Gate"),
                              cex = 1.75)),
       scales = list(tck = c(1, 0),
                     alternating = F,
                     x = list(cex = 1.75,
                              at = seq(from = 2, to = 7, by = 0.5)),
                     y = list(cex = 1.75,
                              at = seq(from = 2, to = 7, by = 0.5))),
       par.settings = list(par.ylab.text = list(cex = 1.75),
                           par.xlab.text = list(cex = 1.75)),
       panel = function(...) {
           panel.xyplot(...)
           panel.points(x = high_list[[3]]$log_GFP,
                        y = high_list[[3]]$log_RFP,
                        pch = 19, col = high_col, cex = 0.25)
           panel.points(x = low_list[[3]]$log_GFP,
                        y = low_list[[3]]$log_RFP,
                        pch = 19, col = low_col, cex = 0.25)
}
))
## code end
dev.off()


## -----
## density plots


## -----
## ODC TFT
pdf(file = paste0(results_dir, "/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_ODC_sort_density.pdf"))

print(
densityplot(dat_sub[[1]]$TFT_ratio,
            plot.points = F,
            ## xlim = c(-2.25, 0.25),
            xlim = c(-2.5, 2.5),
            xlab = expression("UPS Activity (-log"["2"]*" RFP / GFP)"),
            key = list(lines = list(col = c(low_col,
                                            high_col),
                                    lty = 1,
                                    lwd = 4),
                       text = list(labels = c("2% Low UPS activity gate",
                                              "2% High UPS activity gate"),
                                   cex = 1.75)),
            scales = list(tck = c(1, 0),
                          x = list(cex = 1.75),
                          y = list(cex = 1.75)),

            col = gray(0.1),
            par.settings = list(par.xlab.text = list(cex = 1.75),
                                par.ylab.text = list(cex = 1.75)),
            lwd = 2,
            panel = function(...) {
                panel.densityplot(...)
                panel.abline(h = 0,
                             col = gray(0.9),
                             lwd = 2)
                panel.abline(v = c(high_val[[1]],
                                   low_val[[1]]),
                             col = c(high_col,
                                     low_col),
                             lwd = 2.5, lty = 1)

})
)
## code end 

dev.off()


## -----
## rpn4
pdf(file = paste0(results_dir, "/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_Rpn4_sort_density.pdf"))

print(
densityplot(dat_sub[[3]]$TFT_ratio,
            plot.points = F,
            xlim = c(-2.25, 0.25),
            xlab = expression("UPS Activity (-log"["2"]*" RFP / GFP)"),
            key = list(lines = list(col = c(low_col,
                                            high_col),
                                    lty = 1,
                                    lwd = 4),
                       text = list(labels = c("2% Low UPS activity gate",
                                              "2% High UPS activity gate"),
                                   cex = 1.75)),
            scales = list(tck = c(1, 0),
                          x = list(cex = 1.75),
                          y = list(cex = 1.75)),

            col = gray(0.1),
            par.settings = list(par.xlab.text = list(cex = 1.75),
                                par.ylab.text = list(cex = 1.75)),
            lwd = 2,
            panel = function(...) {
                panel.densityplot(...)
                panel.abline(h = 0,
                             col = gray(0.9),
                             lwd = 2)
                panel.abline(v = c(high_val[[3]],
                                   low_val[[3]]),
                             col = c(high_col,
                                     low_col),
                             lwd = 2.5, lty = 1)

})
)
## code end
dev.off()
