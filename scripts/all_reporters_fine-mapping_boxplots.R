## -----
base_dir  <- "~/data/flow/"
individual_dirs <- sapply(X = c("2022.03.09_RPT6_fine-mapping_ODC_TFT_flow",
                                "2022.03.09_RPT6_fine-mapping_rpn4_TFT_flow",
                                "2022.03.10_RPT6_fine-mapping_pro_TFT_flow",
                                "2022.03.10_RPT6_fine-mapping_thr_flow",
                                "2022.03.11_RPT6_fine-mapping_ser_flow",
                                "2022.03.11_RPT6_fine-mapping_trp_flow"),
                                FUN = function(x) {
                                    paste0(base_dir, x)
                          })

## make sure we have the data for each reporter:
sapply(X = individual_dirs,
       FUN = function(x,...) {
           stopifnot(file.exists(paste0(x, "/dataframes/all_medians.R")))
})



## -----
## load packages
source("~/emacs/R/functions/load_flow_packages.R")
source("~/emacs/R/functions/color_setup_BY_RM.R")

## -----
## create lists for storing datasets

## 'all' lists contain all data points
all_genotype_list <- list()
all_strain_list   <- list()


## -----
## load in the data for each reporter
i <- 1
for(i in 1:length(individual_dirs)) {

    ## for each reporter, load the object 'all_medians',
    ## a dataset containing the median of independent
    ## biological replicates for 1e4 cells analyzed by flow
    load(paste0(individual_dirs[i], "/dataframes/all_medians.R"))

    
    
    ## use gated cells for analysis
    all_gated <- all_medians[all_medians$gating == "gated", ]
    all_gated <- all_gated[all_gated$strain != "BY_strain", ]
    ## 2 strains * 2 gRNAs * 24 bio. replicates = 96

    ## convert TFT ratios to z-scores and scale rel. to BY full median
    all_gated$genotype <- ifelse(test = grepl(pattern = "BY.*",
                                              x = all_gated$strain),
                                 yes = "BY",
                                 no = "RM")

    scale_val <- median(all_gated$TFT_loess[all_gated$genotype == "BY"])
    all_gated$TFT_scaled <- as.vector(scale(x = all_gated$TFT_loess,
                                            center = scale_val,
                                            scale = T))
    
    all_genotype_list[[i]] <- all_gated

}
 
## (6 reporters * 4 strains * 24 replicates) = 576 obs
all_genotype_data <- do.call("rbind", all_genotype_list)
all_genotype_data$genotype_factor <- as.factor(all_genotype_data$genotype)

## levels(all_genotype_data$reporter_factor)
col_by <- "#2166ACFF"
col_rm <- "#BF3232FF"

reporter_titles <- gsub(pattern = "ODC_TFT",
                        replacement = "ODC TFT",
                        x = all_genotype_data$reporter_factor)

reporter_titles <- gsub(pattern = "rpn4_TFT",
                        replacement = "Rpn4 TFT",
                        x = reporter_titles)

reporter_titles <- gsub(pattern = "pro_TFT",
                        replacement = "Pro TFT",
                        x = reporter_titles)

reporter_titles <- gsub(pattern = "thr_TFT",
                        replacement = "Thr TFT",
                        x = reporter_titles)

reporter_titles <- gsub(pattern = "ser_TFT",
                        replacement = "Ser TFT",
                        x = reporter_titles)

reporter_titles <- gsub(pattern = "trp_TFT",
                        replacement = "Trp TFT",
                        x = reporter_titles) 

all_genotype_data$reporter_titles <- as.factor(reporter_titles)

## 4 x 12 layout
## save.session(file = "~/Desktop/2022.06.05_RPT6_combined_boxplot.Rsession")

## use this w/ 'index.cond' below: 
## levels(all_genotype_data$reporter_titles)[c(1, 3, 2, 4, 5, 6)]

pdf(file = "~/emacs/ubi_QTL_paper/results/RPT6_all_reporters_boxplot.pdf",
    height = 4, width = 12)
print(
xyplot(TFT_scaled ~ genotype_factor | reporter_titles,
       groups = genotype_factor,
       index.cond = list(c(1, 3, 2, 4, 5, 6)),
       data = all_genotype_data,
       ylim = c(-2.5, 3.6),
       xlab = "",
       layout = c(6, 1),
       ylab = "",
       scales = list(tck = c(1, 0),
                     alternating = F,
                     x = list(labels = c("BY", "RM"),
                              cex = 1.45,
                              rot = 0),
                     y = list(cex = 1.45)),
       par.strip.text = list(cex = 1.35),
       par.settings = list(box.dot = list(pch = "|"),
                           box.umbrella = list(lty = 1,
                                               lwd = 1,
                                               col = gray(0.4)),
                           box.rectangle = list(lty = 1,
                                                lwd = 1,
                                                col = gray(0.4),
                                                fill = gray(0.97)),
                           strip.background = list(col = gray(0.9)),
                           par.ylab.text = list(cex = 1.45),
                           par.xlab.text = list(cex = 1.45),
                           layout.widths = list(left.padding = 9),
                           clip = list(panel = F)),
       panel = function(subscripts, ...) {
           ## 0 line
           panel.abline(h = 0, lty = 2, col = gray(0.7))
           
           ## boxplot
           panel.bwplot(...,
                        subscripts = subscripts,
                        ## pch = "|",
                        box.width = 0.6,
                        horizontal = F,
                        do.out = F)

           ## BY median lines
           panel.segments(x0 = 0.7,
                          x1 = 1.3,
                          y0 = tapply(X = all_genotype_data$TFT_scaled[all_genotype_data$genotype == "BY"],
                                      INDEX = all_genotype_data$reporter_titles[all_genotype_data$genotype == "BY"],
                                      FUN = median)[which.packet()],
                          y1 = tapply(X = all_genotype_data$TFT_scaled[all_genotype_data$genotype == "BY"],
                                      INDEX = all_genotype_data$reporter_titles[all_genotype_data$genotype == "BY"],
                                      FUN = median)[which.packet()],
                          lwd = 2.5,
                          col = "black")
           
           ## RM median lines
           panel.segments(x0 = 1.7,
                          x1 = 2.3,
                          y0 = tapply(X = all_genotype_data$TFT_scaled[all_genotype_data$genotype == "RM"],
                                      INDEX = all_genotype_data$reporter_titles[all_genotype_data$genotype == "RM"],
                                      FUN = median)[which.packet()],
                          y1 = tapply(X = all_genotype_data$TFT_scaled[all_genotype_data$genotype == "RM"],
                                      INDEX = all_genotype_data$reporter_titles[all_genotype_data$genotype == "RM"],
                                      FUN = median)[which.packet()],
                          lwd = 2.5,
                          col = "black")
           
           ## all data points
           panel.stripplot(...,
                           subscripts = subscripts,
                           pch = 21,
                           cex = 0.8,
                           fill = c(col_by, col_rm),
                           col = "black",
                           jitter.data = T,
                           amount = 0.14,
                           horizontal = F)
           ## pvals for BY_full vs. RM_full lmer:
           ## ODC = 2.84e-06
           ## rpn4 = 0.415
           ## pro = 0.0209
           ## thr = 4.99e-05
           ## ser = 0.001053
           ## trp = 0.06220
           panel.text(x = 1.5,
                      y = 3.2,
                      fontface = "plain",
                      srt = 0,
                      cex = 1.3,
                      labels = c(expression(italic(p)*" = 2.8e-6"),
                                 expression(italic(p)*" = 0.021"),
                                 expression(italic(p)*" = 0.42"),
                                 expression(italic(p)*" = 1.1e-3"),
                                 expression(italic(p)*" = 4.9e-5"),
                                 expression(italic(p)*" = 0.062"))[packet.number()])
           ## genotype indicator 
           panel.text(x = -0.45,
                      y = -3.05,
                      fontface = "plain",
                      srt = 0,
                      cex = 1.45,
                      labels = c(expression(italic("RPT6")*" -175"),
                                 rep("", 5))[packet.number()])
           ## y axis label
           panel.text(x = c(-0.8, rep(-10, 5)),
                      y = 1,
                      fontface = "plain",
                      srt = 90,
                      cex = 1.45,
                      labels = c("Proteasome Activity Relative",
                                 rep("", 5))[packet.number()])
           panel.text(x = c(-0.4, rep(-10, 5)),
                      y = 1,
                      fontface = "plain",
                      srt = 90,
                      cex = 1.45,
                      labels = c("to BY median (SD units)",
                                 rep("", 5))[packet.number()])
       })
## code end
)
dev.off()
