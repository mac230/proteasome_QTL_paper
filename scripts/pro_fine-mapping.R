## -----
## START USER INPUT
## no trailing '/' at the end!
base_dir  <- "~/data/flow/2022.03.10_RPT6_fine-mapping_pro_TFT_flow"
frame_dir <- paste0(base_dir, "/dataframes")
stats_log <- paste0(base_dir, "/results/stats_log")
## END USER INPUT 


## -----
## load packages and process fcs files 
source("~/emacs/R/functions/load_flow_packages.R")
source("~/emacs/R/functions/color_setup_BY_RM.R")
source("~/emacs/R/functions/auto_process_fcs_files.R")


## -----
## load the raw and processed data:
load(file = paste0(frame_dir, "/ungated_final.R"))
## load the directory structure for saving output:
load(file = paste0(frame_dir, "/dir_structure.R"))
## load the extracted median for ea. parameter and bio. replicate: 
load(file = paste0(frame_dir, "/all_medians.R"))


## -----
## inferential statistics using 'aov' and 'TukeyHSD':
source("~/emacs/R/functions/univariate_fcs_analysis.R")


## -----
## lmer regression analysis:
source("~/emacs/R/functions/lmer_RPT6_flow_analysis.R")


## -----
## by genotype plot

## addn. variable: genotype including edit 
all_medians$genotype <- ifelse(grepl(pattern = ".*RM.*",
                                     x = all_medians$strain),
                               "edit_RM",
                        ifelse(grepl(pattern = ".*gRNA.*",
                                     x = all_medians$strain),
                               "edit_BY", "BY"))
all_medians$genotype_factor <- as.factor(all_medians$genotype)

## guide RNA
all_medians$guide <- ifelse(!grepl(pattern = ".*gRNA.*",
                                   x = all_medians$strain),
                            "wild-type",
                     ifelse(grepl(pattern = ".*gRNA_02.*",
                                  x = all_medians$strain),
                            "gRNA_02", "gRNA_03"))
all_medians$guide_factor <- as.factor(all_medians$guide)

## relevel to make BY / wild-ytpe the reference level
all_medians$strain_factor <- relevel(x = all_medians$strain_factor,
                                     ref = "BY_strain")

genotype_cols <- ifelse(test = grepl(pattern = ".*BY",
                                     x = levels(all_medians$genotype_factor)),
                        yes = col_by,
                        no = col_rm)

## y axis ranges
TFT_l_range  <- range(all_medians$TFT_loess[all_medians$gating == "gated"])
l_down_range <- TFT_l_range[1] - (0.02 * TFT_l_range[1])
l_up_range   <- TFT_l_range[2] + (0.02 * TFT_l_range[2])
yl_range     <- c(l_down_range, 0.71)

per_genotype_pro_plot <- xyplot(TFT_loess ~ genotype_factor,
                                 groups = genotype_factor,
                                 subset = all_medians$gating == "gated",
                                 data = all_medians,
                                 ylim = yl_range,
                                 xlab = "",
                                 ylab = expression("Proteasome Activity (-log"["2"]*" RFP / GFP)"),
                                 scales = list(tck = c(1, 0),
                                               alternating = F,
                                               x = list(labels = c("BY wild-type",
                                                                   "BY edited",
                                                                   "BY edited"),
                                                        cex = 1.25),
                                               y = list(cex = 1.25)),
                                 ## at = seq(from = 0.31, to = 0.42, by = 0.02))),
                                 par.settings = list(box.dot = list(pch = "|"),
                                                     box.umbrella = list(lty = 1,
                                                                         lwd = 1,
                                                                         col = gray(0.4)),
                                                     box.rectangle = list(lty = 1,
                                                                          lwd = 1,
                                                                          col = gray(0.4),
                                                                          fill = gray(0.97)),
                                                     strip.background = list(col = gray(0.9)),
                                                     par.ylab.text = list(cex = 1.25),
                                                     par.xlab.text = list(cex = 1.25),
                                                     clip = list(panel = F)),
                                 
                                 panel = function(...) {

                                     ## boxplot
                                     panel.bwplot(...,
                                                  ## pch = "|",
                                                  box.width = 0.4,
                                                  horizontal = F,
                                                  do.out = F)
                                     
                                     ## thicker median lines
                                     dat <- all_medians[all_medians$gating == "gated", ]
                                     panel.segments(x0 = 1:length(levels(dat$genotype_factor)) - 0.2,
                                                    x1 = 0.2 + 1:length(levels(dat$genotype_factor)),
                                                    y0 = tapply(X = dat$TFT_loess,
                                                                INDEX = dat$genotype_factor,
                                                                FUN = median),
                                                    y1 = tapply(X = dat$TFT_loess,
                                                                INDEX = dat$genotype_factor,
                                                                FUN = median),
                                                    col = gray(0.4), lwd = 4)
                                     panel.stripplot(...,
                                                     pch = 21,
                                                     cex = 1.1,
                                                     fill = genotype_cols,
                                                     col = "black",
                                                     jitter.data = T,
                                                     amount = 0.1,
                                                     horizontal = F)
                                     
                                     ## strain name and x axis labels 
                                     panel.text(x = 1:3,
                                                y = rep(0.517, 3),
                                                labels = c("-",
                                                           expression(italic("RPT6")*" -175 BY"),
                                                           expression(italic("RPT6")*" -175 RM")),
                                                cex = 1.25)
                                     panel.text(x = c(0.22, 0.1),
                                                y = c(0.526, 0.517),
                                                labels = c("Strain",
                                                           "Allelic Edit"),
                                                cex = 1.25)
                                     ## [x] pvals and connecting segments
                                     panel.text(x = c(2.5, 2.0),
                                                y = c(0.69, 0.703),
                                                labels = c(expression(italic("p")*" = 2.1e-2"),
                                                           expression(italic("p")*" = 7.0e-3")),
                                                cex = 1.25)
                                     panel.segments(x0 = c(2, 1), x1 = c(3, 3),
                                                    y0 = c(0.685, 0.698), y1 = c(0.685, 0.698),
                                                    lwd = 1.2, col = gray(0.4))
                                     ## [x] reporter name as title
                                     panel.text(x = 2, y = 0.715,
                                                labels = "Pro TFT",
                                                cex = 1.25,
                                                fontface = "bold")
                                 })

pdf(file = "~/emacs/ubi_QTL_paper/figures/per_genotype_plot_pro_TFT_final.pdf")
print(per_genotype_pro_plot)
dev.off()


## -----
## per strain plot
strain_cols <- ifelse(test = grepl(pattern = ".*BY.*",
                                     x = levels(all_medians$strain_factor)),
                        yes = col_by,
                        no = col_rm)

per_strain_pro_plot <- xyplot(TFT_loess ~ strain_factor,
                               groups = strain_factor,
                               subset = all_medians$gating == "gated",
                               data = all_medians,
                               ylim = yl_range,
                               xlab = "",
                               ylab = expression("Proteasome Activity (-log"["2"]*" RFP / GFP)"),
                               scales = list(tck = c(1, 0),
                                             alternating = F,
                                             x = list(labels = c("BY wild-type",
                                                                 "BY edited",
                                                                 "BY edited",
                                                                 "BY edited",
                                                                 "BY edited"),
                                                      cex = 1.2),
                                             y = list(cex = 1.25)),
                               par.settings = list(box.dot = list(pch = "|"),
                                                   box.umbrella = list(lty = 1,
                                                                       lwd = 1,
                                                                       col = gray(0.4)),
                                                   box.rectangle = list(lty = 1,
                                                                        lwd = 1,
                                                                        col = gray(0.4),
                                                                        fill = gray(0.97)),
                                                   strip.background = list(col = gray(0.9)),
                                                   par.ylab.text = list(cex = 1.25),
                                                   par.xlab.text = list(cex = 1.25),
                                                   clip = list(panel = F),
                                                   layout.widths = list(left.padding = 5),
                                                   layout.heights = list(bottom.padding = 5)),
                               panel = function(...) {
                                   ## boxplot
                                   panel.bwplot(...,
                                                ## pch = "|",
                                                box.width = 0.4,
                                                horizontal = F,
                                                do.out = F)
                                   
                                   ## thicker median lines
                                   dat <- all_medians[all_medians$gating == "gated", ]
                                   panel.segments(x0 = 1:length(levels(dat$strain_factor)) - 0.2,
                                                  x1 = 0.2 + 1:length(levels(dat$strain_factor)),
                                                  y0 = tapply(X = dat$TFT_loess,
                                                              INDEX = dat$strain_factor,
                                                              FUN = median),
                                                  y1 = tapply(X = dat$TFT_loess,
                                                              INDEX = dat$strain_factor,
                                                              FUN = median),
                                                  col = gray(0.4), lwd = 4)
                                   panel.stripplot(...,
                                                   pch = 21,
                                                   cex = 1.1,
                                                   fill = strain_cols,
                                                   col = "black",
                                                   jitter.data = T,
                                                   amount = 0.1,
                                                   horizontal = F)
                                   ## strain name and x axis labels 
                                   panel.text(x = c(1, 2.5, 4.5),
                                              y = rep(0.517, 3),
                                              labels = c("-",
                                                         expression(italic("RPT6")*" -175 BY"),
                                                         expression(italic("RPT6")*" -175 RM")),
                                              cex = 1.25)
                                   ## gRNA
                                   panel.text(x = 1:5,
                                              y = rep(0.508, 5),
                                              labels = c("-",
                                                         "1", "2",
                                                         "1", "2"),
                                              cex = 1.25)
                                   ## pvals and connecting segments
                                   panel.text(x = seq(from = 1.5, to = 3.0, by = 0.5),
                                              y = seq(from = 0.681, by = 0.008,
                                                      length.out = 4),
                                              labels = c(expression(italic("p")*" = 0.14"),
                                                         expression(italic("p")*" = 0.85"),
                                                         expression(italic("p")*" = 0.049"),
                                                         expression(italic("p")*" = 2e-3")),
                                              cex = 1.25)
                                   panel.segments(x0 = rep(1, 4),
                                                  x1 = 2:5,
                                                  y0 = seq(from = 0.678,
                                                           by = 0.008,
                                                           length.out = 4),
                                                  y1 = seq(from = 0.678,
                                                           by = 0.008,
                                                           length.out = 4),
                                                  lwd = 1.2, col = gray(0.4))
                                   ## x axis descriptors
                                     panel.text(x = c(0.08, -0.14, 0.07),
                                                y = c(0.526, 0.517, 0.508),
                                                labels = c("Strain",
                                                           "Allelic Edit",
                                                           "gRNA"),
                                                cex = 1.25)
                                   ## reporter name as title
                                   panel.text(x = 3, y = 0.715,
                                              labels = "Pro TFT",
                                              cex = 1.25,
                                              fontface = "bold")
                               })

pdf(file = "~/emacs/ubi_QTL_paper/figures/per_strain_plot_pro_TFT.pdf")
print(per_strain_pro_plot)
dev.off()

## -----
## ANOVA and lmer analysis of:
## BY full vs. RM full
## BY full vs. BY

pro_analysis_log <- paste0(base_dir, "/results/pro_analysis_log")


## -----
## BY full vs. RM full

## ANOVA
subd <- all_medians[all_medians$strain != "BY_strain", ]
subd <- subd[subd$gating == "gated", ]
str(subd)
unique(subd$genotype)

pro_aov <- aov(formula = TFT_loess ~ genotype_factor * guide_factor,
                data = subd)

write(x = paste0("\n\n", Sys.time(),
                 " ANOVA of gated BY_full vs RM_full for pro TFT:\n\n"),
      file = pro_analysis_log,
      append = T)
capture.output(summary(pro_aov),
               file = pro_analysis_log,
               append = T)

## lmer
pro_model <- lmer(TFT_loess ~ genotype_factor + guide_factor + (1|plate),
                       data = subd,
                       REML = F)
write(x = paste0("\n\n", Sys.time(),
                 " lmer of gated BY_full vs RM_full for pro TFT:\n\n"),
      file = pro_analysis_log,
      append = T)
capture.output(summary(pro_model),
               file = pro_analysis_log,
               append = T)


## -----
## BY wild-type vs. BY full 
## ANOVA
subd <- all_medians[all_medians$genotype != "edit_RM", ]
subd <- subd[subd$gating == "gated", ]
str(subd)
unique(subd$genotype)

pro_aov <- aov(formula = TFT_loess ~ genotype_factor,
                data = subd)

write(x = paste0("\n\n", Sys.time(),
                 " ANOVA of gated BY_wild_type vs BY_full for pro TFT:\n\n"),
      file = pro_analysis_log,
      append = T)
capture.output(summary(pro_aov),
               file = pro_analysis_log,
               append = T)

## lmer
pro_model <- lmer(TFT_loess ~ genotype_factor + (1|plate),
                       data = subd,
                       REML = F)
write(x = paste0("\n\n", Sys.time(),
                 " lmer of gated BY_wild_type vs BY_full for pro TFT:\n\n"),
      file = pro_analysis_log,
      append = T)
capture.output(summary(pro_model),
               file = pro_analysis_log,
               append = T)
