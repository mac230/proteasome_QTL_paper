## -----
## set up additional variables: genotype, edit, and guide RNA for lmer analysis
## genotype
all_medians$background <- ifelse(grepl(pattern = ".*BY.*",
                                     x = all_medians$strain),
                               "BY", "RM")
all_medians$background_factor <- as.factor(all_medians$background)

## genotype including edit 
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


## relevel factors to make BY / wild-ytpe the reference level
all_medians$strain_factor <- relevel(x = all_medians$strain_factor,
                                     ref = "BY_strain")

all_medians$guide_factor <- relevel(x = all_medians$guide_factor,
                                     ref = "wild-type")


## -----
## lmer analysis
## write("", stats_log, append = F)
write("\n\n\n------------- lmer analysis", stats_log, append = T)
sapply(X = c("ungated", "gated"),
       FUN = function(x) {

           ## subset the data 
           dat <- all_medians[all_medians$gating == x, ]

           ## write output
           write(paste0(Sys.time(), " analyis of TFT ratio on ",
                        x, " cells"),
                 file = stats_log, append = T)
     
           ## lmer model w/ genetic background and edit as fixed factors
           genotype_model <- lmer(TFT_loess ~ genotype_factor + guide_factor + (1|plate),
                                  data = dat,
                                  REML = F)
           write(paste0("\n\n---------------------\n", x, " genotype lmer\n"),
                 file = stats_log, append = T)
           capture.output(summary(genotype_model),
                          file = stats_log, append = T)

           ## pairwise t-test comparing genetic backgrounds 
           genotype_t_test <- pairwise.t.test(x = dat$TFT_loess,
                                              g = dat$genotype_factor,
                                              p.adjust.method = "none")
           write(paste0("\n\n---------------------\n", x, " genotype pairwise t-test\n"),
                 file = stats_log, append = T)
           capture.output(genotype_t_test,
                          file = stats_log, append = T)

           ## lmer model comparing all strains individually
           per_strain_model <- lmer(TFT_loess ~ strain_factor + (1|plate),
                                    data = dat,
                                    REML = F)
           write(paste0("\n\n---------------------\n", x, " per strain lmer\n"),
                 file = stats_log, append = T)
           capture.output(summary(per_strain_model),
                          file = stats_log, append = T)

           ## pairwise t-test comparing all strains individually
           per_strain_t_test <- pairwise.t.test(x = dat$TFT_loess,
                                                g = dat$strain_factor,
                                                p.adjust.method = "none")
           write(paste0("\n\n---------------------\n", x, " per strain pairwise t-test\n"),
                 file = stats_log, append = T)
           capture.output(per_strain_t_test,
                          file = stats_log, append = T)
       })


## -----
## plots 
## [1] TFT ratio by time and TFT loess by time 
## plot the effect of the adjustment
## get slopes of linear regression of
## time on raw and corrected TFT ratio
## for adding to the plots we'll produce
u_text <- lm(TFT_ratio ~ time_rel,
             subset = all_medians$gating == "ungated",
             data = all_medians)$coefficients["time_rel"]

c_text <- lm(TFT_loess ~ time_rel,
             subset = all_medians$gating == "ungated",
             data = all_medians)$coefficients["time_rel"]

## get y axis range (R default +/- 4% of max / min, respectively)
TFT_r_range  <- range(all_medians$TFT_ratio[all_medians$gating == "ungated"])
r_down_range <- TFT_r_range[1] - (0.04 * TFT_r_range[1])
r_up_range   <- TFT_r_range[2] + (0.04 * TFT_r_range[2])
yr_range    <- c(r_down_range, r_up_range)

TFT_l_range <- range(all_medians$TFT_loess[all_medians$gating == "ungated"])
l_down_range <- TFT_l_range[1] - (0.04 * TFT_l_range[1])
l_up_range   <- TFT_l_range[2] + (0.04 * TFT_l_range[2])
yl_range    <- c(l_down_range, l_up_range)


## -----
## TFT ratio by time_rel before / after loess adjustment
## 'TFT_ratio_time_adjustment.pdf'
## two scatter plots w/ no grouping by genotype
uncorrected_plot <- xyplot(TFT_ratio ~ time_rel,
                           subset = all_medians$gating == "ungated", 
                           data = all_medians,
                           ylim = yr_range,
                           ylab = "Raw TFT Ratio",
                           xlab = "Relative Time (sec)",
                           scales = list(tck = c(1, 0)),
                           panel = function(...) {
                               panel.xyplot(...)
                               panel.abline(lm(TFT_ratio ~ time_rel,
                                               subset = all_medians$gating == "ungated",
                                               data = all_medians),
                                            lty = 2, col = gray(0.7))
                               panel.text(labels = paste0("slope = ", u_text),
                                          x = 500,
                                          y = quantile(x = all_medians$TFT_ratio,
                                                       probs = 0.9))
                           })

corrected_plot <- xyplot(TFT_loess ~ time_rel,
                         subset = all_medians$gating == "ungated", 
                         data = all_medians,
                         ylim = yr_range,
                         ylab = "Adjusted TFT Ratio",
                         xlab = "Relative Time (sec)",
                         scales = list(tck = c(1, 0)),
                         panel = function(...) {
                             panel.xyplot(...)
                             panel.abline(lm(TFT_loess ~ time_rel,
                                             subset = all_medians$gating == "ungated",
                                             data = all_medians),
                                          lty = 2, col = gray(0.7))
                             panel.text(labels = paste0("slope = ", c_text),
                                        x = 500,
                                        y = quantile(x = all_medians$TFT_loess,
                                                     probs = 0.9))
                         })

pdf(file = paste0(results_dir, "/TFT_ratio_time_adjustment.pdf"),
    height = 7, width = 14)
print(uncorrected_plot,
      split = c(1, 1, 2, 1),
      newpage = T)
print(corrected_plot,
      split = c(2, 1, 2, 1),
      newpage = F)
dev.off()


## -----
## TFT ratio by background (BY, BY edit, or RM)
## 'TFT_ratio_adjustment_by_genotype.pdf'
## two scatter plots w/ grouping by genotype
genotype_cols <- ifelse(test = grepl(pattern = ".*BY",
                                     x = levels(all_medians$genotype_factor)),
                        yes = col_by,
                        no = col_rm)

raw_genotype_plot <- xyplot(TFT_ratio ~ genotype_factor,
                            groups = genotype_factor,
                            subset = all_medians$gating == "ungated",
                            data = all_medians,
                            main = unique(all_medians$reporter),                            
                            ylim = yr_range,
                            scales = list(tck = c(1, 0),
                                          alternating = F,
                                          x = list(rot = 45,
                                                   cex = 1.25),
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
                                                clip = list(panel = F)),
                            
                            panel = function(...) {

                                ## boxplot
                                panel.bwplot(...,
                                             ## pch = "|",
                                             box.width = 0.4,
                                             horizontal = F,
                                             do.out = F)
                                
                                ## thicker median lines
                                dat <- all_medians[all_medians$gating == "ungated", ]
                                panel.segments(x0 = 1:length(levels(dat$genotype_factor)) - 0.2,
                                               x1 = 0.2 + 1:length(levels(dat$genotype_factor)),
                                               y0 = tapply(X = dat$TFT_ratio,
                                                           INDEX = dat$genotype_factor,
                                                           FUN = median),
                                               y1 = tapply(X = dat$TFT_ratio,
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
                            })

adj_genotype_plot <- xyplot(TFT_loess ~ genotype_factor,
                            groups = genotype_factor,
                            subset = all_medians$gating == "ungated",
                            data = all_medians,
                            main = unique(all_medians$reporter),
                            ylim = yr_range,
                            scales = list(tck = c(1, 0),
                                          alternating = F,
                                          x = list(rot = 45,
                                                   cex = 1.25),
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
                                                clip = list(panel = F)),
                            
                            panel = function(...) {

                                ## boxplot
                                panel.bwplot(...,
                                             ## pch = "|",
                                             box.width = 0.4,
                                             horizontal = F,
                                             do.out = F)
                                
                                ## thicker median lines
                                dat <- all_medians[all_medians$gating == "ungated", ]
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
                            })

pdf(file = paste0(results_dir, "/TFT_ratio_adjustment_by_genotype.pdf"),
    height = 7, width = 12)
print(raw_genotype_plot, split = c(1, 1, 2, 1), newpage = T)
print(adj_genotype_plot, split = c(2, 1, 2, 1), newpage = F)
dev.off()


## -----
## TFT ratio by background (BY, BY edit, or RM)
## 'TFT_ratio_adjustment_by_strain.pdf'
## two scatter plots w/ grouping by strain
strain_cols <- ifelse(test = grepl(pattern = ".*BY",
                                     x = levels(all_medians$strain_factor)),
                        yes = col_by,
                        no = col_rm)

raw_strain_plot <- xyplot(TFT_ratio ~ strain_factor,
                          groups = strain_factor,
                          subset = all_medians$gating == "ungated",
                          data = all_medians,
                          main = unique(all_medians$reporter),                          
                          ylim = yr_range,
                          scales = list(tck = c(1, 0),
                                        alternating = F,
                                        x = list(rot = 45,
                                                 cex = 1.25),
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
                                              clip = list(panel = F)),
                          
                          panel = function(...) {

                              ## boxplot
                              panel.bwplot(...,
                                           ## pch = "|",
                                           box.width = 0.4,
                                           horizontal = F,
                                           do.out = F)
                              
                              ## thicker median lines
                              dat <- all_medians[all_medians$gating == "ungated", ]
                              panel.segments(x0 = 1:length(levels(dat$strain_factor)) - 0.2,
                                             x1 = 0.2 + 1:length(levels(dat$strain_factor)),
                                             y0 = tapply(X = dat$TFT_ratio,
                                                         INDEX = dat$strain_factor,
                                                         FUN = median),
                                             y1 = tapply(X = dat$TFT_ratio,
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
                          })

adj_strain_plot <- xyplot(TFT_loess ~ strain_factor,
                          groups = strain_factor,
                          subset = all_medians$gating == "ungated",
                          data = all_medians,
                          main = unique(all_medians$reporter),                          
                          ylim = yr_range,
                          scales = list(tck = c(1, 0),
                                        alternating = F,
                                        x = list(rot = 45,
                                                 cex = 1.25),
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
                                              clip = list(panel = F)),
                          
                          panel = function(...) {

                              ## boxplot
                              panel.bwplot(...,
                                           ## pch = "|",
                                           box.width = 0.4,
                                           horizontal = F,
                                           do.out = F)
                              
                              ## thicker median lines
                              dat <- all_medians[all_medians$gating == "ungated", ]
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
                          })

pdf(file = paste0(results_dir, "/TFT_ratio_adjustment_by_strain.pdf"),
    height = 7, width = 12)
print(raw_strain_plot, split = c(1, 1, 2, 1), newpage = T)
print(adj_strain_plot, split = c(2, 1, 2, 1), newpage = F)
dev.off()


## -----
## final plots to separate pdfs
## use 'update' to adjust ylims 
pdf(file = paste0(results_dir, "/TFT_adjusted_ratio_by_genotype.pdf"),
    height = 7, width = 6)
print(update(adj_genotype_plot, ylim = yl_range))
dev.off()

pdf(file = paste0(results_dir, "/TFT_adjusted_ratio_by_strain.pdf"),
    height = 7, width = 6)
print(update(adj_strain_plot, ylim = yl_range))
dev.off()
