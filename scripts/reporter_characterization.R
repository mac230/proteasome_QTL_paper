## -----
## load all the required packages
source("~/emacs/R/functions/load_flow_packages.R")

#############
## USER INPUT
#############
## the only things that should need to be changed for this script to work
## for later when we plot
reporter_names <- c("ODC_TFT", "Rpn4_TFT")
## no trailing '/' at the end of 'base_dir'!
base_dir <- "/home/mahlon/data/flow/2021.10.17_TDH3pr_mCh_TFT_ODC_rpn4_fcs"
needed_dirs <- c("/fcs", "/results", "/tables",
                 "/scripts", "/dataframes", "/sessions",
                 "/dataframes/gated", "/dataframes/ungated")
dir_maker <- function(x){
    ifelse(!dir.exists(paths = paste0(base_dir, x)),
           dir.create(path = paste0(base_dir, x)),
           paste0("dir ", paste0(base_dir, x), " exists."))
}
sapply(X = needed_dirs, FUN = dir_maker)
work_dir       <- paste0(base_dir, "/fcs")
results_dir    <- paste0(base_dir, "/results")
tables_dir     <- paste0(base_dir, "/tables")
sessions_dir   <- paste0(base_dir, "/sessions")
frame_dir      <- paste0(base_dir, "/dataframes")
gated_dir      <- paste0(frame_dir, "/gated/")
ungated_dir    <- paste0(frame_dir, "/ungated/")
#################
## END USER INPUT
#################

## -----
## run the loop to generate dataframes for
## each set of fcs files for a given reporter
## we'll merge these sets of dataframes to
## a single dataframe using 'rbind' below.
source("~/emacs/UPS_QTL_paper/scripts/process_fcs_files.R")

## load colors for strains
source("~/emacs/UPS_QTL_paper/scripts/flow_color_setup.R")

## make sure the strains look good:
all_cols

## now, read each reporter's dataframe in and
## combine into a single dataframe
## generate a list of files in a directory using
## the 'dir' command, e.g.:
## dir(gated_dir)
## dir(ungated_dir)

## '_u' = ungated
out_u <- vector(mode = "list", length = length(dir(ungated_dir)))
for (o in 1:length(dir(gated_dir))) {
    out_u[[o]] <- read.table(file = paste0(ungated_dir, dir(ungated_dir)[o]),
                             header = T, sep = ",")
       }
## 5 strains * 8 replicates * 2 reporters * 10,000 cells = 800000
str(out_u) ## looks good

## '_g' = gated
out_g <- vector(mode = "list", length = length(dir(gated_dir)))
for (o in 1:length(dir(gated_dir))) {
    out_g[[o]] <- read.table(file = paste0(gated_dir, dir(gated_dir)[o]),
                             header = T, sep = ",")
       }

str(out_g)

## for the ungated set, it'll be
## 1e5 cells * 5 strains/reporter * 20 reporters = 1e7 rows
## fsc gating reduces 1e5 to ~2e4, so ~2e6 rows
## nrow(out_all) = 2284942
## if you want to look at ungated cells
## out_all <- do.call("rbind", out_u)
out_all <- do.call("rbind", out_g)

str(out_all)
levels(out_all$reporter)
levels(out_all$strain)
levels(out_all$file)

aa_order <- c(1, 2)
out_all$aa_order <- factor(out_all$reporter,
                           levels = levels(out_all$reporter)[aa_order])

out_all$strain_char <- as.character(out_all$strain)
out_all <- out_all[out_all$strain == "BY_strain"   |
                   out_all$strain == "RM_strain"   |
                   out_all$strain == "rpn4_strain" , ]


## get the strain factor in the desired order
## levels(out_all$strain)[c(2, 1, 6, 5, 4, 7)]
strain_order <- c(1, 3, 4, 5, 2)
levels(out_all$strain)[strain_order]
out_all$strain_order <- factor(out_all$strain,
                               levels = levels(out_all$strain)[strain_order])

strain_paste <- expand.grid(unique(out_all$replicate),
                          levels(out_all$strain_order))

strain_paste <- paste0(strain_paste$Var2, "_", strain_paste$Var1)

out_all$strain_rep <- factor(paste0(out_all$strain_order, "_", out_all$replicate),
                             levels = strain_paste)

params <- colnames(out_all)[unlist(lapply(X = out_all, FUN = is.numeric))]
params[10] <- "log2 TFT Ratio"

## -----
## extract the median from each biological replicate and use this
## value to build a dataframe w/ 8 observations per strain per
## reporter.  This dataframe is what we'll use for stats and for
## creating stripcharts and heatmaps

## 'aggregate' creates a new dataframe from x by applying FUN to
## all unique combinations of the factors supplied to the 'by'
## argument - in this case, grab the mean of numeric data and
## keep everything else a factor
out_agg <- aggregate.data.frame(x = out_all,
                                by = list(out_all$strain_order,
                                          out_all$replicate,
                                          out_all$aa_order),
                                FUN = function(x) {
                                    ifelse(is.numeric(x), median(x), as.character(x))
                                },
                                simplify = T)

## 'aggregate' seems to strip the levels from factors, so add
## these back using the values present in the original dataframe
out_agg$strain_fac <- factor(x = out_agg$strain,
                             levels = unique(out_agg$strain),
                             labels = unique(out_agg$strain))

out_agg$file_fac <- factor(out_agg$file,
                           levels = unique(out_agg$file),
                           labels = unique(out_agg$file))

out_agg$aa_order_fac <- factor(out_agg$aa_order,
                               levels = unique(out_agg$aa_order),
                               labels = unique(out_agg$aa_order))

out_agg$strain_rep_fac <- factor(out_agg$strain_rep,
                                 levels = unique(out_agg$strain_rep),
                                 labels = unique(out_agg$strain_rep))


out_agg$n_TFT_ratio <- -1 * out_agg$TFT_ratio

## -----
## transform the data
## Now we need to transform the data so that
## the TFT ratio scales such that:
## high TFT ratio = high UPS activity
## low TFT ratio = low UPS activity
## We'll do this by:
## 1. multiplying log2 TFT ratios by -1
## 2. subtracting the min. TFT value from each replicate
## this has to be done on a per-reporter basis, so
## we'll use the usual for() loop followed by do.call("rbind")
out_agg_sub <- out_agg[out_agg$strain_order == "BY_strain"   |
                       out_agg$strain_order == "RM_strain"   |
                       out_agg$strain_order == "rpn4_strain", ]

m <- 1
out_medians <- list()
for (m in 1:length(unique(out_agg_sub$aa_order))) {
    reporter <- unique(out_agg_sub$aa_order)[m]
    current_frame <- out_agg_sub[out_agg_sub$aa_order == reporter, ]
    current_mean <- mean(current_frame$n_TFT_ratio)
    current_sd   <- sd(current_frame$n_TFT_ratio)
    current_frame$z_TFT <- (current_frame$n_TFT_ratio - current_mean) / current_sd
    out_medians[[m]] <- current_frame
}

out_medians <- do.call("rbind", out_medians)

save(out_medians,
     file = "~/emacs/ubi_QTL_paper/results/2021.11.20_ODC_rpn4_flow_out.RData")

## -----
## statistical analysis
## statistical analysis of out medians on s_TFT_ratio

## -----
## ODC TFT
ODC_out <- out_medians[out_medians$aa_order == "ODC TFT", ]
ODC_aov <- aov(ODC_out$TFT_ratio ~ ODC_out$strain_order)
ODC_ph  <- TukeyHSD(ODC_aov)
ODC_ph_out <- as.data.frame(ODC_ph$ODC_out)
ODC_ph_out$reporter <- rep("ODC TFT", nrow(ODC_ph_out))

## save output
write.table(x = as.data.frame(ODC_ph$ODC_out),
            file = paste0(tables_dir,
                          "/",
                          gsub(pattern = " .*",
                               replacement = "",
                               x = Sys.time()), "_ODC_TFT_ph_out.csv"),
            append = F, quote = F, sep = ",",
            row.names = F, col.names = T)


## -----
## rpn4 TFT
rpn4_out <- out_medians[out_medians$aa_order == "Rpn4 TFT", ]
rpn4_out <- out_medians[out_medians$aa_order == "Rpn4 TFT", ]
rpn4_aov <- aov(rpn4_out$n_TFT_ratio ~ rpn4_out$strain_order)
rpn4_ph  <- TukeyHSD(rpn4_aov)
rpn4_ph_out <- as.data.frame(rpn4_ph$rpn4_out)
rpn4_ph_out$reporter <- rep("rpn4 TFT", nrow(rpn4_ph_out))

t.test(x = rpn4_out$TFT_ratio[rpn4_out$strain == "BY_strain"],
       y = ODC_out$TFT_ratio[ODC_out$strain == "BY_strain"],
       alternative = "two.sided")


## -----
## interaction
out_subset <- out_medians[out_medians$strain == "BY_strain" |
                          out_medians$strain == "RM_strain", ]

int_aov <- aov(TFT_ratio ~ strain * reporter,
               data = out_subset)

## save output
write.table(x = as.data.frame(rpn4_ph$rpn4_out),
            file = paste0(tables_dir,
                          "/",
                          gsub(pattern = " .*",
                               replacement = "",
                               x = Sys.time()), "_rpn4_TFT_ph_out.csv"),
            append = F, quote = F, sep = ",",
            row.names = F, col.names = T)


save(out_agg, file = paste0(tables_dir, "/",
                            gsub(pattern = " .*",
                                 replacement = "",
                                 x = Sys.time()), "ODC_rpn4_TFT_out_agg.Rdata"))


## color setup prior to plotting
color_subset <- as.logical(sapply(X = names(all_cols),
                                  FUN = function(x) {
                                      max(x == as.character(unique(out_subset_medians$strain_order)))
                                  }))

out_medians$sub_strain <- as.factor(out_medians$strain_char)

all_cols <- c("#4579AC", "#B65656", "#BA9C5E")

## -----
## strip plot

s_plot <- xyplot(n_TFT_ratio ~ sub_strain | aa_order,
                 data = out_medians,
                 type = c("p"),
                 ylim = c(0.8, 3.3),
                 cex = 2,
                 col = "black",
                 ylab = expression("UPS Activty (-log"[2]*" RFP / GFP)"),
                 xlab = "",
                 scales = list(alternating = F,
                               tck = c(1, 0),
                               x = list(labels = c("BY",
                                                   "RM",
                                                   expression(paste("BY "*italic("rpn4")*Delta))),
                                        cex = 1.75),
                               y = list(at = seq(from = 0.0, to = 3.5, by = 0.5),
                                        cex = 1.75)),
                 par.strip.text = list(cex = 1.75),
                 horizontal = F,
                 par.settings = list(strip.background = list(col = gray(0.9)),
                                     par.ylab.text = list(cex = 1.75)),
                 panel = function(...) {
                     panel.stripplot(...,
                                     jitter.data = T,
                                     fill = all_cols[1:3],
                                     pch = 21,
                                     factor = 1.25)
                 })

pdf(file = paste0(results_dir, "/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_ODC_rpn4_stripchart.pdf"),
    width = 12, height = 7)
print(s_plot)
grid.text(label = rep("#", 3),
          x = c(0.45, 0.75, 0.9),
          y = c(0.64, 0.64, 0.64),
          default.units = "npc",
          gp = gpar(cex = 2.0))
dev.off()


## -----
## density plot
out_all_sub <- out_all[out_all$strain == "BY_strain" |
                       out_all$strain == "RM_strain" |
                       out_all$strain == "rpn4_strain", ]

out_all_sub$sub_strain <- as.factor(out_all_sub$strain_char)

out_all_sub$n_TFT_ratio <- -1 * out_all_sub$TFT_ratio

all_t_cols <- c("#2166AC77", "#BF323277", "#CCA14A77")

rep_cols <- unlist(lapply(X = 1:3, FUN = function(x) {
                              rep(all_cols[x],
                                  times = sum(grepl(pattern = names(all_cols[x]),
                                                    x = levels(out_all$strain_rep))))
                          }))

d_plot <- densityplot(~ n_TFT_ratio | aa_order,
                      data = out_all_sub,
                      plot.points = F,
                      xlim = c(-0.25, 5.25),
                      ylim = c(-0.1, 1.7),
                      col = all_t_cols,
                      index.cond = list(aa_order),
                      xlab = expression("UPS Activty (-log"[2]*" RFP / GFP)"),
                      scales = list(tck = c(1, 0),
                                    alternating = F,
                                    x = list(at = seq(from = 0.0, to = 5, by = 1),
                                             cex = 1.75),
                                    y = list(at = seq(from = 0.0, to = 1.5, by = 0.5),
                                             cex = 1.75)),
                      par.strip.text = list(cex = 2),
                      par.settings = list(strip.background = list(col = gray(0.9)),
                                          par.xlab.text = list(cex = 1.75),
                                          par.ylab.text = list(cex = 1.75)),
                      panel = function(x, y, q, subscripts, ...) {
                          panel.densityplot(x,
                                            plot.points = F,
                                            groups = out_all_sub$strain_rep,
                                            subscripts = subscripts,
                                            lty = 1,
                                            lwd = 1.5,
                                            col = rep_cols,
                                            ylim = c(0, 2))
                          panel.densityplot(x,
                                            groups = out_all_sub$sub_strain,
                                            data = out_all_sub,
                                            subscripts = subscripts,
                                            ...,
                                            lwd = 7.5)
                      })
pdf(file = paste0(results_dir, "/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_ODC_rpn4_density_plot.pdf"),
    height = 7, width = 12)
print(d_plot)
grid.lines(x = c(0.09, 0.13),
           y = c(0.85, 0.85),
           default.units = "npc",
           gp = gpar(col = all_cols[1],
                     lwd = 10,
                     lineend = "butt"))
grid.lines(x = c(0.09, 0.13),
           y = c(0.80, 0.80),
           default.units = "npc",
           gp = gpar(col = all_cols[2],
                     lwd = 10,
                     lineend = "butt"))
grid.lines(x = c(0.09, 0.13),
           y = c(0.75, 0.75),
           default.units = "npc",
           gp = gpar(col = all_cols[3],
                     lwd = 10,
                     lineend = "butt"))
grid.text(label = c("BY", "RM",
                    expression(paste("BY "*italic("rpn4")*Delta))),
          x = rep(0.135, 3),
          y = c(0.85, 0.80, 0.745),
          just = "left",
          default.units = "npc",
          gp = gpar(cex = 1.6))
dev.off()
