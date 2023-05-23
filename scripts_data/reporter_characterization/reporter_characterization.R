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
base_dir <- "~/emacs/proteasome_QTL_paper/scripts/reporter_characterization"
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
str(out_u) 

## '_g' = gated
out_g <- vector(mode = "list", length = length(dir(gated_dir)))
for (o in 1:length(dir(gated_dir))) {
    out_g[[o]] <- read.table(file = paste0(gated_dir, dir(gated_dir)[o]),
                             header = T, sep = ",")
       }
str(out_g)


## load colors for strains
source("~/emacs/UPS_QTL_paper/scripts/flow_color_setup.R")


## for the ungated set, it'll be
## 1e5 cells * 5 strains/reporter * 20 reporters = 1e7 rows
## fsc gating reduces 1e5 to ~2e4, so ~2e6 rows
## nrow(out_all) = 2284942
## if you want to look at ungated cells
## out_all <- do.call("rbind", out_u)
out_all <- do.call("rbind", out_g)
str(out_all)

out_all$strain_factor <- as.factor(out_all$strain)
levels(out_all$strain_factor)

out_all$reporter_factor <- as.factor(out_all$reporter)
levels(out_all$reporter_factor)

out_all$file_factor <- as.factor(out_all$file)
levels(out_all$file_factor)

aa_order <- c(1, 2)
out_all$aa_order <- factor(out_all$reporter,
                           levels = levels(out_all$reporter_factor)[aa_order])

out_all$strain_char <- as.character(out_all$strain)
out_all <- out_all[out_all$strain == "BY_strain"   |
                   out_all$strain == "RM_strain"   |
                   out_all$strain == "rpn4_strain" , ]


## get the strain factor in the desired order
strain_order <- c(1, 3, 4)
out_all$strain_order <- factor(out_all$strain,
                               levels = levels(out_all$strain_factor)[strain_order])

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
## high TFT ratio = high proteasome activity
## low TFT ratio = low proteasome activity
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
     file = "~/emacs/proteasome_QTL_paper/scripts/reporter_characterization/plot_analysis_data.RData")
