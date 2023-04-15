stopifnot(file.exists(paste0("~/emacs/proteasome_QTL_paper/",
                                "scripts_data/replicating_peaks_table/",
                                "ODC_Rpn4_merged_replicating_peaks.csv")))

dat <- read.table(file = paste0("~/emacs/proteasome_QTL_paper/",
                                "scripts_data/replicating_peaks_table/",
                                "ODC_Rpn4_merged_replicating_peaks.csv"),
                 header = T, sep = ",", stringsAsFactors = F)

dat$avg_LOD <- (dat$rep_1_LOD + dat$rep_2_LOD) / 2
dat$avg_delta_AF <- (dat$rep_1_delta_AF + dat$rep_2_delta_AF) / 2
dat$avg_right_Index <- (dat$rep_1_right_Index + dat$rep_2_right_Index) / 2
dat$avg_left_Index <- (dat$rep_1_left_Index + dat$rep_2_left_Index) / 2
dat$avg_max_Index <- (dat$rep_1_max_Index + dat$rep_2_max_Index) / 2

data.frame(reporter = dat$reporter,
           chr = dat$chr,
           LOD = dat$avg_LOD,
           delta_AF = dat$avg_delta_AF,
           left_index = dat$avg_left_Index,
           max_index = dat$avg_max_Index,
           right_index = dat$avg_right_Index)
