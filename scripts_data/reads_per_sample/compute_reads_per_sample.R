dat   <- read.csv(file = paste0("~/emacs/proteasome_QTL_paper/",
                                "scripts_data/reads_per_sample/reads.csv"),
                  header = T,
                  sep = ",")

reads <- dat$reads
cover <- dat$median_coverage
paste0("total reads = ", sum(reads))
## "total reads = 153887828"

paste0("median reads = ", median(reads))
## "median reads = 9757090.5"

paste0("min reads = ", min(reads))
"min reads = 5994921"

paste0("max reads = ", max(reads))
"max reads = 14753319"

paste0("coverage median = ", median(cover), 
       "; coverage range = ", min(cover), " to ", max(cover))
"coverage median = 21; coverage range = 16 to 25"
