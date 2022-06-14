## load the dataframe w/ information on whether the RM allele is derived
load(file = "~/emacs/UPS_QTL_paper/variant_analysis/R_RMDerived_withInfo_200718.RData")

str(RMDerived_withInfo)
## CHROM = chromosome
## POS = variant location
## ID = ?
## REF = reference allele
## ALT_allele = alternative allele
## RMDerived = whether RM allele is derived (T = yes, F = no, NA = unknown)

## simplify the name of the dataset
dat <- RMDerived_withInfo

## information on our QTGs for extracting variants from defined genomic ranges 
genes <- "RPT6"
chr   <- "chrVII"
## genomic ranges for QTGs; go 500 bp up/down to get promoter/terminator
gene_lower_value <- 411400
gene_upper_value <- 411600

qtg_table <- data.frame(gene = genes,
                        chr = chr,
                        lower = gene_lower_value,
                        upper = gene_upper_value)

qtg_out <- list()

for(i in 1:nrow(qtg_table)) {

    ## extract the information for each gene's genomic range;
    ## add to a list that we'll then collapse to a dataframe
    q <- dat[dat$CHROM == qtg_table[i, "chr"] &
             dat$POS > qtg_table[i, "lower"] &
             dat$POS < qtg_table[i, "upper"], ]
    q$gene <- rep(qtg_table[i, "gene"], nrow(q))
    qtg_out[[i]] <- q

}

qtg_final <- do.call("rbind", qtg_out)

qtg_final ## RM allele at 411461 is derived
