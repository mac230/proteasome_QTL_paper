## module load R/3.6.0
## R

## -----
## setup
## load the big, or a reduced VCF, call it "af_data"
base_dir <- "~/emacs/UPS_QTL_paper/variant_analysis/"
load(paste0(base_dir, "PeterGenotypesRMAlt_200505.RData"))

## simplify the name for subsequent work
af_data <- PeterGenotypesRMAlt

## load a table with each variant's chromosome and position
variants <- data.frame(Position = 411461, Chr = "chrVII", Gene = "RPT6", Causal = NA)

## here's the information on the fields supplied in the matrix:
## "GT:AD:DP:GQ:PGT:PID:PL"

## GT = The genotype of this sample at this site. For a diploid organism, the GT
## field indicates the two alleles carried by the sample, encoded by a 0 for the
## REF allele, 1 for the first ALT allele, 2 for the second ALT allele,
## etc. When there is a single ALT allele (by far the more common case), GT will
## be either:
## 0/0 : the sample is homozygous reference
## 0/1 : the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
## 1/1 : the sample is homozygous alternate

## AD is the unfiltered allele depth, i.e. the number of reads that support each
## of the reported alleles. All reads at the position (including reads that did
## not pass the variant caller’s filters) are included in this number, except
## reads that were considered uninformative. Reads are considered uninformative
## when they do not provide enough statistical evidence to support one allele
## over another.

## DP is the filtered depth, at the sample level. This gives you the number of
## filtered reads that support each of the reported alleles. You can check the
## variant caller’s documentation to see which filters are applied by
## default. Only reads that passed the variant caller’s filters are included in
## this number. However, unlike the AD calculation, uninformative reads are
## included in DP.
head(af_data)
## GQ = genotype quality (phred score)
ret <- vector(length = nrow(variants))
for (i in 1:nrow(variants)) {
    for (i in 1:nrow(variants)) {

        ## assign a variant to 'x'
        x  <- variants[i, ]

        ## get Yhe row corresponding to the variant from the 1,011 matrix
        gt <- af_data[af_data$CHROM == x[1, "Chr"] &
                      as.numeric(af_data$POS) == as.numeric(x[1, "Position"]), ]

        ## if we unambiguously get a single variant ('if(nrow(gt) == 1)'),
        if(nrow(gt) == 1){

        thisINFO <- gt[,"INFO"]

        ## double 'strsplit' call here to go from, e.g.,
        ## "1/1:1,247:248:99:.:.:10388,717,0" to "1", "1"
        RMAlleles <- strsplit(strsplit(gt$AAA, ":")[[1]][1], "/")[[1]]

        RMAltAllele <- as.numeric(RMAlleles[which(!RMAlleles %in% c("0"))][1])

        ## the matrix already contains the allele frequencies, so
        ## just get the corresponding field, 'strsplit' it, and convert
        ## from chr to numeric; the allele frequency is in the first
        ## element of the list and the 2nd element of the vector
        thisAFField <- strsplit(as.character(thisINFO), ";")[[1]][2]
        if(grepl(pattern = "AF=", x = thisAFField)){
            ret[i] <- as.numeric(strsplit(strsplit(thisAFField, "=")[[1]][2], ",")[[1]][RMAltAllele])
        }
    }
}
}

var_info <- list()

## build the final table and write out
variants$pop_RM_AF <- ret

write.table(x = variants,
            file = "~/emacs/ubi_QTL_paper/results/411461_pop_AF_frequency.csv",
            append = F, quote = F, sep = ",",
            row.names = T, col.names = T)


## -----
## make a table for creating the variant history tree
i <- 1
for (i in 1:nrow(variants)) {

    ## assign a variant to 'x'
    x  <- variants[i, ]

    ## get the row corresponding to the variant from the 1,011 matrix
    gt <- af_data[af_data$CHROM == x[1, "Chr"] &
                  as.numeric(af_data$POS) == as.numeric(x[1, "Position"]), ]

    ## returns a vector of genotypes for the rows of the matrix
    ## we'll later cbind these to create the dataframe that we need
    var_info[[i]] <- sapply(X = gt[, 10:ncol(gt)],
                            FUN = function(strain) {
                                strsplit(x = strain, split = ":")[[1]][1]
                                })
}


var_info_final <- do.call("cbind", var_info)
var_info_final <- as.data.frame(var_info_final)

## these will become the rownames for the final table
strain_list <- colnames(gt[, 10:ncol(gt)])

## these will be the column names for the final table
variant_list <- paste0(variants$Position, "_",
                       variants$Gene)

rownames(var_info_final) <- strain_list
colnames(var_info_final) <- variant_list

write.table(x = var_info_final,
            file = "~/emacs/ubi_QTL_paper/results/QTN_tree_data.csv",
            append = F, quote = F, sep = ",",
            row.names = T, col.names = T)
