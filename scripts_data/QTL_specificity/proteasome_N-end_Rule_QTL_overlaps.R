## -----
## setup
library("lattice")
library("viridis")

dat <- read.table(file = paste0("~/emacs/proteasome_QTL_paper",
                                "/scripts_data/QTL_specificity/",
                                "proteasome_activity_replicating_QTLs.csv"),
                  sep = ",",
                  header = T,
                  stringsAsFactors = F)

                 ## Ac/N-end
amino_acids <- c("Ala", "Cys", "Gly", "Met", "Pro", "Ser", "Thr", "Val",
                 ## Arg/N-end type I
                 "Arg", "Asn", "Asp", "Gln", "Glu", "His", "Lys",
                 ## Arg/N-end type II
                 "Ile", "Leu", "Phe", "Trp", "Tyr")

## shared qtls = V, VIIa, XII
          ## ODC
qtls <- c(paste0("ODC_", c("IIa", "IIb", "IVa", "V", "VIIa", "VIIb",
                           "X", "XII", "XIIIa", "XIIIb", "XIVa")),
          ## rpn4
          paste0("rpn4_", c("IVb", "V", "VIIa",
                            "VIIc", "XII", "XIVb", "XV")))
qtls_rev <- rev(qtls)


## -----
## unordered heatmap:

## create a dataframe for making a grid:
shared_qtls <- as.data.frame(matrix(data = 0,
                                    nrow = 18, ncol = 20))

row.names(shared_qtls) <- qtls_rev
colnames(shared_qtls) <- amino_acids

cols <- c(gray(0.8), magma(n = 13, alpha = 0.8)[2:13])
col_range <- c(0, seq(from = 1, to = 1.0001e5, by = 1e4))


## -----
## ordered heatmap
## this is a 'waffle' plot that shows which
## ubiquitin-independent QTLs overlap w/
## which N-end rule QTLs:

ordered_qtls <- as.data.frame(matrix(data = 0,
                                     nrow = 18, ncol = 11))

row.names(ordered_qtls) <- qtls_rev
colnames(ordered_qtls)  <- as.character(1:11)

ordered_aa_names <- as.data.frame(matrix(data = 0,
                                     nrow = 18, ncol = 11))

aas_ordered <- c("Ala", "Cys", "Gly", "Met",
                 "Pro", "Ser", "Thr", "Val",
                 "Arg", "Asn", "Asp", "Gln",
                 "Glu", "His", "Lys", "Ile",
                 "Leu", "Phe", "Trp", "Tyr")


## rpn4 xv         ## Ala Cys Gly Met   Pro Ser Thr Val Arg   Asn   Asp    Gln Glu   His Lys    Ile Leu Phe   Trp Tyr
shared_qtls[1, ] <- c(0,  0,  0,  4050, 0,  0,  0,  0,  8450, 4300, 10350, 0,  7100, 9550, 300, 0,  0,  7600, 0,  0)
out <- shared_qtls[1, shared_qtls[1, ] > 0][order(shared_qtls[1, shared_qtls[1, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[1, ]     <- out
ordered_aa_names[1, ] <- out_names

## rpn4 xiv        ## Ala Cys Gly Met Pro Ser Thr    Val Arg Asn Asp     Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[2, ] <- c(0,  0,  0,  0,  0,  0,  70200, 0,  0,  0,  76400,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[2, shared_qtls[2, ] > 0][order(shared_qtls[2, shared_qtls[2, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[2, ]     <- out
ordered_aa_names[2, ] <- out_names

## rpn4 xii        ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn Asp     Gln Glu    His     Lys Ile     Leu Phe Trp Tyr
shared_qtls[3, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  19100,  0,  7250,  54650,  0,  23200,  0,  0,  0,  0)
out <- shared_qtls[3, shared_qtls[3, ] > 0][order(shared_qtls[3, shared_qtls[3, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[3, ]     <- out
ordered_aa_names[3, ] <- out_names

## rpn4 VIIc       ## Ala Cys Gly Met Pro Ser Thr Val Arg     Asn    Asp Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[4, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  20550,  8600,  9250,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[4, shared_qtls[4, ] > 0][order(shared_qtls[4, shared_qtls[4, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[4, ]     <- out
ordered_aa_names[4, ] <- out_names

## rpn4 VIIa       ## Ala     Cys    Gly    Met Pro    Ser    Thr    Val    Arg Asn    Asp Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[5, ] <- c(43800,  37300, 60000, 0,  56400, 36000, 22800, 53300, 0,  74350, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[5, shared_qtls[5, ] > 0][order(shared_qtls[5, shared_qtls[5, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[5, ]     <- out
ordered_aa_names[5, ] <- out_names

## rpn4 v          ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn Asp Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[6, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[6, shared_qtls[6, ] > 0][order(shared_qtls[6, shared_qtls[6, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[6, ]     <- out
ordered_aa_names[6, ] <- out_names

## rpn4 iv          ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn Asp Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[7, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[7, shared_qtls[7, ] > 0][order(shared_qtls[7, shared_qtls[7, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[7, ]     <- out
ordered_aa_names[7, ] <- out_names

## odc xiv         ## Ala Cys Gly Met    Pro    Ser Thr Val Arg Asn Asp Gln Glu His    Lys    Ile Leu    Phe    Trp    Tyr
shared_qtls[8, ] <- c(0,  0,  0,  41850, 17450, 0,  0,  0,  0,  0,  0,  0,  0,  17000, 23050, 0,  17200, 22000, 24350, 26250)
out <- shared_qtls[8, shared_qtls[8, ] > 0][order(shared_qtls[8, shared_qtls[8, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[8, ]     <- out
ordered_aa_names[8, ] <- out_names

## odc xiiib       ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn    Asp Gln Glu His    Lys Ile Leu Phe Trp Tyr
shared_qtls[9, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  86050, 0,  0,  0,  88050, 0,  0,  0,  0,  0,  0)
out <- shared_qtls[9, shared_qtls[9, ] > 0][order(shared_qtls[9, shared_qtls[9, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[9, ]     <- out
ordered_aa_names[9, ] <- out_names

## odc xiiia       ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn Asp Gln Glu His    Lys Ile Leu Phe   Trp Tyr
shared_qtls[10, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  27200, 0,  0,  0,  9750, 0,  0)
out <- shared_qtls[10, shared_qtls[10, ] > 0][order(shared_qtls[10, shared_qtls[10, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[10, ]     <- out
ordered_aa_names[10, ] <- out_names

## odc xii          ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn Asp Gln Glu His    Lys Ile Leu Phe Trp Tyr
shared_qtls[11, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  40650, 0,  0,  0,  0,  0,  0)
out <- shared_qtls[11, shared_qtls[11, ] > 0][order(shared_qtls[11, shared_qtls[11, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- c("His", rep("", 10))
ordered_qtls[11, ]     <- out
ordered_aa_names[11, ] <- out_names

## odc x          ## Ala     Cys     Gly Met Pro Ser Thr Val     Arg Asn Asp Gln    Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[12, ] <- c(55550, 40750,  0,  0,  0,  0,  0, 11500,  0,  0,  0,  33700, 0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[12, shared_qtls[12, ] > 0][order(shared_qtls[12, shared_qtls[12, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[12, ]     <- out
ordered_aa_names[12, ] <- out_names

## odc viib        ## Ala    Cys    Gly   Met Pro    Ser   Thr   Val   Arg Asn     Asp Gln Glu His Lys Ile Leu Phe  Trp    Tyr
shared_qtls[13, ] <- c(13950, 25000, 450,  0,  15150, 3350, 4600, 8100, 0,  74050,  0,  0,  0,  0,  0,  0,  0,  700, 12600, 41900)
out <- shared_qtls[13, shared_qtls[13, ] > 0][order(shared_qtls[13, shared_qtls[13, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[13, ]     <- out
ordered_aa_names[13, ] <- out_names

## odc viia         ## Ala Cys Gly Met Pro Ser Thr    Val Arg Asn   Asp Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[14, ] <- c(0,  0,  0,  0,  0,  0,  91350, 0,  0,  5800, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[14, shared_qtls[14, ] > 0][order(shared_qtls[14, shared_qtls[14, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[14, ]     <- out
ordered_aa_names[14, ] <- out_names

## odc v           ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn Asp Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[15, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[15, shared_qtls[15, ] > 0][order(shared_qtls[15, shared_qtls[15, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[15, ]     <- out
ordered_aa_names[15, ] <- out_names

## odc iva         ## Ala Cys Gly    Met Pro Ser Thr Val Arg Asn Asp Gln   Glu    His    Lys Ile Leu Phe Trp Tyr
shared_qtls[16, ] <- c(0,  0,  48400, 0,  0,  0,  0,  0,  0,  0,  0,  8450, 34150, 12650, 0,  0,  0,  0,  0,  0)
out <- shared_qtls[16, shared_qtls[16, ] > 0][order(shared_qtls[16, shared_qtls[16, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[16, ]     <- out
ordered_aa_names[16, ] <- out_names

## odc iib         ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn Asp Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[17, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[17, shared_qtls[17, ] > 0][order(shared_qtls[17, shared_qtls[17, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[17, ]     <- out
ordered_aa_names[17, ] <- out_names

## odc iia         ## Ala Cys Gly Met Pro Ser Thr Val Arg Asn Asp Gln Glu His Lys Ile Leu Phe Trp Tyr
shared_qtls[18, ] <- c(0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0)
out <- shared_qtls[18, shared_qtls[18, ] > 0][order(shared_qtls[18, shared_qtls[18, ] > 0])]
out <- unlist(c(out, rep(0, 11 - length(out))))
out_names <- names(out)
ordered_qtls[18, ]     <- out
ordered_aa_names[18, ] <- out_names

ac_n_ends <- c("Ala", "Cys", "Gly", "Met", "Pro", "Ser", "Thr", "Val")

## make the heatmap
cols <- c(gray(0.95), magma(n = 13, alpha = 0.8)[2:13])

pdf(file = "~/emacs/ubi_QTL_paper/figures/final_figures/figure_05_QTL_sharing.pdf",
    height = 7, width = 8)
levelplot(t(as.matrix(ordered_qtls)),
          xlab = "",
          ylab = "",
          at = col_range,
          col.regions = cols,
          scales = list(tck = c(1, 0),
                        x = list(cex = 1.2,
                                 labels = as.character(1:11)),
                        y = list(cex = 1.2,
                                 labels = gsub(pattern = ".+_",
                                               replacement = "",
                                               x = row.names(ordered_qtls)))),
          colorkey = list(space = "right",
                          labels = as.character(c(seq(from = 0, to = 100, by = 10))),
                          at = c(seq(from = 0, to = 1e5, by = 1e4)),
                          col = cols[2:length(cols)],
                          tck = 1.2,
                          axis.text = list(cex = 1.1)),
          par.settings = list(par.xlab.text = list(cex = 1.4),
                              clip = list(panel = F),
                              layout.widths = list(left.padding = 3,
                                                   right.padding = 0.5),
                              layout.heights = list(bottom.padding = 3,
                                                    top.padding = 0)),
          panel = function(...) {
              panel.levelplot(...,
                              border = gray(0.7),
                              border.lty = 1,
                              border.lwd = 0.5)
              ## reporter separation lines
              panel.segments(x0 = c(-1.6, -1.6),
                             x1 = c(-1.6, -1.6),
                             y0 = c(1, 8),
                             y1 = c(7, 18),
                             lwd = 1.75,
                             col = gray(0.2))
              ## reporter text
              panel.text(x = c(-3.0, -3.0),
                         y = c(4, 13),
                         labels = c("Rpn4 TFT\nQTLs",
                                    "ODC TFT\nQTLs"))
              ## colorkey title
              panel.text(x = c(14.5, 15.2),
                         y = c(9.5, 9.5),
                         labels = c("Distance between QTL Peaks Detected with",
                                    "N-Degron and Proteasome Activity Reporters (kb)"),
                         srt = 90,
                         cex = 1.2)
              for(i in 1:nrow(t(ordered_aa_names))) {
                  for(j in 1:ncol(t(ordered_aa_names))) {
                      panel.text(x = i,
                                 y = j,
                                 labels = t(ordered_aa_names)[i, j],
                                 col = ifelse(max(t(ordered_aa_names)[i, j] == ac_n_ends) >= 1,
                                               yes = "skyblue",
                                               no = "mintcream"),
                                 font = ifelse(max(t(ordered_aa_names)[i, j] == ac_n_ends) >= 1,
                                               yes = "plain",
                                               no = "plain"))
                  }}
              panel.text(x = c(6, 6),
                         y = c(-1, -1.7),
                         labels = c("N-Degrons with a QTL Overlapping",
                                    "Indicated Proteasome Activity QTL"),
                         cex = 1.2)
              ## N-end rule pathway legend
              panel.rect(xleft = 7.6,
                         xright = 11.4,
                         ytop = 18.4,
                         ybottom = 16.6,
                         col = gray(0.4),
                         border = gray(0.2))
              ## legend text 
              panel.text(x = c(9.45, 9.55),
                         y = c(17.9, 17.1),
                         labels = c("Ac/N-end Rule", "Arg/N-end Rule"),
                         col = c("skyblue", "mintcream"),
                         cex = 1.0)
})
## code end 
dev.off()
