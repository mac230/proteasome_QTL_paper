require("ggseqlogo")
require("ggplot2")

## YAP1 position weight matrix from YeTFaSCo
A <- c(0.833333333333333, 0.833333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.166666666666667)
T <- c(0.0, 0.0, 0.0, 0.166666666666667, 1.0, 1.0, 1.0, 0.833333333333333, 0.833333333333333)
G <- c(0.166666666666667, 0.166666666666667, 0.833333333333333, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
C <- c(0.0, 0.0, 0.166666666666667, 0.833333333333333, 0.0, 0.0, 0.0, 0.166666666666667, 0.0)

## bind it as a matrix
yap1 <- as.matrix(rbind(A, T, G, C))

pdf(file = "~/emacs/proteasome_QTL_paper/results/yap1_binding_motif.pdf",
    width = 12, height = 6)
    ggseqlogo(yap1, method = "bits")
dev.off()
