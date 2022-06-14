library("lattice")
vals <- c(0.3, 0.6, -log2(0.3 / 0.6), 0.8, 0.8, -log2(1))
fac  <- as.factor(1:length(vals))

low_col <- "#ec6ab7ff"
high_col <-  "#5AC168ff"
TFT_col <- gray(0.9)

## 'high_TFT' = high UPS activity, so low TFT ratio
low_col  <- "#ec6ab7"
high_col <- "#88d392"

## plot
TFT_bar <- barchart(vals ~ fac,
                    type = c("g"),
                    xlab = "Degradation Rate",
                    ylab = "Intensity (a.u.)",
                    ylim = c(-0.1, 1.1),
                    xlim = c(0.7, 5.7),
                    ## border col and lwd for the border around the bars
                    border = gray(0.25),
                    lwd = 2.5,
                    ## Ratio of bar width to inter-bar space
                    ## default seems to be > 1
                    box.ratio = 1.2,
                    col = c(low_col, high_col, TFT_col),
                    scales = list(x = list(cex = 3,
                                           at = c(2, 4.4),
                                           labels = c("High", "Low")),
                                  y = list(at = seq(from = 0, to = 1, by = 0.2),
                                           cex = 3),
                                  tck = c(1, 0)),
                    par.settings = list(par.ylab.text = list(cex = 3),
                                        par.xlab.text = list(cex = 3)),
                    key = list(rectangles = list(col = c(low_col,
                                                         high_col,
                                                         TFT_col),
                                                 border = gray(0.25),
                                                 lwd = 10,
                                                 size = 6,
                                                 height = 0.7),
                               text = list(labels = c("RFP",
                                                      "GFP",
                                                      expression("-log"[2]*" RFP / GFP")),
                                           cex = 3),
                               ## corner = c(0.99, 0.99),
                               columns = 3,
                               padding.text = 5,
                               cex.border = 10,
                               between = 1),
                               ## background = "white",
                               ## border = gray(0.4)),
                    panel = function(x, y, ...) {
                        panel.abline(h = seq(from = 0, to = 1, by = 0.2),
                                     v = c(0.8, 2, 3.2, 4.4, 5.6),
                                     col = gray(0.9))
                        ## taken from the 'panel.barchart' help
                        ## specify the exact location of the bars
                        ## x: Extent of Bars. By default, bars start at left of panel,
                        ## unless ‘origin’ is specified, in which case they start there.
                        ## y: Horizontal location of bars. Possibly a factor.
                        ## but, flip them, since 'horizontal = F'
                        panel.barchart(x = c(1.4, 2, 2.6, 3.8, 4.4, 5),
                                       y = vals,
                                       ...)
                    })


pdf(file = paste0("~/emacs/ubi_QTL_paper/results/",
                  gsub(pattern = " .*",
                       replacement = "",
                       x = Sys.time()),
                  "_TFT_schematic_barplot.pdf"
                  ),
    height = 9, width = 9)
print(TFT_bar)
dev.off()
