## load in the summarized flow data
load(paste0("~/emacs/proteasome_QTL_paper/",
            "scripts_data/reporter_characterization/",
            "reporter_characterization_data.RData"))


## -----
## ODC TFT
ODC_out <- reporter_data[reporter_data$reporter == "ODC TFT", ]
str(ODC_out)

## BY vs. RM; p = 0.0001986
t.test(x = ODC_out$TFT_ratio[ODC_out$strain == "BY_strain"],
       y = ODC_out$TFT_ratio[ODC_out$strain == "RM_strain"],
       alternative = "two.sided")

## BY vs. rpn4; p = 1.422e-06
t.test(x = ODC_out$TFT_ratio[ODC_out$strain == "BY_strain"],
       y = ODC_out$TFT_ratio[ODC_out$strain == "rpn4_strain"],
       alternative = "two.sided")

## BY vs. rpn4; p = 1.218e-06
t.test(x = ODC_out$TFT_ratio[ODC_out$strain == "RM_strain"],
       y = ODC_out$TFT_ratio[ODC_out$strain == "rpn4_strain"],
       alternative = "two.sided")


## -----
## rpn4 TFT
rpn4_out <- reporter_data[reporter_data$reporter == "Rpn4 TFT", ]
str(rpn4_out)

## BY vs. RM; p = 1.184e-08
t.test(x = rpn4_out$TFT_ratio[rpn4_out$strain == "BY_strain"],
       y = rpn4_out$TFT_ratio[rpn4_out$strain == "RM_strain"],
       alternative = "two.sided")

## BY vs. rpn4; p = 1.647e-13
t.test(x = rpn4_out$TFT_ratio[rpn4_out$strain == "BY_strain"],
       y = rpn4_out$TFT_ratio[rpn4_out$strain == "rpn4_strain"],
       alternative = "two.sided")

## BY vs. rpn4; p = 7.42e-13
t.test(x = rpn4_out$TFT_ratio[rpn4_out$strain == "RM_strain"],
       y = rpn4_out$TFT_ratio[rpn4_out$strain == "rpn4_strain"],
       alternative = "two.sided")


## -----
## BY per reporter; p = 6.914e-10
by_out <- reporter_data[reporter_data$strain == "BY_strain", ]

t.test(x = by_out$TFT_ratio[by_out$reporter == "ODC TFT"],
       y = by_out$TFT_ratio[by_out$reporter == "Rpn4 TFT"],
       alternative = "two.sided")


## -----
## RM per reporter; p = 5.29e-14
rm_out <- reporter_data[reporter_data$strain == "RM_strain", ]

t.test(x = rm_out$TFT_ratio[rm_out$reporter == "ODC TFT"],
       y = rm_out$TFT_ratio[rm_out$reporter == "Rpn4 TFT"],
       alternative = "two.sided")


## -----
## rpn4 per reporter; p = 0.3604
rpn4_out <- reporter_data[reporter_data$strain == "rpn4_strain", ]

t.test(x = rpn4_out$TFT_ratio[rpn4_out$reporter == "ODC TFT"],
       y = rpn4_out$TFT_ratio[rpn4_out$reporter == "Rpn4 TFT"],
       alternative = "two.sided")


## -----
## strain by reporter interaction; p = 0.0125
out_subset <- reporter_data[reporter_data$strain == "BY_strain" |
                            reporter_data$strain == "RM_strain", ]

int_aov <- aov(TFT_ratio ~ strain * reporter,
               data = out_subset)

