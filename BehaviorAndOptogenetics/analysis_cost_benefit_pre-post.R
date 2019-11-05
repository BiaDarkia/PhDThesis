# This script provides a tool to further analyse data that were previously extracted using 
# custom Python scripts from the Med Associates result files that record choice behavior and 
# reaction times during behavioral training and testing on the benefit-benefit, cost-cost and cost-benefit 
# decision-making paradimg. This program compares changes in choice behavior and reaction times
# that are recorded during testing and analyses differences between choice behavior and reaction times
# pre- and post-surgery. This program runs statistical tests and generates figures to visualize the data.

# Import required libraries
library(plyr)
library(data.table)
library(car)
library(ez)
library(sjstats)
library(RColorBrewer)

# Read in the lever pressing, training and testing data
tidy_data_training <- read.table("./tidy_data_training.csv", header=TRUE, sep = ",")
tidy_data_training_cb <- read.table("./tidy_data_training_cb.csv", header=TRUE, sep = ",")
tidy_data <- read.table("./tidy_data.csv", header=TRUE, sep = ",")
tidy_data_day <- read.table("./tidy_data_day.csv", header=TRUE, sep = ",")
arch <- read.table("./arch.csv", header=TRUE, sep = ",")

# Change data into a format to calculate basic statistics and to run statistical tests
colnames(tidy_data_training)[colnames(tidy_data_training)=="block"] <- "opto_block"
colnames(tidy_data_training_cb)[colnames(tidy_data_training_cb)=="block"] <- "opto_block"
colnames(tidy_data)[colnames(tidy_data)=="opto"] <- "opto_block"
colnames(tidy_data_day)[colnames(tidy_data_day)=="opto"] <- "opto_block"

# Change the format to include archaerhodopsin vs. control as a factor
data_training <- join(arch, tidy_data_training, by="subject", type="inner")
data_training_cb <- join(arch, tidy_data_training_cb, by="subject", type="inner")
data <- join(arch, tidy_data, by="subject", type="inner")
data_day <- join(arch, tidy_data_day, by="subject", type="inner")

# Calculate the number of total trials
data_training$total <- data_training$high + data_training$low + data_training$omit
data_training_cb$total <- data_training_cb$high + data_training_cb$low + data_training_cb$omit

# Calculate the percentage of omitted trials
data_training$omit_pct <- (data_training$omitted/data_training$total)*100
data_training_cb$omit_pct <- (data_training_cb$omitted/data_training_cb$total)*100
data$omit_pct <- (data$omitted/data$total)*100
data_day$omit_pct <- (data_day$omitted/data_day$total)*100

# Change the data into wide format
data_noopto <- data.table(subset(tidy_data, opto_block=='0'))
data_opto <- data.table(subset(tidy_data, opto_block=='1'))

setnames(data_noopto, 5:13, paste0(names(data_noopto)[5:13], '_noopto'))
setnames(data_opto, 5:13, paste0(names(data_opto)[5:13], '_opto'))

data_wide <- join(data_noopto, data_opto[, c("opto_block", "group"):=NULL], 
    by=c("subject", "condition"), type="inner")
data_wide <- join(arch, data_wide, by="subject", type="inner")

# Calculate the changes/difference in choice behavior and reaction time
# without and with delivery of the 590 nm light stimulus to prelimbic
# cortical layer 1
data_wide['diff'] <- data_wide['performance_opto'] -
    data_wide['performance_noopto']
data_wide['diff_rt_diff'] <- data_wide['rt_diff_opto'] -
    data_wide['rt_diff_noopto']

labels_pct_diff <- c("-40","","-20","","0","","20","","40")

# Stripchart and boxplot of the change in the percentage of high benefit, high cost or 
# high benefit-high cost trials out of the total number of trials that were not omitted 
# for archaerhodopsin and control animals for all three decision-making paradigms
tiff("figures/difference-opto.tiff", width=3600, height=6000, res=300)
par(mfrow=c(3,1), mar=c(5.1, 8.6, 10.6, 2.1), mgp=c(6, 1.5, 0))
boxplot(diff~arch, data = data_wide, subset = data_wide$condition == '3', ylim=c(-40,40), 
    axes=FALSE, frame.plot=FALSE, lwd=4, range=1)
stripchart(diff~arch, data = data_wide, subset = data_wide$condition == '3',
    method = "jitter", col=c("blue", "red"), pch=16, cex=3, vertical=TRUE, add=TRUE)
title(ylab="Difference light ON-OFF (%)", cex.lab=2.6, font.lab=2)
axis(side = 1, at = 1:2, labels=c("Control", "Archaerhodopsin"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-40,40,10), labels=labels_pct_diff, las=1, lwd=5, cex.axis=2.2, font=2)
mtext(text="A- Cost-Benefit Decision-Making", side=3, line=3, at=0.95, cex=2.6, font=2)
boxplot(diff~arch, data = data_wide, subset = data_wide$condition == '1', ylim=c(-40,40), 
    axes=FALSE, frame.plot=FALSE, lwd=4, range=1)
stripchart(diff~arch, data = data_wide, subset = data_wide$condition == '1',
    method = "jitter", col=c("blue", "red"), pch=16, cex=3, vertical=TRUE, add=TRUE)
title(ylab="Difference light ON-OFF (%)", cex.lab=2.6, font.lab=2)
axis(side = 1, at = 1:2, labels=c("Control", "Archaerhodopsin"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-40,40,10), labels=labels_pct_diff, las=1, lwd=5, cex.axis=2.2, font=2)
mtext(text="B- Benefit-Benefit Decision-Making", side=3, line=3, at=1, cex=2.6, font=2)
boxplot(diff~arch, data = data_wide, subset = data_wide$condition == '2', ylim=c(-40,40), 
    axes=FALSE, frame.plot=FALSE, lwd=4, range=1)
stripchart(diff~arch, data = data_wide, subset = data_wide$condition == '2',
    method = "jitter", col=c("blue", "red"), pch=16, cex=3, vertical=TRUE, add=TRUE)
title(ylab="Difference light ON-OFF (%)", cex.lab=2.6, font.lab=2)
axis(side = 1, at = 1:2, labels=c("Control", "Archaerhodopsin"),
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-40,40,10), labels=labels_pct_diff, las=1, lwd=5, cex.axis=2.2, font=2)
mtext(text="C- Cost-Cost Decision-Making", side=3, line=3, at=0.88, cex=2.6, font=2)
dev.off()

labels_rt_diff <- c("-4","","-2","","0","","2","","4")

# Stripchart and boxplot of the change in the differences in reaction time on trials, in which animals choose the 
# high benefit, high cost or high benefit-high cost option, as compared to trials, in which animals
# choose the low benefit, low cost or low benefit-low cost option for archaerhodopsin and 
# each control animals for all three decision-making paradigms
tiff("figures/difference-rt-opto.tiff", width=3600, height=6000, res=300)
par(mfrow=c(3,1), mar=c(5.1, 8.6, 10.6, 2.1), mgp=c(6, 1.5, 0))
boxplot(diff_rt_diff~arch, data = data_wide, subset = data_wide$condition == '3', ylim=c(-4,4), 
    axes=FALSE, frame.plot=FALSE, lwd=4, range=1)
stripchart(diff_rt_diff~arch, data = data_wide, subset = data_wide$condition == '3',
    method = "jitter", col=c("blue", "red"), pch=16, cex=3, vertical=TRUE, add=TRUE)
title(ylab="Difference light ON-OFF (s)", cex.lab=2.6, font.lab=2)
axis(side = 1, at = 1:2, labels=c("Control", "Archaerhodopsin"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
mtext(text="A- Cost-Benefit Decision-Making", side=3, line=3, at=0.95, cex=2.6, font=2)
boxplot(diff_rt_diff~arch, data = data_wide, subset = data_wide$condition == '1', ylim=c(-4,4), 
    axes=FALSE, frame.plot=FALSE, lwd=4, range=1)
stripchart(diff_rt_diff~arch, data = data_wide, subset = data_wide$condition == '1',
    method = "jitter", col=c("blue", "red"), pch=16, cex=3, vertical=TRUE, add=TRUE)
title(ylab="Difference light ON-OFF (s)", cex.lab=2.6, font.lab=2)
axis(side = 1, at = 1:2, labels=c("Control", "Archaerhodopsin"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
mtext(text="B- Benefit-Benefit Decision-Making", side=3, line=3, at=1, cex=2.6, font=2)
boxplot(diff_rt_diff~arch, data = data_wide, subset = data_wide$condition == '2', ylim=c(-4,4), 
    axes=FALSE, frame.plot=FALSE, lwd=4, range=1)
stripchart(diff_rt_diff~arch, data = data_wide, subset = data_wide$condition == '2',
    method = "jitter", col=c("blue", "red"), pch=16, cex=3, vertical=TRUE, add=TRUE)
title(ylab="Difference light ON-OFF (s)", cex.lab=2.6, font.lab=2)
axis(side = 1, at = 1:2, labels=c("Control", "Archaerhodopsin"),
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
mtext(text="C- Cost-Cost Decision-Making", side=3, line=3, at=0.88, cex=2.6, font=2)
dev.off()

# Initialize variables to store outputs of statistical tests that compare changes in choice behavior
# and reaction times between archaerhodopsin and control animals
diff_tests_output <- data.frame(condition = numeric(), variable = character(),
    levene = numeric(), ks = numeric(), p_value = numeric(), effect_size = numeric(), stringsAsFactors=FALSE)

levenes_diff = list()
ks_diff = list()
ttests_diff = list()
mw_diff = list()

levenes_diff_rt_diff = list()
ks_diff_rt_diff = list()
ttests_diff_rt_diff = list()
mw_diff_rt_diff = list()

# Calculate Levene's test, Kolmogorov-Smirnov's test, t-tests and Mann-Whitney U tests to statistically
# compare changes in the choice behavior and in reaction times between archaerhodopsin and control animals
for (c in c(1, 2, 3)) {
    data_cond <- subset(data_wide, condition==c)
  
    levenes_diff = append(levenes_diff, list(leveneTest(diff~as.factor(arch), data_cond)))
    ks_diff = append(ks_diff, list(ks.test(data_cond$diff[data_cond$arch=='0'], data_cond$diff[data_cond$arch=='1'])))
    ttests_diff = append(ttests_diff, list(t.test(data_cond$diff[data_cond$arch=='1'], 
        data_cond$diff[data_cond$arch=='0'], paired=FALSE)))
    mw_diff = append(mw_diff, list(wilcox.test(data_cond$diff[data_cond$arch=='1'], 
        data_cond$diff[data_cond$arch=='0'], paired=FALSE)))
  
    levenes_diff_rt_diff = append(levenes_diff_rt_diff, list(leveneTest(diff_rt_diff~as.factor(arch), data_cond)))
    ks_diff_rt_diff = append(ks_diff_rt_diff, list(ks.test(data_cond$diff_rt_diff[data_cond$arch=='0'], data_cond$diff_rt_diff[data_cond$arch=='1'])))
    ttests_diff_rt_diff = append(ttests_diff_rt_diff, list(t.test(data_cond$diff_rt_diff[data_cond$arch=='1'], 
        data_cond$diff_rt_diff[data_cond$arch=='0'], paired=FALSE)))
    mw_diff_rt_diff = append(mw_diff_rt_diff, list(wilcox.test(data_cond$diff_rt_diff[data_cond$arch=='1'], 
        data_cond$diff_rt_diff[data_cond$arch=='0'], paired=FALSE)))
  
    p = mw_diff[[length(mw_diff)]]$p.value
    eff = abs(qnorm(p/2)/sqrt(nrow(data_cond)))
  
    diff_tests_output[nrow(diff_tests_output) + 1,] = list(c, 'p',
        levenes_diff[[length(levenes_diff)]]$`Pr(>F)`, ks_diff[[length(ks_diff)]]$p.value, p, eff)
  
    p_diff_rt = mw_diff_rt_diff[[length(mw_diff_rt_diff)]]$p.value
    eff_diff_rt = abs(qnorm(p_diff_rt/2)/sqrt(nrow(data_cond)))
  
    diff_tests_output[nrow(diff_tests_output) + 1,] = list(c, 'diff_rt', 
        levenes_diff_rt_diff[[length(levenes_diff_rt_diff)]]$`Pr(>F)`, ks_diff_rt_diff[[length(ks_diff_rt_diff)]]$p.value, p_diff_rt, eff_diff_rt)
}

# Initialize arrays to store subject IDs of any "outliers" (1st/3rd quartile +/- 0.5*IQR)
outliers <- c()
outlier_arch <- c()
outlier_condition <- c()

# Detect "outliers"
for (a in c(0, 1)) {
  for (c in c(1, 2, 3)) {
    data_wide_out <- subset(data_wide, arch == a & condition == c)
    out_lower <- quantile(data_wide_out$diff, 0.25) - 0.5*IQR(data_wide_out$diff)
    out_upper <- quantile(data_wide_out$diff, 0.75) + 0.5*IQR(data_wide_out$diff)
    outlier <- subset(data_wide_out, diff <= out_lower | diff >= out_upper)
    for (i in outlier$subject) {
        outliers <- c(outliers, as.character(i))
        outlier_arch <- c(outlier_arch, a)
        outlier_condition <- c(outlier_condition, c)
    }
  }
}

# Store outliers in a data frame
data_outliers <- data.frame(outliers, outlier_arch, outlier_condition)
setnames(data_outliers, c("outliers", "outlier_arch"), c("subject", "arch"))

# Split control and archaerhodopsin animal outliers
data_outliers_sham <- subset(data_outliers, arch==0)
data_outliers_arch <- subset(data_outliers, arch==1)

# Retrieve a list of subject IDs
data_outliers_sham_subject <- unique(data_outliers_sham$subject)
data_outliers_arch_subject <- unique(data_outliers_arch$subject)

# Assign a color to each subject
data_outliers_sham_color <- brewer.pal(n = 9, name = "Blues")[(10-length(data_outliers_sham_subject)):9]
data_outliers_arch_color <- brewer.pal(n = 9, name = "Reds")[(10-length(data_outliers_arch_subject)):9]

data_outliers_color <- rbindlist(list(data.frame(data_outliers_sham_subject, data_outliers_sham_color),
    data.frame(data_outliers_arch_subject, data_outliers_arch_color)), use.names=FALSE)
setnames(data_outliers_color, c("data_outliers_sham_subject", "data_outliers_sham_color"), c("subject", "color"))

data_outliers <- join(data_outliers, data_outliers_color, by="subject", type="inner")

# Initiate a new column to combine information about day and block
data_day$day_block <- data_day$day

data_day$day_block[data_day$condition==1 & data_day$day==1 & data_day$opto_block==1] <- 2
data_day$day_block[data_day$condition==1 & data_day$day==2 & data_day$opto_block==0] <- 3
data_day$day_block[data_day$condition==1 & data_day$day==2 & data_day$opto_block==1] <- 4
data_day$day_block[data_day$condition==1 & data_day$day==3 & data_day$opto_block==0] <- 5
data_day$day_block[data_day$condition==1 & data_day$day==3 & data_day$opto_block==1] <- 6
data_day$day_block[data_day$condition==2 & data_day$day==1 & data_day$opto_block==0] <- 7
data_day$day_block[data_day$condition==2 & data_day$day==1 & data_day$opto_block==1] <- 8
data_day$day_block[data_day$condition==2 & data_day$day==2 & data_day$opto_block==0] <- 9
data_day$day_block[data_day$condition==2 & data_day$day==2 & data_day$opto_block==1] <- 10
data_day$day_block[data_day$condition==2 & data_day$day==3 & data_day$opto_block==0] <- 11
data_day$day_block[data_day$condition==2 & data_day$day==3 & data_day$opto_block==1] <- 12
data_day$day_block[data_day$condition==3 & data_day$day==1 & data_day$opto_block==0] <- 13
data_day$day_block[data_day$condition==3 & data_day$day==1 & data_day$opto_block==1] <- 14
data_day$day_block[data_day$condition==3 & data_day$day==2 & data_day$opto_block==0] <- 15
data_day$day_block[data_day$condition==3 & data_day$day==2 & data_day$opto_block==1] <- 16
data_day$day_block[data_day$condition==3 & data_day$day==3 & data_day$opto_block==0] <- 17
data_day$day_block[data_day$condition==3 & data_day$day==3 & data_day$opto_block==1] <- 18

data_day_outliers <- join(data_day, data_outliers, by=c("subject", "arch") , type="inner")
data_day_outliers <- transform(data_day_outliers, subject = as.character(subject), color = as.character(color))

# Split data between control and archaerhodopsin animals
data_day_sham <- subset(data_day, arch==0)
data_day_arch <- subset(data_day, arch==1)

# Split data between outliers in each decision-making paradigm
data_day_outliers_bb_sham <- subset(data_day_outliers, outlier_condition==1 & arch==0)
data_day_outliers_bb_arch <- subset(data_day_outliers, outlier_condition==1 & arch==1)
data_day_outliers_cc_sham <- subset(data_day_outliers, outlier_condition==2 & arch==0)
data_day_outliers_cc_arch <- subset(data_day_outliers, outlier_condition==2 & arch==1)
data_day_outliers_cb_sham <- subset(data_day_outliers, outlier_condition==3 & arch==0)
data_day_outliers_cb_arch <- subset(data_day_outliers, outlier_condition==3 & arch==1)

labels_day_blocks <- rep(c("1,OFF", "1,ON", "2,OFF", "2,ON", "3,OFF", "3,ON"), times=3)

# Superimpose choice behavior of control animal outliers in the benefit-benefit decision-making paradigm
# in shades of blue on top of the choice behavior of all control animals plotted in gray
tiff("figures/outliers-bb-sham.tiff", width=6000, height=3000, res=300)
par(mar=c(9.1, 12.1, 6.1, 0.1), mgp=c(6.5, 1.5, 0))
interaction.plot(data_day_sham$day_block, data_day_sham$subject, 
    data_day_sham$performance, ylim=c(0,100), axes=FALSE, 
    xlab="Benefit-Benefit                       Cost-Cost                       Cost-Benefit\nDay, Light OFF/ON", 
    ylab="High Benefit, High Cost or\nHigh Benefit-Cost Choices (%)", cex.lab=2.4, font.lab=2, col="gray70",
    type="b", lty=1, lwd=6, pch=16, cex=3, frame.plot=FALSE, legend=FALSE)
interaction.plot(data_day_outliers_bb_sham$day_block, data_day_outliers_bb_sham$subject, 
    data_day_outliers_bb_sham$performance, col=c(unique(data_day_outliers_bb_sham$color)),
    type="b", lty=1, lwd=6, pch=16, cex=3, axes=FALSE, legend=FALSE, add=TRUE)
title(main="A-Control Animal Outliers in the Benefit-Benefit Decision-Making Paradigm", cex.main=2.6, font.main=2, adj=0.3)
axis(side=1, at=1:18, labels=labels_day_blocks, las=1, lwd=5, cex.axis=1.6, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.0, font=2)
dev.off()

# Superimpose choice behavior of archaerhodopsin animal outliers in the benefit-benefit decision-making paradigm
# in shades of blue on top of the choice behavior of all archaerhodopsin animals plotted in gray
tiff("figures/outliers-bb-arch.tiff", width=6000, height=3000, res=300)
par(mar=c(9.1, 12.1, 6.1, 0.1), mgp=c(6.5, 1.5, 0))
interaction.plot(data_day_arch$day_block, data_day_arch$subject, 
    data_day_arch$performance, ylim=c(0,100), axes=FALSE, 
    xlab="Benefit-Benefit                       Cost-Cost                       Cost-Benefit\nDay, Light OFF/ON", 
    ylab="High Benefit, High Cost or\nHigh Benefit-Cost Choices (%)", cex.lab=2.4, font.lab=2, col="gray70",
    type="b", lty=1, lwd=6, pch=16, cex=3, frame.plot=FALSE, legend=FALSE)
interaction.plot(data_day_outliers_bb_arch$day_block, data_day_outliers_bb_arch$subject, 
    data_day_outliers_bb_arch$performance, col=c(unique(data_day_outliers_bb_arch$color)),
    type="b", lty=1, lwd=6, pch=16, cex=3, axes=FALSE, legend=FALSE, add=TRUE)
title(main="B-Archaerhodopsin Animal Outliers in the Benefit-Benefit Decision-Making Paradigm", cex.main=2.6, font.main=2, adj=0.3)
axis(side=1, at=1:18, labels=labels_day_blocks, las=1, lwd=5, cex.axis=1.6, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.0, font=2)
dev.off()

# Superimpose choice behavior of control animal outliers in the cost-cost decision-making paradigm
# in shades of blue on top of the choice behavior of all control animals plotted in gray
tiff("figures/outliers-cc-sham.tiff", width=6000, height=3000, res=300)
par(mar=c(9.1, 12.1, 6.1, 0.1), mgp=c(6.5, 1.5, 0))
interaction.plot(data_day_sham$day_block, data_day_sham$subject, 
    data_day_sham$performance, ylim=c(0,100), axes=FALSE, 
    xlab="Benefit-Benefit                       Cost-Cost                       Cost-Benefit\nDay, Light OFF/ON", 
    ylab="High Benefit, High Cost or\nHigh Benefit-Cost Choices (%)", cex.lab=2.4, font.lab=2, col="gray70",
    type="b", lty=1, lwd=6, pch=16, cex=3, frame.plot=FALSE, legend=FALSE)
interaction.plot(data_day_outliers_cc_sham$day_block, data_day_outliers_cc_sham$subject, 
    data_day_outliers_cc_sham$performance, col=c(unique(data_day_outliers_cc_sham$color)),
    type="b", lty=1, lwd=6, pch=16, cex=3, axes=FALSE, legend=FALSE, add=TRUE)
title(main="A-Control Animal Outliers in the Cost-Cost Decision-Making Paradigm", cex.main=2.6, font.main=2, adj=0.3)
axis(side=1, at=1:18, labels=labels_day_blocks, las=1, lwd=5, cex.axis=1.6, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.0, font=2)
dev.off()

# Superimpose choice behavior of archaerhodopsin animal outliers in the cost-cost decision-making paradigm
# in shades of blue on top of the choice behavior of all archaerhodopsin animals plotted in gray
tiff("figures/outliers-cc-arch.tiff", width=6000, height=3000, res=300)
par(mar=c(9.1, 12.1, 6.1, 0.1), mgp=c(6.5, 1.5, 0))
interaction.plot(data_day_arch$day_block, data_day_arch$subject, 
    data_day_arch$performance, ylim=c(0,100), axes=FALSE, 
    xlab="Benefit-Benefit                       Cost-Cost                       Cost-Benefit\nDay, Light OFF/ON", 
    ylab="High Benefit, High Cost or\nHigh Benefit-Cost Choices (%)", cex.lab=2.4, font.lab=2, col="gray70",
    type="b", lty=1, lwd=6, pch=16, cex=3, frame.plot=FALSE, legend=FALSE)
interaction.plot(data_day_outliers_cc_arch$day_block, data_day_outliers_cc_arch$subject, 
    data_day_outliers_cc_arch$performance, col=c(unique(data_day_outliers_cc_arch$color)),
    type="b", lty=1, lwd=6, pch=16, cex=3, axes=FALSE, legend=FALSE, add=TRUE)
title(main="B-Archaerhodopsin Animal Outliers in the Cost-Cost Decision-Making Paradigm", cex.main=2.6, font.main=2, adj=0.3)
axis(side=1, at=1:18, labels=labels_day_blocks, las=1, lwd=5, cex.axis=1.6, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.0, font=2)
dev.off()

# Superimpose choice behavior of control animal outliers in the cost-benefit decision-making paradigm
# in shades of blue on top of the choice behavior of all control animals plotted in gray
tiff("figures/outliers-cb-sham.tiff", width=6000, height=3000, res=300)
par(mar=c(9.1, 12.1, 6.1, 0.1), mgp=c(6.5, 1.5, 0))
interaction.plot(data_day_sham$day_block, data_day_sham$subject, 
    data_day_sham$performance, ylim=c(0,100), axes=FALSE, 
    xlab="Benefit-Benefit                       Cost-Cost                       Cost-Benefit\nDay, Light OFF/ON", 
    ylab="High Benefit, High Cost or\nHigh Benefit-Cost Choices (%)", cex.lab=2.4, font.lab=2, col="gray70",
    type="b", lty=1, lwd=6, pch=16, cex=3, frame.plot=FALSE, legend=FALSE)
interaction.plot(data_day_outliers_cb_sham$day_block, data_day_outliers_cb_sham$subject, 
    data_day_outliers_cb_sham$performance, col=c(unique(data_day_outliers_cb_sham$color)),
    type="b", lty=1, lwd=6, pch=16, cex=3, axes=FALSE, legend=FALSE, add=TRUE)
title(main="A-Control Animal Outliers in the Cost-Benefit Decision-Making Paradigm", cex.main=2.6, font.main=2, adj=0.3)
axis(side=1, at=1:18, labels=labels_day_blocks, las=1, lwd=5, cex.axis=1.6, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.0, font=2)
dev.off()

# Superimpose choice behavior of archaerhodopsin animal outliers in the cost-benefit decision-making paradigm
# in shades of blue on top of the choice behavior of all archaerhodopsin animals plotted in gray
tiff("figures/outliers-cb-arch.tiff", width=6000, height=3000, res=300)
par(mar=c(9.1, 12.1, 6.1, 0.1), mgp=c(6.5, 1.5, 0))
interaction.plot(data_day_arch$day_block, data_day_arch$subject, 
    data_day_arch$performance, ylim=c(0,100), axes=FALSE, 
    xlab="Benefit-Benefit                       Cost-Cost                       Cost-Benefit\nDay, Light OFF/ON", 
    ylab="High Benefit, High Cost or\nHigh Benefit-Cost Choices (%)", cex.lab=2.4, font.lab=2, col="gray70",
    type="b", lty=1, lwd=6, pch=16, cex=3, frame.plot=FALSE, legend=FALSE)
interaction.plot(data_day_outliers_cb_arch$day_block, data_day_outliers_cb_arch$subject, 
    data_day_outliers_cb_arch$performance, col=c(unique(data_day_outliers_cb_arch$color)),
    type="b", lty=1, lwd=6, pch=16, cex=3, axes=FALSE, legend=FALSE, add=TRUE)
title(main="B-Archaerhodopsin Animal Outliers in the Cost-Benefit Decision-Making Paradigm", cex.main=2.6, font.main=2, adj=0.3)
axis(side=1, at=1:18, labels=labels_day_blocks, las=1, lwd=5, cex.axis=1.6, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.0, font=2)
dev.off()

setnames(data_outliers, c("outlier_condition"), c("condition"))

data_training_outliers <- join(data_training, data_outliers, by=c("subject", "arch", "condition"), type="inner")
data_testing_outliers <- join(data_day, data_outliers, by=c("subject", "arch", "condition"), type="inner")

data_training_lastday_bb_sham <- subset(data_training, day==-1 & opto_block==0 & arch==0 & condition==1)
data_training_lastday_bb_arch <- subset(data_training, day==-1 & opto_block==0 & arch==1 & condition==1)
data_training_outliers_lastday_sham <- subset(data_training_outliers, day==-1 & opto_block==0 & arch==0 & condition==1)
data_training_outliers_lastday_arch <- subset(data_training_outliers, day==-1 & opto_block==0 & arch==1 & condition==1)
data_testing_lastday_bb_sham <- subset(data_day, opto_block==0 & arch==0 & condition==1 & day==3)
data_testing_lastday_bb_arch <- subset(data_day, opto_block==0 & arch==1 & condition==1 & day==3)
data_testing_outliers_lastday_sham <- subset(data_testing_outliers, opto_block==0 & arch==0 & condition==1 & day==3)
data_testing_outliers_lastday_arch <- subset(data_testing_outliers, opto_block==0 & arch==1 & condition==1 & day==3)

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_bb_sham$performance, data_testing_lastday_bb_sham$performance)
#r_arch <- cor(data_training_lastday_bb_arch$performance, data_testing_lastday_bb_arch$performance)

# Plot pre-surgery versus post-surgery choice behavior of control and archaerhodopsin animals on the
# benefit-benefit decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-bb-performance.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0), las=1)
plot(data_training_lastday_bb_sham$performance, data_testing_lastday_bb_sham$performance, xlim=c(0,100), ylim=c(0,100), 
     axes=FALSE, xlab="Pre-Surgery High Benefit Choices (%)", ylab="Post-Surgery High Benefit Choices (%)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_bb_arch$performance, data_testing_lastday_bb_arch$performance, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$performance, data_testing_outliers_lastday_sham$performance, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$performance, data_testing_outliers_lastday_arch$performance, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_bb_sham$performance ~ data_training_lastday_bb_sham$performance), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_bb_arch$performance ~ data_training_lastday_bb_arch$performance), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
#text(90, 8, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(90, 4, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_bb_sham$rt, data_testing_lastday_bb_sham$rt)
#r_arch <- cor(data_training_lastday_bb_arch$rt, data_testing_lastday_bb_arch$rt)

# Plot pre-surgery versus post-surgery reaction time of control and archaerhodopsin animals on the
# benefit-benefit decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-bb-rt.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_training_lastday_bb_sham$rt, data_testing_lastday_bb_sham$rt, xlim=c(0,10), ylim=c(0,10), 
     axes=FALSE, xlab="Pre-Surgery Reaction Time (s)", ylab="Post-Surgery Reaction Time (s)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_bb_arch$rt, data_testing_lastday_bb_arch$rt, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$rt, data_testing_outliers_lastday_sham$rt, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$rt, data_testing_outliers_lastday_arch$rt, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_bb_sham$rt ~ data_training_lastday_bb_sham$rt), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_bb_arch$rt ~ data_training_lastday_bb_arch$rt), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,10,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,10,1), las=1, lwd=5, cex.axis=2.2, font=2)
#text(8, 0.6, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(8, 0.2, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_bb_sham$omit_pct, data_testing_lastday_bb_sham$omit_pct)
#r_arch <- cor(data_training_lastday_bb_arch$omit_pct, data_testing_lastday_bb_arch$omit_pct)

# Plot pre-surgery versus post-surgery percentage of omitted trials of control and archaerhodopsin animals on the
# benefit-benefit decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-bb-omit.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_training_lastday_bb_sham$omit_pct, data_testing_lastday_bb_sham$omit_pct, xlim=c(0,100), ylim=c(0,100), 
     axes=FALSE, xlab="Pre-Surgery Reaction Time (s)", ylab="Post-Surgery Reaction Time (s)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_bb_arch$omit_pct, data_testing_lastday_bb_arch$omit_pct, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$omit_pct, data_testing_outliers_lastday_sham$omit_pct, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$omit_pct, data_testing_outliers_lastday_arch$omit_pct, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_bb_sham$omit_pct ~ data_training_lastday_bb_sham$omit_pct), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_bb_arch$omit_pct ~ data_training_lastday_bb_arch$omit_pct), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
#text(90, 8, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(90, 4, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

data_training_lastday_cc_sham <- subset(data_training, day==-1 & opto_block==0 & arch==0 & condition==2)
data_training_lastday_cc_arch <- subset(data_training, day==-1 & opto_block==0 & arch==1 & condition==2)
data_training_outliers_lastday_sham <- subset(data_training_outliers, day==-1 & opto_block==0 & arch==0 & condition==2)
data_training_outliers_lastday_arch <- subset(data_training_outliers, day==-1 & opto_block==0 & arch==1 & condition==2)
data_testing_lastday_cc_sham <- subset(data_day, opto_block==0 & arch==0 & condition==2 & day==3)
data_testing_lastday_cc_arch <- subset(data_day, opto_block==0 & arch==1 & condition==2 & day==3)
data_testing_outliers_lastday_sham <- subset(data_testing_outliers, opto_block==0 & arch==0 & condition==2 & day==3)
data_testing_outliers_lastday_arch <- subset(data_testing_outliers, opto_block==0 & arch==1 & condition==2 & day==3)

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_cc_sham$performance, data_testing_lastday_cc_sham$performance)
#r_arch <- cor(data_training_lastday_cc_arch$performance, data_testing_lastday_cc_arch$performance)

# Plot pre-surgery versus post-surgery choice behavior of control and archaerhodopsin animals on the
# cost-cost decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cc-performance.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_training_lastday_cc_sham$performance, data_testing_lastday_cc_sham$performance, xlim=c(0,100), ylim=c(0,100), 
     axes=FALSE, xlab="Pre-Surgery High Cost Choices (%)", ylab="Post-Surgery High Cost Choices (%)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_cc_arch$performance, data_testing_lastday_cc_arch$performance, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$performance, data_testing_outliers_lastday_sham$performance, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$performance, data_testing_outliers_lastday_arch$performance, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_cc_sham$performance ~ data_training_lastday_cc_sham$performance), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_cc_arch$performance ~ data_training_lastday_cc_arch$performance), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
#text(90, 8, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(90, 4, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_cc_sham$rt, data_testing_lastday_cc_sham$rt)
#r_arch <- cor(data_training_lastday_cc_arch$rt, data_testing_lastday_cc_arch$rt)

# Plot pre-surgery versus post-surgery reaction time of control and archaerhodopsin animals on the
# cost-cost decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cc-rt.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_training_lastday_cc_sham$rt, data_testing_lastday_cc_sham$rt, xlim=c(0,10), ylim=c(0,10), 
     axes=FALSE, xlab="Pre-Surgery Reaction Time (s)", ylab="Post-Surgery Reaction Time (s)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_cc_arch$rt, data_testing_lastday_cc_arch$rt, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$rt, data_testing_outliers_lastday_sham$rt, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$rt, data_testing_outliers_lastday_arch$rt, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_cc_sham$rt ~ data_training_lastday_cc_sham$rt), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_cc_arch$rt ~ data_training_lastday_cc_arch$rt), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,10,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,10,1), las=1, lwd=5, cex.axis=2.2, font=2)
#text(8, 0.6, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(8, 0.2, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_cc_sham$omit_pct, data_testing_lastday_cc_sham$omit_pct)
#r_arch <- cor(data_training_lastday_cc_arch$omit_pct, data_testing_lastday_cc_arch$omit_pct)

# Plot pre-surgery versus post-surgery percentage of omitted trials of control and archaerhodopsin animals on the
# cost-cost decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cc-omit.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_training_lastday_cc_sham$omit_pct, data_testing_lastday_cc_sham$omit_pct, xlim=c(0,100), ylim=c(0,100), 
     axes=FALSE, xlab="Pre-Surgery High Cost Choices (%)", ylab="Post-Surgery High Cost Choices (%)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_cc_arch$omit_pct, data_testing_lastday_cc_arch$omit_pct, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$omit_pct, data_testing_outliers_lastday_sham$omit_pct, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$omit_pct, data_testing_outliers_lastday_arch$omit_pct, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_cc_sham$omit_pct ~ data_training_lastday_cc_sham$omit_pct), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_cc_arch$omit_pct ~ data_training_lastday_cc_arch$omit_pct), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
#text(90, 8, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(90, 4, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

data_training_cb_lastdays <- as.data.table(subset(data_training_cb, day>=4 & opto_block==0))
data_training_cb_lastday <- data_training_cb_lastdays[data_training_cb_lastdays[, .I[day == max(day)], by=subject]$V1]
data_testing_cb_lastday <- subset(data_day, subject %in% data_training_cb_lastday$subject & opto_block==0 & condition==3 & day==3)

data_training_cb_lastday_outliers <- join(data_training_cb_lastday, data_outliers, by=c("subject", "arch", "condition"), type="inner")
data_testing_cb_lastday_outliers <- join(data_testing_cb_lastday, data_outliers, by=c("subject", "arch", "condition"), type="inner")

data_training_lastday_cb_sham <- subset(data_training_cb_lastday, arch==0)
data_training_lastday_cb_arch <- subset(data_training_cb_lastday, arch==1)
data_training_outliers_lastday_sham <- subset(data_training_cb_lastday_outliers, arch==0)
data_training_outliers_lastday_arch <- subset(data_training_cb_lastday_outliers, arch==1)
data_testing_lastday_cb_sham <- subset(data_testing_cb_lastday, arch==0)
data_testing_lastday_cb_arch <- subset(data_testing_cb_lastday, arch==1)
data_testing_outliers_lastday_sham <- subset(data_testing_cb_lastday_outliers, arch==0)
data_testing_outliers_lastday_arch <- subset(data_testing_cb_lastday_outliers, arch==1)

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_cb_sham$performance, data_testing_lastday_cb_sham$performance)
#r_arch <- cor(data_training_lastday_cb_arch$performance, data_testing_lastday_cb_arch$performance)

# Plot pre-surgery versus post-surgery choice behavior of control and archaerhodopsin animals on the
# cost-benefit decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cb-performance.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_training_lastday_cb_sham$performance, data_testing_lastday_cb_sham$performance, xlim=c(0,100), ylim=c(0,100), 
     axes=FALSE, xlab="Pre-Surgery High Benefit-Cost Choices (%)", ylab="Post-Surgery High Benefit-Cost Choices (%)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_cb_arch$performance, data_testing_lastday_cb_arch$performance, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$performance, data_testing_outliers_lastday_sham$performance, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$performance, data_testing_outliers_lastday_arch$performance, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_cb_sham$performance ~ data_training_lastday_cb_sham$performance), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_cb_arch$performance ~ data_training_lastday_cb_arch$performance), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
#text(90, 8, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(90, 4, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_cb_sham$rt, data_testing_lastday_cb_sham$rt)
#r_arch <- cor(data_training_lastday_cb_arch$rt, data_testing_lastday_cb_arch$rt)

# Plot pre-surgery versus post-surgery reaction time of control and archaerhodopsin animals on the
# cost-benefit decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cb-rt.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_training_lastday_cb_sham$rt, data_testing_lastday_cb_sham$rt, xlim=c(0,10), ylim=c(0,10), 
     axes=FALSE, xlab="Pre-Surgery Reaction Time (s)", ylab="Post-Surgery Reaction Time (s)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_cb_arch$rt, data_testing_lastday_cb_arch$rt, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$rt, data_testing_outliers_lastday_sham$rt, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$rt, data_testing_outliers_lastday_arch$rt, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_cb_sham$rt ~ data_training_lastday_cb_sham$rt), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_cb_arch$rt ~ data_training_lastday_cb_arch$rt), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,10,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,10,1), las=1, lwd=5, cex.axis=2.2, font=2)
#text(8, 0.6, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(8, 0.2, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Calculate Pearson r (optional)
#r_sham <- cor(data_training_lastday_cb_sham$omit_pct, data_testing_lastday_cb_sham$omit_pct)
#r_arch <- cor(data_training_lastday_cb_arch$omit_pct, data_testing_lastday_cb_arch$omit_pct)

# Plot pre-surgery versus post-surgery percentage of omitted trials of control and archaerhodopsin animals on the
# cost-cost decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cb-omit.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_training_lastday_cb_sham$omit_pct, data_testing_lastday_cb_sham$omit_pct, xlim=c(0,100), ylim=c(0,100), 
     axes=FALSE, xlab="Pre-Surgery High Cost-Benefit Choices (%)", ylab="Post-Surgery High Cost-Benefit Choices (%)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_training_lastday_cb_arch$omit_pct, data_testing_lastday_cb_arch$omit_pct, col="red", pch=19, cex=3)
points(data_training_outliers_lastday_sham$omit_pct, data_testing_outliers_lastday_sham$omit_pct, col="blue4", pch=19, cex=3)
points(data_training_outliers_lastday_arch$omit_pct, data_testing_outliers_lastday_arch$omit_pct, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_lastday_cb_sham$omit_pct ~ data_training_lastday_cb_sham$omit_pct), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_lastday_cb_arch$omit_pct ~ data_training_lastday_cb_arch$omit_pct), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
#text(90, 8, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(90, 4, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

data_testing_sham <- subset(data_day, opto_block==0 & arch==0 & condition==3 & day==3)
data_testing_arch <- subset(data_day, opto_block==0 & arch==1 & condition==3 & day==3)
data_testing_outliers_sham <- subset(data_testing_outliers, opto_block==0 & arch==0 & condition==3 & day==3)
data_testing_outliers_arch <- subset(data_testing_outliers, opto_block==0 & arch==1 & condition==3 & day==3)

# Calculate Pearson r (optional)
#r_sham <- cor(data_testing_sham$conc, data_testing_sham$performance)
#r_arch <- cor(data_testing_arch$conc, data_testing_arch$performance)

# Plot estimated dilution of sweetened condensed milk for control and archaerhodopsin animals versus post-surgery
# choice behavior on the cost-benefit decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cb-concentration-performance.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_testing_sham$conc, data_testing_sham$performance, xlim=c(0,20), ylim=c(0,100), 
     axes=FALSE, xlab="Dilution of Sweetened Condensed Milk (%)", ylab="Post-Surgery High Benefit-Cost Choices (%)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_testing_arch$conc, data_testing_arch$performance, col="red", pch=19, cex=3)
points(data_testing_outliers_sham$conc, data_testing_outliers_sham$performance, col="blue4", pch=19, cex=3)
points(data_testing_outliers_arch$conc, data_testing_outliers_arch$performance, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_sham$performance ~ data_testing_sham$conc), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_arch$performance ~ data_testing_arch$conc), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,20,2), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
#text(18, 8, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(18, 4, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Calculate Pearson r (optional)
#r_sham <- cor(data_testing_sham$conc, data_testing_sham$rt)
#r_arch <- cor(data_testing_arch$conc, data_testing_arch$rt)

# Plot estimated dilution of sweetened condensed milk for control and archaerhodopsin animals versus post-surgery
# reaction time on the cost-benefit decision-making paradigm and superimpose "outliers" in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cb-concentration-rt.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_testing_sham$conc, data_testing_sham$rt, xlim=c(0,20), ylim=c(0,10), 
     axes=FALSE, xlab="Dilution of Sweetened Condensed Milk (%)", ylab="Post-Surgery Reaction Time (s)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_testing_arch$conc, data_testing_arch$rt, col="red", pch=19, cex=3)
points(data_testing_outliers_sham$conc, data_testing_outliers_sham$rt, col="blue4", pch=19, cex=3)
points(data_testing_outliers_arch$conc, data_testing_outliers_arch$rt, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_sham$rt ~ data_testing_sham$conc), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_arch$rt ~ data_testing_arch$conc), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,20,2), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,10,1), las=1, lwd=5, cex.axis=2.2, font=2)
#text(18, 0.6, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(18, 0.2, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Calculate Pearson r (optional)
#r_sham <- cor(data_testing_sham$conc, data_testing_sham$omit_pct)
#r_arch <- cor(data_testing_arch$conc, data_testing_arch$omit_pct)

# Plot estimated dilution of sweetened condensed milk for control and archaerhodopsin animals versus post-surgery
# percentage of omitted trials on the cost-benefit decision-making paradigm and superimpose "outliers" 
# in darker shades of blue/red.
# Pearson r and a linear fit may also be plotted (see code in comments above and below)
tiff("figures/pre-post-cb-concentration-omit.tiff", width=3600, height=3600, res=300)
par(mar=c(9.1, 12.1, 4.1, 0.1), mgp=c(5, 1.5, 0))
plot(data_testing_sham$conc, data_testing_sham$omit_pct, xlim=c(0,20), ylim=c(0,100), 
     axes=FALSE, xlab="Dilution of Sweetened Condensed Milk (%)", ylab="Post-Surgery High Cost-Benefit Choices (%)", 
     cex.lab=2.6, font.lab=2, col="blue", lty=1, lwd=6, pch=19, cex=2.4, frame.plot=FALSE)
points(data_testing_arch$conc, data_testing_arch$omit_pct, col="red", pch=19, cex=3)
points(data_testing_outliers_sham$conc, data_testing_outliers_sham$omit_pct, col="blue4", pch=19, cex=3)
points(data_testing_outliers_arch$conc, data_testing_outliers_arch$omit_pct, col="violetred3", pch=19, cex=3)
#abline(lm(data_testing_sham$omit_pct ~ data_testing_sham$conc), col="blue", lty=1, lwd=6)
#abline(lm(data_testing_arch$omit_pct ~ data_testing_arch$conc), col="red", lty=1, lwd=6)
axis(side=1, at=seq(0,20,2), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
#text(18, 8, paste0("r = ", round(r_sham, digits = 2)), col="blue", cex=1.6, font=2)
#text(18, 4, paste0("r = ", round(r_arch, digits = 2)), col="red", cex=1.6, font=2)
dev.off()

# Combine pre- and post-surgery data for each decision-making paradigm in one data frame
data_pre <- rbind(data_training_lastday_bb_sham, data_training_lastday_bb_arch,
    data_training_lastday_cc_sham, data_training_lastday_cc_arch,
    data_training_lastday_cb_sham, data_training_lastday_cb_arch)
data_pre$pre_post <- replicate(nrow(data_pre), 0)

data_post <- rbind(data_testing_lastday_bb_sham, data_testing_lastday_bb_arch,
    data_testing_lastday_cc_sham, data_testing_lastday_cc_arch,
    data_testing_lastday_cb_sham, data_testing_lastday_cb_arch)
data_post <- subset(data_post, select = -c(day_block))
data_post$pre_post <- replicate(nrow(data_post), 1)

data <- data.frame(rbind(data_pre, data_post))

data$pre_post <- as.factor(data$pre_post)
data$arch <- as.factor(data$arch)

# Initialize variables to store outputs from statistical tests
posthoc_tests_output <- data.frame(arch_pre_post = numeric(), condition = numeric(), variable = character(),
    levene = numeric(), shapiro = numeric(), p_value_within = numeric(), effect_size_within = numeric(), 
    wilcox = numeric(), wilcox_effect_size = numeric(), levene_between = numeric(), ks_between = numeric(), 
    p_value_between = numeric(), effect_size_between = numeric(), mann_whitney = numeric(), 
    mann_whitney_effect_size = numeric(),stringsAsFactors=FALSE)

levenes = list()
shapiro = list()
levenes_between = list()
ks_between = list()

anovas = list()
ttests_within = list()
wilcox_within = list()
ttests_between = list()
mw_between = list()

levenes_rt = list()
shapiro_rt = list()
levenes_rt_between = list()
ks_rt_between = list()

anovas_rt = list()
ttests_rt_within = list()
wilcox_rt_within = list()
ttests_rt_between = list()
mw_rt_between = list()

levenes_omit = list()
shapiro_omit = list()
levenes_omit_between = list()
ks_omit_between = list()

anovas_omit = list()
ttests_omit_within = list()
wilcox_omit_within = list()
ttests_omit_between = list()
mw_omit_between = list()


# Calculate ANOVAs using pre- versus post-surgery as within and injected virus as between animal factor
for (c in c(1, 2, 3)) {
    data_cond <- subset(data, condition==c)
  
    anova <- ezANOVA(data = data_cond, dv = performance, wid = subject, within = pre_post, 
        between = arch, type = 2, return_aov = TRUE)
    anovas <- append(anovas, list(anova_stats(anova$aov)))
  
    anova_rt <- ezANOVA(data = data_cond, dv = rt, wid = subject, within = pre_post, 
        between = arch, type = 2, return_aov = TRUE)
    anovas_rt <- append(anovas_rt, list(anova_stats(anova_rt$aov)))
    
    anova_omit <- ezANOVA(data = data_cond, dv = omit_pct, wid = subject, within = pre_post, 
                     between = arch, type = 2, return_aov = TRUE)
    anovas_omit <- append(anovas, list(anova_stats(anova_omit$aov)))
  
    # Calculate Levene's test, Shapiro's test and posthoc tests to compare pre- and post-surgery choice behavior,
    # reaction times and percentage of omitted trials for archaerhodopsin and control animals 
    # on each decision-making paradigm. In addition, calculate Levene's test, Kolmogorov-Smirnov's test,
    # t-tests and Mann-Whitney U tests to statistically compare choice behavior, reaction times and the
    # percentage of omitted trials between archaerhodopsin and control animals before and after the surgery
    for (ap in c(0, 1)) {
        data_arch <- subset(data_cond, arch==ap)
    
        levenes = append(levenes, list(leveneTest(performance~pre_post, data_arch)))
        shapiro = append(shapiro, list(shapiro.test(data_arch$performance[data_arch$pre_post=='0']-data_arch$performance[data_arch$pre_post=='1'])))
    
        ttests_within = append(ttests_within, list(t.test(data_arch$performance[data_arch$pre_post=='0'], 
            data_arch$performance[data_arch$pre_post=='1'], paired=TRUE)))
        wilcox_within = append(wilcox_within, list(wilcox.test(data_arch$performance[data_arch$pre_post=='0'], 
            data_arch$performance[data_arch$pre_post=='1'], paired=TRUE)))

        levenes_rt = append(levenes_rt, list(leveneTest(rt~pre_post, data_arch)))
        shapiro_rt = append(shapiro_rt, list(shapiro.test(data_arch$rt[data_arch$pre_post=='0']-data_arch$rt[data_arch$pre_post=='1'])))
    
        ttests_rt_within = append(ttests_rt_within, list(t.test(data_arch$rt[data_arch$pre_post=='0'], 
            data_arch$rt[data_arch$pre_post=='1'], paired=TRUE)))
        wilcox_rt_within = append(wilcox_rt_within, list(wilcox.test(data_arch$rt[data_arch$pre_post=='0'],
            data_arch$rt[data_arch$pre_post=='1'], paired=TRUE)))
        
        levenes_omit = append(levenes_omit, list(leveneTest(omit_pct~pre_post, data_arch)))
        shapiro_omit = append(shapiro_omit, list(shapiro.test(data_arch$omit_pct[data_arch$pre_post=='0']-data_arch$omit_pct[data_arch$pre_post=='1'])))
        
        ttests_omit_within = append(ttests_omit_within, list(t.test(data_arch$omit_pct[data_arch$pre_post=='0'], 
            data_arch$omit_pct[data_arch$pre_post=='1'], paired=TRUE)))
        wilcox_omit_within = append(wilcox_omit_within, list(wilcox.test(data_arch$omit_pct[data_arch$pre_post=='0'],
            data_arch$omit_pct[data_arch$pre_post=='1'], paired=TRUE)))
    
        data_pre_post <- subset(data_cond, pre_post==ap)
    
        levenes_between = append(levenes_between, list(leveneTest(performance~as.factor(arch), data_pre_post)))
        ks_between = append(ks_between, list(ks.test(data_pre_post$performance[data_pre_post$arch=='0'], data_pre_post$performance[data_pre_post$arch=='1'])))
        ttests_between = append(ttests_between, list(t.test(data_pre_post$performance[data_pre_post$arch=='1'], 
            data_pre_post$performance[data_pre_post$arch=='0'], paired=FALSE)))
        mw_between = append(mw_between, list(wilcox.test(data_pre_post$performance[data_pre_post$arch=='1'], 
            data_pre_post$performance[data_pre_post$arch=='0'], paired=FALSE)))
    
        levenes_rt_between = append(levenes_rt_between, list(leveneTest(rt~as.factor(arch), data_pre_post)))
        ks_rt_between = append(ks_rt_between, list(ks.test(data_pre_post$rt[data_pre_post$arch=='0'], data_pre_post$rt[data_pre_post$arch=='1'])))
        ttests_rt_between = append(ttests_rt_between, list(t.test(data_pre_post$rt[data_pre_post$arch=='1'], 
            data_pre_post$rt[data_pre_post$arch=='0'], paired=FALSE)))
        mw_rt_between = append(mw_rt_between, list(wilcox.test(data_pre_post$rt[data_pre_post$arch=='1'], 
            data_pre_post$rt[data_pre_post$arch=='0'], paired=FALSE)))
        
        levenes_omit_between = append(levenes_omit_between, list(leveneTest(omit_pct~as.factor(arch), data_pre_post)))
        ks_omit_between = append(ks_omit_between, list(ks.test(data_pre_post$omit_pct[data_pre_post$arch=='0'], data_pre_post$omit_pct[data_pre_post$arch=='1'])))
        ttests_omit_between = append(ttests_omit_between, list(t.test(data_pre_post$omit_pct[data_pre_post$arch=='1'], 
            data_pre_post$omit_pct[data_pre_post$arch=='0'], paired=FALSE)))
        mw_omit_between = append(mw_omit_between, list(wilcox.test(data_pre_post$omit_pct[data_pre_post$arch=='1'], 
            data_pre_post$omit_pct[data_pre_post$arch=='0'], paired=FALSE)))

        t = ttests_within[[length(ttests_within)]]$statistic
        df = ttests_within[[length(ttests_within)]]$parameter
        wp = wilcox_within[[length(wilcox_within)]]$p.value
        eff = abs(qnorm(wp/2)/sqrt(nrow(data_pre_post)))
        t_between = ttests_between[[length(ttests_between)]]$statistic
        df_between = ttests_between[[length(ttests_between)]]$parameter
        mwp = mw_between[[length(mw_between)]]$p.value
        mw_eff = abs(qnorm(mwp/2)/sqrt(nrow(data_pre_post)))
    
        posthoc_tests_output[nrow(posthoc_tests_output) + 1,] = list(ap, c, 'p', 
            levenes[[length(levenes)]]$`Pr(>F)`, shapiro[[length(shapiro)]]$p.value,
            ttests_within[[length(ttests_within)]]$p.value, sqrt(((t)^2)/(((t)^2)+df)), wp, eff,
            levenes_between[[length(levenes_between)]]$`Pr(>F)`, ks_between[[length(ks_between)]]$p.value, 
            ttests_between[[length(ttests_between)]]$p.value, sqrt(((t_between)^2)/(((t_between)^2)+df_between)), 
            mwp, mw_eff)

        t_rt = ttests_rt_within[[length(ttests_rt_within)]]$statistic
        df_rt = ttests_rt_within[[length(ttests_rt_within)]]$parameter
        wp_rt = wilcox_rt_within[[length(wilcox_rt_within)]]$p.value
        eff_rt = abs(qnorm(wp_rt/2)/sqrt(nrow(data_pre_post)))
        t_rt_between = ttests_rt_between[[length(ttests_rt_between)]]$statistic
        df_rt_between = ttests_rt_between[[length(ttests_rt_between)]]$parameter
        mwp_rt = mw_rt_between[[length(mw_rt_between)]]$p.value
        mw_eff_rt = abs(qnorm(mwp_rt/2)/sqrt(nrow(data_pre_post)))
    
        posthoc_tests_output[nrow(posthoc_tests_output) + 1,] = list(ap, c, 'rt', 
            levenes_rt[[length(levenes_rt)]]$`Pr(>F)`, shapiro_rt[[length(shapiro_rt)]]$p.value,
            ttests_rt_within[[length(ttests_rt_within)]]$p.value, sqrt(((t_rt)^2)/(((t_rt)^2)+df_rt)), wp_rt, eff_rt,
            levenes_rt_between[[length(levenes_rt_between)]]$`Pr(>F)`, ks_rt_between[[length(ks_rt_between)]]$p.value, 
            ttests_rt_between[[length(ttests_rt_between)]]$p.value, sqrt(((t_rt_between)^2)/(((t_rt_between)^2)+df_rt_between)), 
            mwp_rt, mw_eff_rt)
        
        t_omit = ttests_omit_within[[length(ttests_omit_within)]]$statistic
        df_omit = ttests_omit_within[[length(ttests_omit_within)]]$parameter
        wp_omit = wilcox_omit_within[[length(wilcox_omit_within)]]$p.value
        eff_omit = abs(qnorm(wp_omit/2)/sqrt(nrow(data_pre_post)))
        t_omit_between = ttests_omit_between[[length(ttests_omit_between)]]$statistic
        df_omit_between = ttests_omit_between[[length(ttests_omit_between)]]$parameter
        mwp_omit = mw_omit_between[[length(mw_omit_between)]]$p.value
        mw_eff_omit = abs(qnorm(mwp_omit/2)/sqrt(nrow(data_pre_post)))
        
        posthoc_tests_output[nrow(posthoc_tests_output) + 1,] = list(ap, c, 'o', 
             levenes_omit[[length(levenes_omit)]]$`Pr(>F)`, shapiro_omit[[length(shapiro_omit)]]$p.value,
             ttests_omit_within[[length(ttests_omit_within)]]$p.value, sqrt(((t_omit)^2)/(((t_omit)^2)+df_omit)), wp_omit, eff_omit,
             levenes_omit_between[[length(levenes_omit_between)]]$`Pr(>F)`, ks_omit_between[[length(ks_omit_between)]]$p.value, 
             ttests_omit_between[[length(ttests_omit_between)]]$p.value, sqrt(((t_omit_between)^2)/(((t_omit_between)^2)+df_omit_between)), 
             mwp_omit, mw_eff_omit)

    }
}

# Export the results of all statistical tests
write.csv(diff_tests_output, './analysis/diff_tests_output.csv', row.names = FALSE)
write.csv(posthoc_tests_output, './analysis/pre-post_posthoc_tests_output.csv', row.names = FALSE)

sink('./analysis/pre-post-anovas-output.txt')
print('Benefit-Benefit Decision-Making (Performance):')
print(anovas[[1]])
print(('Cost-Cost Decision-Making (Performance):'))
print(anovas[[2]])
print(('Cost-Benefit Decision-Making (Performance):'))
print(anovas[[3]])
print(('Benefit-Benefit Decision-Making (Reaction Time):'))
print(anovas_rt[[1]])
print(('Cost-Cost Decision-Making (Reaction Time):'))
print(anovas_rt[[2]])
print(('Cost-Benefit Decision-Making (Reaction Time):'))
print(anovas_rt[[3]])
print(('Benefit-Benefit Decision-Making (Omitted Trials):'))
print(anovas_omit[[1]])
print(('Cost-Cost Decision-Making (Omitted Trials):'))
print(anovas_omit[[2]])
print(('Cost-Benefit Decision-Making (Omitted Trials):'))
print(anovas_omit[[3]])
sink()
