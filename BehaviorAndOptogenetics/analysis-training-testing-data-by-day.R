# This script implements statistical tests and generates plots to explore if animals show a within session retention
# effect on the benefit-benefit and cost-cost decision-making paradigm. In addition, the script gives basic statistics
# on the time that animals required to learn lever pressing, the benefit-benefit, cost-cost and cost-benefit 
# decision-making paradigm. It analyses data that were previously extracted using custom Python scripts 
# from the Med Associates result files that record choice behavior and reaction times during the 
# last three days of behavioral traiing and during behavioral testing on the benefit-benefit, 
# cost-cost and cost-benefit decision-making paradimg.

# Import required libraries
library(plyr)
library(data.table)
library(car)

# Read in the lever pressing, training and testing data 
tidy_data_training <- read.table("./tidy_data_training.csv", header=TRUE, sep = ",")
tidy_data_day <- read.table("./tidy_data_day.csv", header=TRUE, sep = ",")
arch <- read.table("./arch.csv", header=TRUE, sep = ",")
data_lp <- read.table("./tidy_data_training_lp.csv", header=TRUE, sep = ",")
data_training_cb <- read.table("./tidy_data_training_cb.csv", header=TRUE, sep = ",")

# Change data into a format to calculate basic statistics and to run statistical tests
colnames(tidy_data_training)[colnames(tidy_data_training)=="block"] <- "opto_block"
colnames(tidy_data_day)[colnames(tidy_data_day)=="opto"] <- "opto_block"

data_training <- join(arch, tidy_data_training, by="subject", type="inner")
data_day <- join(arch, tidy_data_day, by="subject", type="inner")

# Define a function to calculate mean, standard error of the mean and the standard deviation
f_agg <- function(x) c(mean = mean(x), sem = sqrt(var(x)/length(x)), sd = sd(x))

# Calculate basic statistics for lever pressing training
data_lp <- join(arch, data_lp, by="subject", type="inner")
data_lp <- aggregate(date ~ subject, data = data_lp, FUN = length)
data_lp_stats <- f_agg(data_lp$date)

# Calculate basic statistics for training of the cost-benefit decsion-making paradigm
data_training_cb <- join(arch, data_training_cb, by="subject", type="inner")
data_training_cb <- aggregate(date ~ subject, data = unique(data_training_cb[,-c(2, 4:9)]), FUN = length)
data_training_cb_stats <- f_agg(data_training_cb$date)

# Calculate basic statistics for training of the benefit-benefit and cost-cost decsion-making paradigms
data_training_stats <- aggregate(day ~ subject + condition, data = data_training[,-c(2, 4:12)], FUN = min)
data_training_stats$day <- abs(data_training_stats$day)
data_training_bb_stats <- f_agg(data_training_stats$day[data_training_stats$condition == 1])
data_training_cc_stats <- f_agg(data_training_stats$day[data_training_stats$condition == 2])

# Save the calculated basic statistics in a text file
sink('./analysis/output-training.txt')
print(cat('Statistics for Lever Pressing Training: ', data_lp_stats))
print(cat('Statistics for Benefit-Benefit Decision-Making Training: ', data_training_bb_stats))
print(cat('Statistics for Cost-Cost Decision-Making Training: ', data_training_cc_stats))
print(cat('Statistics for Cost-Benefit Decision-Making Training: ', data_training_cb_stats))
sink()

# For statistical testing exclude any training data that is not from the last three days of
# training on the benefit-benefit, cost-cost or cost-benefit decision-making paradigm
data_training <- subset(data_training, day >= -3)

# Combine data from behavioral training and testing into one table
data_long <- rbind.fill(data_training, data_day)

# Initialize variables to save the output from statistical tests
day_tests_output <- data.frame(day = character(), arch = numeric(), condition = numeric(), variable = character(),
    shapiro = numeric(), levene = numeric(), ttest = numeric(), effect_size_ttest = numeric(), 
    wilcox = numeric(), effect_size_wilcox = numeric(), stringsAsFactors = FALSE)

levenes_training = list()
levenes_day = list()

shapiro_training = list()
shapiro_day = list()

ttests_training_within = list()
ttests_day_within = list()

wilcox_training_within = list()
wilcox_day_within = list()

levenes_rt_training = list()
levenes_rt_day = list()

shapiro_rt_training = list()
shapiro_rt_day = list()

ttests_rt_training_within = list()
ttests_rt_day_within = list()

wilcox_rt_training_within = list()
wilcox_rt_day_within = list()

# Calculate Levene's test, Shapiro's test, t-tests and Mann-Whitney U tests to compare the 
# choice behavior and reaction time on the last three days of behavioral training or
# on all days of behavioral testing between the first 20 (for behavioral testing these are the light OFF)
# and last 20 trials (for behavioral testing these are the light ON trials) that were recorded on
# the benefit-benefit and cost-cost decision-making paradigm for archaerhodopsin and control animals
for (c in c(1, 2)) {
    
    data_cond_training <- subset(data_training, condition==c)
    data_cond_day <- subset(data_day, condition==c)
  
    for (a in c(0, 1)) {
        data_arch <- subset(data_cond_training, arch==a)
        
        levenes_training = append(levenes_training, list(leveneTest(performance~as.factor(opto_block), data_arch)))
        shapiro_training = append(shapiro_training, list(shapiro.test(data_arch$performance[data_arch$opto_block=='0']-data_arch$performance[data_arch$opto_block=='1'])))
    
        ttests_training_within = append(ttests_training_within, 
            list(t.test(data_arch$performance[data_arch$opto_block=='0'], 
            data_arch$performance[data_arch$opto_block=='1'], paired=TRUE)))
        wilcox_training_within = append(wilcox_training_within, 
            list(wilcox.test(data_arch$performance[data_arch$opto_block=='0'], 
            data_arch$performance[data_arch$opto_block=='1'], paired=TRUE)))
        
        levenes_rt_training = append(levenes_rt_training, list(leveneTest(rt_diff~as.factor(opto_block), data_arch)))
        shapiro_rt_training = append(shapiro_rt_training, list(shapiro.test(data_arch$rt_diff[data_arch$opto_block=='0']-data_arch$rt_diff[data_arch$opto_block=='1'])))
        
        ttests_rt_training_within = append(ttests_rt_training_within, 
            list(t.test(data_arch$rt_diff[data_arch$opto_block=='0'], 
            data_arch$rt_diff[data_arch$opto_block=='1'], paired=TRUE)))
        wilcox_rt_training_within = append(wilcox_rt_training_within, 
            list(wilcox.test(data_arch$rt_diff[data_arch$opto_block=='0'], 
            data_arch$rt_diff[data_arch$opto_block=='1'], paired=TRUE)))
        
        t = ttests_training_within[[length(ttests_training_within)]]$statistic
        df = ttests_training_within[[length(ttests_training_within)]]$parameter
        wp = wilcox_training_within[[length(wilcox_training_within)]]$p.value
        eff = abs(qnorm(wp/2)/sqrt(nrow(data_arch)))
        
        day_tests_output[nrow(day_tests_output) + 1,] = list('training', a, c, 'p', 
            levenes_training[[length(levenes_training)]]$`Pr(>F)`, shapiro_training[[length(shapiro_training)]]$p.value,
            ttests_training_within[[length(ttests_training_within)]]$p.value, sqrt(((t)^2)/(((t)^2)+df)), wp, eff)
        
        t_rt = ttests_rt_training_within[[length(ttests_rt_training_within)]]$statistic
        df_rt = ttests_rt_training_within[[length(ttests_rt_training_within)]]$parameter
        wp_rt = wilcox_rt_training_within[[length(wilcox_rt_training_within)]]$p.value
        eff_rt = abs(qnorm(wp_rt/2)/sqrt(nrow(data_arch)))
        
        day_tests_output[nrow(day_tests_output) + 1,] = list('training', a, c, 'rt', 
            levenes_rt_training[[length(levenes_rt_training)]]$`Pr(>F)`,  shapiro_rt_training[[length(shapiro_rt_training)]]$p.value, 
            ttests_rt_training_within[[length(ttests_rt_training_within)]]$p.value, sqrt(((t_rt)^2)/(((t_rt)^2)+df_rt)), wp_rt, eff_rt)
    
        data_arch <- subset(data_cond_day, arch==a)
        
        levenes_day = append(levenes_day, list(leveneTest(performance~as.factor(opto_block), data_arch)))
        shapiro_day = append(shapiro_day, list(shapiro.test(data_arch$performance[data_arch$opto_block=='0']-data_arch$performance[data_arch$opto_block=='1'])))
    
        ttests_day_within = append(ttests_day_within, 
            list(t.test(data_arch$performance[data_arch$opto_block=='0'], 
            data_arch$performance[data_arch$opto_block=='1'], paired=TRUE)))
        wilcox_day_within = append(wilcox_day_within, 
            list(wilcox.test(data_arch$performance[data_arch$opto_block=='0'], 
            data_arch$performance[data_arch$opto_block=='1'], paired=TRUE)))
        
        levenes_rt_day = append(levenes_rt_day, list(leveneTest(rt_diff~as.factor(opto_block), data_arch)))
        shapiro_rt_day = append(shapiro_rt_day, list(shapiro.test(data_arch$rt_diff[data_arch$opto_block=='0']-data_arch$rt_diff[data_arch$opto_block=='1'])))
        
        ttests_rt_day_within = append(ttests_rt_day_within, 
            list(t.test(data_arch$rt_diff[data_arch$opto_block=='0'], 
            data_arch$rt_diff[data_arch$opto_block=='1'], paired=TRUE)))
        wilcox_rt_day_within = append(wilcox_rt_day_within, 
            list(wilcox.test(data_arch$rt_diff[data_arch$opto_block=='0'], 
            data_arch$rt_diff[data_arch$opto_block=='1'], paired=TRUE)))
        
        t = ttests_day_within[[length(ttests_day_within)]]$statistic
        df = ttests_day_within[[length(ttests_day_within)]]$parameter
        wp = wilcox_day_within[[length(wilcox_day_within)]]$p.value
        eff = abs(qnorm(wp/2)/sqrt(nrow(data_arch)))
        
        day_tests_output[nrow(day_tests_output) + 1,] = list('testing', a, c, 'p', 
            levenes_day[[length(levenes_day)]]$`Pr(>F)`, shapiro_day[[length(shapiro_day)]]$p.value,
            ttests_day_within[[length(ttests_day_within)]]$p.value, sqrt(((t)^2)/(((t)^2)+df)), wp, eff)
        
        t_rt = ttests_rt_day_within[[length(ttests_rt_day_within)]]$statistic
        df_rt = ttests_rt_day_within[[length(ttests_rt_day_within)]]$parameter
        wp_rt = wilcox_rt_day_within[[length(wilcox_rt_day_within)]]$p.value
        eff_rt = abs(qnorm(wp_rt/2)/sqrt(nrow(data_arch)))
        
        day_tests_output[nrow(day_tests_output) + 1,] = list('testing', a, c, 'rt', 
            levenes_rt_day[[length(levenes_rt_day)]]$`Pr(>F)`,  shapiro_rt_day[[length(shapiro_rt_day)]]$p.value, 
            ttests_rt_day_within[[length(ttests_rt_day_within)]]$p.value, sqrt(((t_rt)^2)/(((t_rt)^2)+df_rt)), wp_rt, eff_rt)
  }
}

# Save the output of statistical tests in a csv file
write.csv(day_tests_output, './analysis/day_tests_output.csv', row.names = FALSE)

# subset data from the benefit-benefit, cost-cost and cost-benefit decision-making paradigm
data_bb <- subset(data_long, condition == 1)
data_cc <- subset(data_long, condition == 2)
data_cb <- subset(data_long, condition == 3)

# Plot the mean percentage of high benefit, high cost or high benefit-high cost trials
# out of the total number of trials that were not omitted for the last three days of behavioral
# training and all three days of behavioral testing for archaerhodopsin and control animals
# separated for the first 20 (for behavioral testing these are the light OFF) and last 20 trials 
# (for behavioral testing these are the light ON trials) for the benefit-benefit and cost-cost
# decision-making paradigm, and only for all three days of behavioral testing for the 
# cost-benefit decision-making paradigm

p <- 1

tiff("figures/benefit-benefit-day.tiff", width=4800, height=2400, res=300)
par(mar=c(8.6, 9.6, 5.1, 4.6), mgp=c(5.5, 1.5, 0), family="Cambria")
plot(performance~day, data=data_bb, type='n', xlim=c(-3,3), ylim=c(20,100), axes=FALSE, xlab="", ylab="", frame.plot=FALSE)

for (i in c(0, 1)) {
  for (j in c(0, 1)) {
      data_train <- do.call(data.frame, aggregate(performance ~ day, data = subset(data_bb, arch == i & opto_block == j & day < 0), f_agg))
      data_test <- do.call(data.frame, aggregate(performance ~ day, data = subset(data_bb, arch == i & opto_block == j & day > 0), f_agg))
      lines(performance.mean~day, data=data_train, type="o", pch=c(15,17,0,2)[p], cex=3)
      lines(performance.mean~day, data=data_test, type="o", pch=c(15,17,0,2)[p], cex=3)
      arrows(data_train$day, data_train$performance.mean-data_train$performance.sem, data_train$day,
          data_train$performance.mean+data_train$performance.sem, code=3, angle=90, length=0.1)
      arrows(data_test$day, data_test$performance.mean-data_test$performance.sem, data_test$day,
          data_test$performance.mean+data_test$performance.sem, code=3, angle=90, length=0.1)
      p = p+1
  }
}

title(ylab="High benefit choices (%)", xlab="Days", cex.lab=2.3, font.lab=2)
mtext(text="A- Benefit-Benefit Decision-Making", side=3, line=2, at=-2.21, cex=2.6, font=2)
axis(side = 1, at=seq(-3,3,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(20,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
dev.off()

p <- 1

tiff("figures/cost-cost-day.tiff", width=4800, height=2400, res=300)
par(mar=c(8.6, 9.6, 5.1, 4.6), mgp=c(5.5, 1.5, 0), family="Cambria")
plot(performance~day, data=data_cc, type='n', xlim=c(-3,3), ylim=c(0,80), axes=FALSE, xlab="", ylab="", frame.plot=FALSE)

for (i in c(0, 1)) {
  for (j in c(0, 1)) {
    data_train <- do.call(data.frame, aggregate(performance ~ day, data = subset(data_cc, arch == i & opto_block == j & day < 0), f_agg))
    data_test <- do.call(data.frame, aggregate(performance ~ day, data = subset(data_cc, arch == i & opto_block == j & day > 0), f_agg))
    lines(performance.mean~day, data=data_train, type="o", pch=c(15,17,0,2)[p], cex=3)
    lines(performance.mean~day, data=data_test, type="o", pch=c(15,17,0,2)[p], cex=3)
    arrows(data_train$day, data_train$performance.mean-data_train$performance.sem, data_train$day,
        data_train$performance.mean+data_train$performance.sem, code=3, angle=90, length=0.1)
    arrows(data_test$day, data_test$performance.mean-data_test$performance.sem, data_test$day,
        data_test$performance.mean+data_test$performance.sem, code=3, angle=90, length=0.1)
    p = p+1
  }
}

title(ylab="High cost choices (%)", xlab="Days", cex.lab=2.3, font.lab=2)
mtext(text="B- Cost-Cost Decision-Making", side=3, line=2, at=-2.47, cex=2.6, font=2)
axis(side = 1, at=seq(-3,3,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(0,80,10), las=1, lwd=5, cex.axis=2.2, font=2)
dev.off()

p <- 1

tiff("figures/cost-benefit-day.tiff", width=4800, height=2400, res=300)
par(mar=c(8.6, 9.6, 5.1, 25.1), mgp=c(5.5, 1.5, 0), family="Cambria")
plot(performance~day, data=data_cb, type='n', xlim=c(0,3), ylim=c(20,100), axes=FALSE, xlab="", ylab="", frame.plot=FALSE)

for (i in c(0, 1)) {
  for (j in c(0, 1)) {
    data_test <- do.call(data.frame, aggregate(performance ~ day, data = subset(data_cb, arch == i & opto_block == j & day > 0), f_agg))
    lines(performance.mean~day, data=data_test, type="o", pch=c(15,17,0,2)[p], cex=3)
    arrows(data_test$day, data_test$performance.mean-data_test$performance.sem, data_test$day,
        data_test$performance.mean+data_test$performance.sem, code=3, angle=90, length=0.1)
    p = p+1
  }
}

title(ylab="High cost-benefit choices (%)", xlab="Days", cex.lab=2.3, font.lab=2)
mtext(text="C- Cost-Benefit Decision-Making", side=3, line=2, at=0.53, cex=2.6, font=2)
axis(side=1, at=seq(0,3,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(20,100,10), las=1, lwd=5, cex.axis=2.2, font=2)
legend('right', bty='n', inset=c(-0.50,0), xpd = TRUE, legend=c("Control,\nBlock 1/Light OFF", 
    "Control,\nBlock 2/Light ON", "Archaerhodopsin,\nBlock 1/Light OFF", "Archaerhodopsin,\nBlock 2/Light ON"), 
    y.intersp=2, pch=c(15,17,0,2), lwd=2, cex=1.8, text.font=1)
dev.off()

# Plot the mean of the difference in reaction time on trials, in which animals choose the 
# high benefit, high cost or high benefit-high cost option, as compared to trials, in which animals
# choose the low benefit, low cost or low benefit-low cost option for the last three days of behavioral
# training and all three days of behavioral testing for archaerhodopsin and control animals
# separated for the first 20 (for behavioral testing these are the light OFF) and last 20 trials 
# (for behavioral testing these are the light ON trials) for the benefit-benefit and cost-cost
# decision-making paradigm, and only for all three days of behavioral testing for the 
# cost-benefit decision-making paradigm

p <- 1

tiff("figures/benefit-benefit-rt-day.tiff", width=4800, height=2400, res=300)
par(mar=c(8.6, 9.6, 4.6, 4.6), mgp=c(5, 1.5, 0), family="Cambria")
plot(performance~day, data=data_bb, type='n', xlim=c(-3,3), ylim=c(-4,4), axes=FALSE, xlab="", ylab="", frame.plot=FALSE)

for (i in c(0, 1)) {
  for (j in c(0, 1)) {
    data_train <- do.call(data.frame, aggregate(rt_diff ~ day, data = subset(data_bb, arch == i & opto_block == j & day < 0), f_agg))
    data_test <- do.call(data.frame, aggregate(rt_diff ~ day, data = subset(data_bb, arch == i & opto_block == j & day > 0), f_agg))
    lines(rt_diff.mean~day, data=data_train, type="o", pch=c(15,17,0,2)[p], cex=3)
    lines(rt_diff.mean~day, data=data_test, type="o", pch=c(15,17,0,2)[p], cex=3)
    arrows(data_train$day, data_train$rt_diff.mean-data_train$rt_diff.sem, data_train$day,
        data_train$rt_diff.mean+data_train$rt_diff.sem, code=3, angle=90, length=0.1)
    arrows(data_test$day, data_test$rt_diff.mean-data_test$rt_diff.sem, data_test$day,
        data_test$rt_diff.mean+data_test$rt_diff.sem, code=3, angle=90, length=0.1)
    p = p+1
  }
}

title(ylab="Difference in mean reaction time\nfor high-low benefit choices (s)", 
    xlab="Days", cex.lab=1.8, font.lab=2)
mtext(text="A- Benefit-Benefit Decision-Making", side=3, line=2, at=-2.21, cex=2.6, font=2)
axis(side = 1, at=seq(-3,3,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), las=1, lwd=5, cex.axis=2.2, font=2)
dev.off()

p <- 1

tiff("figures/cost-cost-rt_day.tiff", width=4800, height=2400, res=300)
par(mar=c(8.6, 9.6, 4.6, 4.6), mgp=c(5, 1.5, 0), family="Cambria")
plot(performance~day, data=data_cc, type='n', xlim=c(-3,3), ylim=c(-4,4), axes=FALSE, xlab="", ylab="", frame.plot=FALSE)

for (i in c(0, 1)) {
  for (j in c(0, 1)) {
    data_train <- do.call(data.frame, aggregate(rt_diff ~ day, data = subset(data_cc, arch == i & opto_block == j & day < 0), f_agg))
    data_test <- do.call(data.frame, aggregate(rt_diff ~ day, data = subset(data_cc, arch == i & opto_block == j & day > 0), f_agg))
    lines(rt_diff.mean~day, data=data_train, type="o", pch=c(15,17,0,2)[p], cex=3)
    lines(rt_diff.mean~day, data=data_test, type="o", pch=c(15,17,0,2)[p], cex=3)
    arrows(data_train$day, data_train$rt_diff.mean-data_train$rt_diff.sem, data_train$day,
        data_train$rt_diff.mean+data_train$rt_diff.sem, code=3, angle=90, length=0.1)
    arrows(data_test$day, data_test$rt_diff.mean-data_test$rt_diff.sem, data_test$day,
        data_test$rt_diff.mean+data_test$rt_diff.sem, code=3, angle=90, length=0.1)
    p = p+1
  }
}

title(ylab="Difference in mean reaction time\nfor high-low cost choices (s)", 
    xlab="Days", cex.lab=1.8, font.lab=2)
mtext(text="B- Cost-Cost Decision-Making", side=3, line=2, at=-2.47, cex=2.6, font=2)
axis(side = 1, at=seq(-3,3,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), las=1, lwd=5, cex.axis=2.2, font=2)
dev.off()

p <- 1

tiff("figures/cost-benefit-rt_day.tiff", width=4800, height=2400, res=300)
par(mar=c(8.6, 9.6, 4.6, 25.1), mgp=c(5, 1.5, 0), family="Cambria")
plot(performance~day, data=data_cb, type='n', xlim=c(0,3), ylim=c(-4,4), axes=FALSE, xlab="", ylab="", frame.plot=FALSE)

for (i in c(0, 1)) {
  for (j in c(0, 1)) {
    data_test <- do.call(data.frame, aggregate(rt_diff ~ day, data = subset(data_cb, arch == i & opto_block == j & day > 0), f_agg))
    lines(rt_diff.mean~day, data=data_test, type="o", pch=c(15,17,0,2)[p], cex=3)
    arrows(data_test$day, data_test$rt_diff.mean-data_test$rt_diff.sem, data_test$day,
        data_test$rt_diff.mean+data_test$rt_diff.sem, code=3, angle=90, length=0.1)
    p = p+1
  }
}

title(ylab="Difference in mean reaction time\nfor high-low cost-benefit choices (s)",
    xlab="Days", cex.lab=1.8, font.lab=2)
mtext(text="C- Cost-Benefit Decision-Making", side=3, line=2, at=0.53, cex=2.6, font=2)
axis(side=1, at=seq(0,3,1), las=1, lwd=5, cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), las=1, lwd=5, cex.axis=2.2, font=2)
legend('right', bty='n', inset=c(-0.50,0), xpd = TRUE, legend=c("Control,\nBlock 1/Light OFF", 
    "Control,\nBlock 2/Light ON", "Archaerhodopsin,\nBlock 1/Light OFF", "Archaerhodopsin,\nBlock 2/Light ON"), 
    y.intersp=2, pch=c(15,17,0,2), lwd=2, cex=1.8, text.font=1)
dev.off()
