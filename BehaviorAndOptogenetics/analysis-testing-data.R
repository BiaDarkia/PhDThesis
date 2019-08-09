# This script provides a tool to further analyse data that were previously extracted using 
# custom Python scripts from the Med Associates result files that record choice behavior and 
# reaction times during behavioral testing on the benefit-benefit, cost-cost and cost-benefit 
# decision-making paradimg. This program runs statistical tests and generates figures to visualize the data.

# Import required libraries
library(plyr)
library(data.table)
library(car)
library(ez)
library(sjstats)
library(RColorBrewer)

# Read the data and change the format to easily run statistical tests
tidy_data <- read.table("./tidy_data.csv", header=TRUE, sep = ",")
arch <- read.table("./arch.csv", header=TRUE, sep = ",")

data_long <- join(arch, tidy_data, by="subject", type="inner")
data_long$opto <- as.factor(data_long$opto)
data_long$arch <- as.factor(data_long$arch)
data_long$omit_pct <- (data_long$omitted/data_long$total)*100

# Initialize variables to store outputs from statistical tests
posthoc_tests_output <- data.frame(arch = numeric(), condition = numeric(), variable = character(),
    levene = numeric(), shapiro = numeric(), p_value = numeric(), effect_size = numeric(), stringsAsFactors=FALSE)

levenes = list()
shapiro = list()

anovas = list()
ttests_within = list()

levenes_rt_diff = list()
shapiro_rt_diff = list()

anovas_rt_diff = list()
ttests_rt_diff_within = list()

levenes_rt = list()
shapiro_rt = list()

anovas_rt = list()
ttests_rt_within = list()

levenes_omit = list()
shapiro_omit = list()

anovas_omit = list()
ttests_omit_within = list()

# Calculate ANOVAs using light ON versus OFF as within and injected virus as between animal factor
for (c in c(1, 2, 3)) {
        data_cond <- subset(data_long, condition==c)
    
        anova <- ezANOVA(data = data_cond, dv = performance, wid = subject, within = opto, 
            between = arch, type = 2, return_aov = TRUE)
        anovas <- append(anovas, list(anova_stats(anova$aov)))
        
        anova_rt_diff <- ezANOVA(data = data_cond, dv = rt_diff, wid = subject, within = opto, 
            between = arch, type = 2, return_aov = TRUE)
        anovas_rt_diff <- append(anovas_rt_diff, list(anova_stats(anova_rt_diff$aov)))
        
        anova_rt <- ezANOVA(data = data_cond, dv = rt, wid = subject, within = opto, 
                                 between = arch, type = 2, return_aov = TRUE)
        anovas_rt <- append(anovas_rt, list(anova_stats(anova_rt$aov)))
        
        anova_omit <- ezANOVA(data = data_cond, dv = omit_pct, wid = subject, within = opto, 
            between = arch, type = 2, return_aov = TRUE)
        anovas_omit <- append(anovas_omit, list(anova_stats(anova_omit$aov)))

        # Calculate Levene's test, Shapiro's test and posthoc tests to compare light ON versus OFF condition
        # for archaerhodopsin and control animals for each decision-making paradigm
        for (a in c(1, 0)) {
                data_arch <- subset(data_cond, arch==a)
                
                levenes = append(levenes, list(leveneTest(performance~as.factor(opto), data_arch)))
                shapiro = append(shapiro, list(shapiro.test(data_arch$performance[data_arch$opto=='0']-data_arch$performance[data_arch$opto=='1'])))
                
                ttests_within = append(ttests_within, list(t.test(data_arch$performance[data_arch$opto=='0'], 
                    data_arch$performance[data_arch$opto=='1'], paired=TRUE)))
                
                levenes_rt_diff = append(levenes_rt_diff, list(leveneTest(rt_diff~as.factor(opto), data_arch)))
                shapiro_rt_diff = append(shapiro_rt_diff, list(shapiro.test(data_arch$rt_diff[data_arch$opto=='0']-data_arch$rt_diff[data_arch$opto=='1'])))
                
                ttests_rt_diff_within = append(ttests_rt_diff_within, list(t.test(data_arch$rt_diff[data_arch$opto=='0'], 
                    data_arch$rt_diff[data_arch$opto=='1'], paired=TRUE)))
                
                levenes_rt = append(levenes_rt, list(leveneTest(rt~as.factor(opto), data_arch)))
                shapiro_rt = append(shapiro_rt, list(shapiro.test(data_arch$rt[data_arch$opto=='0']-data_arch$rt[data_arch$opto=='1'])))
                
                ttests_rt_within = append(ttests_rt_within, list(t.test(data_arch$rt[data_arch$opto=='0'], 
                    data_arch$rt[data_arch$opto=='1'], paired=TRUE)))
                
                levenes_omit = append(levenes_omit, list(leveneTest(omit_pct~as.factor(opto), data_arch)))
                shapiro_omit = append(shapiro_omit, list(shapiro.test(data_arch$omit_pct[data_arch$opto=='0']-data_arch$omit_pct[data_arch$opto=='1'])))
                
                ttests_omit_within = append(ttests_omit_within, list(t.test(data_arch$omit_pct[data_arch$opto=='0'], 
                    data_arch$omit_pct[data_arch$opto=='1'], paired=TRUE)))
                
                t = ttests_within[[length(ttests_within)]]$statistic
                df = ttests_within[[length(ttests_within)]]$parameter
                
                posthoc_tests_output[nrow(posthoc_tests_output) + 1,] = list(a, c, 'p', 
                    levenes[[length(levenes)]]$`Pr(>F)`, shapiro[[length(shapiro)]]$p.value,
                    ttests_within[[length(ttests_within)]]$p.value, sqrt(((t)^2)/(((t)^2)+df)))
                
                t_rt_diff = ttests_rt_diff_within[[length(ttests_rt_diff_within)]]$statistic
                df_rt_diff = ttests_rt_diff_within[[length(ttests_rt_diff_within)]]$parameter
                
                posthoc_tests_output[nrow(posthoc_tests_output) + 1,] = list(a, c, 'rt_diff', 
                    levenes_rt_diff[[length(levenes_rt_diff)]]$`Pr(>F)`, shapiro_rt_diff[[length(shapiro_rt_diff)]]$p.value,
                    ttests_rt_diff_within[[length(ttests_rt_diff_within)]]$p.value, sqrt(((t_rt_diff)^2)/(((t_rt_diff)^2)+df_rt_diff)))
                
                t_rt = ttests_rt_within[[length(ttests_rt_within)]]$statistic
                df_rt = ttests_rt_within[[length(ttests_rt_within)]]$parameter
                
                posthoc_tests_output[nrow(posthoc_tests_output) + 1,] = list(a, c, 'rt', 
                     levenes_rt[[length(levenes_rt)]]$`Pr(>F)`, shapiro_rt[[length(shapiro_rt)]]$p.value,
                     ttests_rt_within[[length(ttests_rt_within)]]$p.value, sqrt(((t_rt)^2)/(((t_rt)^2)+df_rt)))
                
                t_omit = ttests_omit_within[[length(ttests_omit_within)]]$statistic
                df_omit = ttests_omit_within[[length(ttests_omit_within)]]$parameter
                
                posthoc_tests_output[nrow(posthoc_tests_output) + 1,] = list(a, c, 'o', 
                    levenes_omit[[length(levenes_omit)]]$`Pr(>F)`, shapiro_omit[[length(shapiro_omit)]]$p.value,
                    ttests_omit_within[[length(ttests_omit_within)]]$p.value, sqrt(((t_omit)^2)/(((t_omit)^2)+df_omit)))
                
       }
}

data_bb_arch <- subset(data_long, condition==1 & arch==1)
data_bb_sham <- subset(data_long, condition==1 & arch==0)
data_cc_arch <- subset(data_long, condition==2 & arch==1)
data_cc_sham <- subset(data_long, condition==2 & arch==0)
data_cb_arch <- subset(data_long, condition==3 & arch==1)
data_cb_sham <- subset(data_long, condition==3 & arch==0)

labels_pct <- c("0","","","","","50","","","","","100")

# Interaction plot of the omitted trials for each archaerhodopsin and 
# each control animal for all three decision-making paradigms
tiff("figures/interaction-omit.tiff", width=5600, height=7200, res=300)
par(mfrow=c(3,2), mar=c(5.1, 8.6, 10.6, 0.1), mgp=c(6, 1.5, 0))
boxplot(omit_pct~opto, data = data_cb_arch, ylim=c(0,100), axes=FALSE, xlab="", 
     ylab="High cost-benefit choices (%)", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cb_arch$opto, data_cb_arch$subject, data_cb_arch$omit_pct, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 85, labels="n=10", cex=2.2)
mtext(text="A- Cost-Benefit Decision-Making", side=3, line=3, at=1.3, cex=2.8, font=2)
boxplot(omit_pct~opto, data = data_cb_sham, ylim=c(0,100), axes=FALSE, xlab="", 
     ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cb_sham$opto, data_cb_sham$subject, data_cb_sham$omit_pct, type="b", 
     legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 85, labels="n=12", cex=2.2)
boxplot(omit_pct~opto, data = data_bb_arch, ylim=c(0,100), axes=FALSE, xlab="", 
     ylab="High benefit choices (%)", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_arch$opto, data_bb_arch$subject, data_bb_arch$omit_pct, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main = 2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 85, labels="n=10", cex=2.2)
mtext(text="B- Benefit-Benefit Decision-Making", side=3, line=3, at=1.34, cex=2.8, font=2)
boxplot(omit_pct~opto, data = data_bb_sham, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_sham$opto, data_bb_sham$subject, data_bb_sham$omit_pct, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 85, labels="n=12", cex=2.2)
boxplot(omit_pct~opto, data = data_cc_arch, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="High cost choices (%)", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cc_arch$opto, data_cc_arch$subject, data_cc_arch$omit_pct, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 85, labels="n=10", cex=2.2)
mtext(text="C- Cost-Cost Decision-Making", side=3, line=3, at=1.26, cex=2.8, font=2)
boxplot(omit_pct~opto, data = data_cc_sham, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cc_sham$opto, data_cc_sham$subject, data_cc_sham$omit_pct, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 85, labels="n=12", cex=2.2)
dev.off()

# Interaction plot of the percentage of high benefit, high cost or high benefit-high cost trials
# out of the total number of trials that were not omitted for each archaerhodopsin and 
# each control animal for all three decision-making paradigms
tiff("figures/interaction-opto.tiff", width=5600, height=8400, res=300)
par(mfrow=c(3,2), mar=c(5.1, 8.6, 15.1, 0.1), mgp=c(6, 1.5, 0))
boxplot(performance~opto, data = data_cb_arch, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="High cost-benefit choices (%)", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cb_arch$opto, data_cb_arch$subject, data_cb_arch$performance, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main=2.2, line=5, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 0.5, labels="n=10", cex=2.2)
mtext(text="A- Cost-Benefit Decision-Making", side=3, line=9, at=1.3, cex=2.8, font=2)
boxplot(performance~opto, data = data_cb_sham, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cb_sham$opto, data_cb_sham$subject, data_cb_sham$performance, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=5, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 0.5, labels="n=12", cex=2.2)
boxplot(performance~opto, data = data_bb_arch, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="High benefit choices (%)", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_arch$opto, data_bb_arch$subject, data_bb_arch$performance, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main = 2.2, line=5, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 0.5, labels="n=10", cex=2.2)
mtext(text="B- Benefit-Benefit Decision-Making", side=3, line=9, at=1.34, cex=2.8, font=2)
boxplot(performance~opto, data = data_bb_sham, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_sham$opto, data_bb_sham$subject, data_bb_sham$performance, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=5, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 0.5, labels="n=12", cex=2.2)
boxplot(performance~opto, data = data_cc_arch, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="High cost choices (%)", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cc_arch$opto, data_cc_arch$subject, data_cc_arch$performance, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main=2.2, line=5, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 0.5, labels="n=10", cex=2.2)
mtext(text="C- Cost-Cost Decision-Making", side=3, line=9, at=1.26, cex=2.8, font=2)
boxplot(performance~opto, data = data_cc_sham, ylim=c(0,100), axes=FALSE, xlab="", 
    ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cc_sham$opto, data_cc_sham$subject, data_cc_sham$performance, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=5, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,100,10), labels=labels_pct, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, 0.5, labels="n=12", cex=2.2)
dev.off()

labels_rt <- c("0","","2","","6","","2","","8")

# Interaction plot of the overall mean reaction time on all trials that were not omitted
# for each archaerhodopsin and  each control animal for all three decision-making paradigms
tiff("figures/interaction-rt-opto.tiff", width=5600, height=7200, res=300)
par(mfrow=c(3,2), mar=c(5.1, 10.6, 10.6, 0.1), mgp=c(6, 1.5, 0))
boxplot(rt~opto, data = data_cb_arch, ylim=c(0,8), axes=FALSE, xlab="", 
     ylab="Mean reaction time (s)",
     cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cb_arch$opto, data_cb_arch$subject, data_cb_arch$rt, type="b", 
     legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,8,1), labels=labels_rt, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=10", cex=2.2)
mtext(text="A- Cost-Benefit Decision-Making", side=3, line=3, at=1.26, cex=2.8, font=2)
boxplot(rt~opto, data = data_cb_sham, ylim=c(0,8), axes=FALSE, xlab="", 
     ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cb_sham$opto, data_cb_sham$subject, data_cb_sham$rt, type="b", 
     legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,8,1), labels=labels_rt, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=12", cex=2.2)
boxplot(rt~opto, data = data_bb_arch, ylim=c(0,8), axes=FALSE, xlab="", 
     ylab="Mean reaction time (s)",
     cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_arch$opto, data_bb_arch$subject, data_bb_arch$rt, type="b", 
     legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main = 2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,8,1), labels=labels_rt, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=10", cex=2.2)
mtext(text="B- Benefit-Benefit Decision-Making", side=3, line=3, at=1.32, cex=2.8, font=2)
boxplot(rt~opto, data = data_bb_sham, ylim=c(0,8), axes=FALSE, xlab="", 
     ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_sham$opto, data_bb_sham$subject, data_bb_sham$rt, type="b", 
     legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,8,1), labels=labels_rt, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=12", cex=2.2)
boxplot(rt~opto, data = data_cc_arch, ylim=c(0,8), axes=FALSE, xlab="", 
     ylab="Mean reaction time (s)",
     cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cc_arch$opto, data_cc_arch$subject, data_cc_arch$rt, type="b", 
     legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,8,1), labels=labels_rt, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=10", cex=2.2)
mtext(text="C- Cost-Cost Decision-Making", side=3, line=3, at=1.2, cex=2.8, font=2)
boxplot(rt~opto, data = data_bb_sham, ylim=c(0,8), axes=FALSE, xlab="", 
     ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_sham$opto, data_bb_sham$subject, data_bb_sham$rt, type="b", 
     legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(0,8,1), labels=labels_rt, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=12", cex=2.2)
dev.off()

labels_rt_diff <- c("-4","","-2","","0","","2","","4")

# Interaction plot of the difference in reaction time on trials, in which animals choose the 
# high benefit, high cost or high benefit-high cost option, as compared to trials, in which animals
# choose the low benefit, low cost or low benefit-low cost option for each archaerhodopsin and 
# each control animal for all three decision-making paradigms
tiff("figures/interaction-rt_diff-opto.tiff", width=5600, height=7200, res=300)
par(mfrow=c(3,2), mar=c(5.1, 10.6, 10.6, 0.1), mgp=c(6, 1.5, 0))
boxplot(rt_diff~opto, data = data_cb_arch, ylim=c(-4,4), axes=FALSE, xlab="", 
    ylab="Difference in mean reaction time\nfor high-low cost-benefit choices (s)",
    cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cb_arch$opto, data_cb_arch$subject, data_cb_arch$rt_diff, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=10", cex=2.2)
mtext(text="A- Cost-Benefit Decision-Making", side=3, line=3, at=1.3, cex=2.8, font=2)
boxplot(rt_diff~opto, data = data_cb_sham, ylim=c(-4,4), axes=FALSE, xlab="", 
    ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cb_sham$opto, data_cb_sham$subject, data_cb_sham$rt_diff, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=12", cex=2.2)
boxplot(rt_diff~opto, data = data_bb_arch, ylim=c(-4,4), axes=FALSE, xlab="", 
    ylab="Difference in mean reaction time\nfor high-low benefit choices (s)",
    cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_arch$opto, data_bb_arch$subject, data_bb_arch$rt_diff, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main = 2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=10", cex=2.2)
mtext(text="B- Benefit-Benefit Decision-Making", side=3, line=3, at=1.34, cex=2.8, font=2)
boxplot(rt_diff~opto, data = data_bb_sham, ylim=c(-4,4), axes=FALSE, xlab="", 
    ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_sham$opto, data_bb_sham$subject, data_bb_sham$rt_diff, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=12", cex=2.2)
boxplot(rt_diff~opto, data = data_cc_arch, ylim=c(-4,4), axes=FALSE, xlab="", 
    ylab="Difference in mean reaction time\nfor high-low cost choices (s)",
    cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_cc_arch$opto, data_cc_arch$subject, data_cc_arch$rt_diff, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Archaerhodopsin", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
    lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=10", cex=2.2)
mtext(text="C- Cost-Cost Decision-Making", side=3, line=3, at=1.26, cex=2.8, font=2)
boxplot(rt_diff~opto, data = data_bb_sham, ylim=c(-4,4), axes=FALSE, xlab="", 
    ylab="", cex.lab=2.6, font.lab=2, frame.plot=FALSE, legend=FALSE, lty=1, lwd=6)
interaction.plot(data_bb_sham$opto, data_bb_sham$subject, data_bb_sham$rt_diff, type="b", 
    legend=FALSE, lty=1, lwd=6, pch=16, cex=3, col=brewer.pal(9, "Set1"), add=TRUE)
title(main="Control Virus", cex.main=2.2, line=-2, font.main=1)
axis(side = 1, at = 1:2, labels=c("Light OFF", "Light ON"), 
     lwd=5, mgp=c(3, 2, 0), cex.axis=2.2, font=2)
axis(side=2, at=seq(-4,4,1), labels=labels_rt_diff, las=1, lwd=5, cex.axis=2.2, font=2)
text(2, -3.5, labels="n=12", cex=2.2)
dev.off()

# Change the data into wide format
data_noopto <- data.table(subset(tidy_data, opto=='0'))
data_opto <- data.table(subset(tidy_data, opto=='1'))

setnames(data_noopto, 5:13, paste0(names(data_noopto)[5:13], '_noopto'))
setnames(data_opto, 5:13, paste0(names(data_opto)[5:13], '_opto'))

data_wide <- join(data_noopto, data_opto[, c("opto", "group"):=NULL], 
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

# Export the results of all statistical tests
write.csv(posthoc_tests_output, './analysis/posthoc_tests_output.csv', row.names = FALSE)
write.csv(diff_tests_output, './analysis/diff_tests_output.csv', row.names = FALSE)

sink('./analysis/anovas-output.txt')
print('Benefit-Benefit Decision-Making (Performance):')
print(anovas[[1]])
print(('Cost-Cost Decision-Making (Performance):'))
print(anovas[[2]])
print(('Cost-Benefit Decision-Making (Performance):'))
print(anovas[[3]])
print(('Benefit-Benefit Decision-Making (Reaction Time Difference):'))
print(anovas_rt_diff[[1]])
print(('Cost-Cost Decision-Making (Reaction Time Difference):'))
print(anovas_rt_diff[[2]])
print(('Cost-Benefit Decision-Making (Reaction Time Difference):'))
print(anovas_rt_diff[[3]])
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
