# This scripts generates a figure that summarizes the total number of contacts between 
# VM axon terminals in the uppermost 40 mum of prelimbic cortex and dendrites of pyramidal neurons 
# marked with MAP2, dendrites of corticostriatal neurons marked with RGS14 and layer 1 inhibitory 
# interneurons marked with 5HT3aR estimated using systematic random sampling.
# The estimated number of contacts is plotted as black dots, and the mean and standard error 
# of the mean are plotted in red. In addition, this script calculates and outputs the means
# and standard errors of the mean that were plotted in the figure

# Import required libraries
library(psych)
library(dplyr)
library(ggplot2)

# Format the data for plotting
data_long <- read.csv("stereology_summary_long.csv")
data_wide <- read.csv("stereology_summary_wide.csv")

stats = describe(data_wide)

# Generate the figure
png('./estimated_total.png', width = 800, height = 600)
plot(ggplot(data_long, aes(marker, estimate)) + geom_point(size = 6) + 
    theme(axis.title.x = element_blank(), axis.title.y = element_text(colour = "black",
    size = "28", face = "bold", margin = margin(t = 0, r = 20, b = 0, l = 0)),
    axis.text = element_text(colour = "black", size = "24")) + 
    ylab("Estimated Total") + scale_x_discrete(breaks=c("MAP2", "RGS14", "HT3aR"), 
    limits=c("MAP2", "RGS14", "HT3aR"), labels = c("Pyramidal\nNeurons (MAP2)", 
    "IT Neurons\n(RGS14)", "Layer 1 Inhibitory\nInterneurons (5HT3aR)")) + 
    scale_y_continuous(limits = c(0, 32000)) +
    stat_summary(fun.data_long=mean_se, geom="errorbar", color="red", width=0.3, size=2) + 
    stat_summary(fun.y=mean, geom="point", color="red", size=6))
dev.off()

# Output means and standard errors of the mean that were plotted
print(stats)
