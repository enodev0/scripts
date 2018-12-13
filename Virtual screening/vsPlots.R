# Generate plots of virtual screening results
#
#
#

library(tidyverse)
#library(ggsignif)
library(reshape2)

df <- read_csv(file = "ReplicationStats-3-repeats--2018-09-21.csv")
fd <- read_csv(file="RawReplicationData-3-repeats--2018-09-21.csv", col_names=FALSE)

jpeg("plot.jpg")
ggplot(df, aes(x=df$'Ligands', y=df$'mean-dG')) +
geom_line(color="black") +
geom_point(size=4, colour="black") +
geom_errorbar(aes(ymin=df$'mean-dG'-df$'stdev-dG', ymax=df$'mean-dG'+df$'stdev-dG', width=0.2)) +
labs(title = "Binding Affinity Comparison",
    subtitle="#replicates=5") +
xlab("Ligands") + ylab("dG (kCal/mol)") + theme_light() + theme(text = element_text(size=20))
# geom_hline(yintercept = -6.0)
dev.off()

dfa <- fd %>% melt
dfb <- dfa[,-2]
colnames(dfb) <- c("Ligands","dG")

jpeg("plot4.jpg")
# Explicitly defining the factor levels allows us to re-arrange the x-axis at will for easier comparisons
# dfb$Ligands <- factor(dfb$Ligands, levels = c("CDB2", "R439", "R134", "R287", "R419", "R83"))
ggplot(data=dfb, aes(x=Ligands, y=dG)) + geom_boxplot(fill="#0eb7f4") +
geom_jitter(width=0.09, size=0.75) +
stat_summary(fun.y=mean, colour="black", geom="point", shape=3, size=5,show.legend = TRUE) +
# geom_signif(comparisons=list(c("CDB2", "R439"), c("R134", "R439"), c("R134", "R287"), c("R287", "R419"),
#                             c("R419", "R83")),
#             test="t.test", textsize = 4,
#             map_signif_level=FALSE, step_increase=c(0.0)) +
# geom_hline(yintercept = -6.0, linetype="dashed") +
labs(
    title="Binding affinity comparisons",
    subtitle="No. of replicates per ligand = 5"
) +
xlab("Ligand") + ylab("Predicted binding affinity (kcal/mol)") + theme_light() +
theme(text = element_text(size=20))

dev.off()

colnames(fd) <- c()

jpeg("file3.jpg")
par(mfrow=c(3,2))

for (i in 1:nrow(fd))
{
    fd[i,-1] %>% as.numeric %>% density %>% plot(main=fd[i,1], col="blue",
                  xlab = "delta G (kcal/mole)", xlim=c(-7.9, -5.0)) +
                  abline(v=-5.3, col="red") +
                  rug(jitter(as.numeric(fd[i,-1])), col="blue")
}
dev.off()

sr <- read_csv("screenres-isoform1-1-repeats--2018-09-21.csv", col_names = TRUE)

jpeg("file1.jpg")
qplot(sr$'mean-dG', geom="histogram", binwidth=0.1) + geom_vline(xintercept = -5.3, col="red") +
geom_vline(xintercept = -7.0, col="blue")+
ggtitle("Isoform 1 binding energy histogram")+
xlab("Predicted binding affinity (kcal/mole)")+
ylab("Frequency")
dev.off()

jpeg("file2.jpg")
hist(sr$'mean-dG',col="#0AC8DF", border="white", breaks=30, main="Binding affinity", xlab="dG (kCal/mole)")
abline(v=c(-5.25), col="red", lwd=3, lty=2)
dev.off()
