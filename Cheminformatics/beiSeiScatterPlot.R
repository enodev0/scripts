
library(tidyverse)

df <- read_csv("screenlib-bei-sei.csv", col_names=TRUE)

head(df)

dfa <- filter(df, Isoform == 1)
dfb <- filter(df, Isoform == 2)

# ISOFORM-1
ggplot(dfa,aes(x=SEI, y=BEI)) + geom_point(alpha=0.3) +
coord_cartesian(xlim = c(0, 50), ylim = c(0,50)) +
geom_abline(intercept=0, slope=1) + theme_linedraw() + theme(text = element_text(size=30))

filter(dfa, BEI > 15 & SEI > 16 & SEI < 25)

# ISOFORM-2
ggplot(dfb,aes(x=SEI, y=BEI)) + geom_point(alpha=0.2) +
coord_cartesian(xlim = c(0, 50), ylim = c(0,50)) +
geom_abline(intercept=0, slope=1) + theme_linedraw() + theme(text = element_text(size=30))

filter(dfb, BEI > 16 & BEI < 22 & SEI > 17.5 & SEI < 22)

head(dfa)

# ISOFORM-1
ab1 <- data.frame(filter(dfa, Ligand %in% c("reduce-134", "reduce-439", "chembl-dbank-min-red_2", "triacsinC")))
ab1[1,1] <- "R134"
ab1[2,1] <- "R439"
ab1[3,1] <- "Triacsin C"
ab1[4,1] <- "CDB2"
ggplot(NULL,aes(x=SEI, y=BEI)) + geom_point(data=dfa,alpha=0.5, color="#f0c454") +
coord_cartesian(xlim = c(0, 50), ylim = c(0,50)) +
geom_point(data=ab1, aes(x=SEI,y=BEI, color=Ligand)) +
geom_abline(intercept=0, slope=1)+
theme_linedraw() +
theme(text = element_text(size=15))

# ISOFORM-2
ab2 <- data.frame(filter(dfb, Ligand %in% c("reduce-134", "reduce-38", "D1Tric", "triacsinC")))
ab2[1,1] <- "R134"
ab2[2,1] <- "MCD1"
ab2[3,1] <- "R38"
ab2[4,1] <- "Triacsin C"
ggplot(NULL,aes(x=SEI, y=BEI)) + geom_point(data=dfb,alpha=0.5, color="#f0c454") +
coord_cartesian(xlim = c(0, 50), ylim = c(0,50)) +
geom_point(data=ab2, aes(x=SEI,y=BEI, color=Ligand)) +
geom_abline(intercept=0, slope=1)+
theme_linedraw() +
theme(text = element_text(size=15))

# Isoform1
ab1 <- data.frame(filter(dfa, Ligand %in% c("reduce-134", "reduce-439", "chembl-dbank-min-red_2", "triacsinC")))
ab1[1,1] <- "R134"
ab1[2,1] <- "R439"
ab1[3,1] <- "Triacsin C"
ab1[4,1] <- "CDB2"
ggplot(NULL,aes(x=SEI, y=BEI)) + geom_point(data=dfa, alpha=0.3) + geom_point(data=ab1, color="red",size=2) +
coord_cartesian(xlim = c(0, 50), ylim = c(0,50)) +
geom_abline(intercept=0, slope=1) + theme_linedraw() + theme(text = element_text(size=30))

# ISOFORM-2
ab2 <- data.frame(filter(dfb, Ligand %in% c("reduce-134", "reduce-38", "D1Tric", "triacsinC")))
ab2[1,1] <- "R134"
ab2[2,1] <- "MCD1"
ab2[3,1] <- "R38"
ab2[4,1] <- "Triacsin C"
ggplot(NULL,aes(x=SEI, y=BEI)) + geom_point(data=dfb,alpha=0.2) + geom_point(data=ab2, color="red",size=2) +
coord_cartesian(xlim = c(0, 50), ylim = c(0,50)) +
geom_abline(intercept=0, slope=1) + theme_linedraw() + theme(text = element_text(size=30))

ab2

dfa_ratio <- (dfa$BEI/dfa$SEI)
dfb_ratio <- (dfb$BEI/dfb$SEI)

head(dfa_ratio)

efficiency_ratio <- (df$BEI/df$SEI)

head(efficiency_ratio)

df <- cbind(df, efficiency_ratio)

colnames(df)[9] <- "BSRatio"

head(df)

write_csv(df, path = "screenlib-bei-sei-withRatio.csv", col_names=TRUE)

filter(df, BSRatio <= 1.1 & BSRatio >= 0.9)

a <- filter(df, Isoform == 1)

head(a)

filter(a, grepl("chembl", Ligand))
