# Script to plot the molecular property distributions of a database of molecules.
#
#
#
library(tidyverse)

df <- read_csv("LigandProperties.csv", col_names=TRUE)

colnames(df)

g1 <- ggplot(df, aes(x=df$'Heavy.Atoms.Count')) + geom_histogram(fill="#0dc2dd", color="black", position="identity", binwidth=1) +
  labs(
    x="Heavy Atoms Count",
    y="Frequency") +
theme_light() + theme(text = element_text(size=10))

g2 <- ggplot(df, aes(x=df$'XLogP')) + geom_histogram(fill="#0dc2dd", color="black", position="identity", binwidth=0.1) +
  labs(
    x="XLogP",
    y="Frequency") +
theme_light() + theme(text = element_text(size=10))

g3 <- ggplot(df, aes(x=df$'Topological.Polar.Surface.Area')) + geom_histogram(fill="#0dc2dd", color="black", position="identity", binwidth=4) +
  labs(
    x="TPSA",
    y="Frequency") +
theme_light() + theme(text = element_text(size=10))

g4 <- ggplot(df, aes(x=df$'Rotatable.Bonds.Count')) + geom_histogram(fill="#0dc2dd", color="black", position="identity", binwidth=1) +
  labs(
    x="Rotatable Bonds",
    y="Frequency") +
theme_light() + theme(text = element_text(size=10))

g5 <- ggplot(df, aes(x=df$'Mannhold.LogP')) + geom_histogram(fill="#0dc2dd", color="black", position="identity", binwidth=0.15) +
  labs(
    x="Mannhold LogP",
    y="Frequency") +
theme_light() + theme(text = element_text(size=10))

g6 <- ggplot(df, aes(x=df$'Molar.Mass')) + geom_histogram(fill="#0dc2dd", color="black", position="identity", binwidth=7) +
  labs(
    x="Molar Mass",
    y="Frequency") +
theme_light() + theme(text = element_text(size=10))

par(mfrow=c(3,3))
g1
g2
g3
g4
g5
g6
