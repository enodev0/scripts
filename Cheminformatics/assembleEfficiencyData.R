
library(tidyverse)

iso1_dG <- read_csv(file = "screenres-isoform1-Modified.csv", col_names=TRUE)
iso2_dG <- read_csv(file = "screenres-isoform2-Modified.csv", col_names=TRUE)

head(iso1_dG)
head(iso2_dG)

tpsa_mw <- read_csv(file = "tpsa-mw-screenlib.csv", col_names=TRUE)

colnames(tpsa_mw)[1] <- "Ligands"

iso1_dG <- iso1_dG[,-c(3,4)]

iso2_dG <- iso2_dG[,-c(3,4)]

i1class <- replicate(length(iso1_dG$Ligands), "1")
i2class <- replicate(length(iso2_dG$Ligands), "2")

df1 <- inner_join(iso1_dG, tpsa_mw)
df2 <- inner_join(iso2_dG, tpsa_mw)

colnames(df1)[2] <- "dG"
colnames(df2)[2] <- "dG"

df1 <- cbind(df1, i1class)
df2 <- cbind(df2, i2class)

head(df1)
head(df2)

colnames(df1)[5] <- "Isoform"
colnames(df2)[5] <- "Isoform"

head(df1)
head(df2)

led <- rbind(df1,df2)

write_csv(led, path = "LigandEfficiencyRawData.csv", col_names=TRUE)

sessionInfo()
