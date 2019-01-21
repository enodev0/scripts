
library(tidyverse)

df <- read_csv("LigandEfficiencyRawData.csv", col_names=TRUE)

Isoform <- df$Isoform
Ligand <- df$Ligands
TPSA <- df$TPSA
MW <- df$MW
deltaG <- df$dG

# Ki in micromolars
Ki <- function(dG)
{
    return(exp((dG*1000)/(1.986*298.15))*(1/(10**(-6))))
}

# pKi
pKi <- function(Ki)
{
    # Converting Ki from uM/L
    # to M/L
    return(-log10(Ki/1000000))
}

# Ligand Efficiency
LE <- function()
{
    return(-(deltaG/HAC))
}

# Binding Efficiency Index
BEI <- function()
{
    # Convert MW to kDa
    kDa <- MW/1000
    return(pki/kDa)
}

# Surface Binding Efficiency Index
SEI <- function()
{
    # normalisze
    norm <- TPSA/100
    return(pki/norm)
}

ki <- round(Ki(deltaG), digits=3)

pki <- round(pKi(ki), digits=3)

pki %>% head

bei <- round(BEI(), digits=3)
sei <- round(SEI(), digits=3)

ki %>% head

pki %>% head

bei %>% head

sei %>% head

bindTable <- cbind(Ligand, Isoform, deltaG,MW,ki,pki,bei,sei)

colnames(bindTable) <- c("Ligand", "Isoform", "deltaG", "MW", "Ki", "pKi", "BEI", "SEI")

df <- data.frame(bindTable)

df %>% head

write_csv(df, path = "screenlib-bei-sei.csv", col_names=TRUE)

sessionInfo()
