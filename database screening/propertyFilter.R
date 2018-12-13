# Filter a database of molecules based on the mean and SD of their descriptors.
#
#
#
#
library(tidyverse)
# library(doMC)
library(caret)
# library(pROC)
library(FSelector)
library(corrplot)
# library(e1071)

# Load CDK descriptors of Triacsin C analogs obtained from literature, and bioisosteres from scaffolod hopping
ana <- read_csv(file="SCAFFOLD_PLUS_ANALOGS.csv")
# Load CDK descriptors of all the molecules in ChEMBL 23 and DrugBank 5.0.1
cdk <- read_csv(file="cdk-trimmed-allProperties.csv")
cdk <- cdk[complete.cases(cdk), ] # ditch rows with NAs

head(ana)
colnames(ana)

head(cdk)
colnames(cdk)

# How many molecules are being contributed by each database?
table(cdk$Database)

# Keep columns which are common to both the dataframes

# Trim non-matching columns
kdc <- cdk[,c(colnames(ana[,-28]))]
dfa <- kdc[,-(ncol(kdc))]             # From databases
dfb <- ana[,-(c(28, ncol(ana)))]      # From literature + scaffold hop

colnames(dfa)

colnames(dfa) == colnames(dfb)

head(dfa)

head(dfb)

skewAnalogs <- apply(dfb, 2,skewness, na.rm=TRUE)

skewCDK <- apply(dfa, 2, skewness, na.rm=TRUE)

skewAnalogs

skewCDK

# The actual filter function to screen molecules from ChEMBL and DrugBank based on molecular properties

threshold_filter <- function(toBeFiltered, suppliesUandS, Z=1) { # dataframe without nominal columns

    # Check that the dataframe columns match, else halt execution
    if(is.element(FALSE, (colnames(toBeFiltered) == colnames(suppliesUandS)))) {
        stop("Dataframe columns mismatch")
    }
    avg <- c()
    std <- c()
    nam <- colnames(suppliesUandS)  # save column headers for later access
    suppliesUandS <- unname(suppliesUandS) # strip column headers
    for (m in 1:(ncol(suppliesUandS))) {
            #print(m)
            avg <- append(avg, lapply(as.list(suppliesUandS[,m]), mean)) # get avg of mth column
            std <- append(std, lapply(as.list(suppliesUandS[,m]), sd))   # get stdev of mth column
    }

    statusDF <- c("ColumnIndex", "#Moleules", "ColumnName", "Status")

    for (i in 1:(ncol(toBeFiltered))) {
            statusDF <- rbind(statusDF, c(i, nrow(toBeFiltered), colnames(toBeFiltered[i]), "before"))
            toBeFiltered <- filter(toBeFiltered, (toBeFiltered[i]) <= (avg[[i]]+(Z*std[[i]])) & (toBeFiltered[i]) >= (avg[[i]]-(Z*std[[i]])))
            statusDF <- rbind(statusDF, c(i, nrow(toBeFiltered), colnames(toBeFiltered[i]), "after"))
    }

    return(list(toBeFiltered, statusDF))
}

# Remove columns with zero variance
nearZeroVar(dfa)

colnames(dfa[,40])

gc()

nearZeroVar(dfb)

colnames(dfb[,c(40,41,42)])

colnames(dfa[,c(40,41,42)]) == colnames(dfb[,c(40,41,42)])

dfa <- dfa[, -(c(40,41,42))] # CDK
dfb <- dfb[, -(c(40,41,42))] # Analogs

colnames(dfa) == colnames(dfb)

colnames(dfa)

apply(dfa, 2, class)

apply(dfb, 2, class)

nrow(dfa)
ncol(dfa)

dim(dfb)

# Pearson's correlation matrix of the molecules from the database and analogs.

cor(dfb, method="pearson", use="complete.obs") %>% corrplot(method="color", order="hclust", tl.cex=0.5)

# CDK correlation

cor(dfa, method="pearson", use="complete.obs") %>% corrplot(method="color", order="hclust", tl.cex=0.5)

gc()

# Apply the filter on the database
filtered <- threshold_filter(toBeFiltered = dfa, suppliesUandS = dfb, Z=1.7)

nrow(filtered[[1]])
ncol(filtered[[1]])

colnames(filtered[[1]]) == colnames(dfb)

nrow(dfb)
ncol(dfb)

# Obtain the SMILES strings of the molecules which have been accepted by the threshold filter
extract <- semi_join(kdc, filtered[[1]]) # return all rows from kdc where there are matching values in filtered
# keeping just columns from kdc.

extract

nrow(extract)

# Write out the SMILES strings of the molecules which have passed the threshold filter
write_csv(extract[,"SMILES"], path="filteredSMILES-Z17.csv")
write_csv(extract, path="filteredSamples-Z17.csv")

sessionInfo()


