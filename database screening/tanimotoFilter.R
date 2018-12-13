# This script goes through the .tar.gz files of PubChem DB
# and compares the fingerprints of molecules therein with that
# of a molecule of interest.
#
# Deletes processed archives in order to save HDD space. Make sure
# to have a backup of the same.
#
library(R.utils)
library(iterators)
library(rcdk)
library(rcdklibs)
library(fingerprint)


database <- "ChEMBL"
logfile <- "/home/subhashree/Documents/Somdeb/ChEMBL/LOG3.txt"

write(sprintf("Screening database: %s\n\n", database), logfile, append = TRUE)
screen_start <- Sys.time()



csv_loc <- "/home/subhashree/Documents/Somdeb/ChEMBL/SCREENDATA-CHEMBL-3.csv"
# Set the below two variables according to the environment
root_path <- "/home/subhashree/Documents/Somdeb/ChEMBL/DB/Part3/"
files <- Sys.glob("/home/subhashree/Documents/Somdeb/ChEMBL/DB/Part3/*.gz")
total_archives <- length(files)
# Log the details of skipped molecules
#skipped_loc <- "/home/bioinfo/Documents/Somdeb/Dissertation/DrugBank/Part3/SKIPPED2.csv"
#skipped <- cbind("SMILES", "CHEMBL_ID") # dataframe to contain all skipped molecules
#write.table(skipped, skipped_loc, sep = ",", col.names = T, append = T, row.names = F)

# Keep track of number of archives processed
archives_processed <<- 0 # Globally updatable, hence the <<-
# Keep track of the total number of molecules screened
total_screened <<- 0
# Keep track of the number of molecules skipped
total_skipped <<- 0
# Keep track of fatal exceptions
fatals_ignored <<- 0


# gades = Get All DEScriptors
gades <- function(molecule)
{
    dc <- get.desc.categories()

    dn_hybrid         <- get.desc.names(dc[1])
    dn_constitutional <- get.desc.names(dc[2])
    dn_topological    <- get.desc.names(dc[3])
    dn_electronic     <- get.desc.names(dc[4])
    dn_geometrical    <- get.desc.names(dc[5])

    eval_hyb <- eval.desc(molecule, dn_hybrid)
    eval_con <- eval.desc(molecule, dn_constitutional)
    eval_top <- eval.desc(molecule, dn_topological)
    eval_ele <- eval.desc(molecule, dn_electronic)
    eval_geo <- eval.desc(molecule, dn_geometrical)

    desc <- cbind(eval_hyb, eval_con, eval_top, eval_ele, eval_geo);

    return(desc)
}



# Load 2D structure of Triasin C for fingerprint comparison
tc2 <- load.molecules(c("/home/subhashree/Documents/Somdeb/ChEMBL/Tric2D.sdf"))
tcf <- get.fingerprint(tc2[[1]], type="extended")
# tanimoto distance between a supplied SDF and Triacsin C
deltatan <- function(mol)
{
    query <- get.fingerprint(mol, type="extended")
    tanimoto <- distance(tcf, query, method="tanimoto")

    return (tanimoto)
}
##### Ensure this chunk of code is not executed more than once per screening ###########
triacsin_desc <- gades(tc2[[1]])
#SMILES  <- c(get.property(tc2[[1]], "similarity"))
SMILES  <- get.smiles(tc2[[1]], flavor = smiles.flavors(c('Canonical')), smigen = NULL)
# IUPAC <-   c(get.property(tc2[[1]], "JCHEM_IUPAC"))
Tanimoto <- c(deltatan(tc2[[1]]))
triacsin_desc <- cbind(triacsin_desc,Tanimoto,SMILES)
# All future results will be rbinded to this one.
tric_final <- triacsin_desc[colSums(!is.na(triacsin_desc)) > 0] # keep only those descriptors for 2D SDFs
write.table(tric_final, csv_loc, sep = ",", col.names = T, append = T, row.names = F)
########################################################################################



# Worker function to process SDFs
process_file <- function(f, archives_processed) {
  # cat(sprintf("ENTERED FUNCTION\n\n"))
    # Current number of molecules screened
    current_screened <- 0
    # Current number of molecules skipped
    current_skipped <- 0

    dump_file_loc <- sprintf("/home/subhashree/Documents/Somdeb/ChEMBL/SkippedChemBl3/skippedFromArchive-%d", archives_processed)
    dump_mols <- c()

    # cat(sprintf("dump-loc: %s\n\n", dump_file_loc))
    # cat(sprintf("dump-size: %d\n\n", length(dump_mols)))

    # cat(sprintf("LOADING MOLECULES ... "))
    pubmols <- iload.molecules(f, type="sdf", skip=TRUE) # continue even if invalid molecule
    df <- c() # details of all processed molecules
    # skipped_df <- c() # details of all skipped molecules (IUPAC, SMILES, CID)
    # cat(sprintf(" DONE. \n\nENTERING LOOP ...\n\n"))
    while (hasNext(pubmols))
    {
      current_screened <- current_screened + 1
      total_screened  <<- total_screened + 1
      # Main processing loop
      start <- Sys.time()
      mol <- nextElem(pubmols)
      # cat(sprintf("GENERATING DESCRIPTORS ..."))
      descriptors <- try(gades(mol))

      # Some SDFs leads to errors in rcdk leading to crash of the program when running gades() above.
      # To avoid seeing glorious Java exception vomits, the above function is wrapped under a try block
      # and promptly ignored.

      if(inherits(descriptors, "try-error"))
      {
            cat(sprintf("FATAL: Recovering and continuing"))
            #error handling code, maybe just skip this iteration using next
            fatals_ignored <<- fatals_ignored + 1
            next # skip it.
      }

      # cat(sprintf(" DONE\n\n"))
      # IUPAC <-   c(get.property(tc2[[1]], "JCHEM_IUPAC"))
      # print(IUPAC)
      SMILES  <- get.smiles(mol, flavor = smiles.flavors(c('Canonical')), smigen = NULL)
      # print(SMILES)
      Tanimoto <- c(deltatan(mol))
      # print(Tanimoto)
      descriptors <- cbind(descriptors,Tanimoto,SMILES)
      final <- descriptors[colSums(!is.na(descriptors)) > 0] # keep only those descriptors for 2D SDFs

      if(ncol(tric_final) != ncol(final)) {
        # cat(sprintf("IN INNER SEGMENT\n\n"))
        current_skipped <- current_skipped + 1
        total_skipped <<- total_skipped + 1
        # RANK_ID <- c(get.property(mol, "result-rank"))
        # skipped_df <- rbind(skipped_df, cbind(SMILES, IUPAC)) # log details of skipped molecules
        dump_mols <- append(dump_mols, mol)
        next # skip this molecule.
      }
        # cat(sprintf("IN OUTER SEGMENT\n\n"))
        df <- rbind(df, final)
        stop <- Sys.time()
        latency <- stop-start # time to screen one molecule (string)
        # convert to float
        dlat <- as.double(latency)

        # cat(sprintf("PASS\n\n"))

        cat(sprintf("TMS: %d, TSK: %d, CMS: %d, CSK: %d, FATALS: %d\n",
                    total_screened,
                    total_skipped,
                    current_screened,
                    current_skipped,
                    fatals_ignored))
        cat(sprintf("Archive: %d of %d, latency: %0.2f seconds/molecule\n\n",
                    archives_processed,
                    total_archives,
                    dlat))
    }
    # cat(sprintf("dump_size(after screen): %d\n\n", length(dump_mols)))
    #write.table(skipped_df, skipped_loc, sep = ",", col.names = F, append = T, row.names = F)
    if (length(dump_mols) > 0) {
        write.molecules(dump_mols, dump_file_loc, together=FALSE, write.props=TRUE)
        dump_mols <- NULL
    }
    write.table(df, csv_loc, sep = ",", col.names = F, append = T, row.names = F)
}


# Main database processing loop to process archives
for (file in files)
{
    archives_processed <- archives_processed + 1
    cat(sprintf("Processing file: %d (%s)\n\n", archives_processed, basename(file)))
    cat(sprintf("Running gunzip ...\n\n "))
    gunzip(file, remove=FALSE) # don't remove parent archive

    working_file <- paste(root_path,tools::file_path_sans_ext(basename(file)),sep="")

    if (file.exists(working_file)) {
        cat(sprintf("%s exists\n\n", basename(file)))
    } else {
        cat(sprintf("WARNING: %s doesn't exist, continuing anyway.\n\n", working_file))
    }

    # cat(sprintf("HANDING OVER CONTROL\n\n"))
    cat(sprintf("Processing file ...\n\n"))

    process_file(working_file, archives_processed) # hand over to SDF processor function



    cat(sprintf("Deleting %s\n\n", tools::file_path_sans_ext(basename(file))))
    file.remove(working_file)
}

write(sprintf("\n\nTotal number of archives processed: %d, Total molecules screened: %d, Total molecules skipped: %d\n\n", total_archives, total_screened, total_skipped), logfile, append = TRUE)

screen_stop <- Sys.time()
screen_latency <- as.double(screen_stop - screen_start)
write(sprintf("\n\nTime taken to screen %d molecules from %d archives: %0.2f seconds\n\n",
              total_screened-total_skipped,
              total_archives,
              screen_latency), logfile, append = TRUE)
