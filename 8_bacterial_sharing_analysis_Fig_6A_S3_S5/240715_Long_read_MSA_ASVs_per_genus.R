################################################################################
##### Per-genus multiple sequence alignments of ASVs 
### Author(s): Asier Fernández / Johanne E. Spreckels
### Last updated: 15th July, 2024
################################################################################

#****************
# Load libraries
#****************
library(stringr)
library(seqinr)
library(msa)

#********************
# Defining functions
#********************
#A. MSA_result: generate a MSA of the ASVs of selected genera using ClustalOmega
# It outputs a FASTA file with ASVs sequences, a FASTA file with the alignment and a PDF with the alignment (including logos)
# output_dir : full path to directory where one folder will be stored for each genus with its results
# description: information about sequencing technology and pairs of samples to compare
MSA_result <- function(ASVs, names_ASVs, output_dir, description) {
  for (i in 1:length(names_ASVs)) {
    genus <- names_ASVs [i] #name of the genus
    dir <- paste0(output_dir ,genus)
    dir.create(dir) #Create a directory for the genus
    setwd (dir)  #Move to the directory
    genus_ASVs <- ASVs[grep (genus, names(ASVs))] #select ASVs of the specific genus
    filename <- paste0(paste(genus,description, sep = "_"), ".fa")
    write.fasta(genus_ASVs, names(genus_ASVs), filename) #generate a FASTA formatted file
    MSA <-  msa(filename, method = "ClustalOmega", type = "DNA", order = "input") #do the MSA
    ASV_names_for_plot <- names(genus_ASVs) #get names of the ASVs sequences to be displayed later in the alignment and the tree
    names(MSA@unmasked) <- ASV_names_for_plot
    try({msaPrettyPrint(MSA, output="pdf", file = paste(genus, description, "MSA.pdf", sep= "_"), #get a PDF with the alignment
                   showNames ="left", shadingMode="similar",
                   shadingColors="blues", showLogo="top", logoColors="rasmol",
                   showLegend=T, askForOverwrite=FALSE)})
    msa_align <- msaConvert(MSA, type="seqinr::alignment") #convert the alignment for later processing with seqinr
    msa_align$nam <- ASV_names_for_plot #add names for the plot
    filename_aln <- paste0(paste(genus,description, sep = "_"), "_aln.fa")
    write.fasta(as.list(msa_align$seq), names(genus_ASVs), filename_aln) #write alignment fasta (as ALN file is empty)
    }
}

# Set working directory
setwd("~/Desktop/PhD/Projects/Johanne Milk paper/")

# Read metadata, ASV counts and taxonomy
metadata <- read.delim("230321_sample_information_full-length_16S_test_n94.csv", sep = ",") 
ASV_counts <- read.delim("23Jan316_ASV_counts.tsv") 
ASV_taxonomy <- read.delim("23Jan316_taxonomy_GTDB.tsv") 

######################################################################
# Genus-level MSAs of ASVs in maternal milk -baby fecal samples
######################################################################

#################################################
#A. 16S - Mother Milk vs Baby Feces 
#################################################
setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/Mother_Milk-Baby_Feces/")

#*****************
# A1.Subset samples
#*****************
metadata_milk_feces_infant <- metadata[grepl("pilot", metadata$sample_set) &
                       grepl("PSK", metadata$sample_type) &
                       (grepl("feces", metadata$sample_type) | grepl("milk", metadata$sample_type)) &
                       !grepl("mother-feces", metadata$sample_type), ]

ASV_counts <- ASV_counts[rownames(ASV_counts) %in% metadata_milk_feces_infant$sample_ID, ]
ASV_counts <- ASV_counts[match(metadata_milk_feces_infant$sample_ID, rownames(ASV_counts)), ]

# Merge metadata and counts tables
ASV_counts$sample_ID <- rownames(ASV_counts)

# Merge metadata with ASV_counts by sample_ID
metadata_milk_feces_infant <- merge(metadata_milk_feces_infant, ASV_counts, by = "sample_ID", all.x = TRUE)

# Add origin and type to metadata
metadata_milk_feces_infant$origin <- ifelse(grepl("milk", metadata_milk_feces_infant$sample_type), "milk", 
                          ifelse(grepl("feces", metadata_milk_feces_infant$sample_type), "feces", NA))
metadata_milk_feces_infant$type <- ifelse(grepl("baby", metadata_milk_feces_infant$sample_type), "infant", "mother")
metadata_milk_feces_infant <- metadata_milk_feces_infant[, c("sample_ID", "origin", "type", setdiff(names(metadata_milk_feces_infant), c("sample_ID", "origin", "type")))]


#***********************************************************************
# A2.Get the genera (and their ASVs) present in at least 2 paired samples
#***********************************************************************
# For each sample select those ASVs that are present
# Use the sample_ID + lowest-level taxonomic assignment as ASV name
# Concatenate all sequences in a list
sample_ASVs <- metadata_milk_feces_infant[,17:ncol(metadata_milk_feces_infant)]
rownames(sample_ASVs) <- paste(metadata_milk_feces_infant$origin, metadata_milk_feces_infant$type, metadata_milk_feces_infant$pair_ID, sep = "_")
ASV_final_list <- list()

for (i in 1:nrow(sample_ASVs)) {
  sample_name <- rownames(sample_ASVs)[i]
  ASV <- colnames(sample_ASVs)[which(sample_ASVs[i,] != 0)] #get ASV name
  ASV_tax <- ASV_taxonomy[ASV_taxonomy$name %in% ASV, "Genus"]
  ASV_name_full <- paste(sample_name, ASV_tax, sep = "_")
  ASV_seq <- ASV_taxonomy[ASV_taxonomy$name %in% ASV, "seq"]
  ASV_list <- ASV_seq
  names(ASV_list) <- ASV_name_full
  ASV_final_list <- append(ASV_final_list, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list)

# Create a dataframe with the counts of each genus for each participant 
names_genera <- gsub(".*g__", "", grep("g__", names(ASV_final_list), value = T))
names_genera <- unique(names_genera) # 70 unique genera
genera_counts_milk_feces_infant <- data.frame(matrix(ncol = length(rownames(sample_ASVs)) + 1, nrow = length(names_genera)))
colnames(genera_counts_milk_feces_infant) <- c("genera", rownames(sample_ASVs))
genera_counts_milk_feces_infant$genera <- names_genera

for (i in 2:ncol(genera_counts_milk_feces_infant)) {
  microbes <- names(ASV_final_list)[grep(colnames(genera_counts_milk_feces_infant)[i], names(ASV_final_list))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  genera_counts_milk_feces_infant[i] <- counts$Freq[match(genera_counts_milk_feces_infant$genera,counts$genera)]
}
genera_counts_milk_feces_infant[is.na(genera_counts_milk_feces_infant)] <- 0 #set NAs to 0

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
genera_in_pairs_milk_feces_infant <- list()
IDs <- unique(metadata_milk_feces_infant$pair_ID)
for (i in IDs) {
  index <- which(IDs == i)
  pattern <- paste0("_",i, "$")
  pair_counts <- genera_counts_milk_feces_infant[, grep(pattern, colnames(genera_counts_milk_feces_infant))]
  genera_in_pairs_milk_feces_infant [[index]] <- genera_counts_milk_feces_infant$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair 
}

genera_in_pairs_milk_feces_infant <- genera_in_pairs_milk_feces_infant[!sapply(genera_in_pairs_milk_feces_infant,is.null)] # removes null elements
names(genera_in_pairs_milk_feces_infant) <- IDs

# Get names of the shared genera present in 2 or more pairs
shared_genera_milk_feces_infant <- names(which(table(unlist(genera_in_pairs_milk_feces_infant)) > 1))

# Get all ASVs of selected genera (212)
names(ASV_final_list) <- make.unique(names(ASV_final_list), sep = '_') # distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matching_ASVs_milk_feces_infant <- ASV_final_list[unique (grep(paste(shared_genera_milk_feces_infant,collapse="|"), # get all ASVs belonging to shared genera
                                             names(ASV_final_list)))]

#***********************************************************************
# A3.Get the MSA for each genus that is detected among paired individuals
#***********************************************************************
MSA_result(matching_ASVs_milk_feces_infant,shared_genera_milk_feces_infant, "~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/Mother_Milk-Baby_Feces/", "16S_mother_milk_baby_feces")



#################################################
#B. 16S Mother Milk vs Mother Feces 
#################################################
setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/Mother_Milk-Mother_Feces/")

#*****************
# B1.Subset samples
#*****************
metadata_milk_feces_mother <- metadata[grepl("pilot", metadata$sample_set) &
                                         grepl("PSK", metadata$sample_type) &
                                         (grepl("feces", metadata$sample_type) | grepl("milk", metadata$sample_type)) &
                                         !grepl("baby-feces", metadata$sample_type), ]

ASV_counts <- ASV_counts[rownames(ASV_counts) %in% metadata_milk_feces_mother$sample_ID, ]
ASV_counts <- ASV_counts[match(metadata_milk_feces_mother$sample_ID, rownames(ASV_counts)), ]

# Merge metadata and counts tables
ASV_counts$sample_ID <- rownames(ASV_counts)

# Merge metadata with ASV_counts by sample_ID
metadata_milk_feces_mother <- merge(metadata_milk_feces_mother, ASV_counts, by = "sample_ID", all.x = TRUE)

# Add origin and type to metadata
metadata_milk_feces_mother$origin <- ifelse(grepl("milk", metadata_milk_feces_mother$sample_type), "milk", 
                          ifelse(grepl("feces", metadata_milk_feces_mother$sample_type), "feces", NA))
metadata_milk_feces_mother$type <- ifelse(grepl("baby", metadata_milk_feces_mother$sample_type), "infant", "mother")
metadata_milk_feces_mother <- metadata_milk_feces_mother[, c("sample_ID", "origin", "type", setdiff(names(metadata_milk_feces_mother), c("sample_ID", "origin", "type")))]


#***********************************************************************
# A2.Get the genera (and their ASVs) present in at least 2 paired samples
#***********************************************************************
# For each sample select those ASVs that are present
# Use the sample_ID + lowest-level taxonomic assignment as ASV name
# Concatenate all sequences in a list
sample_ASVs <- metadata_milk_feces_mother[,17:ncol(metadata_milk_feces_mother)]
rownames(sample_ASVs) <- paste(metadata_milk_feces_mother$origin, metadata_milk_feces_mother$type, metadata_milk_feces_mother$pair_ID, sep = "_")
ASV_final_list <- list()

for (i in 1:nrow(sample_ASVs)) {
  sample_name <- rownames(sample_ASVs)[i]
  ASV <- colnames(sample_ASVs)[which(sample_ASVs[i,] != 0)] #get ASV name
  ASV_tax <- ASV_taxonomy[ASV_taxonomy$name %in% ASV, "Genus"]
  ASV_name_full <- paste(sample_name, ASV_tax, sep = "_")
  ASV_seq <- ASV_taxonomy[ASV_taxonomy$name %in% ASV, "seq"]
  ASV_list <- ASV_seq
  names(ASV_list) <- ASV_name_full
  ASV_final_list <- append(ASV_final_list, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list)

# Create a dataframe with the counts of each genus for each participant 
names_genera <- gsub(".*g__", "", grep("g__", names(ASV_final_list), value = T))
names_genera <- unique(names_genera) # 81 unique genera
genera_counts_milk_feces_mother <- data.frame(matrix(ncol = length(rownames(sample_ASVs)) + 1, nrow = length(names_genera)))
colnames(genera_counts_milk_feces_mother) <- c("genera", rownames(sample_ASVs))
genera_counts_milk_feces_mother$genera <- names_genera

for (i in 2:ncol(genera_counts_milk_feces_mother)) {
  microbes <- names(ASV_final_list)[grep(colnames(genera_counts_milk_feces_mother)[i], names(ASV_final_list))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  genera_counts_milk_feces_mother[i] <- counts$Freq[match(genera_counts_milk_feces_mother$genera,counts$genera)]
}
genera_counts_milk_feces_mother[is.na(genera_counts_milk_feces_mother)] <- 0 #set NAs to 0

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
genera_in_pairs_milk_feces_mother <- list()
IDs <- unique(metadata_milk_feces_mother$pair_ID)
for (i in IDs) {
  index <- which(IDs == i)
  pattern <- paste0("_",i, "$")
  pair_counts <- genera_counts_milk_feces_mother[, grep(pattern, colnames(genera_counts_milk_feces_mother))]
  genera_in_pairs_milk_feces_mother [[index]] <- genera_counts_milk_feces_mother$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair 
}

genera_in_pairs_milk_feces_mother <- genera_in_pairs_milk_feces_mother[!sapply(genera_in_pairs_milk_feces_mother,is.null)] # removes null elements
names(genera_in_pairs_milk_feces_mother) <- IDs

# Get names of the shared genera present in 2 or more pairs
shared_genera_milk_feces_mother <- names(which(table(unlist(genera_in_pairs_milk_feces_mother)) > 1))

# Get all ASVs of selected genera (120)
names(ASV_final_list) <- make.unique(names(ASV_final_list), sep = '_') # distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matching_ASVs_milk_feces_mother <- ASV_final_list[unique (grep(paste(shared_genera_milk_feces_mother,collapse="|"), # get all ASVs belonging to shared genera
                                             names(ASV_final_list)))]

#***********************************************************************
# A3.Get the MSA for each genus that is detected among paired individuals
#***********************************************************************
MSA_result(matching_ASVs_milk_feces_mother,shared_genera_milk_feces_mother, "~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/Mother_Milk-Mother_Feces//", "16S_mother_milk_mother_feces")


########################
#C. Save results 
########################

setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/Mother_Milk-Baby_Feces/")

# 16S Milk-Infant Feces results
write.table(metadata_milk_feces_infant,"Metadata_milk_feces_infant.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(genera_counts_milk_feces_infant,"Genera_counts.txt", sep = "\t", row.names = F, quote = FALSE) 
saveRDS(genera_in_pairs_milk_feces_infant, "Genera_counts_pairs.rds")
saveRDS(matching_ASVs_milk_feces_infant, "ASVs_shared_genera.rds")

setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/Mother_Milk-Mother_Feces/")

# 16S Milk- Maternal Feces results
write.table(metadata_milk_feces_mother,"Metadata_milk_feces_mother.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(genera_counts_milk_feces_mother,"Genera_counts.txt", sep = "\t", row.names = F, quote = FALSE) 
saveRDS(genera_in_pairs_milk_feces_mother, "Genera_counts_pairs.rds")
saveRDS(matching_ASVs_milk_feces_mother, "ASVs_shared_genera.rds")

