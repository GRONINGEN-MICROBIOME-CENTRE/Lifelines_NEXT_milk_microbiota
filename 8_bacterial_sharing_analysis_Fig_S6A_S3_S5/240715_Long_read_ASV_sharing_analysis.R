################################################################################
##### ASVs sharing analysis 
### Author(s): Asier Fernández / Johanne E. Spreckels
### Last updated: 15th July, 2024
################################################################################

#****************
# Load libraries
#****************
library(stringr)

#********************
# Defining functions
#********************
# Function to estimate the number of related and unrelated distances == or != 0.
# It also estimates the number of pairs in which PSE are detected
# SNP_distances: Dataframe with SNP distances of the ASVs from the shared genera 
# ASVs: list of sequences of ASVs belonging to the specific genus (e.g. matching_ASVs_genus_16S23S)
exact_PSE_analysis <- function(SNP_distances, ASVs) {
  #Define vectors to store results
  counts_related_below <- 0
  counts_related_above <-0
  counts_unrelated_below<-0
  counts_unrelated_above <-0
  pairs_with_PSE<-c()
  unrelated_comp_with_PSE <- list()
  unrelated_comp_without_PSE <- list()
  matching_ASVs_genus <- ASVs[names(ASVs) %in% rownames(SNP_distances)] #select the ASVs from the specific genus
  # Iterate over the distance matrix
  for (i in 1:nrow(SNP_distances)) {
    for (j in 1:ncol(SNP_distances)) {
      if(!is.na(SNP_distances[i,j])) {
        ASV1_fullname <- rownames(SNP_distances)[i] #get full ASV same
        ASV1_sample_origin <- str_extract(ASV1_fullname, "mother|infant") #will be mother/infant
        ASV1_pair <- str_extract(ASV1_fullname, "(?<=pair_)\\d+") #get pair_ID
        ASV1_seq_length <- str_length(matching_ASVs_genus[ASV1_fullname]) #get sequence length
        ASV2_fullname <- colnames(SNP_distances) [j] # get full ASV name
        ASV2_sample_origin <- str_extract(ASV2_fullname, "mother|infant") 
        ASV2_pair <- str_extract(ASV2_fullname, "(?<=pair_)\\d+")
        ASV2_seq_length <- str_length(matching_ASVs_genus[ASV2_fullname])
        if (SNP_distances[i,j] == 0) { #if the distance is below our threshold (=0)
          if (ASV1_sample_origin != ASV2_sample_origin) { #exclude comparisons from ASVs from the same sample origin
            if (ASV1_fullname != ASV2_fullname) { # exclude the same ASVs (obviously their distance is 0)
              if (ASV1_pair == ASV2_pair) { # if they are related (same pair_ID)
                pairs_with_PSE<-c(pairs_with_PSE, ASV1_pair) #add name of the pair
                counts_related_below <- counts_related_below + 1 #add 1 count to the related group
              } else { # if they are unrelated (different pair_ID)
                counts_unrelated_below <- counts_unrelated_below + 1 #add 1 count to the unrelated group
                unrelated_comp_with_PSE <-c(unrelated_comp_with_PSE, list(paste(paste(ASV1_sample_origin,ASV1_pair, sep ="_"), 
                                                                                paste(ASV2_sample_origin,ASV2_pair, sep ="_") , sep="_"))) #keep the unrelated pairs with PSE
                unrelated_comp_with_PSE <-c(unrelated_comp_with_PSE, list(paste(paste(ASV2_sample_origin,ASV2_pair, sep ="_"), 
                                                                                paste(ASV1_sample_origin,ASV1_pair, sep ="_") , sep="_"))) # all combinations
              }
            }
          }
        } else { #if the distance is above our threshold
          if (ASV1_sample_origin != ASV2_sample_origin) { #exclude comparisons from ASVs from the same sample origin
            if (ASV1_fullname != ASV2_fullname) { # exclude the same ASVs (obviously their distance is 0)
              if (ASV1_pair == ASV2_pair) { # if they are related (same pair_ID)
                counts_related_above <- counts_related_above + 1 #add 1 count to the related group
              } else { # if they are unrelated (different pair_ID)
                counts_unrelated_above <- counts_unrelated_above + 1 #add 1 count to the unrelated group
                unrelated_comp_without_PSE <-c(unrelated_comp_without_PSE, list(paste(paste(ASV1_sample_origin,ASV1_pair, sep ="_"), 
                                                                                      paste(ASV2_sample_origin,ASV2_pair, sep ="_") , sep="_"))) #keep the unrelated pairs without PSE
                unrelated_comp_without_PSE <-c(unrelated_comp_without_PSE, list(paste(paste(ASV2_sample_origin,ASV2_pair, sep ="_"), 
                                                                                      paste(ASV1_sample_origin,ASV1_pair, sep ="_") , sep="_"))) # all combinations
              }
            }
          }
        }
      }
    }
  }
  # Create contingency tables 
  cont_table <- data.frame(rbind(c(counts_related_above, counts_related_below),
                                 c(counts_unrelated_above, counts_unrelated_below)))
  colnames(cont_table) <- c("Above_threshold", "Below_threshold")
  rownames(cont_table) <- c("Related_distances", "Unrelated_distances")
  
  return(list("Contingency_table"= cont_table, "Pairs_with_PSEs" = pairs_with_PSE, 
              "Unrelated_pairs_with_PSE" = unrelated_comp_with_PSE, 
              "Unrelated_pairs_without_PSE" = unrelated_comp_without_PSE))
}


exact_PSE_analysis_same_origin <- function(SNP_distances, ASVs) {
  #Define vectors to store results
  counts_related_below <- 0
  counts_related_above <-0
  counts_unrelated_below<-0
  counts_unrelated_above <-0
  pairs_with_PSE<-c()
  unrelated_comp_with_PSE <- list()
  unrelated_comp_without_PSE <- list()
  matching_ASVs_genus <- ASVs[names(ASVs) %in% rownames(SNP_distances)] #select the ASVs from the specific genus
  # Iterate over the distance matrix
  for (i in 1:nrow(SNP_distances)) {
    for (j in 1:ncol(SNP_distances)) {
      if(!is.na(SNP_distances[i,j])) {
        ASV1_fullname <- rownames(SNP_distances)[i] #get full ASV same
        ASV1_sample_type <- str_extract(ASV1_fullname, "milk|feces") 
        ASV1_pair <- str_extract(ASV1_fullname, "(?<=pair_)\\d+") #get pair_ID
        ASV1_seq_length <- str_length(matching_ASVs_genus[ASV1_fullname]) #get sequence length
        ASV2_fullname <- colnames(SNP_distances) [j] # get full ASV name
        ASV2_sample_type <- str_extract(ASV2_fullname, "milk|feces") 
        ASV2_pair <- str_extract(ASV2_fullname, "(?<=pair_)\\d+")
        ASV2_seq_length <- str_length(matching_ASVs_genus[ASV2_fullname])
        if (SNP_distances[i,j] == 0) { #if the distance is below our threshold (=0)
          if (ASV1_sample_type != ASV2_sample_type) { #exclude comparisons from ASVs from the same sample type
            if (ASV1_fullname != ASV2_fullname) { # exclude the same ASVs (obviously their distance is 0)
              if (ASV1_pair == ASV2_pair) { # if they are related (same pair_ID)
                pairs_with_PSE<-c(pairs_with_PSE, ASV1_pair) #add name of the pair
                counts_related_below <- counts_related_below + 1 #add 1 count to the related group
              } else { # if they are unrelated (different pair_ID)
                counts_unrelated_below <- counts_unrelated_below + 1 #add 1 count to the unrelated group
                unrelated_comp_with_PSE <-c(unrelated_comp_with_PSE, list(paste(paste(ASV1_sample_type,ASV1_pair, sep ="_"), 
                                                                                paste(ASV2_sample_type,ASV2_pair, sep ="_") , sep="_"))) #keep the unrelated pairs with PSE
                unrelated_comp_with_PSE <-c(unrelated_comp_with_PSE, list(paste(paste(ASV2_sample_type,ASV2_pair, sep ="_"), 
                                                                                paste(ASV1_sample_type,ASV1_pair, sep ="_") , sep="_"))) # all combinations
              }
            }
          }
        } else { #if the distance is above our threshold
          if (ASV1_sample_type != ASV2_sample_type) { #exclude comparisons from ASVs from the same sample type
            if (ASV1_fullname != ASV2_fullname) { # exclude the same ASVs (obviously their distance is 0)
              if (ASV1_pair == ASV2_pair) { # if they are related (same pair_ID)
                counts_related_above <- counts_related_above + 1 #add 1 count to the related group
              } else { # if they are unrelated (different pair_ID)
                counts_unrelated_above <- counts_unrelated_above + 1 #add 1 count to the unrelated group
                unrelated_comp_without_PSE <-c(unrelated_comp_without_PSE, list(paste(paste(ASV1_sample_type,ASV1_pair, sep ="_"), 
                                                                                      paste(ASV2_sample_type,ASV2_pair, sep ="_") , sep="_"))) #keep the unrelated pairs without PSE
                unrelated_comp_without_PSE <-c(unrelated_comp_without_PSE, list(paste(paste(ASV2_sample_type,ASV2_pair, sep ="_"), 
                                                                                      paste(ASV1_sample_type,ASV1_pair, sep ="_") , sep="_"))) # all combinations
              }
            }
          }
        }
      }
    }
  }
  # Create contingency tables 
  cont_table <- data.frame(rbind(c(counts_related_above, counts_related_below),
                                 c(counts_unrelated_above, counts_unrelated_below)))
  colnames(cont_table) <- c("Above_threshold", "Below_threshold")
  rownames(cont_table) <- c("Related_distances", "Unrelated_distances")
  
  return(list("Contingency_table"= cont_table, "Pairs_with_PSEs" = pairs_with_PSE, 
              "Unrelated_pairs_with_PSE" = unrelated_comp_with_PSE, 
              "Unrelated_pairs_without_PSE" = unrelated_comp_without_PSE))
}



# Set working directory
setwd("~/Desktop/PhD/Projects/Johanne Milk paper/")

# Read metadata and previous results files
metadata_milk_feces_infant <- read.delim("RESULTS/Mother_Milk-Baby_Feces/Metadata_milk_feces_infant.txt") 
metadata_milk_feces_mother <- read.delim("RESULTS/Mother_Milk-Mother_Feces/Metadata_milk_feces_mother.txt") 

# Milk-Infant feces results
genera_counts_milk_feces_infant <- read.delim("RESULTS/Mother_Milk-Baby_Feces/Genera_counts.txt", header = T, check.names = F) 
genera_in_pairs_milk_feces_infant <- readRDS("RESULTS/Mother_Milk-Baby_Feces/Genera_counts_pairs.rds")
matching_ASVs_milk_feces_infant <- readRDS("RESULTS/Mother_Milk-Baby_Feces/ASVs_shared_genera.rds")

# Milk- Maternal feces results
genera_counts_milk_feces_mother <- read.delim("RESULTS/Mother_Milk-Mother_Feces/Genera_counts.txt", header = T, check.names = F) 
genera_in_pairs_milk_feces_mother <- readRDS("RESULTS/Mother_Milk-Mother_Feces/Genera_counts_pairs.rds")
matching_ASVs_milk_feces_mother <- readRDS("RESULTS/Mother_Milk-Mother_Feces/ASVs_shared_genera.rds")

#################################################################################
# SNP distances analysis in ASVs of shared genera in maternal-baby oral samples
#################################################################################

# Input: Distance matrices created by **snp-dists** using the MSAs generated in the script MSA_ASVs_per_genus.R

#################################################
#A. 16S - Mother Milk vs Baby Feces 
#################################################

setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/Mother_Milk-Baby_Feces/")

#***************************************
# Load the SNP distance matrix
#***************************************
# List the distance matrices
filenames_milk_feces_infant <- list.files(pattern="*.tab", recursive = TRUE)

# Create list of data frame names excluding extension
names_shared_genera <-sub("/.*", "", filenames_milk_feces_infant)
names_milk_feces_infant <-paste(names_shared_genera , "milk_feces_infant", sep = "_")

# Load all files (add 1st column as rownames)
SNP_distances_milk_feces_infant <- lapply(filenames_milk_feces_infant, function(file) read.delim(file, check.names = FALSE))
names(SNP_distances_milk_feces_infant) <- names_milk_feces_infant
for (i in 1:length(SNP_distances_milk_feces_infant)){
  rownames(SNP_distances_milk_feces_infant[[i]]) <- SNP_distances_milk_feces_infant[[i]][,1]
  SNP_distances_milk_feces_infant[[i]][,1] <- NULL
}

# Get half of the distance matrices (avoid duplicated values)
for (i in 1:length(SNP_distances_milk_feces_infant)) {
  SNP_distances_milk_feces_infant[[i]][upper.tri(SNP_distances_milk_feces_infant[[i]])] <- NA
}


#***************************************************************
# Estimate the PSEs in related vs unrelated mother-baby pairs
#***************************************************************
# Generate all combinations of pair IDs 
mothers <- paste("mother",1:14, sep="_")
babies <- paste("infant", 1:14, sep="_")

all_pairs <- expand.grid(mothers, babies)
all_pairs_combined <- paste(all_pairs$Var1, all_pairs$Var2, sep="_")
related_pairs <- paste("mother", 1:14, "infant",1:14, sep="_")
unrelated_pairs <- all_pairs_combined[!all_pairs_combined %in% related_pairs]
related_pairs <- sub(".+_", "", related_pairs)
related_pairs <- paste0("pair_", related_pairs)

# Generate the empty final data table
final_PSE_data_milk_feces_infant <- data.frame(matrix(nrow=length(all_pairs_combined), ncol=length(names_milk_feces_infant)))
colnames(final_PSE_data_milk_feces_infant) <- names_milk_feces_infant
rownames(final_PSE_data_milk_feces_infant) <- c(related_pairs, unrelated_pairs)
number_related_pairs <- length(unique(metadata_milk_feces_infant$pair_ID))
final_PSE_data_milk_feces_infant$Relatedness <- c(rep("Related", number_related_pairs), 
                                rep("Unrelated", length(all_pairs_combined) - number_related_pairs))

# Number of pairs in which each genus is detected
all_pairs <- table(unlist(genera_in_pairs_milk_feces_infant))
total_pairs <- all_pairs[names(all_pairs) %in% names_shared_genera]

# Pairs in which each genus is detected
Streptococcus_in_pairs <- names(grep("Streptococcus", genera_in_pairs_milk_feces_infant, value = T))
Staphylococcus_in_pairs <- names(grep("Staphylococcus", genera_in_pairs_milk_feces_infant, value = T))
Bifidobacterium_in_pairs <- names(grep("Bifidobacterium", genera_in_pairs_milk_feces_infant, value = T))
Bacteroides_in_pairs <- names(grep("Bacteroides", genera_in_pairs_milk_feces_infant, value = T))
Veillonella_in_pairs <- names(grep("Veillonella", genera_in_pairs_milk_feces_infant, value = T))
Rothia_in_pairs <- names(grep("Rothia", genera_in_pairs_milk_feces_infant, value = T))

# Get the contigency table and the pairs in which we observe PSEs. Add results to final_PSE_data_milk_feces_infant
PSE_Streptococcus_milk_feces_infant <- exact_PSE_analysis(SNP_distances_milk_feces_infant$Streptococcus_milk_feces_infant, matching_ASVs_milk_feces_infant)
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% Streptococcus_in_pairs,"Streptococcus_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(paste0("pair_",PSE_Streptococcus_milk_feces_infant$Pairs_with_PSEs)),"Streptococcus_milk_feces_infant" ] <- "PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Streptococcus_milk_feces_infant$Unrelated_pairs_without_PSE),"Streptococcus_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Streptococcus_milk_feces_infant$Unrelated_pairs_with_PSE),"Streptococcus_milk_feces_infant" ] <- "PSE"

PSE_Staphylococcus_milk_feces_infant <- exact_PSE_analysis(SNP_distances_milk_feces_infant$Staphylococcus_milk_feces_infant, matching_ASVs_milk_feces_infant)
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% Staphylococcus_in_pairs,"Staphylococcus_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(paste0("pair_",PSE_Staphylococcus_milk_feces_infant$Pairs_with_PSEs)),"Staphylococcus_milk_feces_infant" ] <- "PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Staphylococcus_milk_feces_infant$Unrelated_pairs_without_PSE),"Staphylococcus_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Staphylococcus_milk_feces_infant$Unrelated_pairs_with_PSE),"Staphylococcus_milk_feces_infant" ] <- "PSE"

PSE_Bifidobacterium_milk_feces_infant <- exact_PSE_analysis(SNP_distances_milk_feces_infant$Bifidobacterium_milk_feces_infant, matching_ASVs_milk_feces_infant)
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% Bifidobacterium_in_pairs,"Bifidobacterium_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(paste0("pair_",PSE_Bifidobacterium_milk_feces_infant$Pairs_with_PSEs)),"Bifidobacterium_milk_feces_infant" ] <- "PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Bifidobacterium_milk_feces_infant$Unrelated_pairs_without_PSE),"Bifidobacterium_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Bifidobacterium_milk_feces_infant$Unrelated_pairs_with_PSE),"Bifidobacterium_milk_feces_infant" ] <- "PSE"

PSE_Bacteroides_milk_feces_infant <- exact_PSE_analysis(SNP_distances_milk_feces_infant$Bacteroides_milk_feces_infant, matching_ASVs_milk_feces_infant)
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% Bacteroides_in_pairs,"Bacteroides_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(paste0("pair_",PSE_Bacteroides_milk_feces_infant$Pairs_with_PSEs)),"Bacteroides_milk_feces_infant" ] <- "PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Bacteroides_milk_feces_infant$Unrelated_pairs_without_PSE),"Bacteroides_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Bacteroides_milk_feces_infant$Unrelated_pairs_with_PSE),"Bacteroides_milk_feces_infant" ] <- "PSE"

PSE_Veillonella_milk_feces_infant <- exact_PSE_analysis(SNP_distances_milk_feces_infant$Veillonella_milk_feces_infant, matching_ASVs_milk_feces_infant)
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% Veillonella_in_pairs,"Veillonella_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(paste0("pair_",PSE_Veillonella_milk_feces_infant$Pairs_with_PSEs)),"Veillonella_milk_feces_infant" ] <- "PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Veillonella_milk_feces_infant$Unrelated_pairs_without_PSE),"Veillonella_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Veillonella_milk_feces_infant$Unrelated_pairs_with_PSE),"Veillonella_milk_feces_infant" ] <- "PSE"

PSE_Rothia_milk_feces_infant <- exact_PSE_analysis(SNP_distances_milk_feces_infant$Rothia_milk_feces_infant, matching_ASVs_milk_feces_infant)
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% Rothia_in_pairs,"Rothia_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(paste0("pair_",PSE_Rothia_milk_feces_infant$Pairs_with_PSEs)),"Rothia_milk_feces_infant" ] <- "PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Rothia_milk_feces_infant$Unrelated_pairs_without_PSE),"Rothia_milk_feces_infant" ] <- "no_PSE"
final_PSE_data_milk_feces_infant[rownames(final_PSE_data_milk_feces_infant) %in% unique(PSE_Rothia_milk_feces_infant$Unrelated_pairs_with_PSE),"Rothia_milk_feces_infant" ] <- "PSE"

# Save contingency tables in list
cont_tables_milk_feces_infant <- list(PSE_Streptococcus_milk_feces_infant$Contingency_table,
                         PSE_Staphylococcus_milk_feces_infant$Contingency_table,
                         PSE_Bifidobacterium_milk_feces_infant$Contingency_table,
                         PSE_Bacteroides_milk_feces_infant$Contingency_table,
                         PSE_Veillonella_milk_feces_infant$Contingency_table,
                         PSE_Rothia_milk_feces_infant$Contingency_table)


names (cont_tables_milk_feces_infant) <- c("Streptococcus", "Staphylococcus", "Bifidobacterium",
                              "Bacteroides", "Veillonella", "Rothia")



#################################################
#B. 16S - Mother Milk vs Mother Feces
#################################################
setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/Mother_Milk-Mother_Feces/")

#***************************************
# Load the SNP distance matrix
#***************************************
# List the distance matrices
filenames_milk_feces_mother <- list.files(pattern="*.tab", recursive = TRUE)

# Create list of data frame names excluding extension
names_shared_genera <-sub("/.*", "", filenames_milk_feces_mother)
names_milk_feces_mother <-paste(names_shared_genera , "milk_feces_mother", sep = "_")

# Load all files (add 1st column as rownames)
SNP_distances_milk_feces_mother <- lapply(filenames_milk_feces_mother, function(file) read.delim(file, check.names = FALSE))
names(SNP_distances_milk_feces_mother) <- names_milk_feces_mother
for (i in 1:length(SNP_distances_milk_feces_mother)){
  rownames(SNP_distances_milk_feces_mother[[i]]) <- SNP_distances_milk_feces_mother[[i]][,1]
  SNP_distances_milk_feces_mother[[i]][,1] <- NULL
}

# Get half of the distance matrices (avoid duplicated values)
for (i in 1:length(SNP_distances_milk_feces_mother)) {
  SNP_distances_milk_feces_mother[[i]][upper.tri(SNP_distances_milk_feces_mother[[i]])] <- NA
}


#***************************************************************
# Estimate the PSEs in related vs unrelated mother-baby pairs
#***************************************************************
# Generate all combinations of pair IDs 
mothers <- paste("milk",1:14, sep="_")
babies <- paste("feces", 1:14, sep="_")

all_pairs <- expand.grid(mothers, babies)
all_pairs_combined <- paste(all_pairs$Var1, all_pairs$Var2, sep="_")
related_pairs <- paste("milk", 1:14, "feces",1:14, sep="_")
unrelated_pairs <- all_pairs_combined[!all_pairs_combined %in% related_pairs]
related_pairs <- sub(".+_", "", related_pairs)
related_pairs <- paste0("pair_", related_pairs)

# Generate the empty final data table
final_PSE_data_milk_feces_mother <- data.frame(matrix(nrow=length(all_pairs_combined), ncol=length(names_milk_feces_mother)))
colnames(final_PSE_data_milk_feces_mother) <- names_milk_feces_mother
rownames(final_PSE_data_milk_feces_mother) <- c(related_pairs, unrelated_pairs)
number_related_pairs <- length(unique(metadata_milk_feces_mother$pair_ID))
final_PSE_data_milk_feces_mother$Relatedness <- c(rep("Related", number_related_pairs), 
                                rep("Unrelated", length(all_pairs_combined) - number_related_pairs))

# Number of pairs in which each genus is detected
all_pairs <- table(unlist(genera_in_pairs_milk_feces_mother))
total_pairs <- all_pairs[names(all_pairs) %in% names_shared_genera]

# Pairs in which each genus is detected
Phocaeicola_in_pairs <- names(grep("Phocaeicola", genera_in_pairs_milk_feces_mother, value = T))
Bifidobacterium_in_pairs <- names(grep("Bifidobacterium", genera_in_pairs_milk_feces_mother, value = T))
Bacteroides_in_pairs <- names(grep("Bacteroides", genera_in_pairs_milk_feces_mother, value = T))
Streptococcus_in_pairs <- names(grep("Streptococcus", genera_in_pairs_milk_feces_mother, value = T))

# Get the contigency table and the pairs in which we observe PSEs. Add results to final_PSE_data_milk_feces_mother
PSE_Phocaeicola_milk_feces_mother <- exact_PSE_analysis_same_origin(SNP_distances_milk_feces_mother$Phocaeicola_milk_feces_mother, matching_ASVs_milk_feces_mother)
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% Phocaeicola_in_pairs,"Phocaeicola_milk_feces_mother" ] <- "no_PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(paste0("pair_",PSE_Phocaeicola_milk_feces_mother$Pairs_with_PSEs)),"Phocaeicola_milk_feces_mother" ] <- "PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(PSE_Phocaeicola_milk_feces_mother$Unrelated_pairs_without_PSE),"Phocaeicola_milk_feces_mother" ] <- "no_PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(PSE_Phocaeicola_milk_feces_mother$Unrelated_pairs_with_PSE),"Phocaeicola_milk_feces_mother" ] <- "PSE"

PSE_Bifidobacterium_milk_feces_mother <- exact_PSE_analysis_same_origin(SNP_distances_milk_feces_mother$Bifidobacterium_milk_feces_mother, matching_ASVs_milk_feces_mother)
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% Bifidobacterium_in_pairs,"Bifidobacterium_milk_feces_mother" ] <- "no_PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(paste0("pair_",PSE_Bifidobacterium_milk_feces_mother$Pairs_with_PSEs)),"Bifidobacterium_milk_feces_mother" ] <- "PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(PSE_Bifidobacterium_milk_feces_mother$Unrelated_pairs_without_PSE),"Bifidobacterium_milk_feces_mother" ] <- "no_PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(PSE_Bifidobacterium_milk_feces_mother$Unrelated_pairs_with_PSE),"Bifidobacterium_milk_feces_mother" ] <- "PSE"

PSE_Bacteroides_milk_feces_mother <- exact_PSE_analysis_same_origin(SNP_distances_milk_feces_mother$Bacteroides_milk_feces_mother, matching_ASVs_milk_feces_mother)
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% Bacteroides_in_pairs,"Bacteroides_milk_feces_mother" ] <- "no_PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(paste0("pair_",PSE_Bacteroides_milk_feces_mother$Pairs_with_PSEs)),"Bacteroides_milk_feces_mother" ] <- "PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(PSE_Bacteroides_milk_feces_mother$Unrelated_pairs_without_PSE),"Bacteroides_milk_feces_mother" ] <- "no_PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(PSE_Bacteroides_milk_feces_mother$Unrelated_pairs_with_PSE),"Bacteroides_milk_feces_mother" ] <- "PSE"

PSE_Streptococcus_milk_feces_mother <- exact_PSE_analysis_same_origin(SNP_distances_milk_feces_mother$Streptococcus_milk_feces_mother, matching_ASVs_milk_feces_mother)
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% Streptococcus_in_pairs,"Streptococcus_milk_feces_mother" ] <- "no_PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(paste0("pair_",PSE_Streptococcus_milk_feces_mother$Pairs_with_PSEs)),"Streptococcus_milk_feces_mother" ] <- "PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(PSE_Streptococcus_milk_feces_mother$Unrelated_pairs_without_PSE),"Streptococcus_milk_feces_mother" ] <- "no_PSE"
final_PSE_data_milk_feces_mother[rownames(final_PSE_data_milk_feces_mother) %in% unique(PSE_Streptococcus_milk_feces_mother$Unrelated_pairs_with_PSE),"Streptococcus_milk_feces_mother" ] <- "PSE"

# Save contingency tables in list
cont_tables_milk_feces_mother <- list(PSE_Phocaeicola_milk_feces_mother$Contingency_table,
                         PSE_Bifidobacterium_milk_feces_mother$Contingency_table,
                         PSE_Bacteroides_milk_feces_mother$Contingency_table,
                         PSE_Streptococcus_milk_feces_mother$Contingency_table)


names (cont_tables_milk_feces_mother) <- c("Phocaeicola","Bifidobacterium","Bacteroides",
                               "Streptococcus")


########################
#C. Save results 
########################
setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/")

# Milk-Infant Feces results
write.table(final_PSE_data_milk_feces_infant,"Mother_Milk-Baby_Feces/PSE_table_milk_infant_feces.txt", sep = "\t", row.names = T, quote = FALSE) 
saveRDS(cont_tables_milk_feces_infant, "Mother_Milk-Baby_Feces//Contingency_tables_milk_infant_feces.rds")

# Maternal-Infant feces results
write.table(final_PSE_data_milk_feces_mother,"Mother_Milk-Mother_Feces/PSE_table_milk_maternal_feces.txt", sep = "\t", row.names = T, quote = FALSE) 
saveRDS(cont_tables_milk_feces_mother, "Mother_Milk-Mother_Feces/Contingency_tables_milk_maternal_feces.rds")

