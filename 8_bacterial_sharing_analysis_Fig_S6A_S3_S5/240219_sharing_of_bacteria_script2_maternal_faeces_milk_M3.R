################################################################################################################
### INVESTIGATE SHARING OF BACTERIA BETWEEN MATERNAL FAECES AND MILK - SCRIPT 2 - M3 ###########################
################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 1. FUNCTION FOR ASV SHARING ANALYSIS BETWEEN SAMPLE TYPES/INDIVIDUALS PER TIME POINT
# 
# 2. LOAD FILES AND PREPARE DATA
# 2.1 IMPORT FILES FROM SCRIPT 1 AND MSA
# 2.2 ESTIMATE POTENTIAL SHARING EVENTS (PSEs) IN RELATED VS UNRELATED MOTHER-INFANT PAIRS
# 2.3 SAVE RESULTS


##### =========================== 1. FUNCTION FOR ASV SHARING ANALYSIS BETWEEN SAMPLE TYPES/INDIVIDUALS PER TIME POINT =========================== #####

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
  print("Creating empty vectors")
  counts_related_below <- 0
  counts_related_above <-0
  counts_unrelated_below<-0
  counts_unrelated_above <-0
  pairs_with_PSE<-c()
  unrelated_comp_with_PSE <- list()
  unrelated_comp_without_PSE <- list()
  matching_ASVs_genus <- ASVs[names(ASVs) %in% rownames(SNP_distances)] #select the ASVs from the specific genus
  
  # Iterate over the distance matrix
  print("Start iterating over distance matrix")
  for (i in 1:nrow(SNP_distances)) {
    for (j in 1:ncol(SNP_distances)) {
      print(paste0("Starting analysis for ", i, " ", rownames(SNP_distances)[i], " and ", j, " ", colnames(SNP_distances)[j]))
      
      if(!is.na(SNP_distances[i,j])) {
        print("Extracting information from rownames for ASV1")
        ASV1_fullname <- rownames(SNP_distances)[i] #get full ASV same
        ASV1_sampleorigin <- gsub("_FAM.*", "",  ASV1_fullname) #will be mother_human_milk/mother_faeces
        ASV1_pair <- gsub(".*_FAM", "FAM", gsub("\\.g__.*","",  ASV1_fullname)) #i #get pair_ID
        # ASV1_seq_length <- str_length(matching_ASVs_M3[ASV1_fullname]) #get sequence length

        print("Extracting information from colnames for ASV2")
        ASV2_fullname <- colnames(SNP_distances) [j] # get full ASV name
        ASV2_pair <- gsub(".*_FAM", "FAM", gsub("\\.g__.*","",  ASV2_fullname)) # get pair ID
        ASV2_sampleorigin <- gsub("_FAM.*", "",  ASV2_fullname)
        # ASV2_seq_length <- str_length(matching_ASVs_M3[ASV2_fullname]) #get sequence length
        # mean_seq_length <- mean(c(ASV1_seq_length, ASV2_seq_length)) #estimate the mean length of both ASVs
        
        print("Comparing ASVs")
        if (SNP_distances[i,j] == 0) { #if the distance is below our threshold
          print("A - ASVs have SNP distance = 0")
          
          if (ASV1_sampleorigin != ASV2_sampleorigin) { #exclude comparisons from ASVs from the same sample origin
            print("A - The ASVs that were detected in the different sample types")
            if (ASV1_fullname != ASV2_fullname) { # exclude the same ASVs (obviously their distance is 0)
              print("A - The ASVs have different names")
              
              if (ASV1_pair == ASV2_pair) { # if they are related (same pair_ID)
                print("A - Identical ASVs were detected in a related pair")
                pairs_with_PSE <- c(pairs_with_PSE, ASV1_pair) #add name of the pair
                counts_related_below <- counts_related_below + 1 #add 1 count to the related group
                
              } else { # if they are unrelated (different pair_ID)
                print("A - Identical ASVs were detected in an unrelated pair")
                counts_unrelated_below <- counts_unrelated_below + 1 #add 1 count to the unrelated group
                unrelated_comp_with_PSE <-c(unrelated_comp_with_PSE, list(paste(paste(ASV1_sampleorigin,ASV1_pair, sep ="_"), 
                                                                                paste(ASV2_sampleorigin,ASV2_pair, sep ="_") , sep="_"))) #keep the unrelated pairs with PSE
                unrelated_comp_with_PSE <-c(unrelated_comp_with_PSE, list(paste(paste(ASV2_sampleorigin,ASV2_pair, sep ="_"), 
                                                                                paste(ASV1_sampleorigin,ASV1_pair, sep ="_") , sep="_"))) # all combinations
              }
            }else{print("A - Excluded the ASVs because they had the same name")}
          }else{print("A - Excluded the ASVs because they were detected in the same type of samples")}
        } else { #if the distance is above our threshold
          print("B - ASVs have SNP distance > 0")
          
          if (ASV1_sampleorigin != ASV2_sampleorigin) { #exclude comparisons from ASVs from the same sample origin
            print("B - The ASVs that were detected in the different sample types")
            if (ASV1_fullname != ASV2_fullname) { # exclude the same ASVs (obviously their distance is 0)
              print("B - The ASVs have different names")
              
              if (ASV1_pair == ASV2_pair) { # if they are related (same pair_ID)
                print("B - ASVs were detected in a related pair")
                counts_related_above <- counts_related_above + 1 #add 1 count to the related group
                
              } else { # if they are unrelated (different pair_ID)
                print("B - ASVs were detected in an unrelated pair")
                counts_unrelated_above <- counts_unrelated_above + 1 #add 1 count to the unrelated group
                unrelated_comp_without_PSE <-c(unrelated_comp_without_PSE, list(paste(paste(ASV1_sampleorigin,ASV1_pair, sep ="_"), 
                                                                                      paste(ASV2_sampleorigin,ASV2_pair, sep ="_") , sep="_"))) #keep the unrelated pairs without PSE
                unrelated_comp_without_PSE <-c(unrelated_comp_without_PSE, list(paste(paste(ASV2_sampleorigin,ASV2_pair, sep ="_"), 
                                                                                      paste(ASV1_sampleorigin,ASV1_pair, sep ="_") , sep="_"))) # all combinations
              }
            }else{print("B - Excluded the ASVs because they had the same name")}
          }else{print("B - Excluded the ASVs because they were detected in the same type of samples")}
        }
      }
    }
  }
  
  # Create contingency tables
  print("Creating contingency tables")
  cont_table <- data.frame(rbind(c(counts_related_above, counts_related_below),
                                 c(counts_unrelated_above, counts_unrelated_below)))
  colnames(cont_table) <- c("Above_threshold", "Below_threshold")
  rownames(cont_table) <- c("Related_distances", "Unrelated_distances")
  
  # unrelated_comp_with_PSE <- gsub("mother", "mother_human_milk", unrelated_comp_with_PSE)
  # unrelated_comp_with_PSE <- gsub("infant", "infant_faeces", unrelated_comp_with_PSE)
  # 
  # unrelated_comp_without_PSE <- gsub("mother", "mother_human_milk", unrelated_comp_without_PSE)
  # unrelated_comp_without_PSE <- gsub("infant", "infant_faeces", unrelated_comp_without_PSE)
  
  print("Returning output")
  return(list("Contingency_table"= cont_table, "Pairs_with_PSEs" = pairs_with_PSE, 
              "Unrelated_pairs_with_PSE" = unrelated_comp_with_PSE, 
              "Unrelated_pairs_without_PSE" = unrelated_comp_without_PSE))
}


##### =========================== 2. LOAD FILES AND PREPARE DATA =========================== #####

### ===== 2.1 IMPORT FILES FROM SCRIPT 1 AND MSA ===== ###

# Set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/genera_in_at_least_5_pairs/allfa/out/")

# Read metadata and previous results files
## import files, phenotypes+microbiota data on genus level
data <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/231121_milk_composition_paper_sel_phenotypes_sel_microbiota_data_ASVs_n1515.rds")
#1515x7282

## exclude ASVs from unclassified genera
ASVcf <- data[,-grep("g__NA", colnames(data))] #1515x5917

# M3
genera_counts_M3 <- read.delim("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/M3_maternal_faeces_milk_genera_counts.txt", header = T, check.names = F) 
genera_in_pairs_M3 <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/M3_maternal_faeces_milk_genera_counts_pairs.rds")
matching_ASVs_M3 <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/M3_maternal_faeces_milk_ASVs_shared_genera.rds")


# Import: Distance matrices created by **snp-dists** using the MSAs generated in the script MSA_mother_milk_baby_oral.R
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/genera_in_at_least_5_pairs/allfa/out/")

# Load and List the SNP distance matrices
filenames_M3 <- list.files(pattern="*.tab", recursive = TRUE)

# Create list of data frame names excluding extension
names_shared_genera_M3 <- sub("_M3.*", "", filenames_M3)
names_M3 <- paste(names_shared_genera_M3 , "M3", sep = "_")

# Load all files (add 1st column as rownames)
SNP_distances_M3 <- lapply(filenames_M3, function(file) read.delim(file, check.names = FALSE))
names(SNP_distances_M3) <- names_M3
for (i in 1:length(SNP_distances_M3)){
  rownames(SNP_distances_M3[[i]]) <- SNP_distances_M3[[i]][,1]
  SNP_distances_M3[[i]][,1] <- NULL
}

# Get half of the distance matrices (avoid duplicated values)
for (i in 1:length(SNP_distances_M3)) {
  SNP_distances_M3[[i]][upper.tri(SNP_distances_M3[[i]])] <- NA
}


### ===== 2.2 ESTIMATE POTENTIAL SHARING EVENTS (PSEs) IN THE SAME VS DIFFERENT MOTHER-MOTHER PAIRS ===== ###

# Generate all combinations of pair IDs 
mothers_milk_M3 <- colnames(genera_counts_M3)[grep("mother_human_milk", colnames(genera_counts_M3))]
mothers_faeces_M3 <- colnames(genera_counts_M3)[grep("mother_faeces", colnames(genera_counts_M3))]

all_pairs_M3 <- expand.grid(mothers_milk_M3, mothers_faeces_M3)
all_pairs_combined_M3 <- paste(all_pairs_M3$Var1, all_pairs_M3$Var2, sep="_")

all_pairs_tmp_M3 <- all_pairs_M3
all_pairs_tmp_M3$Var1 <- gsub("mother_human_milk_", "", all_pairs_tmp_M3$Var1)
all_pairs_tmp_M3$Var2 <- gsub("mother_faeces_", "", all_pairs_tmp_M3$Var2)

related_famIDs_M3 <- c()
unrelated_famIDs_M3 <- c()

for (i in 1:nrow(all_pairs_tmp_M3)){
  
  if (all_pairs_tmp_M3$Var1[i] == all_pairs_tmp_M3$Var2[i]){
    related_famIDs_M3 <- c(related_famIDs_M3, paste0("mother_human_milk_", all_pairs_tmp_M3$Var1[i], "_mother_faeces_", all_pairs_tmp_M3$Var2[i]))
    
  }else{
    unrelated_famIDs_M3 <- c(unrelated_famIDs_M3, paste0("mother_human_milk_", all_pairs_tmp_M3$Var1[i], "_mother_faeces_", all_pairs_tmp_M3$Var2[i]))
    
  }
}

related_pairs_M3 <- unique(related_famIDs_M3)
unrelated_pairs_M3 <- unique(unrelated_famIDs_M3)
related_pairs_M3 <- gsub("_mother_faeces.*","", gsub("mother_human_milk_", "", related_pairs_M3))


# Generate the empty final data table
final_PSE_data_M3 <- data.frame(matrix(nrow=length(all_pairs_combined_M3), ncol=length(names_M3)))
colnames(final_PSE_data_M3) <- names_M3
rownames(final_PSE_data_M3) <- c(related_pairs_M3, unrelated_pairs_M3)
number_related_pairs_M3 <- length(unique(all_pairs_tmp_M3$Var1))
final_PSE_data_M3$Relatedness <- c(rep("Related", number_related_pairs_M3),
                                   rep("Unrelated", length(all_pairs_combined_M3) - number_related_pairs_M3))

# Number of pairs in which each genus is detected
all_pairs_M3 <- table(unlist(genera_in_pairs_M3))
total_pairs_M3 <- all_pairs_M3[names(all_pairs_M3) %in% names_shared_genera_M3]

# Pairs in which each genus is detected
Bacteroides_in_pairs_M3 <- names(grep("Bacteroides", genera_in_pairs_M3, value = T))
Bifidobacterium_in_pairs_M3 <- names(grep("Bifidobacterium", genera_in_pairs_M3, value = T))
Blautia_in_pairs_M3 <- names(grep("Blautia", genera_in_pairs_M3, value = T))
Clostridium_sensu_stricto_in_pairs_M3 <- names(grep("Clostridium_sensu_stricto", genera_in_pairs_M3, value = T))
Escherichia_Shigella_in_pairs_M3 <- names(grep("Escherichia_Shigella", genera_in_pairs_M3, value = T))
Faecalibacterium_in_pairs_M3 <- names(grep("Faecalibacterium", genera_in_pairs_M3, value = T))
Haemophilus_in_pairs_M3 <- names(grep("Haemophilus", genera_in_pairs_M3, value = T))
Prevotella_in_pairs_M3 <- names(grep("Prevotella", genera_in_pairs_M3, value = T))
Streptococcus_in_pairs_M3 <- names(grep("Streptococcus", genera_in_pairs_M3, value = T))
Veillonella_in_pairs_M3 <- names(grep("Veillonella", genera_in_pairs_M3, value = T))


## Get the contingency table and the pairs in which we observe PSEs. Add results to final_PSE_data_M3
PSE_Bacteroides_M3 <- exact_PSE_analysis(SNP_distances_M3$Bacteroides_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Bacteroides_in_pairs_M3,"Bacteroides_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Bacteroides_M3$Pairs_with_PSEs),"Bacteroides_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Bacteroides_M3$Unrelated_pairs_without_PSE),"Bacteroides_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Bacteroides_M3$Unrelated_pairs_with_PSE),"Bacteroides_M3" ] <- "PSE"

PSE_Bifidobacterium_M3 <- exact_PSE_analysis(SNP_distances_M3$Bifidobacterium_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Bifidobacterium_in_pairs_M3,"Bifidobacterium_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Bifidobacterium_M3$Pairs_with_PSEs),"Bifidobacterium_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Bifidobacterium_M3$Unrelated_pairs_without_PSE),"Bifidobacterium_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Bifidobacterium_M3$Unrelated_pairs_with_PSE),"Bifidobacterium_M3" ] <- "PSE"

PSE_Blautia_M3 <- exact_PSE_analysis(SNP_distances_M3$Blautia_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Blautia_in_pairs_M3,"Blautia_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Blautia_M3$Pairs_with_PSEs),"Blautia_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Blautia_M3$Unrelated_pairs_without_PSE),"Blautia_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Blautia_M3$Unrelated_pairs_with_PSE),"Blautia_M3" ] <- "PSE"

PSE_Clostridium_sensu_stricto_M3 <- exact_PSE_analysis(SNP_distances_M3$Clostridium_sensu_stricto_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Clostridium_sensu_stricto_in_pairs_M3,"Clostridium_sensu_stricto_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Clostridium_sensu_stricto_M3$Pairs_with_PSEs),"Clostridium_sensu_stricto_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Clostridium_sensu_stricto_M3$Unrelated_pairs_without_PSE),"Clostridium_sensu_stricto_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Clostridium_sensu_stricto_M3$Unrelated_pairs_with_PSE),"Clostridium_sensu_stricto_M3" ] <- "PSE"

PSE_Escherichia_Shigella_M3 <- exact_PSE_analysis(SNP_distances_M3$Escherichia_Shigella_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Escherichia_Shigella_in_pairs_M3,"Escherichia_Shigella_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Escherichia_Shigella_M3$Pairs_with_PSEs),"Escherichia_Shigella_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Escherichia_Shigella_M3$Unrelated_pairs_without_PSE),"Escherichia_Shigella_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Escherichia_Shigella_M3$Unrelated_pairs_with_PSE),"Escherichia_Shigella_M3" ] <- "PSE"

PSE_Faecalibacterium_M3 <- exact_PSE_analysis(SNP_distances_M3$Faecalibacterium_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Faecalibacterium_in_pairs_M3,"Faecalibacterium_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Faecalibacterium_M3$Pairs_with_PSEs),"Faecalibacterium_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Faecalibacterium_M3$Unrelated_pairs_without_PSE),"Faecalibacterium_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Faecalibacterium_M3$Unrelated_pairs_with_PSE),"Faecalibacterium_M3" ] <- "PSE"

PSE_Haemophilus_M3 <- exact_PSE_analysis(SNP_distances_M3$Haemophilus_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Haemophilus_in_pairs_M3,"Haemophilus_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Haemophilus_M3$Pairs_with_PSEs),"Haemophilus_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Haemophilus_M3$Unrelated_pairs_without_PSE),"Haemophilus_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Haemophilus_M3$Unrelated_pairs_with_PSE),"Haemophilus_M3" ] <- "PSE"

PSE_Prevotella_M3 <- exact_PSE_analysis(SNP_distances_M3$Prevotella_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Prevotella_in_pairs_M3,"Prevotella_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Prevotella_M3$Pairs_with_PSEs),"Prevotella_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Prevotella_M3$Unrelated_pairs_without_PSE),"Prevotella_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Prevotella_M3$Unrelated_pairs_with_PSE),"Prevotella_M3" ] <- "PSE"

PSE_Streptococcus_M3 <- exact_PSE_analysis(SNP_distances_M3$Streptococcus_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Streptococcus_in_pairs_M3,"Streptococcus_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Streptococcus_M3$Pairs_with_PSEs),"Streptococcus_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Streptococcus_M3$Unrelated_pairs_without_PSE),"Streptococcus_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Streptococcus_M3$Unrelated_pairs_with_PSE),"Streptococcus_M3" ] <- "PSE"

PSE_Veillonella_M3 <- exact_PSE_analysis(SNP_distances_M3$Veillonella_M3, matching_ASVs_M3)
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% Veillonella_in_pairs_M3,"Veillonella_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Veillonella_M3$Pairs_with_PSEs),"Veillonella_M3" ] <- "PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Veillonella_M3$Unrelated_pairs_without_PSE),"Veillonella_M3" ] <- "no_PSE"
final_PSE_data_M3[rownames(final_PSE_data_M3) %in% unique(PSE_Veillonella_M3$Unrelated_pairs_with_PSE),"Veillonella_M3" ] <- "PSE"


# Save contingency tables in list
cont_tables_M3 <- list(PSE_Bacteroides_M3$Contingency_table,
                       PSE_Bifidobacterium_M3$Contingency_table,
                       PSE_Blautia_M3$Contingency_table,
                       PSE_Clostridium_sensu_stricto_M3$Contingency_table,
                       PSE_Escherichia_Shigella_M3$Contingency_table,
                       PSE_Faecalibacterium_M3$Contingency_table,
                       PSE_Haemophilus_M3$Contingency_table,
                       PSE_Prevotella_M3$Contingency_table,
                       PSE_Streptococcus_M3$Contingency_table,
                       PSE_Veillonella_M3$Contingency_table)

names(cont_tables_M3) <- names_shared_genera_M3


# Move relatedness column to the first position
relatedness_column <- final_PSE_data_M3$Relatedness
final_PSE_data_M3 <- final_PSE_data_M3[, !(names(final_PSE_data_M3) == "Relatedness")]
final_PSE_data_M3 <- data.frame(Relatedness = relatedness_column, final_PSE_data_M3)


### ===== 2.3 SAVE RESULTS ===== ###

setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/")

write.table(final_PSE_data_M3,"M3/240219_PSE_table_maternal_faeces_milk_M3.txt", sep = "\t", row.names = T, quote = FALSE) 

# save contingecy table
saveRDS(cont_tables_M3, "M3/240219_maternal_faeces_milk_contingency_tables_M3.rds")



