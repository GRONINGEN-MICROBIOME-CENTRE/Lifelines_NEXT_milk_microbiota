################################################################################################################
### INVESTIGATE SHARING OF BACTERIA BETWEEN MILK AND INFANT FAECES - SCRIPT 2 - M2 ############################
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
        ASV1_sampleorigin <- gsub("_.*", "",  ASV1_fullname) #will be mother/baby, except for oral samples that will be oral (we know its a baby then)
        ASV1_pair <- gsub(".*_FAM", "FAM", gsub("\\.g__.*","",  ASV1_fullname)) #i #get pair_ID
        # ASV1_seq_length <- str_length(matching_ASVs_M2[ASV1_fullname]) #get sequence length

        print("Extracting information from colnames for ASV2")
        ASV2_fullname <- colnames(SNP_distances) [j] # get full ASV name
        ASV2_pair <- gsub(".*_FAM", "FAM", gsub("\\.g__.*","",  ASV2_fullname)) # get pair ID
        ASV2_sampleorigin <- gsub("_.*", "",  ASV2_fullname)
        # ASV2_seq_length <- str_length(matching_ASVs_M2[ASV2_fullname]) #get sequence length
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
  
  unrelated_comp_with_PSE <- gsub("mother", "mother_human_milk", unrelated_comp_with_PSE)
  unrelated_comp_with_PSE <- gsub("infant", "infant_faeces", unrelated_comp_with_PSE)
  
  unrelated_comp_without_PSE <- gsub("mother", "mother_human_milk", unrelated_comp_without_PSE)
  unrelated_comp_without_PSE <- gsub("infant", "infant_faeces", unrelated_comp_without_PSE)
  
  print("Returning output")
  return(list("Contingency_table"= cont_table, "Pairs_with_PSEs" = pairs_with_PSE, 
              "Unrelated_pairs_with_PSE" = unrelated_comp_with_PSE, 
              "Unrelated_pairs_without_PSE" = unrelated_comp_without_PSE))
}


##### =========================== 2. LOAD FILES AND PREPARE DATA =========================== #####

### ===== 2.1 IMPORT FILES FROM SCRIPT 1 AND MSA ===== ###

# Set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/")

# Read metadata and previous results files
## import files, phenotypes+microbiota data on genus level
data <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/231121_milk_composition_paper_sel_phenotypes_sel_microbiota_data_ASVs_n1515.rds")
#1515x7282

## exclude ASVs from unclassified genera
ASVcf <- data[,-grep("g__NA", colnames(data))] #1515x5917

# M2
genera_counts_M2 <- read.delim("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_genera_counts.txt", header = T, check.names = F) 
genera_in_pairs_M2 <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_genera_counts_pairs.rds")
matching_ASVs_M2 <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_ASVs_shared_genera.rds")


# Import: Distance matrices created by **snp-dists** using the MSAs generated in the script MSA_mother_milk_baby_oral.R
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/SNP_distance_tables_2024_02_16/milk_infant_feces/selected_genera/M2/")

# Load and List the SNPdistance matrices
filenames_M2 <- list.files(pattern="*.tab", recursive = TRUE)

# Create list of data frame names excluding extension
names_shared_genera_M2 <- sub("_M2.*", "", filenames_M2)
names_M2 <- paste(names_shared_genera_M2 , "M2", sep = "_")

# Load all files (add 1st column as rownames)
SNP_distances_M2 <- lapply(filenames_M2, function(file) read.delim(file, check.names = FALSE))
names(SNP_distances_M2) <- names_M2
for (i in 1:length(SNP_distances_M2)){
  rownames(SNP_distances_M2[[i]]) <- SNP_distances_M2[[i]][,1]
  SNP_distances_M2[[i]][,1] <- NULL
}

# Get half of the distance matrices (avoid duplicated values)
for (i in 1:length(SNP_distances_M2)) {
  SNP_distances_M2[[i]][upper.tri(SNP_distances_M2[[i]])] <- NA
}


### ===== 2.2 ESTIMATE POTENTIAL SHARING EVENTS (PSEs) IN RELATED VS UNRELATED MOTHER-INFANT PAIRS ===== ###

# Generate all combinations of pair IDs 
mothers_M2 <- colnames(genera_counts_M2)[grep("mother_human_milk", colnames(genera_counts_M2))]
babies_M2 <- colnames(genera_counts_M2)[grep("infant_faeces", colnames(genera_counts_M2))]

all_pairs_M2 <- expand.grid(mothers_M2, babies_M2)
all_pairs_combined_M2 <- paste(all_pairs_M2$Var1, all_pairs_M2$Var2, sep="_")

all_pairs_tmp_M2 <- all_pairs_M2
all_pairs_tmp_M2$Var1 <- gsub("mother_human_milk_", "", all_pairs_tmp_M2$Var1)
all_pairs_tmp_M2$Var2 <- gsub("infant_faeces_", "", all_pairs_tmp_M2$Var2)

related_famIDs_M2 <- c()
unrelated_famIDs_M2 <- c()

for (i in 1:nrow(all_pairs_tmp_M2)){
  
  if (all_pairs_tmp_M2$Var1[i] == all_pairs_tmp_M2$Var2[i]){
    related_famIDs_M2 <- c(related_famIDs_M2, paste0("mother_human_milk_", all_pairs_tmp_M2$Var1[i], "_infant_faeces_", all_pairs_tmp_M2$Var2[i]))
    
  }else{
    unrelated_famIDs_M2 <- c(unrelated_famIDs_M2, paste0("mother_human_milk_", all_pairs_tmp_M2$Var1[i], "_infant_faeces_", all_pairs_tmp_M2$Var2[i]))
    
  }
}

related_pairs_M2 <- unique(related_famIDs_M2)
unrelated_pairs_M2 <- unique(unrelated_famIDs_M2)
related_pairs_M2 <- gsub("_infant_faeces.*","", gsub("mother_human_milk_", "", related_pairs_M2))


# Generate the empty final data table
final_PSE_data_M2 <- data.frame(matrix(nrow=length(all_pairs_combined_M2), ncol=length(names_M2)))
colnames(final_PSE_data_M2) <- names_M2
rownames(final_PSE_data_M2) <- c(related_pairs_M2, unrelated_pairs_M2)
number_related_pairs_M2 <- length(unique(all_pairs_tmp_M2$Var1))
final_PSE_data_M2$Relatedness <- c(rep("Related", number_related_pairs_M2),
                                   rep("Unrelated", length(all_pairs_combined_M2) - number_related_pairs_M2))

# Number of pairs in which each genus is detected
all_pairs_M2 <- table(unlist(genera_in_pairs_M2))
total_pairs_M2 <- all_pairs_M2[names(all_pairs_M2) %in% names_shared_genera_M2]

# Pairs in which each genus is detected
Acinetobacter_in_pairs_M2 <- names(grep("Acinetobacter", genera_in_pairs_M2, value = T))
Actinomyces_in_pairs_M2 <- names(grep("Actinomyces", genera_in_pairs_M2, value = T))
Bacteroides_in_pairs_M2 <- names(grep("Bacteroides", genera_in_pairs_M2, value = T))
Bifidobacterium_in_pairs_M2 <- names(grep("Bifidobacterium", genera_in_pairs_M2, value = T))
Cutibacterium_in_pairs_M2 <- names(grep("Cutibacterium", genera_in_pairs_M2, value = T))
Enterobacter_in_pairs_M2 <- names(grep("Enterobacter", genera_in_pairs_M2, value = T))
Enterococcus_in_pairs_M2 <- names(grep("Enterococcus", genera_in_pairs_M2, value = T))
Escherichia_Shigella_in_pairs_M2 <- names(grep("Escherichia_Shigella", genera_in_pairs_M2, value = T))
Gemella_in_pairs_M2 <- names(grep("Gemella", genera_in_pairs_M2, value = T))
Haemophilus_in_pairs_M2 <- names(grep("Haemophilus", genera_in_pairs_M2, value = T))
Lactobacillus_in_pairs_M2 <- names(grep("Lactobacillus", genera_in_pairs_M2, value = T))
Prevotella_in_pairs_M2 <- names(grep("Prevotella", genera_in_pairs_M2, value = T))
Rothia_in_pairs_M2 <- names(grep("Rothia", genera_in_pairs_M2, value = T))
Staphylococcus_in_pairs_M2 <- names(grep("Staphylococcus", genera_in_pairs_M2, value = T))
Streptococcus_in_pairs_M2 <- names(grep("Streptococcus", genera_in_pairs_M2, value = T))
Veillonella_in_pairs_M2 <- names(grep("Veillonella", genera_in_pairs_M2, value = T))


## Get the contingency table and the pairs in which we observe PSEs. Add results to final_PSE_data_M2
PSE_Acinetobacter_M2 <- exact_PSE_analysis(SNP_distances_M2$Acinetobacter_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Acinetobacter_in_pairs_M2,"Acinetobacter_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Acinetobacter_M2$Pairs_with_PSEs),"Acinetobacter_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Acinetobacter_M2$Unrelated_pairs_without_PSE),"Acinetobacter_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Acinetobacter_M2$Unrelated_pairs_with_PSE),"Acinetobacter_M2" ] <- "PSE"

PSE_Actinomyces_M2 <- exact_PSE_analysis(SNP_distances_M2$Actinomyces_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Actinomyces_in_pairs_M2,"Actinomyces_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Actinomyces_M2$Pairs_with_PSEs),"Actinomyces_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Actinomyces_M2$Unrelated_pairs_without_PSE),"Actinomyces_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Actinomyces_M2$Unrelated_pairs_with_PSE),"Actinomyces_M2" ] <- "PSE"

PSE_Bacteroides_M2 <- exact_PSE_analysis(SNP_distances_M2$Bacteroides_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Bacteroides_in_pairs_M2,"Bacteroides_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Bacteroides_M2$Pairs_with_PSEs),"Bacteroides_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Bacteroides_M2$Unrelated_pairs_without_PSE),"Bacteroides_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Bacteroides_M2$Unrelated_pairs_with_PSE),"Bacteroides_M2" ] <- "PSE"

PSE_Bifidobacterium_M2 <- exact_PSE_analysis(SNP_distances_M2$Bifidobacterium_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Bifidobacterium_in_pairs_M2,"Bifidobacterium_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Bifidobacterium_M2$Pairs_with_PSEs),"Bifidobacterium_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Bifidobacterium_M2$Unrelated_pairs_without_PSE),"Bifidobacterium_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Bifidobacterium_M2$Unrelated_pairs_with_PSE),"Bifidobacterium_M2" ] <- "PSE"

PSE_Cutibacterium_M2 <- exact_PSE_analysis(SNP_distances_M2$Cutibacterium_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Cutibacterium_in_pairs_M2,"Cutibacterium_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Cutibacterium_M2$Pairs_with_PSEs),"Cutibacterium_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Cutibacterium_M2$Unrelated_pairs_without_PSE),"Cutibacterium_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Cutibacterium_M2$Unrelated_pairs_with_PSE),"Cutibacterium_M2" ] <- "PSE"

PSE_Enterobacter_M2 <- exact_PSE_analysis(SNP_distances_M2$Enterobacter_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Enterobacter_in_pairs_M2,"Enterobacter_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Enterobacter_M2$Pairs_with_PSEs),"Enterobacter_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Enterobacter_M2$Unrelated_pairs_without_PSE),"Enterobacter_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Enterobacter_M2$Unrelated_pairs_with_PSE),"Enterobacter_M2" ] <- "PSE"

PSE_Enterococcus_M2 <- exact_PSE_analysis(SNP_distances_M2$Enterococcus_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Enterococcus_in_pairs_M2,"Enterococcus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Enterococcus_M2$Pairs_with_PSEs),"Enterococcus_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Enterococcus_M2$Unrelated_pairs_without_PSE),"Enterococcus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Enterococcus_M2$Unrelated_pairs_with_PSE),"Enterococcus_M2" ] <- "PSE"

PSE_Escherichia_Shigella_M2 <- exact_PSE_analysis(SNP_distances_M2$Escherichia_Shigella_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Escherichia_Shigella_in_pairs_M2,"Escherichia_Shigella_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Escherichia_Shigella_M2$Pairs_with_PSEs),"Escherichia_Shigella_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Escherichia_Shigella_M2$Unrelated_pairs_without_PSE),"Escherichia_Shigella_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Escherichia_Shigella_M2$Unrelated_pairs_with_PSE),"Escherichia_Shigella_M2" ] <- "PSE"

PSE_Gemella_M2 <- exact_PSE_analysis(SNP_distances_M2$Gemella_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Gemella_in_pairs_M2,"Gemella_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Gemella_M2$Pairs_with_PSEs),"Gemella_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Gemella_M2$Unrelated_pairs_without_PSE),"Gemella_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Gemella_M2$Unrelated_pairs_with_PSE),"Gemella_M2" ] <- "PSE"

PSE_Haemophilus_M2 <- exact_PSE_analysis(SNP_distances_M2$Haemophilus_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Haemophilus_in_pairs_M2,"Haemophilus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Haemophilus_M2$Pairs_with_PSEs),"Haemophilus_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Haemophilus_M2$Unrelated_pairs_without_PSE),"Haemophilus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Haemophilus_M2$Unrelated_pairs_with_PSE),"Haemophilus_M2" ] <- "PSE"

PSE_Lactobacillus_M2 <- exact_PSE_analysis(SNP_distances_M2$Lactobacillus_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Lactobacillus_in_pairs_M2,"Lactobacillus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Lactobacillus_M2$Pairs_with_PSEs),"Lactobacillus_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Lactobacillus_M2$Unrelated_pairs_without_PSE),"Lactobacillus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Lactobacillus_M2$Unrelated_pairs_with_PSE),"Lactobacillus_M2" ] <- "PSE"

PSE_Prevotella_M2 <- exact_PSE_analysis(SNP_distances_M2$Prevotella_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Prevotella_in_pairs_M2,"Prevotella_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Prevotella_M2$Pairs_with_PSEs),"Prevotella_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Prevotella_M2$Unrelated_pairs_without_PSE),"Prevotella_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Prevotella_M2$Unrelated_pairs_with_PSE),"Prevotella_M2" ] <- "PSE"

PSE_Rothia_M2 <- exact_PSE_analysis(SNP_distances_M2$Rothia_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Rothia_in_pairs_M2,"Rothia_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Rothia_M2$Pairs_with_PSEs),"Rothia_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Rothia_M2$Unrelated_pairs_without_PSE),"Rothia_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Rothia_M2$Unrelated_pairs_with_PSE),"Rothia_M2" ] <- "PSE"

PSE_Staphylococcus_M2 <- exact_PSE_analysis(SNP_distances_M2$Staphylococcus_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Staphylococcus_in_pairs_M2,"Staphylococcus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Staphylococcus_M2$Pairs_with_PSEs),"Staphylococcus_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Staphylococcus_M2$Unrelated_pairs_without_PSE),"Staphylococcus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Staphylococcus_M2$Unrelated_pairs_with_PSE),"Staphylococcus_M2" ] <- "PSE"

PSE_Streptococcus_M2 <- exact_PSE_analysis(SNP_distances_M2$Streptococcus_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Streptococcus_in_pairs_M2,"Streptococcus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Streptococcus_M2$Pairs_with_PSEs),"Streptococcus_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Streptococcus_M2$Unrelated_pairs_without_PSE),"Streptococcus_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Streptococcus_M2$Unrelated_pairs_with_PSE),"Streptococcus_M2" ] <- "PSE"

PSE_Veillonella_M2 <- exact_PSE_analysis(SNP_distances_M2$Veillonella_M2, matching_ASVs_M2)
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% Veillonella_in_pairs_M2,"Veillonella_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Veillonella_M2$Pairs_with_PSEs),"Veillonella_M2" ] <- "PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Veillonella_M2$Unrelated_pairs_without_PSE),"Veillonella_M2" ] <- "no_PSE"
final_PSE_data_M2[rownames(final_PSE_data_M2) %in% unique(PSE_Veillonella_M2$Unrelated_pairs_with_PSE),"Veillonella_M2" ] <- "PSE"


# Save contingency tables in list
cont_tables_M2 <- list(PSE_Acinetobacter_M2$Contingency_table,
                       PSE_Actinomyces_M2$Contingency_table,
                       PSE_Bacteroides_M2$Contingency_table,
                       PSE_Bifidobacterium_M2$Contingency_table,
                       PSE_Cutibacterium_M2$Contingency_table,
                       PSE_Enterobacter_M2$Contingency_table,
                       PSE_Enterococcus_M2$Contingency_table,
                       PSE_Escherichia_Shigella_M2$Contingency_table,
                       PSE_Gemella_M2$Contingency_table,
                       PSE_Haemophilus_M2$Contingency_table,
                       PSE_Lactobacillus_M2$Contingency_table,
                       PSE_Prevotella_M2$Contingency_table,
                       PSE_Rothia_M2$Contingency_table,
                       PSE_Staphylococcus_M2$Contingency_table,
                       PSE_Streptococcus_M2$Contingency_table,
                       PSE_Veillonella_M2$Contingency_table)

names(cont_tables_M2) <- names_shared_genera_M2


# Move relatedness column to the first position
relatedness_column <- final_PSE_data_M2$Relatedness
final_PSE_data_M2 <- final_PSE_data_M2[, !(names(final_PSE_data_M2) == "Relatedness")]
final_PSE_data_M2 <- data.frame(Relatedness = relatedness_column, final_PSE_data_M2)


### ===== 2.3 SAVE RESULTS ===== ###

setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/")

write.table(final_PSE_data_M2,"M2/240219_PSE_table_maternal_milk_infant_faeces_M2.txt", sep = "\t", row.names = T, quote = FALSE) 

# save contingecy table
saveRDS(cont_tables_M2, "M2/240219_maternal_milk_infant_faeces_contingency_tables_M2.rds")



