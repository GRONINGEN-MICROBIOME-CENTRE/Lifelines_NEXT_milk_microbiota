################################################################################################################
### INVESTIGATE SHARING OF BACTERIA BETWEEN MATERNAL FAECES-MILK AND BETWEEN MILK-INFANT FAECES - SCRIPT 1 #####
################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 1. DATA IMPORT & LOADING OF LIBRARIES
# 1.1 DATA IMPORT
# 1.2 LOADING OF LIBRARIES
# 
# 2. SELECT SAMPLES AND COLUMNS OF INTEREST
# 2.1 MILK - INFANT FAECES: SELECT SAMPLES AND COLUMNS OF INTEREST
# 2.2 MATERNAL FAECES - MILK: SELECT SAMPLES AND COLUMNS OF INTEREST
# 
# 3. GET THE GENERA (AND THEIR ASVs) PRESENT IN PAIRED SAMPLES
# 3.1 MILK - INFANT FAECES AT M1
# 3.2 MILK - INFANT FAECES AT M2
# 3.3 MILK - INFANT FAECES AT M3
# 3.4 MILK - INFANT FAECES AT M6
# 3.5 MATERNAL FAECES - MILK AT M3
# 3.6 OVERVIEW OF POTENTIALLY SHARED GENERA
# 
# 4. GET THE MULTIPLE SEQUENCE ALIGNMENT (MSA) FOR EACH GENUS DETECTED AMONG PAIRED INDIVIDUALS
# 4.1 FUNCTION FOR MULTIPLE SEQUENCE ALIGNMENT (MSA) WITH ClustalW
# 4.2 PREPARE FILES FOR MSA

## Notes:
## I manually selected the bacterial genera present in ≥5 pairs at ≥1 time point
## Ranko ran ClustalO instead of ClustalW
## Asier used snp-dists to get SNP distances


##### =========================== 1. DATA IMPORT & LOADING OF LIBRARIES =========================== #####

### ===== 1.1 DATA IMPORT ===== ###

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk")

## import files, phenotypes+microbiota data on genus level
data <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/231121_milk_composition_paper_sel_phenotypes_sel_microbiota_data_ASVs_n1515.rds")
#1515x7282

## exclude ASVs from unclassified genera
ASVcf <- data[,-grep("g__NA", colnames(data))]
#1515 x 5917

## Columns   Data
## 1:526     Metadata
## 527:5917  ASVs


### ===== 1.2 LOADING OF LIBRARIES ===== ###

library(stringr)
library(seqinr)
library(msa)


##### =========================== 2. SELECT SAMPLES AND COLUMNS OF INTEREST =========================== #####

### ===== 2.1 MILK - INFANT FAECES: SELECT SAMPLES AND COLUMNS OF INTEREST ===== ###

## select only milk and infant faecal samples
milkinff <- ASVcf[ASVcf$sample_origin_type=="mother_human_milk" | ASVcf$sample_origin_type=="infant_faeces", c(1:14,503,506:508,527:ncol(ASVcf))]
## 1289 x 5409
#737 milk samples
#552 infant faecal samples

## add a column showing IDs of mother-infant pairs
milkinff$pair_ID <- c(paste0(milkinff$family_ID, "_", milkinff$NEXT_participation_number, "_", milkinff$infant_number))

## resort columns
milkinff <- milkinff[,c(5410,1:5409)]

## split the file per time point (because samples should be matched per time point)
milkinffM1 <- milkinff[milkinff$time_point=="1_month",]  #532 x 5410; 303 milk and 229 infant faeces
milkinffM2 <- milkinff[milkinff$time_point=="2_months",] # 90 x 5410;  52 milk and  38 infant faeces
milkinffM3 <- milkinff[milkinff$time_point=="3_months",] #466 x 5410; 271 milk and 195 infant faeces
milkinffM6 <- milkinff[milkinff$time_point=="6_months",] #201 x 5410; 111 milk and  90 infant faeces

## per time point, select only individuals with matched samples
## M1
milkM1 <- milkinffM1[milkinffM1$sample_origin_type=="mother_human_milk",]
inffM1 <- milkinffM1[milkinffM1$sample_origin_type=="infant_faeces",]
milkinffM1IDs <- intersect(milkM1$pair_ID, inffM1$pair_ID) #210 matched samples
selmilkinffM1 <- milkinffM1[milkinffM1$pair_ID %in% milkinffM1IDs,] #420 x 5410

## M2
milkM2 <- milkinffM2[milkinffM2$sample_origin_type=="mother_human_milk",]
inffM2 <- milkinffM2[milkinffM2$sample_origin_type=="infant_faeces",]
milkinffM2IDs <- intersect(milkM2$pair_ID, inffM2$pair_ID) #34 matched samples
selmilkinffM2 <- milkinffM2[milkinffM2$pair_ID %in% milkinffM2IDs,] #68 x 5410

## M3
milkM3 <- milkinffM3[milkinffM3$sample_origin_type=="mother_human_milk",]
inffM3 <- milkinffM3[milkinffM3$sample_origin_type=="infant_faeces",]
milkinffM3IDs <- intersect(milkM3$pair_ID, inffM3$pair_ID) #190 matched samples
selmilkinffM3 <- milkinffM3[milkinffM3$pair_ID %in% milkinffM3IDs,] #380 x 5410

## M6
milkM6 <- milkinffM6[milkinffM6$sample_origin_type=="mother_human_milk",]
inffM6 <- milkinffM6[milkinffM6$sample_origin_type=="infant_faeces",]
milkinffM6IDs <- intersect(milkM6$pair_ID, inffM6$pair_ID) #75 matched samples
selmilkinffM6 <- milkinffM6[milkinffM6$pair_ID %in% milkinffM6IDs,] #150 x 5410


### ===== 2.2 MATERNAL FAECES - MILK: SELECT SAMPLES AND COLUMNS OF INTEREST ===== ###

## select only milk and maternal faecal samples, all from M3
matfmilkM3 <- ASVcf[(ASVcf$sample_origin_type=="mother_human_milk" | ASVcf$sample_origin_type=="mother_faeces") & ASVcf$time_point=="3_months", c(1:14,503,506:508,527:ncol(ASVcf))]
## 443 x 5409
#172 maternal faecal samples
#271 milk samples

## exclude duplicated maternal data
matfmilkM3 <- matfmilkM3[grep(";1", matfmilkM3$mother_ID),]
## 438 x 5409
#170 maternal faecal samples
#268 milk samples

## add a column showing IDs of mother-infant pairs
matfmilkM3$pair_ID <- c(paste0(matfmilkM3$family_ID, "_", matfmilkM3$NEXT_participation_number, "_", matfmilkM3$infant_number))

## resort columns
matfmilkM3 <- matfmilkM3[,c(5410,1:5409)]

## select only individuals with matched samples
milkM3 <- matfmilkM3[matfmilkM3$sample_origin_type=="mother_human_milk",]
matfM3 <- matfmilkM3[matfmilkM3$sample_origin_type=="mother_faeces",]
matfmilkM3IDs <- intersect(milkM3$pair_ID, matfM3$pair_ID) #164 matched samples
selmatfmilkM3 <- matfmilkM3[matfmilkM3$pair_ID %in% matfmilkM3IDs,] #328 x 5410


##### =========================== 3. GET THE GENERA (AND THEIR ASVs) PRESENT IN PAIRED SAMPLES =========================== #####

#************************************************************************
# A2.Get the genera (and their ASVs) present in at least x paired samples
#************************************************************************
# For each sample select those ASVs that are present
# Use the sequencing technology +  sample_ID + lowest-level taxonomic assignment as ASV name
# Concatenate all sequences in a list


### ===== 3.1 MILK - INFANT FAECES AT M1 ===== ###

rownames(selmilkinffM1) <- c(paste0(selmilkinffM1$sample_origin_type, "_", selmilkinffM1$pair_ID))
sample_ASVs_M1 <- selmilkinffM1[,20:ncol(milkinffM1)] #420 x 5391

ASV_final_list_M1 <- list()

for (i in 1:nrow(sample_ASVs_M1)) {
  print(paste0("Starting analysis for sample ", i))
  
  print("Getting sample names")
  sample_name <- rownames(sample_ASVs_M1)[i]
  
  print("Getting sequences")
  ASV <- gsub(".*_", "",  colnames(sample_ASVs_M1)[which(sample_ASVs_M1[i,] != 0)]) #sequence
  
  print("Getting bacterial names")
  name_full <- gsub("_ASV.*", "\\1", colnames(sample_ASVs_M1)[which(sample_ASVs_M1[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everything after genus level
  name <- gsub(".*g__", "g__", name) #remove everything before genus level (get genus name)
  
  print("Getting bacterial names if they were unclassified on genus level (if required)")
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  
  print("Creating final ASV list")
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  ASV_final_list_M1 <- append(ASV_final_list_M1, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list_M1)
# all represent classified genera (expected as I excluded unclassified ones before)

# Create a dataframe with the counts of each genus for each participant 
names_genera_M1 <- gsub(".*g__", "", grep("g__", names(ASV_final_list_M1), value = T))
names_genera_M1 <- unique(names_genera_M1) #326 unique genera
genera_counts_M1 <- data.frame(matrix(ncol = length(rownames(sample_ASVs_M1)) + 1, nrow = length(names_genera_M1)))
colnames(genera_counts_M1) <- c("genera", rownames(sample_ASVs_M1))
genera_counts_M1$genera <- names_genera_M1

for (i in 2:ncol(genera_counts_M1)) {
  microbes <- names(ASV_final_list_M1)[grep(paste0(colnames(genera_counts_M1)[i],"\\."), names(ASV_final_list_M1))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  genera_counts_M1[i] <- counts$Freq[match(genera_counts_M1$genera,counts$genera)]
}
genera_counts_M1[is.na(genera_counts_M1)] <- 0 #set NAs to 0

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
genera_in_pairs_M1 <- list()
IDs <- milkinffM1IDs # IDs <- unique(metadata_M1$pair_ID)

for (i in IDs) {
  print("Getting index")
  index <- which(IDs == i)
  print("Getting pattern")
  pattern <- paste0("_",i, "$")
  print("Getting pair counts")
  pair_counts <- genera_counts_M1[, grep(pattern, colnames(genera_counts_M1))]
  print("Selecting genera deteced in at least 1 pair")
  genera_in_pairs_M1 [[index]] <- genera_counts_M1$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair 
}

genera_in_pairs_M1 <- genera_in_pairs_M1[!sapply(genera_in_pairs_M1,is.null)] # removes null elements
names(genera_in_pairs_M1) <- IDs

# Get names of the shared genera present in 1 or more pairs
shared_genera_M1 <- names(which(table(unlist(genera_in_pairs_M1)) > 0)) #change this to a diff number if I want to be stricter on which genera to include
# shared_genera_5_pairs_M1 <- names(which(table(unlist(genera_in_pairs_M1)) > 4)) 

# Get all ASVs of selected genera
names(ASV_final_list_M1) <- make.unique(names(ASV_final_list_M1), sep = '_') # distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matching_ASVs_M1 <- ASV_final_list_M1[unique (grep(paste(shared_genera_M1,collapse="|"), # get all ASVs belonging to shared genera
                                                   names(ASV_final_list_M1)))]

## save results 
write.table(genera_counts_M1,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M1/M1_milk_infant_faeces_genera_counts.txt", sep = "\t", row.names = T, quote = FALSE)
saveRDS(genera_in_pairs_M1, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M1/M1_milk_infant_faeces_genera_counts_pairs.rds")
saveRDS(matching_ASVs_M1, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M1/M1_milk_infant_faeces_ASVs_shared_genera.rds")
# saveRDS(shared_genera_M1, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M1/M1_milk_infant_faeces_ASVs_shared_genera_at_least_1_pair.rds")
# saveRDS(shared_genera_5_pairs_M1, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M1/M1_milk_infant_faeces_ASVs_shared_genera_at_least_5_pairs.rds")
write.table(shared_genera_M1,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M1/M1_milk_infant_faeces_ASVs_shared_genera_at_least_1_pair.txt", sep = "\t", row.names = F, quote = FALSE)
write.table(shared_genera_5_pairs_M1,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M1/M1_milk_infant_faeces_ASVs_shared_genera_at_least_5_pairs.txt", sep = "\t", row.names = F, quote = FALSE)


### ===== 3.2 MILK - INFANT FAECES AT M2 ===== ###

rownames(selmilkinffM2) <- c(paste0(selmilkinffM2$sample_origin_type, "_", selmilkinffM2$pair_ID))
sample_ASVs_M2 <- selmilkinffM2[,20:ncol(milkinffM2)] #68 x 5391

ASV_final_list_M2 <- list()

for (i in 1:nrow(sample_ASVs_M2)) {
  print(paste0("Starting analysis for sample ", i))
  
  print("Getting sample names")
  sample_name <- rownames(sample_ASVs_M2)[i]
  
  print("Getting sequences")
  ASV <- gsub(".*_", "",  colnames(sample_ASVs_M2)[which(sample_ASVs_M2[i,] != 0)]) #sequence
  
  print("Getting bacterial names")
  name_full <- gsub("_ASV.*", "\\1", colnames(sample_ASVs_M2)[which(sample_ASVs_M2[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everything after genus level
  name <- gsub(".*g__", "g__", name) #remove everything before genus level (get genus name)
  
  print("Getting bacterial names if they were unclassified on genus level (if required)")
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  
  print("Creating final ASV list")
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  ASV_final_list_M2 <- append(ASV_final_list_M2, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list_M2)
# all represent classified genera (expected as I excluded unclassified ones before)

# Create a dataframe with the counts of each genus for each participant 
names_genera_M2 <- gsub(".*g__", "", grep("g__", names(ASV_final_list_M2), value = T))
names_genera_M2 <- unique(names_genera_M2) #326 unique genera
genera_counts_M2 <- data.frame(matrix(ncol = length(rownames(sample_ASVs_M2)) + 1, nrow = length(names_genera_M2)))
colnames(genera_counts_M2) <- c("genera", rownames(sample_ASVs_M2))
genera_counts_M2$genera <- names_genera_M2

for (i in 2:ncol(genera_counts_M2)) {
  microbes <- names(ASV_final_list_M2)[grep(paste0(colnames(genera_counts_M2)[i],"\\."), names(ASV_final_list_M2))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  genera_counts_M2[i] <- counts$Freq[match(genera_counts_M2$genera,counts$genera)]
}
genera_counts_M2[is.na(genera_counts_M2)] <- 0 #set NAs to 0

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
genera_in_pairs_M2 <- list()
IDs <- milkinffM2IDs # IDs <- unique(metadata_M2$pair_ID)

for (i in IDs) {
  print("Getting index")
  index <- which(IDs == i)
  print("Getting pattern")
  pattern <- paste0("_",i, "$")
  print("Getting pair counts")
  pair_counts <- genera_counts_M2[, grep(pattern, colnames(genera_counts_M2))]
  print("Selecting genera deteced in at least 1 pair")
  genera_in_pairs_M2 [[index]] <- genera_counts_M2$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair 
}

genera_in_pairs_M2 <- genera_in_pairs_M2[!sapply(genera_in_pairs_M2,is.null)] # removes null elements
names(genera_in_pairs_M2) <- IDs

# Get names of the shared genera present in 1 or more pairs
shared_genera_M2 <- names(which(table(unlist(genera_in_pairs_M2)) > 0)) #change this to a diff number if I want to be stricter on which genera to include
# shared_genera_5_pairs_M2 <- names(which(table(unlist(genera_in_pairs_M2)) > 4))

# Get all ASVs of selected genera
names(ASV_final_list_M2) <- make.unique(names(ASV_final_list_M2), sep = '_') # distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matching_ASVs_M2 <- ASV_final_list_M2[unique (grep(paste(shared_genera_M2,collapse="|"), # get all ASVs belonging to shared genera
                                                   names(ASV_final_list_M2)))]

## save results 
write.table(genera_counts_M2,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_genera_counts.txt", sep = "\t", row.names = T, quote = FALSE)
saveRDS(genera_in_pairs_M2, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_genera_counts_pairs.rds")
saveRDS(matching_ASVs_M2, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_ASVs_shared_genera.rds")
# saveRDS(shared_genera_M2, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_ASVs_shared_genera_at_least_1_pair.rds")
# saveRDS(shared_genera_5_pairs_M2, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_ASVs_shared_genera_at_least_5_pairs.rds")
write.table(shared_genera_M2,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_ASVs_shared_genera_at_least_1_pair.txt", sep = "\t", row.names = F, quote = FALSE)
write.table(shared_genera_5_pairs_M2,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/M2_milk_infant_faeces_ASVs_shared_genera_at_least_5_pairs.txt", sep = "\t", row.names = F, quote = FALSE)


### ===== 3.3 MILK - INFANT FAECES AT M3 ===== ###

rownames(selmilkinffM3) <- c(paste0(selmilkinffM3$sample_origin_type, "_", selmilkinffM3$pair_ID))
sample_ASVs_M3 <- selmilkinffM3[,20:ncol(milkinffM3)] #380 x 5391

ASV_final_list_M3 <- list()

for (i in 1:nrow(sample_ASVs_M3)) {
  print(paste0("Starting analysis for sample ", i))
  
  print("Getting sample names")
  sample_name <- rownames(sample_ASVs_M3)[i]
  
  print("Getting sequences")
  ASV <- gsub(".*_", "",  colnames(sample_ASVs_M3)[which(sample_ASVs_M3[i,] != 0)]) #sequence
  
  print("Getting bacterial names")
  name_full <- gsub("_ASV.*", "\\1", colnames(sample_ASVs_M3)[which(sample_ASVs_M3[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everything after genus level
  name <- gsub(".*g__", "g__", name) #remove everything before genus level (get genus name)
  
  print("Getting bacterial names if they were unclassified on genus level (if required)")
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  
  print("Creating final ASV list")
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  ASV_final_list_M3 <- append(ASV_final_list_M3, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list_M3)
# all represent classified genera (expected as I excluded unclassified ones before)

# Create a dataframe with the counts of each genus for each participant 
names_genera_M3 <- gsub(".*g__", "", grep("g__", names(ASV_final_list_M3), value = T))
names_genera_M3 <- unique(names_genera_M3) #326 unique genera
genera_counts_M3 <- data.frame(matrix(ncol = length(rownames(sample_ASVs_M3)) + 1, nrow = length(names_genera_M3)))
colnames(genera_counts_M3) <- c("genera", rownames(sample_ASVs_M3))
genera_counts_M3$genera <- names_genera_M3

for (i in 2:ncol(genera_counts_M3)) {
  microbes <- names(ASV_final_list_M3)[grep(paste0(colnames(genera_counts_M3)[i],"\\."), names(ASV_final_list_M3))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  genera_counts_M3[i] <- counts$Freq[match(genera_counts_M3$genera,counts$genera)]
}
genera_counts_M3[is.na(genera_counts_M3)] <- 0 #set NAs to 0

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
genera_in_pairs_M3 <- list()
IDs <- milkinffM3IDs # IDs <- unique(metadata_M3$pair_ID)

for (i in IDs) {
  print("Getting index")
  index <- which(IDs == i)
  print("Getting pattern")
  pattern <- paste0("_",i, "$")
  print("Getting pair counts")
  pair_counts <- genera_counts_M3[, grep(pattern, colnames(genera_counts_M3))]
  print("Selecting genera deteced in at least 1 pair")
  genera_in_pairs_M3 [[index]] <- genera_counts_M3$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair 
}

genera_in_pairs_M3 <- genera_in_pairs_M3[!sapply(genera_in_pairs_M3,is.null)] # removes null elements
names(genera_in_pairs_M3) <- IDs

# Get names of the shared genera present in 1 or more pairs
shared_genera_M3 <- names(which(table(unlist(genera_in_pairs_M3)) > 0)) #change this to a diff number if I want to be stricter on which genera to include
# shared_genera_5_pairs_M3 <- names(which(table(unlist(genera_in_pairs_M3)) > 4))

# Get all ASVs of selected genera
names(ASV_final_list_M3) <- make.unique(names(ASV_final_list_M3), sep = '_') # distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matching_ASVs_M3 <- ASV_final_list_M3[unique (grep(paste(shared_genera_M3,collapse="|"), # get all ASVs belonging to shared genera
                                                   names(ASV_final_list_M3)))]

## save results 
write.table(genera_counts_M3,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M3/M3_milk_infant_faeces_genera_counts.txt", sep = "\t", row.names = T, quote = FALSE)
saveRDS(genera_in_pairs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M3/M3_milk_infant_faeces_genera_counts_pairs.rds")
saveRDS(matching_ASVs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M3/M3_milk_infant_faeces_ASVs_shared_genera.rds")
# saveRDS(shared_genera_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M3/M3_milk_infant_faeces_ASVs_shared_genera_at_least_1_pair.rds")
# saveRDS(shared_genera_5_pairs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M3/M3_milk_infant_faeces_ASVs_shared_genera_at_least_5_pairs.rds")
write.table(shared_genera_M3,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M3/M3_milk_infant_faeces_ASVs_shared_genera_at_least_1_pair.txt", sep = "\t", row.names = F, quote = FALSE)
write.table(shared_genera_5_pairs_M3,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M3/M3_milk_infant_faeces_ASVs_shared_genera_at_least_5_pairs.txt", sep = "\t", row.names = F, quote = FALSE)


### ===== 3.4 MILK - INFANT FAECES AT M6 ===== ###

rownames(selmilkinffM6) <- c(paste0(selmilkinffM6$sample_origin_type, "_", selmilkinffM6$pair_ID))
sample_ASVs_M6 <- selmilkinffM6[,20:ncol(milkinffM6)] #150 x 5391

ASV_final_list_M6 <- list()

for (i in 1:nrow(sample_ASVs_M6)) {
  print(paste0("Starting analysis for sample ", i))
  
  print("Getting sample names")
  sample_name <- rownames(sample_ASVs_M6)[i]
  
  print("Getting sequences")
  ASV <- gsub(".*_", "",  colnames(sample_ASVs_M6)[which(sample_ASVs_M6[i,] != 0)]) #sequence
  
  print("Getting bacterial names")
  name_full <- gsub("_ASV.*", "\\1", colnames(sample_ASVs_M6)[which(sample_ASVs_M6[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everything after genus level
  name <- gsub(".*g__", "g__", name) #remove everything before genus level (get genus name)
  
  print("Getting bacterial names if they were unclassified on genus level (if required)")
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  
  print("Creating final ASV list")
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  ASV_final_list_M6 <- append(ASV_final_list_M6, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list_M6)
# all represent classified genera (expected as I excluded unclassified ones before)

# Create a dataframe with the counts of each genus for each participant 
names_genera_M6 <- gsub(".*g__", "", grep("g__", names(ASV_final_list_M6), value = T))
names_genera_M6 <- unique(names_genera_M6) #326 unique genera
genera_counts_M6 <- data.frame(matrix(ncol = length(rownames(sample_ASVs_M6)) + 1, nrow = length(names_genera_M6)))
colnames(genera_counts_M6) <- c("genera", rownames(sample_ASVs_M6))
genera_counts_M6$genera <- names_genera_M6

for (i in 2:ncol(genera_counts_M6)) {
  microbes <- names(ASV_final_list_M6)[grep(paste0(colnames(genera_counts_M6)[i],"\\."), names(ASV_final_list_M6))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  genera_counts_M6[i] <- counts$Freq[match(genera_counts_M6$genera,counts$genera)]
}
genera_counts_M6[is.na(genera_counts_M6)] <- 0 #set NAs to 0

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
genera_in_pairs_M6 <- list()
IDs <- milkinffM6IDs # IDs <- unique(metadata_M6$pair_ID)

for (i in IDs) {
  print("Getting index")
  index <- which(IDs == i)
  print("Getting pattern")
  pattern <- paste0("_",i, "$")
  print("Getting pair counts")
  pair_counts <- genera_counts_M6[, grep(pattern, colnames(genera_counts_M6))]
  print("Selecting genera deteced in at least 1 pair")
  genera_in_pairs_M6 [[index]] <- genera_counts_M6$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair 
}

genera_in_pairs_M6 <- genera_in_pairs_M6[!sapply(genera_in_pairs_M6,is.null)] # removes null elements
names(genera_in_pairs_M6) <- IDs

# Get names of the shared genera present in 1 or more pairs
shared_genera_M6 <- names(which(table(unlist(genera_in_pairs_M6)) > 0)) #change this to a diff number if I want to be stricter on which genera to include
# shared_genera_5_pairs_M6 <- names(which(table(unlist(genera_in_pairs_M6)) > 4))

# Get all ASVs of selected genera
names(ASV_final_list_M6) <- make.unique(names(ASV_final_list_M6), sep = '_') # distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matching_ASVs_M6 <- ASV_final_list_M6[unique (grep(paste(shared_genera_M6,collapse="|"), # get all ASVs belonging to shared genera
                                                   names(ASV_final_list_M6)))]

## save results 
write.table(genera_counts_M6,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M6/M6_milk_infant_faeces_genera_counts.txt", sep = "\t", row.names = T, quote = FALSE)
saveRDS(genera_in_pairs_M6, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M6/M6_milk_infant_faeces_genera_counts_pairs.rds")
saveRDS(matching_ASVs_M6, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M6/M6_milk_infant_faeces_ASVs_shared_genera.rds")
# saveRDS(shared_genera_M6, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M6/M6_milk_infant_faeces_ASVs_shared_genera_at_least_1_pair.rds")
# saveRDS(shared_genera_5_pairs_M6, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M6/M6_milk_infant_faeces_ASVs_shared_genera_at_least_5_pairs.rds")
write.table(shared_genera_M6,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M6/M6_milk_infant_faeces_ASVs_shared_genera_at_least_1_pair.txt", sep = "\t", row.names = F, quote = FALSE)
write.table(shared_genera_5_pairs_M6,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M6/M6_milk_infant_faeces_ASVs_shared_genera_at_least_5_pairs.txt", sep = "\t", row.names = F, quote = FALSE)


### ===== 3.5 MATERNAL FAECES - MILK AT M3 ===== ###

rownames(selmatfmilkM3) <- c(paste0(selmatfmilkM3$sample_origin_type, "_", selmatfmilkM3$pair_ID))
matf_sample_ASVs_M3 <- selmatfmilkM3[,20:ncol(matfmilkM3)] #328 x 5391

matf_ASV_final_list_M3 <- list()

for (i in 1:nrow(matf_sample_ASVs_M3)) {
  print(paste0("Starting analysis for sample ", i))
  
  print("Getting sample names")
  sample_name <- rownames(matf_sample_ASVs_M3)[i]
  
  print("Getting sequences")
  ASV <- gsub(".*_", "",  colnames(matf_sample_ASVs_M3)[which(matf_sample_ASVs_M3[i,] != 0)]) #sequence
  
  print("Getting bacterial names")
  name_full <- gsub("_ASV.*", "\\1", colnames(matf_sample_ASVs_M3)[which(matf_sample_ASVs_M3[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everything after genus level
  name <- gsub(".*g__", "g__", name) #remove everything before genus level (get genus name)
  
  print("Getting bacterial names if they were unclassified on genus level (if required)")
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  
  print("Creating final ASV list")
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  matf_ASV_final_list_M3 <- append(matf_ASV_final_list_M3, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(matf_ASV_final_list_M3)
# all represent classified genera (expected as I excluded unclassified ones before)

# Create a dataframe with the counts of each genus for each participant 
matf_names_genera_M3 <- gsub(".*g__", "", grep("g__", names(matf_ASV_final_list_M3), value = T))
matf_names_genera_M3 <- unique(matf_names_genera_M3) #410 unique genera
matf_genera_counts_M3 <- data.frame(matrix(ncol = length(rownames(matf_sample_ASVs_M3)) + 1, nrow = length(matf_names_genera_M3)))
colnames(matf_genera_counts_M3) <- c("genera", rownames(matf_sample_ASVs_M3))
matf_genera_counts_M3$genera <- matf_names_genera_M3

for (i in 2:ncol(matf_genera_counts_M3)) {
  microbes <- names(matf_ASV_final_list_M3)[grep(paste0(colnames(matf_genera_counts_M3)[i],"\\."), names(matf_ASV_final_list_M3))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  matf_genera_counts_M3[i] <- counts$Freq[match(matf_genera_counts_M3$genera,counts$genera)]
}
matf_genera_counts_M3[is.na(matf_genera_counts_M3)] <- 0 #set NAs to 0

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
matf_genera_in_pairs_M3 <- list()
IDs <- matfmilkM3IDs # IDs <- unique(metadata_M3$pair_ID)

for (i in IDs) {
  print("Getting index")
  index <- which(IDs == i)
  print("Getting pattern")
  pattern <- paste0("_",i, "$")
  print("Getting pair counts")
  pair_counts <- matf_genera_counts_M3[, grep(pattern, colnames(matf_genera_counts_M3))]
  print("Selecting genera deteced in at least 1 pair")
  matf_genera_in_pairs_M3 [[index]] <- matf_genera_counts_M3$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair 
}

matf_genera_in_pairs_M3 <- matf_genera_in_pairs_M3[!sapply(matf_genera_in_pairs_M3,is.null)] # removes null elements
names(matf_genera_in_pairs_M3) <- IDs

# Get names of the shared genera present in 1 or more pairs
matf_shared_genera_M3 <- names(which(table(unlist(matf_genera_in_pairs_M3)) > 0)) #change this to a diff number if I want to be stricter on which genera to include
# matf_shared_genera_5_pairs_M3 <- names(which(table(unlist(matf_genera_in_pairs_M3)) > 4))

# Get all ASVs of selected genera
names(matf_ASV_final_list_M3) <- make.unique(names(matf_ASV_final_list_M3), sep = '_') # distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matf_matching_ASVs_M3 <- matf_ASV_final_list_M3[unique (grep(paste(matf_shared_genera_M3,collapse="|"), # get all ASVs belonging to shared genera
                                                   names(matf_ASV_final_list_M3)))]

## save results 
write.table(matf_genera_counts_M3,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/M3_maternal_faeces_milk_genera_counts.txt", sep = "\t", row.names = T, quote = FALSE)
saveRDS(matf_genera_in_pairs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/M3_maternal_faeces_milk_genera_counts_pairs.rds")
saveRDS(matf_matching_ASVs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/M3_maternal_faeces_milk_ASVs_shared_genera.rds")
saveRDS(matf_shared_genera_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/M3_maternal_faeces_milk_ASVs_shared_genera_at_least_1_pair.rds")
saveRDS(matf_shared_genera_5_pairs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/M3_maternal_faeces_milk_ASVs_shared_genera_at_least_5_pairs.rds")

# write.table(matf_genera_counts_M3,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/matfaecesM3/M3_maternal_faeces_milk_genera_counts.txt", sep = "\t", row.names = T, quote = FALSE)
# saveRDS(matf_genera_in_pairs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/matfaecesM3/M3_maternal_faeces_milk_genera_counts_pairs.rds")
# saveRDS(matf_matching_ASVs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/matfaecesM3/M3_maternal_faeces_milk_ASVs_shared_genera.rds")
# saveRDS(matf_shared_genera_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/matfaecesM3/M3_maternal_faeces_milk_ASVs_shared_genera_at_least_1_pair.rds")
# saveRDS(matf_shared_genera_5_pairs_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/matfaecesM3/M3_maternal_faeces_milk_ASVs_shared_genera_at_least_5_pairs.rds")
# write.table(matf_shared_genera_M3,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/matfaecesM3/M3_maternal_faeces_milk_ASVs_shared_genera_at_least_1_pair.txt", sep = "\t", row.names = F, quote = FALSE)
# write.table(matf_shared_genera_5_pairs_M3,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/matfaecesM3/M3_maternal_faeces_milk_ASVs_shared_genera_at_least_5_pairs.txt", sep = "\t", row.names = F, quote = FALSE)


### ===== 3.6 OVERVIEW OF POTENTIALLY SHARED GENERA ===== ###

### POTENTIALLY SHARED GENERA IN AT LEAST 1 PAIR ###

## combine potentially shared genera for milk - infant faeces from different time point into 1 vector
shared_genera_milkinff <- unique(c(shared_genera_M1, shared_genera_M2, shared_genera_M3, shared_genera_M6))
# -> 51 genera are potentially shared between milk - infant faeces in ≥1 pair at ≥1 time point

## potentially shared genera for maternal faeces - milk at 3 months
matf_shared_genera_M3
# -> 31 genera are potentially shared between maternal faeces - milk in ≥1 pair at 3 months postpartum

## combine all potentially shared bacterial genera into 1 vector
shared_genera <- unique(c(shared_genera_milkinff, matf_shared_genera_M3)) #69


df_shared_genera_1_pair <- data.frame(matrix(ncol=length(shared_genera), nrow=5))
colnames(df_shared_genera_1_pair) <- shared_genera
rownames(df_shared_genera_1_pair) <- c("milk-infant_faeces_M1" ,"milk-infant_faeces_M2", "milk-infant_faeces_M3", "milk-infant_faeces_M6", "maternal_faeces-milk_M3")

for (i in colnames(df_shared_genera_1_pair)){
  #milkinff M1
  if (i %in% shared_genera_M1){
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="milk-infant_faeces_M1", i] <- 1
  }else{
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="milk-infant_faeces_M1", i] <- 0
  }
  
  #milkinff M2
  if (i %in% shared_genera_M2){
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="milk-infant_faeces_M2", i] <- 1
  }else{
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="milk-infant_faeces_M2", i] <- 0
  }
  
  #milkinff M3
  if (i %in% shared_genera_M3){
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="milk-infant_faeces_M3", i] <- 1
  }else{
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="milk-infant_faeces_M3", i] <- 0
  }
  
  #milkinff M6
  if (i %in% shared_genera_M6){
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="milk-infant_faeces_M6", i] <- 1
  }else{
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="milk-infant_faeces_M6", i] <- 0
  }
  
  #matfmilk M3
  if (i %in% matf_shared_genera_M3){
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="maternal_faeces-milk_M3", i] <- 1
  }else{
    df_shared_genera_1_pair[rownames(df_shared_genera_1_pair)=="maternal_faeces-milk_M3", i] <- 0
  }
}

## transpose table
tdf_shared_genera_1_pair <- as.data.frame(t(df_shared_genera_1_pair))
tdf_shared_genera_1_pair$genus <- rownames(tdf_shared_genera_1_pair)
rownames(tdf_shared_genera_1_pair) <- 1:nrow(tdf_shared_genera_1_pair)
tdf_shared_genera_1_pair <- tdf_shared_genera_1_pair[order(tdf_shared_genera_1_pair$genus),c(6,1:5)]

## save table
write.table(tdf_shared_genera_1_pair,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/potentially_shared_genera_at_least_1_pair_all_sample_types_and_time_points.txt", sep = "\t", row.names = F, quote = FALSE)


### POTENTIALLY SHARED GENERA IN AT LEAST 5 PAIRS ###

shared_genera_M1_5_pairs <- names(which(table(unlist(genera_in_pairs_M1)) > 4)) #change this to a diff number if I want to be stricter on which genera to include
shared_genera_M2_5_pairs <- names(which(table(unlist(genera_in_pairs_M2)) > 4)) #change this to a diff number if I want to be stricter on which genera to include
shared_genera_M3_5_pairs <- names(which(table(unlist(genera_in_pairs_M3)) > 4)) #change this to a diff number if I want to be stricter on which genera to include
shared_genera_M6_5_pairs <- names(which(table(unlist(genera_in_pairs_M6)) > 4)) #change this to a diff number if I want to be stricter on which genera to include

matf_shared_genera_M3_5_pairs <- names(which(table(unlist(matf_genera_in_pairs_M3)) > 4)) #change this to a diff number if I want to be stricter on which genera to include

## combine potentially shared genera for milk - infant faeces from different time point into 1 vector
shared_genera_milkinff_5_pairs <- unique(c(shared_genera_M1_5_pairs, shared_genera_M2_5_pairs, shared_genera_M3_5_pairs, shared_genera_M6_5_pairs))
# -> 20 genera are potentially shared between milk - infant faeces in ≥5 pairs at ≥1 time point
# shared_genera_milkinff_5_pairs
# [1] "Acinetobacter"             "Actinomyces"              
# [3] "Atopobium"                 "Bacteroides"              
# [5] "Bifidobacterium"           "Corynebacterium"          
# [7] "Cutibacterium"             "Enterobacter"             
# [9] "Enterococcus"              "Escherichia_Shigella"     
# [11] "Gemella"                   "Haemophilus"              
# [13] "Lactobacillus"             "Prevotella"               
# [15] "Rothia"                    "Staphylococcus"           
# [17] "Streptococcus"             "Veillonella"              
# [19] "Clostridium_sensu_stricto" "Stenotrophomonas"

# for M1, add:
# setdiff(shared_genera_milkinff_5_pairs, shared_genera_M1_5_pairs)
# [1] "Clostridium_sensu_stricto" "Stenotrophomonas"

# for M2, add:
# setdiff(shared_genera_milkinff_5_pairs, shared_genera_M2_5_pairs)
# [1] "Acinetobacter"             "Actinomyces"              
# [3] "Atopobium"                 "Bacteroides"              
# [5] "Bifidobacterium"           "Corynebacterium"          
# [7] "Cutibacterium"             "Enterobacter"             
# [9] "Enterococcus"              "Escherichia_Shigella"     
# [11] "Gemella"                   "Lactobacillus"            
# [13] "Prevotella"                "Rothia"                   
# [15] "Staphylococcus"            "Clostridium_sensu_stricto"
# [17] "Stenotrophomonas" 
## no Atopobium, Corynebacterium, Clostridium_sensu_strictu and Stenotrophomonas at M2

# for M3, add:
# setdiff(shared_genera_milkinff_5_pairs, shared_genera_M3_5_pairs)
# [1] "Bacteroides"     "Corynebacterium"

# for M6, add:
# setdiff(shared_genera_milkinff_5_pairs, shared_genera_M6_5_pairs)
# [1] "Atopobium"                 "Corynebacterium"          
# [3] "Cutibacterium"             "Enterobacter"             
# [5] "Gemella"                   "Lactobacillus"            
# [7] "Prevotella"                "Rothia"                   
# [9] "Clostridium_sensu_stricto" "Stenotrophomonas"
## no Atopobium, Corynebacterium and Cutibacterium at M6

## potentially shared genera for maternal faeces - milk at 3 months
matf_shared_genera_M3_5_pairs
# -> 10 genera are potentially shared between maternal faeces - milk in ≥5 pairs at 3 months postpartum
# matf_shared_genera_M3_5_pairs
# [1] "Bacteroides"               "Bifidobacterium"          
# [3] "Blautia"                   "Clostridium_sensu_stricto"
# [5] "Escherichia_Shigella"      "Faecalibacterium"         
# [7] "Haemophilus"               "Prevotella"               
# [9] "Streptococcus"             "Veillonella"


## focus on milk - infant faeces
df_shared_genera_milkinff_5_pairs <- data.frame(matrix(ncol=length(shared_genera_milkinff_5_pairs), nrow=4))
colnames(df_shared_genera_milkinff_5_pairs) <- shared_genera_milkinff_5_pairs
rownames(df_shared_genera_milkinff_5_pairs) <- c("milk-infant_faeces_M1" ,"milk-infant_faeces_M2", "milk-infant_faeces_M3", "milk-infant_faeces_M6")

for (i in colnames(df_shared_genera_milkinff_5_pairs)){
  #milkinff M1
  if (i %in% shared_genera_M1){
    df_shared_genera_milkinff_5_pairs[rownames(df_shared_genera_milkinff_5_pairs)=="milk-infant_faeces_M1", i] <- 1
  }else{
    df_shared_genera_milkinff_5_pairs[rownames(df_shared_genera_milkinff_5_pairs)=="milk-infant_faeces_M1", i] <- 0
  }
  
  #milkinff M2
  if (i %in% shared_genera_M2){
    df_shared_genera_milkinff_5_pairs[rownames(df_shared_genera_milkinff_5_pairs)=="milk-infant_faeces_M2", i] <- 1
  }else{
    df_shared_genera_milkinff_5_pairs[rownames(df_shared_genera_milkinff_5_pairs)=="milk-infant_faeces_M2", i] <- 0
  }
  
  #milkinff M3
  if (i %in% shared_genera_M3){
    df_shared_genera_milkinff_5_pairs[rownames(df_shared_genera_milkinff_5_pairs)=="milk-infant_faeces_M3", i] <- 1
  }else{
    df_shared_genera_milkinff_5_pairs[rownames(df_shared_genera_milkinff_5_pairs)=="milk-infant_faeces_M3", i] <- 0
  }
  
  #milkinff M6
  if (i %in% shared_genera_M6){
    df_shared_genera_milkinff_5_pairs[rownames(df_shared_genera_milkinff_5_pairs)=="milk-infant_faeces_M6", i] <- 1
  }else{
    df_shared_genera_milkinff_5_pairs[rownames(df_shared_genera_milkinff_5_pairs)=="milk-infant_faeces_M6", i] <- 0
  }

}

tdf_shared_genera_milkinff_5_pairs <- as.data.frame(t(df_shared_genera_milkinff_5_pairs))
tdf_shared_genera_milkinff_5_pairs$genus <- rownames(tdf_shared_genera_milkinff_5_pairs)
rownames(tdf_shared_genera_milkinff_5_pairs) <- 1:nrow(tdf_shared_genera_milkinff_5_pairs)
tdf_shared_genera_milkinff_5_pairs <- tdf_shared_genera_milkinff_5_pairs[order(tdf_shared_genera_milkinff_5_pairs$genus),c(5,1:4)]

write.table(tdf_shared_genera_milkinff_5_pairs,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/potentially_shared_genera_at_least_5_pairs_milk_infant_faeces_samples_all_time_points.txt", sep = "\t", row.names = F, quote = FALSE)


## combine milk-infant faeces and maternal faeces-milk
## combine all potentially shared bacterial genera into 1 vector
shared_genera_5_pairs <- unique(c(shared_genera_milkinff_5_pairs, matf_shared_genera_M3_5_pairs)) #22


df_shared_genera_5_pairs <- data.frame(matrix(ncol=length(shared_genera_5_pairs), nrow=5))
colnames(df_shared_genera_5_pairs) <- shared_genera_5_pairs
rownames(df_shared_genera_5_pairs) <- c("milk-infant_faeces_M1" ,"milk-infant_faeces_M2", "milk-infant_faeces_M3", "milk-infant_faeces_M6", "maternal_faeces-milk_M3")

for (i in colnames(df_shared_genera_5_pairs)){
  #milkinff M1
  if (i %in% shared_genera_M1){
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="milk-infant_faeces_M1", i] <- 1
  }else{
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="milk-infant_faeces_M1", i] <- 0
  }
  
  #milkinff M2
  if (i %in% shared_genera_M2){
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="milk-infant_faeces_M2", i] <- 1
  }else{
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="milk-infant_faeces_M2", i] <- 0
  }
  
  #milkinff M3
  if (i %in% shared_genera_M3){
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="milk-infant_faeces_M3", i] <- 1
  }else{
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="milk-infant_faeces_M3", i] <- 0
  }
  
  #milkinff M6
  if (i %in% shared_genera_M6){
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="milk-infant_faeces_M6", i] <- 1
  }else{
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="milk-infant_faeces_M6", i] <- 0
  }
  
  #matfmilk M3
  if (i %in% matf_shared_genera_M3){
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="maternal_faeces-milk_M3", i] <- 1
  }else{
    df_shared_genera_5_pairs[rownames(df_shared_genera_5_pairs)=="maternal_faeces-milk_M3", i] <- 0
  }
}

tdf_shared_genera_5_pairs <- as.data.frame(t(df_shared_genera_5_pairs))
tdf_shared_genera_5_pairs$genus <- rownames(tdf_shared_genera_5_pairs)
rownames(tdf_shared_genera_5_pairs) <- 1:nrow(tdf_shared_genera_5_pairs)
tdf_shared_genera_5_pairs <- tdf_shared_genera_5_pairs[order(tdf_shared_genera_5_pairs$genus),c(6,1:5)]

write.table(tdf_shared_genera_5_pairs,"/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/potentially_shared_genera_at_least_5_pairs_all_sample_types_and_time_points.txt", sep = "\t", row.names = F, quote = FALSE)


##### =========================== 4. GET THE MULTIPLE SEQUENCE ALIGNMENT (MSA) FOR EACH GENUS DETECTED AMONG PAIRED INDIVIDUALS =========================== #####

#***********************************************************************
# A3.Get the MSA for each genus that is detected among paired individuals
#***********************************************************************


### ===== 4.1 FUNCTION FOR MULTIPLE SEQUENCE ALIGNMENT (MSA) WITH ClustalW ===== ###

#********************
# Defining functions
#********************
#A. MSA_result: generate a MSA of the ASVs of selected genera using ClustalW -> replaced by ClustalO
# It outputs a FASTA file with ASVs sequences, a FASTA file with the alignment and a PDF with the alignment (including logos)
# output_dir : full path to directory where one folder will be stored for each genus with its results
# description: information about sequencing technology and pairs of samples to compare
MSA_result <- function(ASVs, names_ASVs, output_dir, description) {
  for (i in 1:length(names_ASVs)) {
    print("Setting genus name")
    genus <- names_ASVs [i] #name of the genus
    print(genus)
    
    print("Defining output directory")
    dir <- paste0(output_dir ,genus)
    dir.create(dir) #Create a directory for the genus
    setwd(dir)  #Move to the directory
    print(getwd())
    
    print("Selecting ASVs of the specified genus")
    genus_ASVs <- ASVs[grep (genus, names(ASVs))] #select ASVs of the specific genus
    # print(genus_ASVs)
    
    print("Setting filename for ASVs of the specified genus")
    filename <- paste0(paste(genus,description, sep = "_"), ".fa")
    print(filename)
    
    print("Writing .fasta file for ASVs of the specified genus")
    write.fasta(genus_ASVs, names(genus_ASVs), filename) #generate a FASTA formatted file
    
    # print("Running MSA")
    # # print(paste0(getwd(), "/", filename))
    # MSA <-  msa(filename, method = "ClustalW", type = "DNA", order = "input") #do the MSA
    # 
    # ASV_names_for_plot <- gsub("\\..*", "", names(genus_ASVs)) #get names of the ASVs sequences to be displayed later in the alignment and the tree
    # names(MSA@unmasked) <- ASV_names_for_plot
    # try({msaPrettyPrint(MSA, output="pdf", file = paste(genus, description, "MSA.pdf", sep= "_"), #get a PDF with the alignment
    #                     showNames ="left", shadingMode="similar",
    #                     shadingColors="blues", showLogo="top", logoColors="rasmol",
    #                     showLegend=T, askForOverwrite=FALSE)})
    # msa_align <- msaConvert(MSA, type="seqinr::alignment") #convert the alignment for later processing with seqinr
    # msa_align$nam <- ASV_names_for_plot #add names for the plot
    # filename_aln <- paste0(paste(genus,description, sep = "_"), "_aln.fa")
    # write.fasta(as.list(msa_align$seq), names(genus_ASVs), filename_aln) #write alignment fasta (as ALN file is empty)
  }
}


### ===== 4.2 PREPARE FILES FOR MSA ===== ###

MSA_result(ASVs=matching_ASVs_M1,
           names_ASVs=shared_genera_M1,
           output_dir="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M1/",
           description="M1_milk_infant_feces")


MSA_result(ASVs=matching_ASVs_M2,
           names_ASVs=shared_genera_M2,
           output_dir="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M2/",
           description="M2_milk_infant_feces")


MSA_result(ASVs=matching_ASVs_M3,
           names_ASVs=shared_genera_M3,
           output_dir="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M3/",
           description="M3_milk_infant_feces")


MSA_result(ASVs=matching_ASVs_M6,
           names_ASVs=shared_genera_M6,
           output_dir="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/M6/",
           description="M6_milk_infant_feces")


MSA_result(ASVs=matf_matching_ASVs_M3,
           names_ASVs=matf_shared_genera_M3,
           output_dir="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/M3/",
           description="M3_milk_maternal_feces")


## Notes:
## I manually selected the bacterial genera present in ≥5 pairs at ≥1 time point
## Ranko ran ClustalO instead of ClustalW
## Asier used snp-dists to get SNP distances


