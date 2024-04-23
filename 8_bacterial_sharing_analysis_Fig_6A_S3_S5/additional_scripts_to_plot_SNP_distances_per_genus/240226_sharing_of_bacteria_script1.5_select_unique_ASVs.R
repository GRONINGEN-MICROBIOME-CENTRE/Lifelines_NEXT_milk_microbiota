################################################################################################################
### INVESTIGATE SHARING OF BACTERIA BETWEEN MILK AND INFANT FAECES - SCRIPT 1 ##################################
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
# 2.2 MATERNAL FAECES - MILK : SELECT SAMPLES AND COLUMNS OF INTEREST
# 
# 3. GET THE GENERA (AND THEIR ASVs) PRESENT IN PAIRED SAMPLES
# 3.1 MILK - INFANT FAECES AT M1
# 3.2 MILK - INFANT FAECES AT M2
# 3.3 MILK - INFANT FAECES AT M3
# 3.4 MILK - INFANT FAECES AT M6
# 3.5 MATERNAL FAECES - MILK AT M3
# 3.6 OVERVIEW OF POTENTIALLY SHARED GENERA
# 
# 4. GET THE MULTIPLE SEQUENCE ALIGNMENT (MSA) FOR EACH GENUS DETECED AMONG PAIRED INDIVIDUALS
# 4.1 FUNCTION FOR MULTIPLE SEQUENCE ALIGNMENT (MSA) WITH ClustalW
# 4.2 PREPARE FILES FOR MSA


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

## combine all selected data back into 1 data frame
selmilkinffall <- as.data.frame(rbind(selmilkinffM1, selmilkinffM2, selmilkinffM3, selmilkinffM6)) #1018x5410


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


##### =========================== 3. MILK - INFANT FAECES: GET ALL UNIQUE GENERA (AND THEIR ASVs) PRESENT IN AT LEAST 5 PAIRED SAMPLES =========================== #####

## bacteria that were detected in ≥5 mother milk - infant faeces sample pairs at ≥1 time point
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


## create data frames only with the ASVs of the bacteria of interest
rownames(selmilkinffall) <- selmilkinffall$seq_16S_sample_ID

## function to obtain unique ASVs and save them in a list, per bacterial genus
obtain.unique.ASVs.per.genus.milkinff <- function(inputdata, bugname){
  
  ASV_list <- list()
  my_names <- c()
  
  ## create data frames only with the ASVs of the bacteria of interest
  print("Creating data frame my_df")
  my_df <- inputdata[,grep(bugname, colnames(inputdata))]
  print(dim(my_df))
  
  ## exclude ASVs that were not found in any milk or infant faecal sample
  print("Creating data frame my_df_bug")
  my_df_bug <- my_df[,colSums(my_df)!=0]
  print(dim(my_df_bug))
  
  ## per bacterium, extract all unique ASVs and save them in a list
  print("Creating vector with unique ASVs")
  my_ASVs <- unique(colnames(my_df_bug))
  
  print("Creating list with unique ASVs")
  for (i in 1:length(my_ASVs)){
    # extract ASV sequence
    ASV_list <- c(ASV_list, gsub(".*_", "",  my_ASVs[i]))
    
    # add name for ASV
    my_names <- c(my_names, paste0(bugname, "_ASV_", i))
  }
  
  print("Adding ASV names to list with unique ASVs")
  names(ASV_list) <- my_names
  # print(ASV_list)
  
  ## creating output directory for the genus
  print("Creating output directory")
  dir.create(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/all_milkinff_ASVs/", bugname)) #create directory for the genus
  
  ## save unique ASVs as FASTA formatted file
  print("Saving fasta file")
  write.fasta(ASV_list, names(ASV_list), paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/sharing_with_milk/all_milkinff_ASVs/", bugname, "/", bugname, "_ASVs.fa"))
  
  ## return list
  return(ASV_list)
  
}


Acinetobacter_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Acinetobacter")
Actinomyces_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Actinomyces")
Atopobium_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Atopobium")
Bacteroides_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Bacteroides")
Bifidobacterium_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Bifidobacterium")

Clostridium_sensu_stricto_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Clostridium_sensu_stricto")
Corynebacterium_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Corynebacterium")
Cutibacterium_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Cutibacterium")
Enterobacter_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Enterobacter")
Enterococcus_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Enterococcus")

Escherichia_Shigella_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Escherichia_Shigella")
Gemella_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Gemella")
Haemophilus_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Haemophilus")
Lactobacillus_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Lactobacillus")
Prevotella_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Prevotella")

Rothia_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Rothia")
Staphylococcus_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Staphylococcus")
Stenotrophomonas_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Stenotrophomonas")
Streptococcus_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Streptococcus")
Veillonella_ASV_list <- obtain.unique.ASVs.per.genus.milkinff(inputdata=selmilkinffall, bugname="Veillonella")


##### =========================== 4. MATERNAL FAECES - MILK: GET ALL UNIQUE GENERA (AND THEIR ASVs) PRESENT IN AT LEAST 5 PAIRED SAMPLES =========================== #####

## create data frames only with the ASVs of the bacteria of interest
rownames(selmatfmilkM3) <- selmatfmilkM3$seq_16S_sample_ID

## function to obtain unique ASVs and save them in a list, per bacterial genus
obtain.unique.ASVs.per.genus.matfmilk <- function(inputdata, bugname){
  
  ASV_list <- list()
  my_names <- c()
  
  ## create data frames only with the ASVs of the bacteria of interest
  print("Creating data frame my_df")
  my_df <- inputdata[,grep(bugname, colnames(inputdata))]
  print(dim(my_df))
  
  ## exclude ASVs that were not found in any milk or infant faecal sample
  print("Creating data frame my_df_bug")
  my_df_bug <- my_df[,colSums(my_df)!=0]
  print(dim(my_df_bug))
  
  ## per bacterium, extract all unique ASVs and save them in a list
  print("Creating vector with unique ASVs")
  my_ASVs <- unique(colnames(my_df_bug))
  
  print("Creating list with unique ASVs")
  for (i in 1:length(my_ASVs)){
    # extract ASV sequence
    ASV_list <- c(ASV_list, gsub(".*_", "",  my_ASVs[i]))
    
    # add name for ASV
    my_names <- c(my_names, paste0(bugname, "_ASV_", i))
  }
  
  print("Adding ASV names to list with unique ASVs")
  names(ASV_list) <- my_names
  # print(ASV_list)
  
  ## creating output directory for the genus
  print("Creating output directory")
  dir.create(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/", bugname)) #create directory for the genus
  
  ## save unique ASVs as FASTA formatted file
  print("Saving fasta file")
  write.fasta(ASV_list, names(ASV_list), paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/", bugname, "/", bugname, "_ASVs.fa"))
  
  ## return list
  return(ASV_list)
  
}

## bacteria that were detected in ≥5 maternal faeces - mother milk sample pairs at 3 months
# [1] "Bacteroides"               "Bifidobacterium"          
# [3] "Blautia"                   "Clostridium_sensu_stricto"
# [5] "Escherichia_Shigella"      "Faecalibacterium"         
# [7] "Haemophilus"               "Prevotella"               
# [9] "Streptococcus"             "Veillonella"


Bacteroides_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Bacteroides")
Bifidobacterium_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Bifidobacterium")
Blautia_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Blautia")
Clostridium_sensu_stricto_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Clostridium_sensu_stricto")
Escherichia_Shigella_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Escherichia_Shigella")

Faecalibacterium_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Faecalibacterium")
Haemophilus_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Haemophilus")
Prevotella_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Prevotella")
Streptococcus_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Streptococcus")
Veillonella_ASV_list <- obtain.unique.ASVs.per.genus.matfmilk(inputdata=selmatfmilkM3, bugname="Veillonella")



## Notes:
## I manually selected the bacterial genera present in ≥5 pairs at ≥1 time point
## Ranko ran ClustalO instead of ClustalW
## Asier used snp-dists to get SNP distances


