### ========================== MILK COMPOSITION PAPER - MICROBIOTA BY SAMPLE TYPE ========================== ###

## SCRIPT:           COMPARISON OF MICROBIAL COMPOSITON BETWEEN MILK, INFANT FAECES AND MATERNAL FAECES
## DESCRIPTION:      The data used here is from the dada2 pipeline WITH filtering to only include ASVs with lengths 400-431 bp and, for milk and faecal samples, data was decontaminated using SCRuB.
##                   I also previously sub-selected samples: I only include data from good sequencing runs and if a milk HMO measurement for the same mother-sample-time point was available.
##                   In addition, samples with <6500 clean reads after decontamination were excluded.
##                   Rare taxa (ASVs <30 reads across all samples of a sample type) and unclassified genera were excluded (aggregated their reads into 'Other').
##                   Data was analysed on the classified genus level only.
## AUTHORS:          Johanne Spreckels
## DATE OF CREATION: November 2023

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23):     R

## CONTENTS OF THIS FILE
# 1. DATA IMPORT AND DATA PREPARATION
# 1.1. IMPORT DATA
# 1.2. CHECK AND ENSURE CORRECT DATA STRUCTURE
# 1.3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
# 1.4. SELECT DATA OF INTEREST
# 
# 2. PREPARE DATA SET WITH INVERSE-RANK TRANSFORMED ALPHA DIVERSITY MEASURES
# 2.1. FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 2.2. INVERSE-RANK TRANSFORM ALPHA DIVERSITY MEASURES
# 
# 3. PREPARE DATA SET WITH CLR-TRANSFORMED RELATIVE BACTERIAL ABUNDANCES
# 3.1. CALCULATE RELATIVE ABUNDANCES (BASED ON GENERA)
# 3.2. REMOVE BACTERIAL GENERA NOT DETECTED IN SAMPLES FROM DATA SET
# 3.3. CLR-TRANSFORM RELATIVE ABUNDANCES
# 
# 4. BASIC DATA DESCRIPTION PER SAMPLE TYPE
# [4.1. ALL SAMPLES: MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA]
# 4.2. HUMAN MILK: MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA
# 4.3. MATERNAL FAECES: MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA
# 4.4. INFANT FAECES: MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA
# 4.5. COMBINE MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA FROM DIFFERENT SAMPLE TYPES
# 
# 5. IDENTIFY BACTERIAL GENERA WITH ≥10% PREVALENCE
# 5.1. HUMAN MILK: CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE
# 5.2. MATERNAL FAECES: CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE
# 5.3. INFANT FAECES: CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE
# 5.4. SUPPLEMENTARY TABLES SHOWING CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE PER SAMPLE TYPE (TABLES S9, S11, S13)
# 5.5. HUMAN MILK: PLOT CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE (FIGURE S2A)
# 5.6. INFANT FAECES: PLOT CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE (FIGURE S2B)
# 
# 6. RELATIVE BACTERIAL ABUNDANCES BY SAMPLE TYPE AND TIME POINT
# 6.1. PREPARE DATA FOR PLOTTING
# 6.2. RELATIVE ABUNDANCE PLOTS BY SAMPLE TYPE AND TIME POINT (FIGURE 4D)
# 
# 7. BETA DIVERSITY BY SAMPLE TYPE (AND TIME POINT): ADONIS AND PCoA PLOT (FIGURE 4A)
# 
# 8. PLOT ALPHA DIVERSITY BY SAMPLE TYPES (POOLED OR PER TIME POINT)
# 8.1. PREPARE DATA FOR PLOTTING
# [8.2. PLOT ALPHA DIVERSITY BY SAMPLE TYPE (POOLED TIME POINTS)]
# 8.3. PLOT ALPHA DIVERSITY BY SAMPLE TYPE AND TIME POINT (FIGURE 4C)
# [8.4. PLOT ALPHA DIVERSITY FOR MILK AND INFANT FAECES AND TIME POINT]
# 
# 9. COMPARE ALPHA DIVERSITY BETWEEN SAMPLE TYPES (POOLED OR PER TIME POINT)
# 9.1. FUNCTION TO COMPARE ALPHA DIVERSITY BY SAMPLE TYPE (POOLED TIME POINTS)
# 9.2. COMPARE ALPHA DIVERSITY BETWEEN SAMPLE TYPES (PAIR-WISE COMPARISONS, POOLED TIME POINTS)
# 9.3. COMPARE ALPHA DIVERSITY BETWEEN MILK AND INFANT FAECES PER TIME POINT
# 9.4. COMPARE ALPHA DIVERSITY BETWEEN SAMPLE TYPES (PAIR-WISE COMPARISONS, MATERNAL FAECES M3 VS. HUMAN MILK / INFANT FAECES M6)
# 
# 10. CHECK THE EFFECT OF TIME ON MILK ALPHA DIVERSITY AND MILK RELATIVE ABUNDANCES
# 10.1. MODEL FUNCTION WITH CORRECTION FOR DNA ISOLATION BATCH/METHOD (FIXED EFFECT), CLEAN READ COUNT (FIXED EFFECT) AND INDIVIDUAL ID (RANDOM EFFECT)
# 10.2. COMPARE MILK ALPHA DIVERSITY BY NUMERIC TIME POINTS
# 10.3. COMPARE MILK RELATIVE ABUNDANCES BY NUMERIC TIME POINTS (TABLE S10)
# 
# 11. PLOT EFFECT OF TIME ON MILK RELATIVE ABUNDANCES
# 11.1. FUNCTION TO PLOT SIGNIFICANT EFFECTS OF TIME ON RELATIVE ABUNDANCES
# 11.2. PLOT SIGNIFICANT EFFECTS OF TIME ON MILK RELATIVE ABUNDANCES
# 11.3. PLOT ESTIMATES FOR EACH BACTERIUM (USED FOR FIGURE 4E)
# [11.4. PLOT MEAN RELATIVE ABUNDANCE FOR EACH BACTERIUM BY TIME POINT]
# 
# 12. CHECK THE EFFECT OF TIME ON INFANT FAECAL ALPHA DIVERSITY AND INFANT FAECAL RELATIVE ABUNDANCES
# 12.1. MODEL FUNCTION WITH CORRECTION FOR DNA ISOLATION BATCH/METHOD (FIXED EFFECT), CLEAN READ COUNT (FIXED EFFECT) AND INDIVIDUAL ID (RANDOM EFFECT)
# 12.2. COMPARE INFANT FAECAL ALPHA DIVERSITY BY NUMERIC TIME POINTS
# 12.3. COMPARE INFANT FAECAL RELATIVE ABUNDANCES BY NUMERIC TIME POINTS (TABLE S12)
# 
# 13. PLOT EFFECT OF TIME ON INFANT FAECAL RELATIVE ABUNDANCES
# 13.1. FUNCTION TO PLOT SIGNIFICANT EFFECTS OF TIME ON RELATIVE ABUNDANCES
# 13.2. PLOT SIGNIFICANT EFFECTS OF TIME ON INFANT FAECAL RELATIVE ABUNDANCES
# 13.3. PLOT ESTIMATES FOR EACH BACTERIUM (USED FOR FIGURE 4E)
# 
# 14. COMBINED PLOT FOR EFFECT OF TIME ON MILK AND INFANT FAECAL RELATIVE ABUNDANCE ESTIMATES FOR EACH BACTERIUM (FIGURE 4E)
# 
# 15. HUMAN MILK BETA DIVERSITY BY TIME POINT AND BY INDIVIDUAL ID
# 15.1. RUN ADONIS BY TIME POINT (FACTOR) AND INDIVIDUAL ID WITH CORRECTION FOR DNA ISOLATION BATCH/METHOD AND CLEAN READ COUNT
# 15.2. MILK PLOT PCOA BY TIME POINT (FACTOR)
# 15.3. MILK PLOT PCOA BY INDIVIDUAL ID
# 
# 16. INFANT FAECAL BETA DIVERSITY BY TIME POINT AND BY INDIVIDUAL ID
# 16.1. RUN ADONIS BY TIME POINT (FACTOR) AND INDIVIDUAL ID WITH CORRECTION FOR DNA ISOLATION BATCH/METHOD AND CLEAN READ COUNT
# 16.2. INFANT FAECES PLOT PCOA BY TIME POINT (FACTOR)
# 16.3. INFANT FAECES PLOT PCOA BY INDIVIDUAL ID
# 
# 17. PREPARE DATA FRAME FOR COMPARISON OF AITCHISON DISTANCES THAT SHOULD INCLUDE MATERNAL DATA FOR SECOND TWIN INFANT
# 17.1. SELECT NEW DATA FRAME
# 17.2. CALCULATE RELATIVE ABUNDANCES (BASED ON GENERA)
# 17.3. CLR-TRANSFORM RELATIVE ABUNDANCES
# 17.4. ADD COLUMN COMBINING SAMPLE TYPE AND TIME POINT

# 18. HUMAN MILK - INFANT FAECES AITCHISON DISTANCES BY TIME
# 18.1. CREATE DATA FRAME WITH ONLY CORRESPONDING HUMAN MILK AND INFANT FAECES DATA
# 18.2. CALCULATE AITCHISON DISTANCES
# 18.3. FUNCTION TO EXTRACT CORRESPONDING AITCHISON DISTANCES
# 18.4. PLOT HUMAN MILK - INFANT FAECES AITCHISON DISTANCES PER TIME POINT
# 18.5. HISTOGRAMS FOR AITCHISON DISTANCES
# 18.6. COMPARE HUMAN MILK - INFANT FAECES AITCHISON DISTANCES BY TIME
# 
# 19. INFANT FAECES - MATERNAL FAECES AITCHISON DISTANCES BY TIME
# 19.1. CREATE DATA FRAME WITH ONLY CORRESPONDING HUMAN MILK AND INFANT FAECES DATA
# 19.2. CALCULATE AITCHISON DISTANCES
# 19.3. FUNCTION TO EXTRACT CORRESPONDING AITCHISON DISTANCES
# 19.4. PLOT INFANT FAECES - MATERNAL FAECES AITCHISON DISTANCES PER TIME POINT
# 19.5. HISTOGRAMS FOR AITCHISON DISTANCES
# 19.6. COMPARE INFANT FAECES - MATERNAL FAECES AITCHISON DISTANCES BY TIME
# 
# 20. MATERNAL FAECES - HUMAN MILK AITCHISON DISTANCES BY TIME
# 20.1. CREATE DATA FRAME WITH ONLY CORRESPONDING HUMAN MILK AND MATERNAL FAECES DATA
# 20.2. CALCULATE AITCHISON DISTANCES
# 20.3. FUNCTION TO EXTRACT CORRESPONDING AITCHISON DISTANCES
# 20.4. PLOT MATERNAL FAECES - HUMAN MILK AITCHISON DISTANCES PER TIME POINT
# 20.5. HISTOGRAMS FOR AITCHISON DISTANCES
# 20.6. COMPARE MATERNAL FAECES - HUMAN MILK AITCHISON DISTANCES BY TIME
# 
# 21. CREATE COMBINED PLOT OF AITCHISON DISTANCES BETWEEN CORRESPONDING SAMPLES BY SAMPLE TYPE COMPARISON AND TIME (FIGURE 4B)
# 21.1. MERGE AND PREPARE DATA FRAMES
# 21.2. COMBINED PLOT
# 21.3. COMBINE STATISTICS RESULTS AND ADJUST FOR MULTIPLE TESTING
# 21.4. COMPARE PAIR-WISE DISTANCES BETWEEN CORRESPONDING SAMPLES FROM ALL TIME POINTS


##### =========================== 1. DATA IMPORT AND DATA PREPARATION =========================== #####

### ===== 1.1. IMPORT DATA ===== ###

## import files, phenotypes+microbiota data on genus level
data <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/231129_milk_composition_paper_sel_phenotypes_sel_microbiota_data_genera_n1515.rds")
#1515x1154

# Notes:
# This file contains the following:
# - Phenotypes linked to milk, maternal faecal and infant faecal microbiota data.
#   The file also includes real measured levels of 24 single HMOs and 4 grouped HMOs in µg/ml and HMO-based Le and Se status and milk groups as phenotypes.
#   !! Note that maternal phenotype and HMO data was duplicated for twins. This is indicated by ";1" (single baby/first baby) and ";2" (twin baby/duplicated data) in the mother_ID and mother_sample_ID.
#      -> For association with maternal phenotypes (i.e. non-duplicated data), grep only ";1"
#      -> For association with infant phenotypes, use both ";1" and ";2"


### ===== 1.2. CHECK AND ENSURE CORRECT DATA STRUCTURE ===== ###

## check data structure
str(data[,1:50])
str(data[,51:100])
str(data[,100:150])
str(data[,151:200])
str(data[,200:250])
str(data[,251:300])
str(data[,300:350])
str(data[,351:400])
str(data[,400:450])
str(data[,451:500])
str(data[,500:530])
str(data[,c(531,1153,1154)])
# -> looks good

## check sample numbers per sample type in data frame 'data' -> note that this data frame includes duplicated maternal data if a mother had twins!
dim(data)
table(data$sample_origin_type)
# 737 human milk
# 552 infant faeces
# 172 maternal faeces
# 7 isolation negative controls
# 25 lib prep negative controls
# 1 isolation positive control
# 21 lib prep positive controls


### ===== 1.3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES ===== ###

# change order for milk groups
data$mother_milk_HMO_milk_group <- factor(data$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


### ===== 1.4. SELECT DATA OF INTEREST ===== ###

## select only biological samples, exclude 54 control samples
df <- data[-grep("control", data$sample_origin),] #1461x1154

## split dataframe by sample origin
m <- df[df$sample_origin=="mother",]
i <- df[df$sample_origin=="infant",]

## from data frame m (mothers), select every mother only 1x (i.e. remove duplication of data generated for twin babies) and select IDs and maternal phenotype columns
m2 <- m[grep(";1", m$mother_sample_ID),]

## re-merge data frames from mothers and infants
df2 <- as.data.frame(rbind(m2, i)) #1447x1154


## check sample type distribution
table(df2$sample_origin_type, useNA="ifany")
# 725 human milk
# 552 infant faeces
# 170 maternal faeces

## Not included in this data frame:
# 7 isolation negative controls
# 25 lib prep negative controls
# 1 isolation positive control
# 21 lib prep positive controls


## drop empty levels
df2 <- droplevels(df2)


##### =========================== 2. PREPARE DATA SET WITH INVERSE-RANK TRANSFORMED ALPHA DIVERSITY MEASURES  =========================== #####

### ===== 2.1. FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 2.2. INVERSE-RANK TRANSFORM ALPHA DIVERSITY MEASURES ===== ###

## save data frame under new name
sinvr <- df2

## Infant faecal alpha diversity measures
sinvr[,grep("alpha_div", colnames(sinvr))] <- as.data.frame(apply(sinvr[,grep("alpha_div", colnames(sinvr))], 2, invrank))

## update colnames to show transformation
colnames(sinvr)[528:530] <- c(paste0(colnames(sinvr)[528:530], "_invr"))


##### =========================== 3. PREPARE DATA SET WITH CLR-TRANSFORMED RELATIVE BACTERIAL ABUNDANCES  =========================== #####

### ===== 3.1. CALCULATE RELATIVE ABUNDANCES (BASED ON GENERA) ===== ###

## save only columns with absolute bacterial abundances
rownames(sinvr) <- paste0(sinvr$sample_origin_type, "__", sinvr$time_point, "__", sinvr$seq_16S_sample_ID)
sinvrbact <- sinvr[,531:ncol(sinvr)] #1447x624

## calculate relative bacterial abundances
sinvrbact_relab <- (sinvrbact/rowSums(sinvrbact))
table(rowSums(sinvrbact_relab), useNA="ifany") #all samples have 1

## merge back with metadata
dfRA <- merge(sinvr[,1:530], sinvrbact_relab, by="row.names")
rownames(dfRA) <- dfRA$Row.names
dfRA <- dfRA[,-1]


### ===== 3.2. REMOVE BACTERIAL GENERA NOT DETECTED IN SAMPLES FROM DATA SET ===== ###

## exclude column 'Other'
sinvrbact_relab_v2 <- sinvrbact_relab[,-624]

## exclude columns for bacteria that were not identified in any sample
sinvrbact_relab_v2 <- sinvrbact_relab_v2[,colSums(sinvrbact_relab_v2)>0]
dim(sinvrbact_relab_v2)
#1447x611 -> we identified 611 unique classified bacterial genera in the different maternal and infant samples


### ===== 3.3. CLR-TRANSFORM RELATIVE ABUNDANCES ===== ###

# function for clr transformation
do_clr_externalWeighting = function(interest_matrix, core_matrix) {
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  #estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

dfRAclr <- dfRA[,-1154] #exclude column 'Other'
dfRAclr[,531:ncol(dfRAclr)] <- as.data.frame(do_clr_externalWeighting(dfRAclr[,531:ncol(dfRAclr)], dfRAclr[,531:ncol(dfRAclr)]))

## update colnames to show transformation
colnames(dfRAclr)[531:ncol(dfRAclr)] <- c(paste0(colnames(dfRAclr)[531:ncol(dfRAclr)], "_clr"))


##### =========================== 4. BASIC DATA DESCRIPTION PER SAMPLE TYPE =========================== #####

## Continue working with the sinvrbact_relab_v2 data frame, which only contains the columns with counts for bacteria detected in the samples and which does not include the column 'Other' 


# ### ===== 4.1. ALL SAMPLES: MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA ===== ###
# 
# ## PREVALENCE
# ## 1447 all samples, 10% of 1447 -> 145 samples
# all_prev50 <- sinvrbact_relab_v2[colSums(sinvrbact_relab_v2>0)>=145]
# length(colnames(all_prev50)) # -> 50 genera are present in at least 10% of the samples
# colnames(all_prev50)
# 
# ## ABUNDANCE
# ## using mean, 1%
# meanvector <- c()
# for (i in 1:ncol(all_prev50)){meanvector[i] <- mean(all_prev50[,i])}
# all_prev50[1448,] <- meanvector
# all_prev50_mean1 <- all_prev50[,all_prev50[1448,]>=0.01]
# all_prev50 <- all_prev50[-1448,]
# all_prev50_mean1 <- all_prev50_mean1[-1448,]
# length(colnames(all_prev50_mean1)) # -> 12 classified bacterial genera are present in at least 50% of the samples with at least 1% mean relative abundance
# colnames(all_prev50_mean1)
# # [1] "g__Acinetobacter"        "g__Bacteroides"         
# # [3] "g__Bifidobacterium"      "g__Corynebacterium"     
# # [5] "g__Enterobacter"         "g__Escherichia_Shigella"
# # [7] "g__Haemophilus"          "g__Parabacteroides"     
# # [9] "g__Pseudomonas"          "g__Staphylococcus"      
# # [11] "g__Streptococcus"        "g__Veillonella"


### ===== 4.2. HUMAN MILK: MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA ===== ###

## PREVALENCE
## 725 milk samples, 50% of 725 -> 363 samples
milk_sinvrbact_relab_v2 <- sinvrbact_relab_v2[grep("mother_human_milk", rownames(sinvrbact_relab_v2)),]
milk_prev50 <- milk_sinvrbact_relab_v2[,colSums(milk_sinvrbact_relab_v2>0)>=363]
length(colnames(milk_prev50)) # -> 6 genera are present in at least 50% of the samples
# colnames(milk_prev50)

## ABUNDANCE
## using mean, 1%
meanvector <- c()
for (i in 1:ncol(milk_prev50)){meanvector[i] <- mean(milk_prev50[,i])}
milk_prev50[726,] <- meanvector
milk_prev50_mean1 <- milk_prev50[,milk_prev50[726,]>=0.01]
milk_prev50 <- milk_prev50[-726,]
milk_prev50_mean1 <- milk_prev50_mean1[-726,]
length(colnames(milk_prev50_mean1)) # -> 6 classified bacterial genera are present in at least 50% of the samples with at least 1% mean relative abundance
colnames(milk_prev50_mean1)
# [1] "g__Acinetobacter"   "g__Corynebacterium" "g__Gemella"        
# [4] "g__Staphylococcus"  "g__Streptococcus"   "g__Veillonella"   


### ===== 4.3. MATERNAL FAECES: MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA ===== ###

## PREVALENCE
## 170 maternal faecal samples, 50% of 170 -> 85 samples]
matf_sinvrbact_relab_v2 <- sinvrbact_relab_v2[grep("mother_faeces", rownames(sinvrbact_relab_v2)),]
matf_prev50 <- matf_sinvrbact_relab_v2[,colSums(matf_sinvrbact_relab_v2>0)>=85]
length(colnames(matf_prev50)) # -> 37 genera are present in at least 50% of the samples
# colnames(matf_prev50)

## ABUNDANCE
## using mean, 1%
meanvector <- c()
for (i in 1:ncol(matf_prev50)){meanvector[i] <- mean(matf_prev50[,i])}
matf_prev50[171,] <- meanvector
matf_prev50_mean1 <- matf_prev50[,matf_prev50[171,]>=0.01]
matf_prev50 <- matf_prev50[-171,]
matf_prev50_mean1 <- matf_prev50_mean1[-171,]
length(colnames(matf_prev50_mean1)) # -> 18 classified bacterial genera are present in at least 50% of the samples with at least 1% mean relative abundance
colnames(matf_prev50_mean1)
# [1] "g__Agathobacter"                  "g__Akkermansia"                  
# [3] "g__Alistipes"                     "g__Bacteroides"                  
# [5] "g__Barnesiella"                   "g__Bifidobacterium"              
# [7] "g__Blautia"                       "g__Christensenellaceae_R_7_group"
# [9] "g__Collinsella"                   "g__Coprococcus"                  
# [11] "g__Dialister"                     "g__Faecalibacterium"             
# [13] "g__Fusicatenibacter"              "g__Parabacteroides"              
# [15] "g__Roseburia"                     "g__Ruminococcus"                 
# [17] "g__Subdoligranulum"               "g__UCG_002"  


### ===== 4.4. INFANT FAECES: MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA ===== ###

## PREVALENCE
## 552 infant faecal samples, 50% of 552 -> 276 samples
inff_sinvrbact_relab_v2 <- sinvrbact_relab_v2[grep("infant_faeces", rownames(sinvrbact_relab_v2)),]
inff_prev50 <- inff_sinvrbact_relab_v2[,colSums(inff_sinvrbact_relab_v2>0)>=276]
length(colnames(inff_prev50)) # -> 6 genera are present in at least 50% of the samples
# colnames(inff_prev50)

## ABUNDANCE
## using mean, 1%
meanvector <- c()
for (i in 1:ncol(inff_prev50)){meanvector[i] <- mean(inff_prev50[,i])}
inff_prev50[553,] <- meanvector
inff_prev50_mean1 <- inff_prev50[,inff_prev50[553,]>=0.01]
inff_prev50 <- inff_prev50[-553,]
inff_prev50_mean1 <- inff_prev50_mean1[-553,]
length(colnames(inff_prev50_mean1)) # -> 6 classified bacterial genera are present in at least 50% of the samples with at least 1% mean relative abundance
colnames(inff_prev50_mean1)
# [1] "g__Bacteroides"          "g__Bifidobacterium"     
# [3] "g__Escherichia_Shigella" "g__Haemophilus"         
# [5] "g__Streptococcus"        "g__Veillonella" 


### ===== 4.5. COMBINE MOST PREVALENT & MOST ABUNDANT CLASSIFIED BACTERIAL GENERA FROM DIFFERENT SAMPLE TYPES ===== ###

## combine in 1 vector
comb_prev50_mean1 <- unique(c(colnames(milk_prev50_mean1), colnames(matf_prev50_mean1), colnames(inff_prev50_mean1)))
comb_prev50_mean1
length(comb_prev50_mean1) #26


##### =========================== 5. IDENTIFY BACTERIAL GENERA WITH ≥10% PREVALENCE =========================== #####

### ===== 5.1. HUMAN MILK: CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE ===== ###

## Identify bacteria with ≥10% prevalence
## 725 milk samples, 10% of 725 -> 73 samples
# milk_sinvrbact_relab_v2 <- sinvrbact_relab_v2[grep("mother_human_milk", rownames(sinvrbact_relab_v2)),]
milk_prev10 <- milk_sinvrbact_relab_v2[,colSums(milk_sinvrbact_relab_v2>0)>=73]
length(colnames(milk_prev10)) # -> 36 genera are present in at least 10% of the samples
# colnames(milk_prev10)

## merge back with metadata
milkRAprev10 <- merge(sinvr[grep("mother_human_milk", rownames(sinvrbact_relab_v2)),1:530], milk_prev10, by="row.names")
rownames(milkRAprev10) <- milkRAprev10$Row.names
milkRAprev10 <- milkRAprev10[,-1]

## clr-transform relative abundances
milkRAprev10clr <- milkRAprev10 #column 'Other' was already excluded before
milkRAprev10clr[,531:ncol(milkRAprev10clr)] <- as.data.frame(do_clr_externalWeighting(milkRAprev10clr[,531:ncol(milkRAprev10clr)], milkRAprev10clr[,531:ncol(milkRAprev10clr)]))

## update colnames to show transformation
colnames(milkRAprev10clr)[531:ncol(milkRAprev10clr)] <- c(paste0(colnames(milkRAprev10clr)[531:ncol(milkRAprev10clr)], "_clr"))


### ===== 5.2. MATERNAL FAECES: CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE ===== ###

## Identify bacteria with ≥10% prevalence
## 170 maternal faecal samples, 10% of 170 -> 17 samples
# matf_sinvrbact_relab_v2 <- sinvrbact_relab_v2[grep("infant_faeces", rownames(sinvrbact_relab_v2)),]
matf_prev10 <- matf_sinvrbact_relab_v2[,colSums(matf_sinvrbact_relab_v2>0)>=17]
length(colnames(matf_prev10)) # -> 83 genera are present in at least 10% of the samples
# colnames(matf_prev10)

## merge back with metadata
matfRAprev10 <- merge(sinvr[grep("mother_faeces", rownames(sinvrbact_relab_v2)),1:530], matf_prev10, by="row.names")
rownames(matfRAprev10) <- matfRAprev10$Row.names
matfRAprev10 <- matfRAprev10[,-1]

## clr-transform relative abundances
matfRAprev10clr <- matfRAprev10 #column 'Other' was already excluded before
matfRAprev10clr[,531:ncol(matfRAprev10clr)] <- as.data.frame(do_clr_externalWeighting(matfRAprev10clr[,531:ncol(matfRAprev10clr)], matfRAprev10clr[,531:ncol(matfRAprev10clr)]))

## update colnames to show transformation
colnames(matfRAprev10clr)[531:ncol(matfRAprev10clr)] <- c(paste0(colnames(matfRAprev10clr)[531:ncol(matfRAprev10clr)], "_clr"))


### ===== 5.3. INFANT FAECES: CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE ===== ###

## Identify bacteria with ≥10% prevalence
## 552 infant faecal samples, 10% of 552 -> 56 samples
# inff_sinvrbact_relab_v2 <- sinvrbact_relab_v2[grep("infant_faeces", rownames(sinvrbact_relab_v2)),]
inff_prev10 <- inff_sinvrbact_relab_v2[,colSums(inff_sinvrbact_relab_v2>0)>=56]
length(colnames(inff_prev10)) # -> 30 genera are present in at least 10% of the samples
# colnames(inff_prev10)

## merge back with metadata
inffRAprev10 <- merge(sinvr[grep("infant_faeces", rownames(sinvrbact_relab_v2)),1:530], inff_prev10, by="row.names")
rownames(inffRAprev10) <- inffRAprev10$Row.names
inffRAprev10 <- inffRAprev10[,-1]

## clr-transform relative abundances
inffRAprev10clr <- inffRAprev10 #column 'Other' was already excluded before
inffRAprev10clr[,531:ncol(inffRAprev10clr)] <- as.data.frame(do_clr_externalWeighting(inffRAprev10clr[,531:ncol(inffRAprev10clr)], inffRAprev10clr[,531:ncol(inffRAprev10clr)]))

## update colnames to show transformation
colnames(inffRAprev10clr)[531:ncol(inffRAprev10clr)] <- c(paste0(colnames(inffRAprev10clr)[531:ncol(inffRAprev10clr)], "_clr"))


### ===== 5.4. SUPPLEMENTARY TABLES SHOWING CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE PER SAMPLE TYPE (TABLES S9, S11, S13) ===== ###

## FUNCTION TO RETRIEVE SUMMARY STATISTICS
get.prev10.sum.stats <- function(inputdata, columns, time_point){

  ## summary statistics
  print(paste0("Generating summary statistics for ", time_point))
  SumStats <- as.data.frame(apply(inputdata[,columns], 2, summary))
  
  ## add SD
  print(paste0("Adding SD"))
  SumSD <- c()
  for (i in columns){SumSD <- c(SumSD, sd(inputdata[,i]))}
  SumStats[7,] <- SumSD
  
  ## calculate number of individuals with/without the bacterium present
  print(paste0("Adding number of individuals with/without the bacterium present"))
  npres <- c()
  for (i in columns){npres <- c(npres, nrow(inputdata[inputdata[,i]>0,]))}
  SumStats[8,] <- npres
  
  nabs <- c()
  for (i in columns){nabs <- c(nabs, nrow(inputdata[inputdata[,i]==0,]))}
  SumStats[9,] <- nabs
  
  ## calculate percentage of individuals with/without the bacterium present
  print(paste0("Adding percentage of individuals with/without the bacterium present"))
  percpres <- c()
  for (i in columns){percpres <- c(percpres, nrow(inputdata[inputdata[,i]>0,])/nrow(inputdata))}
  SumStats[10,] <- percpres
  
  percabs <- c()
  for (i in columns){percabs <- c(percabs, nrow(inputdata[inputdata[,i]==0,])/nrow(inputdata))}
  SumStats[11,] <- percabs
  
  ## add n
  print(paste0("Adding n"))
  SumStats[12,] <- nrow(inputdata)
  
  ## combine summary statistics
  print(paste0("Combining summary statistics in 1 table"))
  SumStats$statistic <- c(rownames(SumStats)[1:6], "SD", "n_present", "n_absent", "perc_present", "perc_absent", "n")
  rownames(SumStats) <- 1:12
  
  ## clean summary statistics table
  print(paste0("Polishing summary statistics table"))
  SumStats <- SumStats[c(12,8:9,10:11,4,7,3,2,5,1,6), c(ncol(SumStats), 1:(ncol(SumStats)-1))]
  
  # transpose and add time point
  rownames(SumStats) <- SumStats$statistic
  SumStats <- SumStats[,-1]
  SumStats2 <- as.data.frame(t(SumStats))
  SumStats2$time_point <- as.factor(as.character(c(rep(time_point, nrow(SumStats2)))))
  SumStats2$bacterium <- as.factor(as.character(gsub("g__", "", rownames(SumStats2))))
  rownames(SumStats2) <- 1:nrow(SumStats2)
  SumStats2 <- SumStats2[,c(13:14,1:12)]

  print(paste0("Returning results"))
  return(SumStats2)
}


## RETRIEVE SUMMARY STATISTICS FOR HUMAN MILK

# for all time points together
milk_RAprev10_sum_all <- get.prev10.sum.stats(milkRAprev10, columns=c(531:566), time_point="all")

# per time point
milk_RAprev10_sum_M1 <- get.prev10.sum.stats(milkRAprev10[milkRAprev10$time_point=="1_month",], columns=c(531:566), time_point="1_month")
milk_RAprev10_sum_M2 <- get.prev10.sum.stats(milkRAprev10[milkRAprev10$time_point=="2_months",], columns=c(531:566), time_point="2_months")
milk_RAprev10_sum_M3 <- get.prev10.sum.stats(milkRAprev10[milkRAprev10$time_point=="3_months",], columns=c(531:566), time_point="3_months")
milk_RAprev10_sum_M6 <- get.prev10.sum.stats(milkRAprev10[milkRAprev10$time_point=="6_months",], columns=c(531:566), time_point="6_months")

# combine results into 1 table
milk_RAprev10_sum <- as.data.frame(rbind(milk_RAprev10_sum_all, milk_RAprev10_sum_M1, milk_RAprev10_sum_M2, milk_RAprev10_sum_M3, milk_RAprev10_sum_M6))

# save table (TABLE S9)
write.table(milk_RAprev10_sum, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/231215_milk_RA_prev10_summary_statistics.txt", sep="\t", row.names=F, quote=F)


## RETRIEVE SUMMARY STATISTICS FOR MATERNAL FAECES

# only M3
matf_RAprev10_sum_M3 <- get.prev10.sum.stats(matfRAprev10, columns=c(531:613), time_point="3_months")

# save table (TABLE S13)
write.table(matf_RAprev10_sum_M3, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/relative_abundances/231215_maternal_faeces_RA_prev10_summary_statistics.txt", sep="\t", row.names=F, quote=F)


## RETRIEVE SUMMARY STATISTICS FOR INFANT FAECES

# for all time points together
inff_RAprev10_sum_all <- get.prev10.sum.stats(inffRAprev10, columns=c(531:560), time_point="all")

# per time point
inff_RAprev10_sum_M1 <- get.prev10.sum.stats(inffRAprev10[inffRAprev10$time_point=="1_month",], columns=c(531:560), time_point="1_month")
inff_RAprev10_sum_M2 <- get.prev10.sum.stats(inffRAprev10[inffRAprev10$time_point=="2_months",], columns=c(531:560), time_point="2_months")
inff_RAprev10_sum_M3 <- get.prev10.sum.stats(inffRAprev10[inffRAprev10$time_point=="3_months",], columns=c(531:560), time_point="3_months")
inff_RAprev10_sum_M6 <- get.prev10.sum.stats(inffRAprev10[inffRAprev10$time_point=="6_months",], columns=c(531:560), time_point="6_months")

# combine results into 1 table
inff_RAprev10_sum <- as.data.frame(rbind(inff_RAprev10_sum_all, inff_RAprev10_sum_M1, inff_RAprev10_sum_M2, inff_RAprev10_sum_M3, inff_RAprev10_sum_M6))

# save table (TABLE S11)
write.table(inff_RAprev10_sum, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/231215_infant_faeces_RA_prev10_summary_statistics.txt", sep="\t", row.names=F, quote=F)


### ===== 5.5. HUMAN MILK: PLOT CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE (FIGURE S2A) ===== ###

### PREPARE HUMAN MILK DATA FOR PLOTTING

## Plot only the bacterial genera with ≥10% prevalence
## These genera are present in this data frame: milkRAprev10

## shorten colnames
colnames(milkRAprev10) <- gsub("g__", "", colnames(milkRAprev10))

## sort bacterial columns alphabetically
# milkRAprev10 <- milkRAprev10[,c(order(colnames(milkRAprev10)))]

## save bacterial genera separately
milkRAprev10bact <- milkRAprev10[,531:566]

## sum up the relative abundances of excluded bacteria into 'Other'
milkRAprev10bact$Other <- (1-rowSums(milkRAprev10bact))
table(rowSums(milkRAprev10bact)) #all have rowSums = 1


## calculate means per time point
milkmeansM1_v2 <- c()
for (i in 1:ncol(milkRAprev10bact[grep("mother_human_milk__1_month", rownames(milkRAprev10bact)),])){milkmeansM1_v2[i] <- mean(milkRAprev10bact[grep("mother_human_milk__1_month", rownames(milkRAprev10bact)), i])}
names(milkmeansM1_v2) <- colnames(milkRAprev10bact)

milkmeansM2_v2 <- c()
for (i in 1:ncol(milkRAprev10bact[grep("mother_human_milk__2_months", rownames(milkRAprev10bact)),])){milkmeansM2_v2[i] <- mean(milkRAprev10bact[grep("mother_human_milk__2_months", rownames(milkRAprev10bact)), i])}
names(milkmeansM2_v2) <- colnames(milkRAprev10bact)

milkmeansM3_v2 <- c()
for (i in 1:ncol(milkRAprev10bact[grep("mother_human_milk__3_months", rownames(milkRAprev10bact)),])){milkmeansM3_v2[i] <- mean(milkRAprev10bact[grep("mother_human_milk__3_months", rownames(milkRAprev10bact)), i])}
names(milkmeansM3_v2) <- colnames(milkRAprev10bact)

milkmeansM6_v2 <- c()
for (i in 1:ncol(milkRAprev10bact[grep("mother_human_milk__6_months", rownames(milkRAprev10bact)),])){milkmeansM6_v2[i] <- mean(milkRAprev10bact[grep("mother_human_milk__6_months", rownames(milkRAprev10bact)), i])}
names(milkmeansM6_v2) <- colnames(milkRAprev10bact)

## combine means in 1 data frame
milkRAplotdata <- as.data.frame(rbind(milkmeansM1_v2, milkmeansM2_v2, milkmeansM3_v2, milkmeansM6_v2))

## add metadata columns needed for plotting
milkRAplotdata$time_point <- as.factor(as.character(c("1_month", "2_months", "3_months", "6_months")))


## convert data frame from wide to long format
# install.packages("reshape")
library(reshape)
milkRAplotdatalong <- melt(milkRAplotdata, id.vars=colnames(milkRAplotdata)[c(38)], variable_name="Genera")


### HUMAN MILK RELATIVE ABUNDANCE PLOTS BY TIME POINT ===== ###

library(ggplot2)

### define colors for genera
manualcolors37<-c("mediumpurple4", "turquoise1", "tomato2", "yellow", "springgreen3",
                  "tan1", "lightcyan", "blue", "yellow3","deeppink3", 
                  "darkorchid", "lightpink", "palegreen1", "brown", "slategray2", 
                  "thistle2", "seagreen", "royalblue1", "wheat", "salmon1",
                  "deeppink", "seashell", "darkcyan", "khaki2", "chocolate3", 
                  "cadetblue3", "chocolate4", "plum", "burlywood4", "tan",
                  "lightpink3", "green", "darkblue", "firebrick3", "gold", 
                  "dodgerblue", "grey")

### create stacked barplots by sample type and time point
barplot1 <- ggplot(milkRAplotdatalong, aes(x=time_point, y=value, fill=Genera))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  labs(x="Month(s) postpartum", y="Average relative abundance")+
  scale_y_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_fill_manual(values=manualcolors37)+
  theme(panel.grid=element_blank(),
        strip.background = element_rect(color = "white", fill="white"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/231215_stacked_barplots_milk_RA_by_time_point.pdf", barplot1, dpi=500, width=30, height=20, units="cm")


### ===== 5.6. INFANT FAECES: PLOT CLASSIFIED BACTERIAL GENERA WITH ≥10% PREVALENCE (FIGURE S2B) ===== ###

### PREPARE INFANT FAECAL DATA FOR PLOTTING

## Plot only the bacterial genera with ≥10% prevalence
## These genera are present in this data frame: inffRAprev10

## shorten colnames
colnames(inffRAprev10) <- gsub("g__", "", colnames(inffRAprev10))

## sort bacterial columns alphabetically
# inffRAprev10 <- inffRAprev10[,c(order(colnames(inffRAprev10)))]

## save bacterial genera separately
inffRAprev10bact <- inffRAprev10[,531:560]

## sum up the relative abundances of excluded bacteria into 'Other'
inffRAprev10bact$Other <- (1-rowSums(inffRAprev10bact))
table(rowSums(inffRAprev10bact)) #all have rowSums = 1


## calculate means per time point
inffmeansM1_v2 <- c()
for (i in 1:ncol(inffRAprev10bact[grep("infant_faeces__1_month", rownames(inffRAprev10bact)),])){inffmeansM1_v2[i] <- mean(inffRAprev10bact[grep("infant_faeces__1_month", rownames(inffRAprev10bact)), i])}
names(inffmeansM1_v2) <- colnames(inffRAprev10bact)

inffmeansM2_v2 <- c()
for (i in 1:ncol(inffRAprev10bact[grep("infant_faeces__2_months", rownames(inffRAprev10bact)),])){inffmeansM2_v2[i] <- mean(inffRAprev10bact[grep("infant_faeces__2_months", rownames(inffRAprev10bact)), i])}
names(inffmeansM2_v2) <- colnames(inffRAprev10bact)

inffmeansM3_v2 <- c()
for (i in 1:ncol(inffRAprev10bact[grep("infant_faeces__3_months", rownames(inffRAprev10bact)),])){inffmeansM3_v2[i] <- mean(inffRAprev10bact[grep("infant_faeces__3_months", rownames(inffRAprev10bact)), i])}
names(inffmeansM3_v2) <- colnames(inffRAprev10bact)

inffmeansM6_v2 <- c()
for (i in 1:ncol(inffRAprev10bact[grep("infant_faeces__6_months", rownames(inffRAprev10bact)),])){inffmeansM6_v2[i] <- mean(inffRAprev10bact[grep("infant_faeces__6_months", rownames(inffRAprev10bact)), i])}
names(inffmeansM6_v2) <- colnames(inffRAprev10bact)

## combine means in 1 data frame
inffRAplotdata <- as.data.frame(rbind(inffmeansM1_v2, inffmeansM2_v2, inffmeansM3_v2, inffmeansM6_v2))

## add metadata columns needed for plotting
inffRAplotdata$time_point <- as.factor(as.character(c("1_month", "2_months", "3_months", "6_months")))


## convert data frame from wide to long format
# install.packages("reshape")
library(reshape)
inffRAplotdatalong <- melt(inffRAplotdata, id.vars=colnames(inffRAplotdata)[c(32)], variable_name="Genera")


### INFANT FAECAL RELATIVE ABUNDANCE PLOTS BY TIME POINT ===== ###

library(ggplot2)

### define colors for genera
manualcolors31<-c("yellowgreen" , "turquoise1", "darkmagenta", "yellow2", "blue",
                  "deeppink3", "sandybrown", "lightblue2", "green1", "coral4",
                  "orchid", "khaki", "brown", "skyblue4", "thistle2",
                  "darkorange2", "purple1", "seagreen", "salmon1", "seashell",
                  "lightsalmon4", "darkseagreen1", "hotpink4", "red", "darkcyan",
                  "navajowhite", "darkblue", "gold", "lightpink2", "dodgerblue",
                  "grey")

### create stacked barplots by sample type and time point
barplot2 <- ggplot(inffRAplotdatalong, aes(x=time_point, y=value, fill=Genera))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  labs(x="Month(s) postpartum", y="Average relative abundance")+
  scale_y_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_fill_manual(values=manualcolors31)+
  theme(panel.grid=element_blank(),
        strip.background = element_rect(color = "white", fill="white"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/231215_stacked_barplots_infant_faeces_RA_by_time_point.pdf", barplot2, dpi=500, width=30, height=20, units="cm")


##### =========================== 6. RELATIVE BACTERIAL ABUNDANCES BY SAMPLE TYPE AND TIME POINT =========================== #####

### ===== 6.1. PREPARE DATA FOR PLOTTING ===== ###

## Plot only the bacterial genera with ≥50% prevalence and ≥1% in ≥1 sample type
## These genera are shown in this vector: comb_prev50_mean1

## extract data for plotting from sinvrbact_relab_v2 data frame
RAplotbactdata <- sinvrbact_relab_v2[,colnames(sinvrbact_relab_v2) %in% comb_prev50_mean1]

## shorten colnames
colnames(RAplotbactdata) <- gsub("g__", "", colnames(RAplotbactdata))

## fix colname for UCG_002
# check for the full bacterial name of g__UCG_002
dataASV <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/231121_milk_composition_paper_sel_phenotypes_sel_microbiota_data_ASVs_n1515.rds")
UCG_002 <- colnames(dataASV[,grep("UCG_002", colnames(dataASV))])
UCG_002[1]
# k__Bacteria_p__Firmicutes_c__Clostridia_o__Oscillospirales_f__Oscillospiraceae_g__UCG_002

colnames(RAplotbactdata)[25] <- "Oscillospiraceae_UCG_002"

## sort bacterial columns alphabetically
RAplotbactdata <- RAplotbactdata[,c(order(colnames(RAplotbactdata)))]

## sum up the relative abundances of excldued bacteria into 'Other'
RAplotbactdata$Other <- (1-rowSums(RAplotbactdata))
table(rowSums(RAplotbactdata)) #all have rowSums = 1

## calculate means per sample type and time point
milkmeansM1 <- c()
for (i in 1:ncol(RAplotbactdata[grep("mother_human_milk__1_month", rownames(RAplotbactdata)),])){milkmeansM1[i] <- mean(RAplotbactdata[grep("mother_human_milk__1_month", rownames(RAplotbactdata)), i])}
names(milkmeansM1) <- colnames(RAplotbactdata)

milkmeansM2 <- c()
for (i in 1:ncol(RAplotbactdata[grep("mother_human_milk__2_months", rownames(RAplotbactdata)),])){milkmeansM2[i] <- mean(RAplotbactdata[grep("mother_human_milk__2_months", rownames(RAplotbactdata)), i])}
names(milkmeansM2) <- colnames(RAplotbactdata)

milkmeansM3 <- c()
for (i in 1:ncol(RAplotbactdata[grep("mother_human_milk__3_months", rownames(RAplotbactdata)),])){milkmeansM3[i] <- mean(RAplotbactdata[grep("mother_human_milk__3_months", rownames(RAplotbactdata)), i])}
names(milkmeansM3) <- colnames(RAplotbactdata)

milkmeansM6 <- c()
for (i in 1:ncol(RAplotbactdata[grep("mother_human_milk__6_months", rownames(RAplotbactdata)),])){milkmeansM6[i] <- mean(RAplotbactdata[grep("mother_human_milk__6_months", rownames(RAplotbactdata)), i])}
names(milkmeansM6) <- colnames(RAplotbactdata)


matfmeansM3 <- c()
for (i in 1:ncol(RAplotbactdata[grep("mother_faeces__3_months", rownames(RAplotbactdata)),])){matfmeansM3[i] <- mean(RAplotbactdata[grep("mother_faeces__3_months", rownames(RAplotbactdata)), i])}
names(matfmeansM3) <- colnames(RAplotbactdata)


inffmeansM1 <- c()
for (i in 1:ncol(RAplotbactdata[grep("infant_faeces__1_month", rownames(RAplotbactdata)),])){inffmeansM1[i] <- mean(RAplotbactdata[grep("infant_faeces__1_month", rownames(RAplotbactdata)), i])}
names(inffmeansM1) <- colnames(RAplotbactdata)

inffmeansM2 <- c()
for (i in 1:ncol(RAplotbactdata[grep("infant_faeces__2_months", rownames(RAplotbactdata)),])){inffmeansM2[i] <- mean(RAplotbactdata[grep("infant_faeces__2_months", rownames(RAplotbactdata)), i])}
names(inffmeansM2) <- colnames(RAplotbactdata)

inffmeansM3 <- c()
for (i in 1:ncol(RAplotbactdata[grep("infant_faeces__3_months", rownames(RAplotbactdata)),])){inffmeansM3[i] <- mean(RAplotbactdata[grep("infant_faeces__3_months", rownames(RAplotbactdata)), i])}
names(inffmeansM3) <- colnames(RAplotbactdata)

inffmeansM6 <- c()
for (i in 1:ncol(RAplotbactdata[grep("infant_faeces__6_months", rownames(RAplotbactdata)),])){inffmeansM6[i] <- mean(RAplotbactdata[grep("infant_faeces__6_months", rownames(RAplotbactdata)), i])}
names(inffmeansM6) <- colnames(RAplotbactdata)


## combine means in 1 data frame
RAplotdata <- as.data.frame(rbind(milkmeansM1, milkmeansM2, milkmeansM3, milkmeansM6,
                                  matfmeansM3, inffmeansM1, inffmeansM2, inffmeansM3, inffmeansM6))

## add metadata columns needed for plotting
RAplotdata$sample_origin_type <- as.factor(as.character(c(rep("mother_human_milk",4), "mother_faeces", rep("infant_faeces",4))))
RAplotdata$time_point <- as.factor(as.character(c("1_month", "2_months", "3_months", "6_months", "3_months", "1_month", "2_months", "3_months", "6_months")))

## convert data frame from wide to long format
# install.packages("reshape")
library(reshape)
RAplotdatalong <- melt(RAplotdata, id.vars=colnames(RAplotdata)[c(28:29)], variable_name="Genera")

## resort factors
RAplotdatalong$sample_origin_type <- factor(RAplotdatalong$sample_origin_type, levels = c("mother_human_milk", "mother_faeces", "infant_faeces"))


### ===== 6.2. RELATIVE ABUNDANCE PLOTS BY SAMPLE TYPE AND TIME POINT (FIGURE 4D) ===== ###

library(ggplot2)

### define colors for genera
manualcolors27<-c("mediumpurple4","olivedrab1","darkmagenta","cornflowerblue","blue",
                  "antiquewhite1","deeppink3","lightblue2","orange1","orchid",
                  "snow1","palegreen1","orangered3","seagreen","turquoise",
                  "thistle", "salmon1", "seashell","violetred4", "navajowhite",
                  "cyan", "slategray4", "darkblue", "gold", "limegreen",
                  "dodgerblue", "grey")

### create stacked barplots by sample type and time point
barplot1 <- ggplot(RAplotdatalong, aes(x=time_point, y=value, fill=Genera))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  labs(x="Month(s) postpartum", y="Average relative abundance")+
  scale_y_continuous(breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  scale_fill_manual(values=manualcolors27)+
  facet_grid(.~sample_origin_type, scales="free_x", space="free")+
  theme(panel.grid=element_blank(),
        strip.background = element_rect(color = "white", fill="white"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/relative_abundances/231208_stacked_barplots_RA_by_sample_type_and_time_point.pdf", barplot1, dpi=500, width=30, height=10, units="cm")


##### =========================== 7. BETA DIVERSITY BY SAMPLE TYPE (AND TIME POINT): ADONIS AND PCoA PLOT (FIGURE 4A) =========================== #####

## Goal: Check overall bacterial composition by sample type
## Use data frame dfRAclr, which contains data per sample type on the classified genus level, samples with ≥6500 reads and ASVs with ≥30 total reads and which does not contain duplicated maternal data for twins

library(vegan)
library(ggplot2)

## Add column combining sample type and time point
dfRAclr$sample_origin_type_time_point <- NA

dfRAclr[dfRAclr$sample_origin_type=="mother_human_milk" & dfRAclr$time_point=="1_month", "sample_origin_type_time_point"] <- "human_milk_1_month"
dfRAclr[dfRAclr$sample_origin_type=="mother_human_milk" & dfRAclr$time_point=="2_months", "sample_origin_type_time_point"] <- "human_milk_2_months"
dfRAclr[dfRAclr$sample_origin_type=="mother_human_milk" & dfRAclr$time_point=="3_months", "sample_origin_type_time_point"] <- "human_milk_3_months"
dfRAclr[dfRAclr$sample_origin_type=="mother_human_milk" & dfRAclr$time_point=="6_months", "sample_origin_type_time_point"] <- "human_milk_6_months"

dfRAclr[dfRAclr$sample_origin_type=="mother_faeces" & dfRAclr$time_point=="3_months", "sample_origin_type_time_point"] <- "maternal_faeces_3_months"

dfRAclr[dfRAclr$sample_origin_type=="infant_faeces" & dfRAclr$time_point=="1_month", "sample_origin_type_time_point"] <- "infant_faeces_1_month"
dfRAclr[dfRAclr$sample_origin_type=="infant_faeces" & dfRAclr$time_point=="2_months", "sample_origin_type_time_point"] <- "infant_faeces_2_months"
dfRAclr[dfRAclr$sample_origin_type=="infant_faeces" & dfRAclr$time_point=="3_months", "sample_origin_type_time_point"] <- "infant_faeces_3_months"
dfRAclr[dfRAclr$sample_origin_type=="infant_faeces" & dfRAclr$time_point=="6_months", "sample_origin_type_time_point"] <- "infant_faeces_6_months"

table(dfRAclr$sample_origin_type_time_point, useNA="ifany")


## Create Aitchison distance matrix and add metadata
ait <- vegdist(dfRAclr[,grep("g__", colnames(dfRAclr))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
aitm <- cmdscale(ait, k=3) #create distance matrix for first 3 components
aitd <- data.frame(aitm)

new_df <- merge(dfRAclr, aitd, by="row.names")
new_df <- new_df[,-1]
rownames(new_df) <- new_df$seq_16S_sample_ID

## Calculate variance explained
beta.cmd <- cmdscale(as.matrix(ait), k=2, eig=T)
PCoA1 <- beta.cmd$eig[1]/sum(beta.cmd$eig)*100
PCoA2 <- beta.cmd$eig[2]/sum(beta.cmd$eig)*100


## Run adonis by sample type
# a1 <- adonis2(ait ~ new_df$sample_origin_type, permutations=10000)
# a1
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 10000
# 
# adonis2(formula = ait ~ new_df$sample_origin_type, permutations = 10000)
#                             Df SumOfSqs      R2      F    Pr(>F)    
# new_df$sample_origin_type    2   218541 0.30606 318.43 9.999e-05 ***
# Residual                  1444   495507 0.69394                     
# Total                     1446   714049 1.00000                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## resort levels for color and shape
new_df$sample_origin_type_time_point <- factor(new_df$sample_origin_type_time_point, levels = c("human_milk_1_month", "human_milk_2_months", "human_milk_3_months", "human_milk_6_months",
                                                                                                "maternal_faeces_3_months",
                                                                                                "infant_faeces_1_month", "infant_faeces_2_months", "infant_faeces_3_months" ,"infant_faeces_6_months"))
new_df$sample_origin_type <- factor(new_df$sample_origin_type, levels = c("mother_human_milk", "mother_faeces", "infant_faeces"))

pcoaplot1 <- ggplot(new_df, aes(x=X1, y=X2, color=sample_origin_type_time_point))+
  geom_point(alpha=0.2, size=3)+
  stat_ellipse(aes(group=sample_origin_type_time_point))+
  labs(x=paste0("PC1 (", signif(PCoA1,5), "%)"), y=paste0("PC2 (", signif(PCoA2,5), "%)"),
       title=paste0("Beta diversity by sample type and time point\n",
                    "Human milk, n=", nrow(new_df[new_df$sample_origin_type=="mother_human_milk",]), "; ",
                    "Maternal faeces, n=", nrow(new_df[new_df$sample_origin_type=="mother_faeces",]), "; ",
                    "Infant faeces, n=", nrow(new_df[new_df$sample_origin_type=="infant_faeces",])))+
  scale_color_manual(values = c("human_milk_1_month" = "#F0E442",
                                "human_milk_2_months" = "gold2",
                                "human_milk_3_months" = "#E69F00",
                                "human_milk_6_months" = "sienna3",
                                "maternal_faeces_3_months" = "#0072B2",
                                "infant_faeces_1_month" = "lightgreen",
                                "infant_faeces_2_months" = "limegreen",
                                "infant_faeces_3_months" = "#009E73",
                                "infant_faeces_6_months" = "darkolivegreen"))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  guides(size="none")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/all_sample_types/231208_PCoA_PC1_PC2_beta_diversity_by_sample_type_and_time_point.pdf"), pcoaplot1, dpi=500, width=20, height=15, units="cm")


##### =========================== 8. PLOT ALPHA DIVERSITY BY SAMPLE TYPES (POOLED OR PER TIME POINT) =========================== #####

## For plotting, use data frame df2 with non-transformed alpha diversity measures!


### ===== 8.1. PREPARE DATA FOR PLOTTING ===== ###

## change the order of the sample types
df2$sample_origin_type <- factor(df2$sample_origin_type,
                                 levels = c("mother_faeces", "mother_human_milk", "infant_faeces"))


# ### ===== 8.2. PLOT ALPHA DIVERSITY BY SAMPLE TYPE (POOLED TIME POINTS) ===== ###
# 
# ### GENERA RICHNESS
# plotalpha1 <- ggplot(df2, aes(x=sample_origin_type, y=alpha_div_genera_richness, color=sample_origin_type))+
#   geom_jitter(alpha=0.5)+
#   geom_boxplot(alpha=0, color="black")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90),
#         legend.position = "none",
#         panel.grid = element_blank())+
#   ggtitle(paste0("Classified genera richness; n=", nrow(df2)))+
#   scale_y_continuous(limits = c(0,80), breaks = c(seq(0,80,20)))+
#   scale_color_manual(values = c("infant_faeces" = "#009E73",
#                                 "mother_faeces" = "#0072B2",
#                                 "mother_human_milk" = "#E69F00"))
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_genera_richness_by_sample_type.pdf", plotalpha1, dpi=500, width=10, height=15, units="cm")
# 
# 
# ## SHANNON DIVERSITY INDEX
# plotalpha2 <- ggplot(df2, aes(x=sample_origin_type, y=alpha_div_shannon, color=sample_origin_type))+
#   geom_jitter(alpha=0.5)+
#   geom_boxplot(alpha=0, color="black")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90),
#         legend.position = "none",
#         panel.grid = element_blank())+
#   ggtitle(paste0("Shannon diversity index; n=", nrow(df2)))+
#   scale_y_continuous(limits = c(0,4), breaks = c(seq(0,4,1)))+
#   scale_color_manual(values = c("infant_faeces" = "#009E73",
#                                 "mother_faeces" = "#0072B2",
#                                 "mother_human_milk" = "#E69F00"))
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_shannon_by_sample_type.pdf", plotalpha2, dpi=500, width=10, height=15, units="cm")
# 
# 
# ## SIMPSON DIVERSITY INDEX
# plotalpha3 <- ggplot(df2, aes(x=sample_origin_type, y=alpha_div_simpson, color=sample_origin_type))+
#   geom_jitter(alpha=0.5)+
#   geom_boxplot(alpha=0, color="black")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90),
#         legend.position = "none",
#         panel.grid = element_blank())+
#   ggtitle(paste0("Simpson diversity index; n=", nrow(df2)))+
#   scale_y_continuous(limits = c(0,1.2), breaks = c(seq(0,1.2,0.2)))+
#   scale_color_manual(values = c("infant_faeces" = "#009E73",
#                                 "mother_faeces" = "#0072B2",
#                                 "mother_human_milk" = "#E69F00"))
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_simpson_by_sample_type.pdf", plotalpha3, dpi=500, width=10, height=15, units="cm")


### ===== 8.3. PLOT ALPHA DIVERSITY BY SAMPLE TYPE AND TIME POINT (FIGURE 4C) ===== ###

### GENERA RICHNESS
plotalpha4 <- ggplot(df2, aes(x=time_point, y=alpha_div_genera_richness, color=sample_origin_type))+
  geom_jitter(alpha=0.5)+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  ggtitle(paste0("Classified genera richness; n=", nrow(df2)))+
  labs(x="Month(s) postpartum")+
  scale_y_continuous(limits = c(0,80), breaks = c(seq(0,80,20)))+
  facet_grid(.~sample_origin_type, scales="free_x", space="free")+
  scale_color_manual(values = c("infant_faeces" = "#009E73",
                                "mother_faeces" = "#0072B2",
                                "mother_human_milk" = "#E69F00"))+
  theme(axis.text.x = element_text(angle=90),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(color = "white", fill="white"))
        # axis.text.y = element_blank())
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_genera_richness_by_sample_type_and_time_point.pdf", plotalpha4, dpi=500, width=15, height=15, units="cm")
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_genera_richness_by_sample_type_and_time_point_without_yaxis_labels.pdf", plotalpha4, dpi=500, width=15, height=15, units="cm")


## SHANNON DIVERSITY INDEX
plotalpha5 <- ggplot(df2, aes(x=time_point, y=alpha_div_shannon, color=sample_origin_type))+
  geom_jitter(alpha=0.5)+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  ggtitle(paste0("Shannon diversity index; n=", nrow(df2)))+
  labs(x="Month(s) postpartum")+
  scale_y_continuous(limits = c(0,4), breaks = c(seq(0,4,1)))+
  facet_grid(.~sample_origin_type, scales="free_x", space="free")+
  scale_color_manual(values = c("infant_faeces" = "#009E73",
                                "mother_faeces" = "#0072B2",
                                "mother_human_milk" = "#E69F00"))+
  theme(axis.text.x = element_text(angle=90),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(color = "white", fill="white"))
        # axis.text.y = element_blank())
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_shannon_by_sample_type_and_time_point.pdf", plotalpha5, dpi=500, width=15, height=15, units="cm")
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_shannon_by_sample_type_and_time_point_without_yaxis_labels.pdf", plotalpha5, dpi=500, width=15, height=15, units="cm")


## SIMPSON DIVERSITY INDEX
plotalpha6 <- ggplot(df2, aes(x=time_point, y=alpha_div_simpson, color=sample_origin_type))+
  geom_jitter(alpha=0.5)+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  ggtitle(paste0("Simpson diversity index; n=", nrow(df2)))+
  labs(x="Month(s) postpartum")+
  facet_grid(.~sample_origin_type, scales="free_x", space="free")+
  scale_color_manual(values = c("infant_faeces" = "#009E73",
                                "mother_faeces" = "#0072B2",
                                "mother_human_milk" = "#E69F00"))+
  theme(axis.text.x = element_text(angle=90),
        legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(color = "white", fill="white"))
        # axis.text.y = element_blank())
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_simpson_by_sample_type_and_time_point.pdf", plotalpha6, dpi=500, width=15, height=15, units="cm")
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_simpson_by_sample_type_and_time_point_without_yaxis_labels.pdf", plotalpha6, dpi=500, width=15, height=15, units="cm")


# ### ===== 8.4. PLOT ALPHA DIVERSITY FOR MILK AND INFANT FAECES AND TIME POINT ===== ###
# 
# ### GENERA RICHNESS
# plotalpha7 <- ggplot(df2[(df2$sample_origin_type=="mother_human_milk" | df2$sample_origin_type=="infant_faeces"),], aes(x=sample_origin_type, y=alpha_div_genera_richness, color=sample_origin_type))+
#   geom_jitter(alpha=0.5)+
#   geom_boxplot(alpha=0, color="black")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90),
#         legend.position = "none",
#         panel.grid = element_blank())+
#   ggtitle(paste0("Classified genera richness; n=", nrow(df2)))+
#   labs(x="Month(s) postpartum")+
#   scale_y_continuous(limits = c(0,60), breaks = c(seq(0,60,20)))+
#   facet_grid(.~time_point, scales="free_x")+
#   scale_color_manual(values = c("infant_faeces" = "#009E73",
#                                 "mother_human_milk" = "#E69F00"))
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/milk_and_infant_faeces/231208_genera_richness_milk_and_infant_faeces_per_time_point.pdf", plotalpha7, dpi=500, width=15, height=15, units="cm")
# 
# 
# ## SHANNON DIVERSITY INDEX
# plotalpha8 <- ggplot(df2[(df2$sample_origin_type=="mother_human_milk" | df2$sample_origin_type=="infant_faeces"),], aes(x=sample_origin_type, y=alpha_div_shannon, color=sample_origin_type))+
#   geom_jitter(alpha=0.5)+
#   geom_boxplot(alpha=0, color="black")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90),
#         legend.position = "none",
#         panel.grid = element_blank())+
#   ggtitle(paste0("Shannon diversity index; n=", nrow(df2)))+
#   labs(x="Month(s) postpartum")+
#   scale_y_continuous(limits = c(0,4), breaks = c(seq(0,4,1)))+
#   facet_grid(.~time_point, scales="free_x")+
#   scale_color_manual(values = c("infant_faeces" = "#009E73",
#                                 "mother_human_milk" = "#E69F00"))
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/milk_and_infant_faeces/231208_shannon_milk_and_infant_faeces_per_time_point.pdf", plotalpha8, dpi=500, width=15, height=15, units="cm")
# 
# 
# ## SIMPSON DIVERSITY INDEX
# plotalpha9 <- ggplot(df2[(df2$sample_origin_type=="mother_human_milk" | df2$sample_origin_type=="infant_faeces"),], aes(x=sample_origin_type, y=alpha_div_simpson, color=sample_origin_type))+
#   geom_jitter(alpha=0.5)+
#   geom_boxplot(alpha=0, color="black")+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle=90),
#         legend.position = "none",
#         panel.grid = element_blank())+
#   ggtitle(paste0("Simpson diversity index; n=", nrow(df2)))+
#   labs(x="Month(s) postpartum")+
#   scale_y_continuous(limits = c(0,1.2), breaks = c(seq(0,1.2,0.2)))+
#   facet_grid(.~time_point, scales="free_x")+
#   scale_color_manual(values = c("infant_faeces" = "#009E73",
#                                 "mother_human_milk" = "#E69F00"))
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/milk_and_infant_faeces/231208_simpson_milk_and_infant_faeces_per_time_point.pdf", plotalpha9, dpi=500, width=15, height=15, units="cm")


##### =========================== 9. COMPARE ALPHA DIVERSITY BETWEEN SAMPLE TYPES (POOLED OR PER TIME POINT) =========================== #####

## For statistical analyses, use data frame sinvr with inverse-rank transformed alpha diversity measures!


### ===== 9.1. FUNCTION TO COMPARE ALPHA DIVERSITY BY SAMPLE TYPE (POOLED TIME POINTS) ===== ###

## function for linear model with correction for clean read count
run.lm.cor.reads.time <- function(datasetname, inputdata, xcolumn, ycolumn){
  p_models <- c()
  
  levels <- c()
  n_levels <- c()
  e <- c()
  
  d <- c()
  x <- c()
  y <- c()
  n <- c()
  stat <- c()
  
  for (i in xcolumn){
    
    for (j in ycolumn){
      print(paste0("i: ",i,"j: ",j))
      
      print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #first mixed model with correction for clean read count and without the phenotype of interest (i)
      print(paste0("Run first model"))
      m0 <- lm(my_df[,j] ~ seq_16S_n_reads_clean, data=my_df)
      sum0 <- summary(m0)
      # print(sum0)
      
      #second mixed model with correction for clean read count and WITH the phenotype of interest (i)
      print(paste0("Run second model"))
      m1 <- lm(my_df[,j] ~ seq_16S_n_reads_clean + my_df[,i], data=my_df)
      sum1 <- summary(m1)
      # print(sum1)
      
      #compare models and save p value from model comparison
      print(paste0("Compare models"))
      an1 <- anova(m0, m1)
      # print(an1)
      
      p_models <- c(p_models, an1[2, "Pr(>F)"]) #the 2 picks the results from the tested x as there always is only one row for the effect from the clean reads
      # print(p_models)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      print(paste0("Extract parameters from model results"))
      levels <- c(levels, paste(collapse=";",levels(my_df[,i])))
      n_levels <- c(n_levels, paste(collapse=";",table(my_df[,i])))
      e <- c(e, paste(collapse=";",c(0,signif(digits=3,sum1$coefficients[grep("my_df",rownames(sum1$coefficients)), "Estimate"])))) #first row is from the clean reads -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, paste0("lm_correction_reads"))
      d <- c(d, datasetname)
      x <- c(x, colnames(inputdata)[i])
      y <- c(y, colnames(inputdata)[j])
      n <- c(n, nrow(my_df))
    }
  }
  
  #save results in data frame
  print(paste0("Combine results in data frame and save results"))
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, n_levels=n_levels, estimate=e, p=p_models)
  
  #adjust p values for multiple testing
  res$FDR_BH <- p.adjust(res$p, method="BH")
  
  #sort results by FDR_BH and p
  res <- res[order(res$FDR_BH, res$p),]
  
  return(res)
}


### ===== 9.2. COMPARE ALPHA DIVERSITY BETWEEN SAMPLE TYPES (PAIR-WISE COMPARISONS, POOLED TIME POINTS) ===== ###

## Compare milk vs infant faeces
alpha_div_milk_vs_inff_lm_results <- run.lm.cor.reads.time(datasetname = "milk_vs_infant_faeces",
                                                          inputdata = sinvr[(sinvr$sample_origin_type=="mother_human_milk" | sinvr$sample_origin_type=="infant_faeces"),],
                                                          xcolumn = c(507), #507 sample_origin_type
                                                          ycolumn = c(528:530)) #inverse-rank transformed alpha diversity measures

## Compare milk vs maternal faeces
alpha_div_milk_vs_matf_lm_results <- run.lm.cor.reads.time(datasetname = "milk_vs_maternal_faeces",
                                                           inputdata = sinvr[(sinvr$sample_origin_type=="mother_human_milk" | sinvr$sample_origin_type=="mother_faeces"),],
                                                           xcolumn = c(507), #507 sample_origin_type
                                                           ycolumn = c(528:530)) #inverse-rank transformed alpha diversity measures

## Compare maternal vs infant faeces
alpha_div_matf_vs_inff_lm_results <- run.lm.cor.reads.time(datasetname = "maternal_vs_infant_faeces",
                                                           inputdata = sinvr[(sinvr$sample_origin_type=="mother_faeces" | sinvr$sample_origin_type=="infant_faeces"),],
                                                           xcolumn = c(507), #507 sample_origin_type
                                                           ycolumn = c(528:530)) #inverse-rank transformed alpha diversity measures

## Combine results, correct for multiple testing, resort data frame and save results
alpha_div_between_sample_types_results <- as.data.frame(rbind(alpha_div_milk_vs_inff_lm_results[,-10], alpha_div_milk_vs_matf_lm_results[,-10], alpha_div_matf_vs_inff_lm_results[,-10])) #exclude FDR_BH columns
alpha_div_between_sample_types_results$FDR_BH <- p.adjust(alpha_div_between_sample_types_results$p, method = "BH")
alpha_div_between_sample_types_results <- alpha_div_between_sample_types_results[order(alpha_div_between_sample_types_results$FDR_BH),]
write.table(alpha_div_between_sample_types_results, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_lm_results_alpha_diversity_between_sample_types_genus_level.txt", sep="\t", row.names=F, quote=F)

## Save results without Simpson diversity, re-calculate FDR
alpha_div_between_sample_types_results_no_simpson <- alpha_div_between_sample_types_results[alpha_div_between_sample_types_results$y!="alpha_div_simpson_invr",]
alpha_div_between_sample_types_results_no_simpson$FDR_BH <- p.adjust(alpha_div_between_sample_types_results_no_simpson$p, method = "BH")
alpha_div_between_sample_types_results_no_simpson <- alpha_div_between_sample_types_results_no_simpson[order(alpha_div_between_sample_types_results_no_simpson$FDR_BH),]
write.table(alpha_div_between_sample_types_results_no_simpson, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/240216_lm_results_alpha_diversity_richness_shannon_between_sample_types_genus_level.txt", sep="\t", row.names=F, quote=F)


### ===== 9.3. COMPARE ALPHA DIVERSITY BETWEEN MILK AND INFANT FAECES PER TIME POINT ===== ###

## Compare milk vs infant faeces, 1 month
alpha_div_milk_vs_inff_lm_results_M1 <- run.lm.cor.reads.time(datasetname = "milk_vs_infant_faeces_M1",
                                                           inputdata = sinvr[(sinvr$sample_origin_type=="mother_human_milk" | sinvr$sample_origin_type=="infant_faeces") & sinvr$time_point=="1_month",],
                                                           xcolumn = c(507), ycolumn = c(528:530))

## Compare milk vs infant faeces, 2 month
alpha_div_milk_vs_inff_lm_results_M2 <- run.lm.cor.reads.time(datasetname = "milk_vs_infant_faeces_M2",
                                                              inputdata = sinvr[(sinvr$sample_origin_type=="mother_human_milk" | sinvr$sample_origin_type=="infant_faeces") & sinvr$time_point=="2_months",],
                                                              xcolumn = c(507), ycolumn = c(528:530))

## Compare milk vs infant faeces, 3 month
alpha_div_milk_vs_inff_lm_results_M3 <- run.lm.cor.reads.time(datasetname = "milk_vs_infant_faeces_M3",
                                                              inputdata = sinvr[(sinvr$sample_origin_type=="mother_human_milk" | sinvr$sample_origin_type=="infant_faeces") & sinvr$time_point=="3_months",],
                                                              xcolumn = c(507), ycolumn = c(528:530))

## Compare milk vs infant faeces, 6 month
alpha_div_milk_vs_inff_lm_results_M6 <- run.lm.cor.reads.time(datasetname = "milk_vs_infant_faeces_M6",
                                                              inputdata = sinvr[(sinvr$sample_origin_type=="mother_human_milk" | sinvr$sample_origin_type=="infant_faeces") & sinvr$time_point=="6_months",],
                                                              xcolumn = c(507), ycolumn = c(528:530))


## Combine results, correct for multiple testing, resort data frame and save results
alpha_div_milk_vs_inff_lm_results_all <- as.data.frame(rbind(alpha_div_milk_vs_inff_lm_results_M1[,-10], alpha_div_milk_vs_inff_lm_results_M2[,-10], alpha_div_milk_vs_inff_lm_results_M3[,-10], alpha_div_milk_vs_inff_lm_results_M6[,-10])) #exclude FDR_BH columns
alpha_div_milk_vs_inff_lm_results_all$FDR_BH <- p.adjust(alpha_div_milk_vs_inff_lm_results_all$p, method = "BH")
alpha_div_milk_vs_inff_lm_results_all <- alpha_div_milk_vs_inff_lm_results_all[order(alpha_div_milk_vs_inff_lm_results_all$FDR_BH),]
write.table(alpha_div_milk_vs_inff_lm_results_all, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/milk_and_infant_faeces/231208_lm_results_alpha_diversity_between_milk_and_infant_faeces_per_time_point_genus_level.txt", sep="\t", row.names=F, quote=F)

## Save results without Simpson diversity, re-calculate FDR
alpha_div_milk_vs_inff_lm_results_all_no_simpson <- alpha_div_milk_vs_inff_lm_results_all[alpha_div_milk_vs_inff_lm_results_all$y!="alpha_div_simpson_invr",]
alpha_div_milk_vs_inff_lm_results_all_no_simpson$FDR_BH <- p.adjust(alpha_div_milk_vs_inff_lm_results_all_no_simpson$p, method = "BH")
alpha_div_milk_vs_inff_lm_results_all_no_simpson <- alpha_div_milk_vs_inff_lm_results_all_no_simpson[order(alpha_div_milk_vs_inff_lm_results_all_no_simpson$FDR_BH),]
write.table(alpha_div_milk_vs_inff_lm_results_all_no_simpson, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/milk_and_infant_faeces/240216_lm_results_alpha_diversity_richness_shannon_between_milk_and_infant_faeces_per_time_point_genus_level.txt", sep="\t", row.names=F, quote=F)


### ===== 9.4. COMPARE ALPHA DIVERSITY BETWEEN SAMPLE TYPES (PAIR-WISE COMPARISONS, MATERNAL FAECES M3 VS. HUMAN MILK / INFANT FAECES M6) ===== ###

## Compare milk vs maternal faeces
alpha_div_milkM6_vs_matfM3_lm_results <- run.lm.cor.reads.time(datasetname = "milkM6_vs_maternal_faecesM3",
                                                               inputdata = sinvr[(sinvr$sample_origin_type=="mother_human_milk" & sinvr$time_point=="6_months") | 
                                                                               (sinvr$sample_origin_type=="mother_faeces" & sinvr$time_point=="3_months"),],
                                                               xcolumn = c(507), ycolumn = c(528:530)) #507 sample_origin_type

## Compare maternal vs infant faeces
alpha_div_matfM3_vs_inffM6_lm_results <- run.lm.cor.reads.time(datasetname = "maternalM3_vs_infantM6_faeces",
                                                               inputdata = sinvr[(sinvr$sample_origin_type=="infant_faeces" & sinvr$time_point=="6_months") | 
                                                                               (sinvr$sample_origin_type=="mother_faeces" & sinvr$time_point=="3_months"),],
                                                               xcolumn = c(507), ycolumn = c(528:530)) #507 sample_origin_type

## Combine results, correct for multiple testing, resort data frame and save results
alpha_div_M6_vs_matfM3_results <- as.data.frame(rbind(alpha_div_milkM6_vs_matfM3_lm_results[,-10], alpha_div_matfM3_vs_inffM6_lm_results[,-10])) #exclude FDR_BH columns
alpha_div_M6_vs_matfM3_results$FDR_BH <- p.adjust(alpha_div_M6_vs_matfM3_results$p, method = "BH")
alpha_div_M6_vs_matfM3_results <- alpha_div_M6_vs_matfM3_results[order(alpha_div_M6_vs_matfM3_results$FDR_BH),]
write.table(alpha_div_M6_vs_matfM3_results, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/231208_lm_results_alpha_diversity_maternal_faeces_M3_vs_milk_and_infant_faeces_M6.txt", sep="\t", row.names=F, quote=F)

## Save results without Simpson diversity, re-calculate FDR
alpha_div_M6_vs_matfM3_results_no_simpson <- alpha_div_M6_vs_matfM3_results[alpha_div_M6_vs_matfM3_results$y!="alpha_div_simpson_invr",]
alpha_div_M6_vs_matfM3_results_no_simpson$FDR_BH <- p.adjust(alpha_div_M6_vs_matfM3_results_no_simpson$p, method = "BH")
alpha_div_M6_vs_matfM3_results_no_simpson <- alpha_div_M6_vs_matfM3_results_no_simpson[order(alpha_div_M6_vs_matfM3_results_no_simpson$FDR_BH),]
write.table(alpha_div_M6_vs_matfM3_results_no_simpson, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/alpha_diversity/all_sample_types/240216_lm_results_alpha_diversity_richness_shannon_maternal_faeces_M3_vs_milk_and_infant_faeces_M6.txt", sep="\t", row.names=F, quote=F)


##### =========================== 10. CHECK THE EFFECT OF TIME ON MILK ALPHA DIVERSITY AND MILK RELATIVE ABUNDANCES =========================== #####

# install.packages("lmerTest")
library(lmerTest)


### ===== 10.1. MODEL FUNCTION WITH CORRECTION FOR DNA ISOLATION BATCH/METHOD (FIXED EFFECT), CLEAN READ COUNT (FIXED EFFECT) AND INDIVIDUAL ID (RANDOM EFFECT) ===== ###

## function for mixed model with corrections 
run.lmer.cor.ybatch.yreads.yID <- function(datasetname, inputdata, xcolumn, ycolumn){
  p <- c()
  
  levels <- c()
  n_levels <- c()
  intercept <- c()
  e <- c()
  
  d <- c()
  x <- c()
  y <- c()
  n <- c()
  stat <- c()
  
  for (i in xcolumn){
    
    for (j in ycolumn){
      print(paste0("i: ",i,"j: ",j))
      
      print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #mixed model with confounders and WITH the phenotype of interest (i)
      m1 <- lmerTest::lmer(my_df[,j] ~ DNA_isolation_batch + seq_16S_n_reads_clean + (1|mother_ID) + my_df[,i], REML=F, data=my_df)
      sum1 <- summary(m1)
      print(paste0("Printing sum1"))
      print(sum1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      p <- c(p, sum1$coefficients[rownames(sum1$coefficients)=="my_df[, i]", "Pr(>|t|)"])
      levels <- c(levels, paste(collapse=";",levels(my_df[,i])))
      n_levels <- c(n_levels, paste(collapse=";",table(my_df[,i])))
      intercept <- c(intercept, c(sum1$coefficients[rownames(sum1$coefficients)=="(Intercept)", "Estimate"]))
      e <- c(e, c(sum1$coefficients[rownames(sum1$coefficients)=="my_df[, i]", "Estimate"])) #first rows are for confounders -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, paste0("lmer_correction_yDNAbatch_yreads_yID"))
      d <- c(d, datasetname)
      x <- c(x, colnames(inputdata)[i])
      y <- c(y, colnames(inputdata)[j])
      n <- c(n, nrow(my_df))
    }
  }
  
  #save results in data frame
  print(paste0("Combine results in data frame and save results"))
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, n_levels=n_levels, intercept=intercept, estimate=e, p=p)
  
  # Multiple testing correction (BH)
  res$FDR_BH <- p.adjust(res$p, method="BH")
  
  # Order data frame by p value
  res <- res[order(res$FDR_BH, res$p),]
  
  return(res)
}


### ===== 10.2. COMPARE MILK ALPHA DIVERSITY BY NUMERIC TIME POINTS ===== ###

## For statistical analyses, use data frame sinvr with inverse-rank transformed alpha diversity measures!
  
## Compare milk alpha diversity between time points (numeric)
milk_lmer_results_alpha_div_by_time <- run.lmer.cor.ybatch.yreads.yID(datasetname = "milk",
                                                                    inputdata = sinvr[sinvr$sample_origin_type=="mother_human_milk",],
                                                                    xcolumn = c(10), #numeric time point
                                                                    ycolumn = c(528:530)) #inverse-rank transformed alpha diversity measures
write.table(milk_lmer_results_alpha_div_by_time, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/alpha_diversity/time_point/231212_lmer_results_milk_alpha_diversity_by_numeric_time_point.txt", sep="\t", row.names=F, quote=F)


## Save results without Simpson diversity, re-calculate FDR
milk_lmer_results_alpha_div_by_time_no_simpson <- milk_lmer_results_alpha_div_by_time[milk_lmer_results_alpha_div_by_time$y!="alpha_div_simpson_invr",]
milk_lmer_results_alpha_div_by_time_no_simpson$FDR <- p.adjust(milk_lmer_results_alpha_div_by_time_no_simpson$p, method="BH")
milk_lmer_results_alpha_div_by_time_no_simpson <- milk_lmer_results_alpha_div_by_time_no_simpson[order(milk_lmer_results_alpha_div_by_time_no_simpson$FDR, milk_lmer_results_alpha_div_by_time_no_simpson$p),]
write.table(milk_lmer_results_alpha_div_by_time_no_simpson, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/alpha_diversity/time_point/240216_lmer_results_milk_alpha_diversity_richness_shannon_by_numeric_time_point.txt", sep="\t", row.names=F, quote=F)


### ===== 10.3. COMPARE MILK RELATIVE ABUNDANCES BY NUMERIC TIME POINTS (TABLE S10) ===== ###

## For statistical analyses, use data frame xx with clr-transformed relative abundances of the bacterial genera with ≥10% prevalence in milk samples.

## Compare milk bacterial relative abundances between time points (numeric)
milk_lmer_results_RA_by_time <- run.lmer.cor.ybatch.yreads.yID(datasetname = "milk",
                                                               inputdata = milkRAprev10clr,
                                                               xcolumn = c(10), #numeric time point
                                                               ycolumn = c(531:566)) #clr-transformed relative abundances of bacteria with ≥10% prevalence
write.table(milk_lmer_results_RA_by_time, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/time_point/231212_lmer_results_milk_RA_by_numeric_time_point.txt", sep="\t", row.names=F, quote=F)


##### =========================== 11. PLOT EFFECT OF TIME ON MILK RELATIVE ABUNDANCES =========================== #####

### ===== 11.1. FUNCTION TO PLOT SIGNIFICANT EFFECTS OF TIME ON RELATIVE ABUNDANCES ===== ###

plot.milk.RA.by.time.sign.results <- function(inputdata, signresults){
  for (i in signresults[,"y"]){
    print(paste0("Plotting effect of time on ", i))
    # print(paste0(intercept=signresults[signresults$y==paste0(i),"intercept"]))
    # print(paste0(slope=signresults[signresults$y==paste0(i),"estimate"]))
    
    RAbytimeplot <- ggplot(inputdata, aes(x=time_point_numeric, y=inputdata[,i], color=time_point))+
      geom_point(alpha=0.5)+
      geom_boxplot(alpha=0, color="black", aes(group=time_point))+
      theme_bw()+
      theme(legend.position = "none",
            panel.grid = element_blank())+
      ggtitle(paste0(i, " RA by time; n=", signresults[signresults$y==i, "n_total"]))+
      labs(x="Month(s) postpartum", y=paste0(i))+
      scale_color_manual(values = c("1_month" = "#F0E442",
                                    "2_months" = "gold2",
                                    "3_months" = "#E69F00",
                                    "6_months" = "sienna3"))+
      scale_x_continuous(limits = c(0.5,6.5), breaks = c(1,2,3,4,5,6))+
      geom_abline(intercept=signresults[signresults$y==paste0(i),"intercept"],
                  slope=signresults[signresults$y==paste0(i),"estimate"])
    ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/time_point/231212_milk_", i, "_RA_by_time.pdf"), RAbytimeplot, dpi=500, width=15, height=15, units="cm")
  }
}

### ===== 11.2. PLOT SIGNIFICANT EFFECTS OF TIME ON MILK RELATIVE ABUNDANCES ===== ###

plot.milk.RA.by.time.sign.results(inputdata = milkRAprev10clr,
                                  signresults = milk_lmer_results_RA_by_time[milk_lmer_results_RA_by_time$FDR_BH<0.05,])


### ===== 11.3. PLOT ESTIMATES FOR EACH BACTERIUM (USED FOR FIGURE 4E) ===== ###

# milk_lmer_results_RA_by_time <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/time_point/231212_lmer_results_milk_RA_by_numeric_time_point.txt", header=T, sep="\t", stringsAsFactors=T)

## save data separately to prepare for plotting
milk_RA_by_time_plot_data <- milk_lmer_results_RA_by_time

## ensure correct data structure
for (i in c(1:4,6:7)){milk_RA_by_time_plot_data[,i] <- as.factor(as.character(milk_RA_by_time_plot_data[,i]))}

## shorten bacteria names
milk_RA_by_time_plot_data$y <- as.factor(as.character(gsub("g__","", milk_RA_by_time_plot_data$y)))
milk_RA_by_time_plot_data$y <- as.factor(as.character(gsub("_clr","", milk_RA_by_time_plot_data$y)))

## sort bacteria by estimate from high to low (so they will be displayed from low to high in the plot)
milk_RA_by_time_plot_data <- milk_RA_by_time_plot_data[order(-milk_RA_by_time_plot_data$estimate),]
milk_ordered_genera <- milk_RA_by_time_plot_data$y
milk_RA_by_time_plot_data$y <- factor(milk_RA_by_time_plot_data$y, levels = c(milk_ordered_genera))

## add column to color significant results
milk_RA_by_time_plot_data$mycolor <- NA
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH>=0.05, "mycolor"] <- "FDR>=0.05"
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<0.05, "mycolor"] <- "FDR<0.05" #FDR<0.05 *
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<0.01, "mycolor"] <- "FDR<0.01" #FDR<0.01 **
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<0.001, "mycolor"] <- "FDR<0.001" #FDR<0.001 ***
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<1.0e-5, "mycolor"] <- "FDR<1.0e-5"
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<1.0e-10, "mycolor"] <- "FDR<1.0e-10"
milk_RA_by_time_plot_data$mycolor <- as.factor(as.character(milk_RA_by_time_plot_data$mycolor))
milk_RA_by_time_plot_data$mycolor <- factor(milk_RA_by_time_plot_data$mycolor, levels = c("FDR<1.0e-10", "FDR<1.0e-5", "FDR<0.001", "FDR<0.01", "FDR<0.05", "FDR>=0.05"))

## add column to for size of significant results
milk_RA_by_time_plot_data$mysize <- NA
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH>=0.05, "mysize"] <- "1"
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<0.05, "mysize"] <- "2" #FDR<0.05 *
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<0.01, "mysize"] <- "3" #FDR<0.01 **
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<0.001, "mysize"] <- "4" #FDR<0.001 ***
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<1.0e-5, "mysize"] <- "5"
milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<1.0e-10, "mysize"] <- "6"

## plot estimates by bacterial genus
milkRAplot1 <- ggplot(milk_RA_by_time_plot_data, aes(x=estimate, y=y, color=mycolor))+
  geom_point(aes(size=mysize))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(limits = c(-0.55,0.55), breaks = c(-0.5, -0.25, 0, 0.25, 0.5))+
  scale_color_manual(values = c("FDR<1.0e-10" = "palevioletred4",
                                "FDR<1.0e-5" = "indianred3",
                                "FDR<0.001" = "coral1",
                                "FDR<0.01" = "darkorange1",
                                "FDR<0.05" = "tan2",
                                "FDR>=0.05" = "#999999"))+
  labs(x="Estimate", y="", title="Human milk")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/time_point/231212_milk_RA_estimates_by_genus_all.pdf"), milkRAplot1, dpi=500, width=20, height=20, units="cm")

## plot estimates by bacterial genus, only plot FDR<0.05 results
milkRAplot2 <- ggplot(milk_RA_by_time_plot_data[milk_RA_by_time_plot_data$FDR_BH<0.05,], aes(x=estimate, y=y, color=mycolor))+
  geom_point(aes(size=mysize))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_x_continuous(limits = c(-0.55,0.55), breaks = c(-0.5, -0.25, 0, 0.25, 0.5))+
  scale_color_manual(values = c("FDR<1.0e-10" = "palevioletred4",
                                "FDR<1.0e-5" = "indianred3",
                                "FDR<0.001" = "coral1",
                                "FDR<0.01" = "darkorange1",
                                "FDR<0.05" = "tan2"))+
  labs(x="Estimate", y="", title="Human milk")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/time_point/231212_milk_RA_estimates_by_genus_only_sign_no_ytext.pdf"), milkRAplot2, dpi=500, width=10, height=15, units="cm")


# ### ===== 11.4. PLOT MEAN RELATIVE ABUNDANCE FOR EACH BACTERIUM BY TIME POINT ===== ###
# 
# ## Plot only the bacterial genera with ≥10% prevalence and ≥1% in the samples
# 
# ## extract data for plotting from milkRAprev10 data frame
# milk_RA_by_time_plot_data2 <- milkRAprev10[,grep("g__", colnames(milkRAprev10))]
# 
# ## shorten colnames
# colnames(milk_RA_by_time_plot_data2) <- gsub("g__", "", colnames(milk_RA_by_time_plot_data2))
# 
# ## calculate means per bacterium a nd time point
# milkRAmeansM1 <- c()
# for (i in 1:ncol(milk_RA_by_time_plot_data2[grep("mother_human_milk__1_month", rownames(milk_RA_by_time_plot_data2)),])){milkRAmeansM1[i] <- 100*mean(milk_RA_by_time_plot_data2[grep("mother_human_milk__1_month", rownames(milk_RA_by_time_plot_data2)), i])}
# names(milkRAmeansM1) <- colnames(milk_RA_by_time_plot_data2)
# 
# milkRAmeansM2 <- c()
# for (i in 1:ncol(milk_RA_by_time_plot_data2[grep("mother_human_milk__2_months", rownames(milk_RA_by_time_plot_data2)),])){milkRAmeansM2[i] <- 100*mean(milk_RA_by_time_plot_data2[grep("mother_human_milk__2_months", rownames(milk_RA_by_time_plot_data2)), i])}
# names(milkRAmeansM2) <- colnames(milk_RA_by_time_plot_data2)
# 
# milkRAmeansM3 <- c()
# for (i in 1:ncol(milk_RA_by_time_plot_data2[grep("mother_human_milk__3_months", rownames(milk_RA_by_time_plot_data2)),])){milkRAmeansM3[i] <- 100*mean(milk_RA_by_time_plot_data2[grep("mother_human_milk__3_months", rownames(milk_RA_by_time_plot_data2)), i])}
# names(milkRAmeansM3) <- colnames(milk_RA_by_time_plot_data2)
# 
# milkRAmeansM6 <- c()
# for (i in 1:ncol(milk_RA_by_time_plot_data2[grep("mother_human_milk__6_months", rownames(milk_RA_by_time_plot_data2)),])){milkRAmeansM6[i] <- 100*mean(milk_RA_by_time_plot_data2[grep("mother_human_milk__6_months", rownames(milk_RA_by_time_plot_data2)), i])}
# names(milkRAmeansM6) <- colnames(milk_RA_by_time_plot_data2)
# 
# 
# ## combine means in 1 data frame
# milk_RA_by_time_plot_data3 <- as.data.frame(rbind(milkRAmeansM1, milkRAmeansM2, milkRAmeansM3, milkRAmeansM6))
# 
# ## add metadata columns needed for plotting
# milk_RA_by_time_plot_data3$sample_origin_type <- as.factor(as.character(c(rep("mother_human_milk",4))))
# milk_RA_by_time_plot_data3$time_point <- as.factor(as.character(c("1_month", "2_months", "3_months", "6_months")))
# 
# 
# ## convert data frame from wide to long format
# # install.packages("reshape")
# library(reshape)
# milk_RA_by_time_plot_data3_long <- melt(milk_RA_by_time_plot_data3, id.vars=colnames(milk_RA_by_time_plot_data3)[c(37:38)], variable_name="Genus")
# 
# ## change the order of bacteria to be the same as in the milkRAplot1
# milk_RA_by_time_plot_data3_long$Genus <- factor(milk_RA_by_time_plot_data3_long$Genus, levels = c(rev(milk_ordered_genera)))
# 
# 
# ## plot bars
# milkRAplot3 <- ggplot(milk_RA_by_time_plot_data3_long, aes(x=time_point, y=value, fill=time_point))+
#   geom_col()+
#   theme_bw()+
#   facet_grid(Genus~., scales="free_y")+
#   scale_fill_manual(values = c("1_month" = "#F0E442",
#                                 "2_months" = "gold2",
#                                 "3_months" = "#E69F00",
#                                 "6_months" = "sienna3"))+
#   theme(panel.grid = element_blank(),
#         legend.position = "none",
#         strip.text.y = element_text(angle=0))
# ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/time_point/231212_milk_RA_prev10_by_time_bars.pdf"), milkRAplot3, dpi=500, width=20, height=60, units="cm")
# 
# ## get mean values
# milkRAplot4 <- ggplot(milk_RA_by_time_plot_data3_long, aes(x=time_point, y=value, fill=time_point))+
#   geom_col()+
#   theme_bw()+
#   facet_grid(Genus~.)+
#   geom_text(aes(x=time_point, y=10, label=signif(value,3)))+
#   scale_fill_manual(values = c("1_month" = "#F0E442",
#                                "2_months" = "gold2",
#                                "3_months" = "#E69F00",
#                                "6_months" = "sienna3"))+
#   theme(panel.grid = element_blank(),
#         legend.position = "none",
#         strip.text.y = element_text(angle=0))
# ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/time_point/231212_milk_RA_prev10_by_time_mean_values.pdf"), milkRAplot4, dpi=500, width=20, height=60, units="cm")


##### =========================== 12. CHECK THE EFFECT OF TIME ON INFANT FAECAL ALPHA DIVERSITY AND INFANT FAECAL RELATIVE ABUNDANCES =========================== #####

# install.packages("lmerTest")
library(lmerTest)


### ===== 12.1. MODEL FUNCTION WITH CORRECTION FOR DNA ISOLATION BATCH/METHOD (FIXED EFFECT), CLEAN READ COUNT (FIXED EFFECT) AND INDIVIDUAL ID (RANDOM EFFECT) ===== ###

## function for mixed model with corrections 
run.lmer.cor.ymethod.yreads.yID <- function(datasetname, inputdata, xcolumn, ycolumn){
  p <- c()
  
  levels <- c()
  n_levels <- c()
  intercept <- c()
  e <- c()
  
  d <- c()
  x <- c()
  y <- c()
  n <- c()
  stat <- c()
  
  for (i in xcolumn){
    
    for (j in ycolumn){
      print(paste0("i: ",i,"j: ",j))
      
      print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #mixed model with confounders and WITH the phenotype of interest (i)
      m1 <- lmerTest::lmer(my_df[,j] ~ DNA_isolation_method + seq_16S_n_reads_clean + (1|infant_ID) + my_df[,i], REML=F, data=my_df)
      sum1 <- summary(m1)
      print(paste0("Printing sum1"))
      print(sum1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      p <- c(p, sum1$coefficients[rownames(sum1$coefficients)=="my_df[, i]", "Pr(>|t|)"])
      levels <- c(levels, paste(collapse=";",levels(my_df[,i])))
      n_levels <- c(n_levels, paste(collapse=";",table(my_df[,i])))
      intercept <- c(intercept, c(sum1$coefficients[rownames(sum1$coefficients)=="(Intercept)", "Estimate"]))
      e <- c(e, c(sum1$coefficients[rownames(sum1$coefficients)=="my_df[, i]", "Estimate"])) #first rows are for confounders -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, paste0("lmer_correction_yDNAmethod_yreads_yID"))
      d <- c(d, datasetname)
      x <- c(x, colnames(inputdata)[i])
      y <- c(y, colnames(inputdata)[j])
      n <- c(n, nrow(my_df))
    }
  }
  
  #save results in data frame
  print(paste0("Combine results in data frame and save results"))
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, n_levels=n_levels, intercept=intercept, estimate=e, p=p)
  
  # Multiple testing correction (BH)
  res$FDR_BH <- p.adjust(res$p, method="BH")
  
  # Order data frame by p value
  res <- res[order(res$FDR_BH, res$p),]
  
  return(res)
}


## function for mixed model with corrections: additional correction for infant_birth_delivery mode
run.lmer.cor.ymethod.yreads.yID.ydelmode <- function(datasetname, inputdata, xcolumn, ycolumn){
  p <- c()
  
  levels <- c()
  n_levels <- c()
  intercept <- c()
  e <- c()
  
  d <- c()
  x <- c()
  y <- c()
  n <- c()
  stat <- c()
  
  for (i in xcolumn){
    
    for (j in ycolumn){
      print(paste0("i: ",i,"j: ",j))
      
      print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]) & !is.na(inputdata$infant_birth_delivery_mode),]
      
      #mixed model with confounders and WITH the phenotype of interest (i)
      m1 <- lmerTest::lmer(my_df[,j] ~ DNA_isolation_method + seq_16S_n_reads_clean + (1|infant_ID) + infant_birth_delivery_mode + my_df[,i], REML=F, data=my_df)
      sum1 <- summary(m1)
      print(paste0("Printing sum1"))
      print(sum1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      p <- c(p, sum1$coefficients[rownames(sum1$coefficients)=="my_df[, i]", "Pr(>|t|)"])
      levels <- c(levels, paste(collapse=";",levels(my_df[,i])))
      n_levels <- c(n_levels, paste(collapse=";",table(my_df[,i])))
      intercept <- c(intercept, c(sum1$coefficients[rownames(sum1$coefficients)=="(Intercept)", "Estimate"]))
      e <- c(e, c(sum1$coefficients[rownames(sum1$coefficients)=="my_df[, i]", "Estimate"])) #first rows are for confounders -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, paste0("lmer_correction_yDNAmethod_yreads_yID_ydelmode"))
      d <- c(d, datasetname)
      x <- c(x, colnames(inputdata)[i])
      y <- c(y, colnames(inputdata)[j])
      n <- c(n, nrow(my_df))
    }
  }
  
  #save results in data frame
  print(paste0("Combine results in data frame and save results"))
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, n_levels=n_levels, intercept=intercept, estimate=e, p=p)
  
  # Multiple testing correction (BH)
  res$FDR_BH <- p.adjust(res$p, method="BH")
  
  # Order data frame by p value
  res <- res[order(res$FDR_BH, res$p),]
  
  return(res)
}


### ===== 12.2. COMPARE INFANT FAECAL ALPHA DIVERSITY BY NUMERIC TIME POINTS ===== ###

## For statistical analyses, use data frame sinvr with inverse-rank transformed alpha diversity measures!

## Compare infant faecal alpha diversity between time points (numeric)
inff_lmer_results_alpha_div_by_time <- run.lmer.cor.ymethod.yreads.yID(datasetname = "infant_faeces",
                                                                      inputdata = sinvr[sinvr$sample_origin_type=="infant_faeces",],
                                                                      xcolumn = c(10), #numeric time point
                                                                      ycolumn = c(528:530)) #inverse-rank transformed alpha diversity measures
write.table(inff_lmer_results_alpha_div_by_time, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/time_point/231212_lmer_results_infant_faeces_alpha_diversity_by_numeric_time_point.txt", sep="\t", row.names=F, quote=F)

## Save results without Simpson diversity, re-calculate FDR
inff_lmer_results_alpha_div_by_time_no_simpson <- inff_lmer_results_alpha_div_by_time[inff_lmer_results_alpha_div_by_time$y!="alpha_div_simpson_invr",]
inff_lmer_results_alpha_div_by_time_no_simpson$FDR <- p.adjust(inff_lmer_results_alpha_div_by_time_no_simpson$p, method="BH")
inff_lmer_results_alpha_div_by_time_no_simpson <- inff_lmer_results_alpha_div_by_time_no_simpson[order(inff_lmer_results_alpha_div_by_time_no_simpson$FDR, inff_lmer_results_alpha_div_by_time_no_simpson$p),]
write.table(inff_lmer_results_alpha_div_by_time_no_simpson, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/time_point/240216_lmer_results_infant_faeces_alpha_diversity_richness_shannon_by_numeric_time_point.txt", sep="\t", row.names=F, quote=F)


## Check that results remain when additionally correcting for delivery mode
inff_lmer_results_alpha_div_by_time_cor_delmode <- run.lmer.cor.ymethod.yreads.yID.ydelmode(datasetname = "infant_faeces",
                                                                       inputdata = sinvr[sinvr$sample_origin_type=="infant_faeces",],
                                                                       xcolumn = c(10), #numeric time point
                                                                       ycolumn = c(528:530)) #inverse-rank transformed alpha diversity measures
write.table(inff_lmer_results_alpha_div_by_time_cor_delmode, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/time_point/240216_lmer_results_infant_faeces_alpha_diversity_by_numeric_time_point_cor_delmode.txt", sep="\t", row.names=F, quote=F)

## Save results without Simpson diversity, re-calculate FDR
inff_lmer_results_alpha_div_by_time_cor_delmode_no_simpson <- inff_lmer_results_alpha_div_by_time_cor_delmode[inff_lmer_results_alpha_div_by_time_cor_delmode$y!="alpha_div_simpson_invr",]
inff_lmer_results_alpha_div_by_time_cor_delmode_no_simpson$FDR <- p.adjust(inff_lmer_results_alpha_div_by_time_cor_delmode_no_simpson$p, method="BH")
inff_lmer_results_alpha_div_by_time_cor_delmode_no_simpson <- inff_lmer_results_alpha_div_by_time_cor_delmode_no_simpson[order(inff_lmer_results_alpha_div_by_time_cor_delmode_no_simpson$FDR, inff_lmer_results_alpha_div_by_time_cor_delmode_no_simpson$p),]
write.table(inff_lmer_results_alpha_div_by_time_cor_delmode_no_simpson, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/time_point/240216_lmer_results_infant_faeces_alpha_diversity_richness_shannon_by_numeric_time_point_cor_delmode.txt", sep="\t", row.names=F, quote=F)


### ===== 12.3. COMPARE INFANT FAECAL RELATIVE ABUNDANCES BY NUMERIC TIME POINTS (TABLE S12) ===== ###

## For statistical analyses, use data frame xx with clr-transformed relative abundances of the bacterial genera with ≥10% prevalence in infant faecal samples.

## Compare infant faecal bacterial relative abundances between time points (numeric)
inff_lmer_results_RA_by_time <- run.lmer.cor.ymethod.yreads.yID(datasetname = "infant_faeces",
                                                                inputdata = inffRAprev10clr,
                                                                xcolumn = c(10), #numeric time point
                                                                ycolumn = c(531:560))
write.table(inff_lmer_results_RA_by_time, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/231212_lmer_results_inff_RA_by_numeric_time_point.txt", sep="\t", row.names=F, quote=F)


## Check that results remain when additionally correcting for delivery mode
inff_lmer_results_RA_by_time_cor_delmode <- run.lmer.cor.ymethod.yreads.yID.ydelmode(datasetname = "infant_faeces",
                                                                inputdata = inffRAprev10clr,
                                                                xcolumn = c(10), #numeric time point
                                                                ycolumn = c(531:560))
write.table(inff_lmer_results_RA_by_time_cor_delmode, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/240216_lmer_results_inff_RA_by_numeric_time_point_cor_delmode.txt", sep="\t", row.names=F, quote=F)


##### =========================== 13. PLOT EFFECT OF TIME ON INFANT FAECAL RELATIVE ABUNDANCES =========================== #####

### ===== 13.1. FUNCTION TO PLOT SIGNIFICANT EFFECTS OF TIME ON RELATIVE ABUNDANCES ===== ###

plot.inff.RA.by.time.sign.results <- function(inputdata, signresults, outputname){
  for (i in signresults[,"y"]){
    print(paste0("Plotting effect of time on ", i))
    # print(paste0(intercept=signresults[signresults$y==paste0(i),"intercept"]))
    # print(paste0(slope=signresults[signresults$y==paste0(i),"estimate"]))
    
    RAbytimeplot <- ggplot(inputdata, aes(x=time_point_numeric, y=inputdata[,i], color=time_point))+
      geom_point(alpha=0.5)+
      geom_boxplot(alpha=0, color="black", aes(group=time_point))+
      theme_bw()+
      theme(legend.position = "none",
            panel.grid = element_blank())+
      ggtitle(paste0(i, " RA by time; n=", signresults[signresults$y==i, "n_total"]))+
      labs(x="Month(s) postpartum", y=paste0(i))+
      scale_color_manual(values = c("1_month" = "lightgreen",
                                    "2_months" = "limegreen",
                                    "3_months" = "#009E73",
                                    "6_months" = "darkolivegreen"))+
      scale_x_continuous(limits = c(0.5,6.5), breaks = c(1,2,3,4,5,6))+
      geom_abline(intercept=signresults[signresults$y==paste0(i),"intercept"],
                  slope=signresults[signresults$y==paste0(i),"estimate"])
    ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/240216_infant_faeces_", i, "_RA_by_time", outputname, ".pdf"), RAbytimeplot, dpi=500, width=15, height=15, units="cm")
  }
}


### ===== 13.2. PLOT SIGNIFICANT EFFECTS OF TIME ON INFANT FAECAL RELATIVE ABUNDANCES ===== ###

## results that were corrected for DNA isolation methods, reads and ID
plot.inff.RA.by.time.sign.results(inputdata = inffRAprev10clr,
                                  signresults = inff_lmer_results_RA_by_time[inff_lmer_results_RA_by_time$FDR_BH<0.05,],
                                  outputname = "")

## results that were corrected for DNA isolation methods, reads, ID and delivery mode
plot.inff.RA.by.time.sign.results(inputdata = inffRAprev10clr,
                                  signresults = inff_lmer_results_RA_by_time_cor_delmode[inff_lmer_results_RA_by_time_cor_delmode$FDR_BH<0.05,],
                                  outputname = "_cor_delmode")


### ===== 13.3. PLOT ESTIMATES FOR EACH BACTERIUM (USED FOR FIGURE 4E) ===== ###

## data corrected for DNA isolation method, reads and ID
inff_lmer_results_RA_by_time <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/231212_lmer_results_inff_RA_by_numeric_time_point.txt", header=T, sep="\t", stringsAsFactors=T)

## save data separately to prepare for plotting
inff_RA_by_time_plot_data <- inff_lmer_results_RA_by_time

## ensure correct data structure
for (i in c(1:4,6:7)){inff_RA_by_time_plot_data[,i] <- as.factor(as.character(inff_RA_by_time_plot_data[,i]))}

## shorten bacteria names
inff_RA_by_time_plot_data$y <- as.factor(as.character(gsub("g__","", inff_RA_by_time_plot_data$y)))
inff_RA_by_time_plot_data$y <- as.factor(as.character(gsub("_clr","", inff_RA_by_time_plot_data$y)))
inff_RA_by_time_plot_data$y <- as.factor(as.character(gsub("\\[Ruminococcus\\]","Ruminococcus", inff_RA_by_time_plot_data$y)))

## sort bacteria by estimate from high to low (so they will be displayed from low to high in the plot)
inff_RA_by_time_plot_data <- inff_RA_by_time_plot_data[order(-inff_RA_by_time_plot_data$estimate),]
inff_ordered_genera <- inff_RA_by_time_plot_data$y
inff_RA_by_time_plot_data$y <- factor(inff_RA_by_time_plot_data$y, levels = c(inff_ordered_genera))

## add column to color significant results
inff_RA_by_time_plot_data$mycolor <- NA
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH>=0.05, "mycolor"] <- "FDR>=0.05"
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<0.05, "mycolor"] <- "FDR<0.05" #FDR<0.05 *
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<0.01, "mycolor"] <- "FDR<0.01" #FDR<0.01 **
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<0.001, "mycolor"] <- "FDR<0.001" #FDR<0.001 ***
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<1.0e-5, "mycolor"] <- "FDR<1.0e-5"
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<1.0e-10, "mycolor"] <- "FDR<1.0e-10"
inff_RA_by_time_plot_data$mycolor <- as.factor(as.character(inff_RA_by_time_plot_data$mycolor))
inff_RA_by_time_plot_data$mycolor <- factor(inff_RA_by_time_plot_data$mycolor, levels = c("FDR<1.0e-10", "FDR<1.0e-5", "FDR<0.001", "FDR<0.01", "FDR<0.05", "FDR>=0.05"))

## add column to for size of significant results
inff_RA_by_time_plot_data$mysize <- NA
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH>=0.05, "mysize"] <- "1"
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<0.05, "mysize"] <- "2" #FDR<0.05 *
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<0.01, "mysize"] <- "3" #FDR<0.01 **
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<0.001, "mysize"] <- "4" #FDR<0.001 ***
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<1.0e-5, "mysize"] <- "5"
inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<1.0e-10, "mysize"] <- "6"


## plot estimates by bacterial genus
inffRAplot1 <- ggplot(inff_RA_by_time_plot_data, aes(x=estimate, y=y, color=mycolor))+
  geom_point(aes(size=mysize))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(limits = c(-0.85,0.85), breaks = c(-0.8, -0.4, 0, 0.4, 0.8))+
  scale_color_manual(values = c("FDR<1.0e-10" = "palevioletred4",
                                "FDR<1.0e-5" = "indianred3",
                                "FDR<0.001" = "coral1",
                                "FDR<0.01" = "darkorange1",
                                "FDR<0.05" = "tan2",
                                "FDR>=0.05" = "#999999"))+
  labs(x="Estimate", y="", title="Infant faeces")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/231212_infant_faeces_RA_estimates_by_genus_all.pdf"), inffRAplot1, dpi=500, width=20, height=20, units="cm")

## plot estimates by bacterial genus, only plot FDR<0.05 results
inffRAplot2 <- ggplot(inff_RA_by_time_plot_data[inff_RA_by_time_plot_data$FDR_BH<0.05,], aes(x=estimate, y=y, color=mycolor))+
  geom_point(aes(size=mysize))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_x_continuous(limits = c(-0.85,0.85), breaks = c(-0.8, -0.4, 0, 0.4, 0.8))+
  scale_color_manual(values = c("FDR<1.0e-10" = "palevioletred4",
                                "FDR<1.0e-5" = "indianred3",
                                "FDR<0.001" = "coral1",
                                "FDR<0.01" = "darkorange1",
                                "FDR<0.05" = "tan2"))+
  labs(x="Estimate", y="", title="Infant faeces")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/231212_infant_faeces_RA_estimates_by_genus_only_sign_no_ytext.pdf"), inffRAplot2, dpi=500, width=10, height=15, units="cm")


## data corrected for DNA isolation method, reads, ID and delivery mode
inff_lmer_results_RA_by_time_cor_delmode <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/240216_lmer_results_inff_RA_by_numeric_time_point_cor_delmode.txt", header=T, sep="\t", stringsAsFactors=T)

## save data separately to prepare for plotting
inff_RA_by_time_plot_data_cor_delmode <- inff_lmer_results_RA_by_time_cor_delmode

## ensure correct data structure
for (i in c(1:4,6:7)){inff_RA_by_time_plot_data_cor_delmode[,i] <- as.factor(as.character(inff_RA_by_time_plot_data_cor_delmode[,i]))}

## shorten bacteria names
inff_RA_by_time_plot_data_cor_delmode$y <- as.factor(as.character(gsub("g__","", inff_RA_by_time_plot_data_cor_delmode$y)))
inff_RA_by_time_plot_data_cor_delmode$y <- as.factor(as.character(gsub("_clr","", inff_RA_by_time_plot_data_cor_delmode$y)))
inff_RA_by_time_plot_data_cor_delmode$y <- as.factor(as.character(gsub("\\[Ruminococcus\\]","Ruminococcus", inff_RA_by_time_plot_data_cor_delmode$y)))

## sort bacteria by estimate from high to low (so they will be displayed from low to high in the plot)
inff_RA_by_time_plot_data_cor_delmode <- inff_RA_by_time_plot_data_cor_delmode[order(-inff_RA_by_time_plot_data_cor_delmode$estimate),]
inff_ordered_genera <- inff_RA_by_time_plot_data_cor_delmode$y
inff_RA_by_time_plot_data_cor_delmode$y <- factor(inff_RA_by_time_plot_data_cor_delmode$y, levels = c(inff_ordered_genera))

## add column to color significant results
inff_RA_by_time_plot_data_cor_delmode$mycolor <- NA
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH>=0.05, "mycolor"] <- "FDR>=0.05"
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<0.05, "mycolor"] <- "FDR<0.05" #FDR<0.05 *
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<0.01, "mycolor"] <- "FDR<0.01" #FDR<0.01 **
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<0.001, "mycolor"] <- "FDR<0.001" #FDR<0.001 ***
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<1.0e-5, "mycolor"] <- "FDR<1.0e-5"
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<1.0e-10, "mycolor"] <- "FDR<1.0e-10"
inff_RA_by_time_plot_data_cor_delmode$mycolor <- as.factor(as.character(inff_RA_by_time_plot_data_cor_delmode$mycolor))
inff_RA_by_time_plot_data_cor_delmode$mycolor <- factor(inff_RA_by_time_plot_data_cor_delmode$mycolor, levels = c("FDR<1.0e-10", "FDR<1.0e-5", "FDR<0.001", "FDR<0.01", "FDR<0.05", "FDR>=0.05"))

## add column to for size of significant results
inff_RA_by_time_plot_data_cor_delmode$mysize <- NA
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH>=0.05, "mysize"] <- "1"
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<0.05, "mysize"] <- "2" #FDR<0.05 *
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<0.01, "mysize"] <- "3" #FDR<0.01 **
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<0.001, "mysize"] <- "4" #FDR<0.001 ***
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<1.0e-5, "mysize"] <- "5"
inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<1.0e-10, "mysize"] <- "6"


## plot estimates by bacterial genus
inffRAplot1cordelmode <- ggplot(inff_RA_by_time_plot_data_cor_delmode, aes(x=estimate, y=y, color=mycolor))+
  geom_point(aes(size=mysize))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank())+
  scale_x_continuous(limits = c(-0.85,0.85), breaks = c(-0.8, -0.4, 0, 0.4, 0.8))+
  scale_color_manual(values = c("FDR<1.0e-10" = "palevioletred4",
                                "FDR<1.0e-5" = "indianred3",
                                "FDR<0.001" = "coral1",
                                "FDR<0.01" = "darkorange1",
                                "FDR<0.05" = "tan2",
                                "FDR>=0.05" = "#999999"))+
  labs(x="Estimate", y="", title="Infant faeces")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/240216_infant_faeces_RA_estimates_by_genus_all_cor_delmode.pdf"), inffRAplot1cordelmode, dpi=500, width=20, height=20, units="cm")

## plot estimates by bacterial genus, only plot FDR<0.05 results
inffRAplot2cordelmode <- ggplot(inff_RA_by_time_plot_data_cor_delmode[inff_RA_by_time_plot_data_cor_delmode$FDR_BH<0.05,], aes(x=estimate, y=y, color=mycolor))+
  geom_point(aes(size=mysize))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  scale_x_continuous(limits = c(-0.85,0.85), breaks = c(-0.8, -0.4, 0, 0.4, 0.8))+
  scale_color_manual(values = c("FDR<1.0e-10" = "palevioletred4",
                                "FDR<1.0e-5" = "indianred3",
                                "FDR<0.001" = "coral1",
                                "FDR<0.01" = "darkorange1",
                                "FDR<0.05" = "tan2"))+
  labs(x="Estimate", y="", title="Infant faeces")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/time_point/240216_infant_faeces_RA_estimates_by_genus_only_sign_no_ytext_cor_delmode.pdf"), inffRAplot2cordelmode, dpi=500, width=10, height=15, units="cm")


##### =========================== 14. COMBINED PLOT FOR EFFECT OF TIME ON MILK AND INFANT FAECAL RELATIVE ABUNDANCE ESTIMATES FOR EACH BACTERIUM (FIGURE 4E) =========================== #####

## corrected for DNA batch or method and reads and ID
RA_by_time_plot_data <- as.data.frame(rbind(milk_RA_by_time_plot_data, inff_RA_by_time_plot_data))
RA_by_time_plot_data$mysize <- as.factor(as.character(RA_by_time_plot_data$mysize))

## ensure desired sortening from low to high in each panel
RA_by_time_plot_data <- RA_by_time_plot_data[order(RA_by_time_plot_data$estimate),]
RA_by_time_plot_data$labels <- c(paste0(RA_by_time_plot_data$dataset, "_", RA_by_time_plot_data$y)) #make labels unique
RA_by_time_plot_data$labels <- as.factor(as.character(RA_by_time_plot_data$labels))
RA_by_time_plot_data$labels <- factor(RA_by_time_plot_data$labels, levels = c("infant_faeces_Staphylococcus","milk_Staphylococcus","infant_faeces_Lactobacillus","milk_Corynebacterium","infant_faeces_Parabacteroides",
                                                                                      "milk_Gemella","milk_Lactobacillus","milk_Streptococcus","milk_Bifidobacterium","infant_faeces_Haemophilus","infant_faeces_Collinsella",
                                                                                      "infant_faeces_Streptococcus","infant_faeces_Gemella","infant_faeces_Cutibacterium","infant_faeces_Bifidobacterium","infant_faeces_Conservatibacter",
                                                                                      "milk_Anaerococcus","milk_Paracoccus","infant_faeces_Citrobacter","infant_faeces_Enterobacter","milk_Finegoldia","infant_faeces_Bilophila",
                                                                                      "milk_Cutibacterium","milk_Atopobium","infant_faeces_Sutterella","infant_faeces_Actinomyces","milk_Allorhizobium_Neorhizobium_Pararhizobium_Rhizobium",
                                                                                      "milk_Bacillus","milk_Bergeyella","infant_faeces_Bacteroides","milk_Veillonella","milk_Enhydrobacter","infant_faeces_Incertae_Sedis",
                                                                                      "milk_Stenotrophomonas","milk_Enterobacter","milk_Moraxella","milk_Chryseobacterium","milk_Sphingomonas","milk_Acinetobacter",
                                                                                      "infant_faeces_Lachnoclostridium","infant_faeces_Blautia","milk_Prevotella_7","milk_Escherichia_Shigella","milk_Bacteroides","milk_Brevundimonas",
                                                                                      "milk_Prevotella","milk_Fusobacterium","infant_faeces_Akkermansia","infant_faeces_Ruminococcus_gnavus_group","milk_Haemophilus",
                                                                                      "milk_Porphyromonas","infant_faeces_Klebsiella","milk_Alloprevotella","infant_faeces_Erysipelatoclostridium","infant_faeces_Lacticaseibacillus",
                                                                                      "infant_faeces_Clostridium_sensu_stricto_1","infant_faeces_Veillonella","infant_faeces_Eggerthella","milk_Pseudomonas","milk_Neisseria",
                                                                                      "milk_Actinomyces","milk_Rothia","milk_Granulicatella","infant_faeces_Aquamonas","infant_faeces_Enterococcus","infant_faeces_Escherichia_Shigella"))


## plot estimates by bacterial genus, only plot FDR<0.05 results
milkinffRAplot2 <- ggplot(RA_by_time_plot_data[RA_by_time_plot_data$FDR_BH<0.05,], aes(x=estimate, y=labels, color=mycolor))+
  coord_flip()+
  geom_point(aes(size=mysize))+
  theme_bw()+
  scale_x_continuous(limits = c(-0.85,0.85), breaks = c(-0.8, 0, 0.8))+
  scale_color_manual(values = c("FDR<1.0e-10" = "palevioletred4",
                                "FDR<1.0e-5" = "indianred3",
                                "FDR<0.001" = "coral1",
                                "FDR<0.01" = "darkorange1",
                                "FDR<0.05" = "tan2"))+
  labs(x="Estimate", y="")+
  facet_grid(.~dataset, scales = "free_x", space = "free_x")+
  theme(panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background = element_rect(color = "white", fill="white"))
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/relative_abundances/231221_milk_and_infant_faeces_RA_estimates_by_genus_only_sign_no_ytext.pdf"), milkinffRAplot2, dpi=500, width=45, height=13, units="cm")


## corrected for DNA batch or method and reads and ID and for delivery mode (for faecal RA)
RA_by_time_plot_data_2 <- as.data.frame(rbind(milk_RA_by_time_plot_data, inff_RA_by_time_plot_data_cor_delmode))
RA_by_time_plot_data_2$mysize <- as.factor(as.character(RA_by_time_plot_data_2$mysize))

## ensure desired sortening from low to high in each panel
RA_by_time_plot_data_2 <- RA_by_time_plot_data_2[order(RA_by_time_plot_data_2$estimate),]
RA_by_time_plot_data_2$labels <- c(paste0(RA_by_time_plot_data_2$dataset, "_", RA_by_time_plot_data_2$y)) #make labels unique
RA_by_time_plot_data_2$labels <- as.factor(as.character(RA_by_time_plot_data_2$labels))
RA_by_time_plot_data_2$labels <- factor(RA_by_time_plot_data_2$labels, levels = c("infant_faeces_Staphylococcus","milk_Staphylococcus","infant_faeces_Lactobacillus","milk_Corynebacterium","infant_faeces_Parabacteroides",
                                                                              "milk_Gemella","milk_Lactobacillus","milk_Streptococcus","milk_Bifidobacterium","infant_faeces_Haemophilus","infant_faeces_Collinsella",
                                                                              "infant_faeces_Streptococcus","infant_faeces_Gemella","infant_faeces_Cutibacterium","infant_faeces_Bifidobacterium","infant_faeces_Conservatibacter",
                                                                              "milk_Anaerococcus","milk_Paracoccus","infant_faeces_Citrobacter","infant_faeces_Enterobacter","milk_Finegoldia","infant_faeces_Bilophila",
                                                                              "milk_Cutibacterium","milk_Atopobium","infant_faeces_Sutterella","infant_faeces_Actinomyces","milk_Allorhizobium_Neorhizobium_Pararhizobium_Rhizobium",
                                                                              "milk_Bacillus","milk_Bergeyella","infant_faeces_Bacteroides","milk_Veillonella","milk_Enhydrobacter","infant_faeces_Incertae_Sedis",
                                                                              "milk_Stenotrophomonas","milk_Enterobacter","milk_Moraxella","milk_Chryseobacterium","milk_Sphingomonas","milk_Acinetobacter",
                                                                              "infant_faeces_Lachnoclostridium","infant_faeces_Blautia","milk_Prevotella_7","milk_Escherichia_Shigella","milk_Bacteroides","milk_Brevundimonas",
                                                                              "milk_Prevotella","milk_Fusobacterium","infant_faeces_Akkermansia","infant_faeces_Ruminococcus_gnavus_group","milk_Haemophilus",
                                                                              "milk_Porphyromonas","infant_faeces_Klebsiella","milk_Alloprevotella","infant_faeces_Erysipelatoclostridium","infant_faeces_Lacticaseibacillus",
                                                                              "infant_faeces_Clostridium_sensu_stricto_1","infant_faeces_Veillonella","infant_faeces_Eggerthella","milk_Pseudomonas","milk_Neisseria",
                                                                              "milk_Actinomyces","milk_Rothia","milk_Granulicatella","infant_faeces_Aquamonas","infant_faeces_Enterococcus","infant_faeces_Escherichia_Shigella"))


## plot estimates by bacterial genus, only plot FDR<0.05 results
milkinffRAplot2_cor <- ggplot(RA_by_time_plot_data_2[RA_by_time_plot_data_2$FDR_BH<0.05,], aes(x=estimate, y=labels, color=mycolor))+
  coord_flip()+
  geom_point(aes(size=mysize))+
  theme_bw()+
  scale_x_continuous(limits = c(-0.85,0.85), breaks = c(-0.8, 0, 0.8))+
  scale_color_manual(values = c("FDR<1.0e-10" = "palevioletred4",
                                "FDR<1.0e-5" = "indianred3",
                                "FDR<0.001" = "coral1",
                                "FDR<0.01" = "darkorange1",
                                "FDR<0.05" = "tan2"))+
  labs(x="Estimate", y="")+
  facet_grid(.~dataset, scales = "free_x", space = "free_x")+
  theme(panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        strip.background = element_rect(color = "white", fill="white"))
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/relative_abundances/240216_milk_and_infant_faeces_RA_estimates_by_genus_only_sign_no_ytext_cor_delmode.pdf"), milkinffRAplot2_cor, dpi=500, width=45, height=13, units="cm")


##### =========================== 15. HUMAN MILK BETA DIVERSITY BY TIME POINT AND BY INDIVIDUAL ID =========================== #####

## Use data frame dfRAclr, which contains clr-transformed relative abundances of all bacteria detected in the samples.


### ===== 15.1. RUN ADONIS BY TIME POINT (FACTOR) AND INDIVIDUAL ID WITH CORRECTION FOR DNA ISOLATION BATCH/METHOD AND CLEAN READ COUNT ===== ###

## extract milk samples from data frame dfRAclr
milkRAclr <- dfRAclr[grep("mother_human_milk", rownames(dfRAclr)),]

## check if there are NAs in time_point and ID columns (this should not be the case!)
nrow(milkRAclr[is.na(milkRAclr$time_point),]) #0
nrow(milkRAclr[is.na(milkRAclr$mother_ID),]) #0

## create Aitchison distance matrix and add metadata
milkait <- vegdist(milkRAclr[,grep("g__", colnames(milkRAclr))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
milkaitm <- cmdscale(milkait, k=3) #create distance matrix for first 3 components
milkaitd <- data.frame(milkaitm)

milkdf <- merge(milkRAclr, milkaitd, by="row.names")
milkdf <- milkdf[,-1]
rownames(milkdf) <- milkdf$seq_16S_sample_ID

milk.beta.cmd <- cmdscale(as.matrix(milkait), k=2, eig=T)
milkPCoA1 <- milk.beta.cmd$eig[1]/sum(milk.beta.cmd$eig)*100
milkPCoA2 <- milk.beta.cmd$eig[2]/sum(milk.beta.cmd$eig)*100

## run adonis
# milk_a1 <- adonis2(milkait ~ milkdf$DNA_isolation_batch + milkdf$seq_16S_n_reads_clean + milkdf$mother_ID + milkdf$time_point, by="margin") #default: 999 permutations
milk_a2 <- adonis2(milkait ~ milkdf$DNA_isolation_batch + milkdf$seq_16S_n_reads_clean + milkdf$mother_ID + milkdf$time_point, by="margin", permutations=10000) # -> here the p values for both time point and mother ID lower to 9.999e-05
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 10000
# 
# adonis2(formula = milkait ~ milkdf$DNA_isolation_batch + milkdf$seq_16S_n_reads_clean + milkdf$mother_ID + milkdf$time_point, permutations = 10000, by = "margin")
# Df SumOfSqs      R2      F    Pr(>F)    
# milkdf$DNA_isolation_batch     1      463 0.00201 1.8903    0.0117 *  
# milkdf$seq_16S_n_reads_clean   1      251 0.00109 1.0237    0.3824    
# milkdf$mother_ID             335   127815 0.55541 1.5586 9.999e-05 ***
# milkdf$time_point              3     1978 0.00860 2.6939 9.999e-05 ***
# Residual                     384    94004 0.40848                     
# Total                        724   230129 1.00000                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## store adonis results in data frame
milk_R2_ID <- milk_a2$R2[3] #extract the R2 value for mother_ID
milk_R2_tp <- milk_a2$R2[4] #extract the R2 value for time_point (factor)
milk_p_ID <- milk_a2$'Pr(>F)'[3] #extract the p value for mother_ID
milk_p_tp <- milk_a2$'Pr(>F)'[4] #extract the p value for time_point (factor)

adonis2.results.milk <- data.frame(dataset=c(rep("human_milk",2)), statistic=c("adonis2_correction_DNAbatch_reads_ID", "adonis2_correction_DNAbatch_reads_timepoint"),
                                   column=c("time_point", "mother_ID"), n=c(rep(nrow(milkdf),2)), R2=c(milk_R2_tp, milk_R2_ID), p=c(milk_p_tp, milk_p_ID))
adonis2.results.milk

# save results
write.table(adonis2.results.milk, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/beta_diversity/231212_milk_adonis_results_corr_isolbatch_reads_by_time_point_and_ID.txt", sep="\t", row.names=F, quote=F)


### ===== 15.2. MILK PLOT PCOA BY TIME POINT (FACTOR) ===== ###

milk_pcoa_tp <- ggplot(milkdf, aes(x=X1, y=X2, color=time_point))+
  geom_point(alpha=0.5)+
  stat_ellipse(aes(group=time_point))+
  labs(x=paste0("PC1 (", signif(milkPCoA1,5), "%)"), y=paste0("PC2 (", signif(milkPCoA2,5), "%)"), color="Time point",
       title=paste0("Human milk, n=", nrow(milkdf), "\n", "Beta diversity by time point (factor)", "\n",
                    "Adonis with correction: R2=", signif(milk_R2_tp,3), " and p=", signif(milk_p_tp,3)))+
  scale_color_manual(values = c("1_month" = "#F0E442",
                                "2_months" = "gold2",
                                "3_months" = "#E69F00",
                                "6_months" = "sienna3"))+
  theme_bw()
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/beta_diversity/231212_milk_PCoA_PC1_PC2_by_time_point.pdf"), milk_pcoa_tp, dpi=500, width=20, height=15, units="cm")


### ===== 15.3. MILK PLOT PCOA BY INDIVIDUAL ID ===== ###

milk_pcoa_ID <- ggplot(milkdf, aes(x=X1, y=X2, color=mother_ID))+
  geom_point(alpha=0.5)+
  labs(x=paste0("PC1 (", signif(milkPCoA1,5), "%)"), y=paste0("PC2 (", signif(milkPCoA2,5), "%)"),
       title=paste0("Human milk, n=", nrow(milkdf), "\n", "Beta diversity by mother ID", "\n",
                    "Adonis with correction: R2=", signif(milk_R2_ID,3), " and p=", signif(milk_p_ID,3)))+
  theme_bw()+
  theme(legend.position = "none")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/beta_diversity/231212_milk_PCoA_PC1_PC2_by_individual_ID.pdf"), milk_pcoa_ID, dpi=500, width=20, height=15, units="cm")


##### =========================== 16. INFANT FAECAL BETA DIVERSITY BY TIME POINT AND BY INDIVIDUAL ID =========================== #####

## Use data frame dfRAclr, which contains clr-transformed relative abundances of all bacteria detected in the samples.


### ===== 16.1. RUN ADONIS BY TIME POINT (FACTOR) AND INDIVIDUAL ID WITH CORRECTION FOR DNA ISOLATION BATCH/METHOD AND CLEAN READ COUNT ===== ###

## extract infant faecal samples from data frame dfRAclr
inffRAclr <- dfRAclr[grep("infant_faeces", rownames(dfRAclr)),]

## check if there are NAs in time_point and ID columns (this should not be the case!)
nrow(inffRAclr[is.na(inffRAclr$time_point),]) #0
nrow(inffRAclr[is.na(inffRAclr$infant_ID),]) #0

## create Aitchison distance matrix and add metadata
inffait <- vegdist(inffRAclr[,grep("g__", colnames(inffRAclr))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
inffaitm <- cmdscale(inffait, k=3) #create distance matrix for first 3 components
inffaitd <- data.frame(inffaitm)

inffdf <- merge(inffRAclr, inffaitd, by="row.names")
inffdf <- inffdf[,-1]
rownames(inffdf) <- inffdf$seq_16S_sample_ID

inff.beta.cmd <- cmdscale(as.matrix(inffait), k=2, eig=T)
inffPCoA1 <- inff.beta.cmd$eig[1]/sum(inff.beta.cmd$eig)*100
inffPCoA2 <- inff.beta.cmd$eig[2]/sum(inff.beta.cmd$eig)*100

## run adonis
# inff_a1 <- adonis2(inffait ~ inffdf$DNA_isolation_method + inffdf$seq_16S_n_reads_clean + inffdf$infant_ID + inffdf$time_point, by="margin") #default: 999 permutations
inff_a2 <- adonis2(inffait ~ inffdf$DNA_isolation_method + inffdf$seq_16S_n_reads_clean + inffdf$infant_ID + inffdf$time_point, by="margin", permutations=10000) # -> here the p values for both time point and infant ID lower to 9.999e-05
# Permutation test for adonis under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 10000
# 
# adonis2(formula = inffait ~ inffdf$DNA_isolation_method + inffdf$seq_16S_n_reads_clean + inffdf$infant_ID + inffdf$time_point, permutations = 10000, by = "margin")
#                               Df SumOfSqs      R2      F    Pr(>F)    
# inffdf$DNA_isolation_method    1      383 0.00217 1.8077    0.0112 *  
# inffdf$seq_16S_n_reads_clean   1      265 0.00150 1.2494    0.1567    
# inffdf$infant_ID             266   108685 0.61492 1.9267 9.999e-05 ***
# inffdf$time_point              3     4568 0.02585 7.1805 9.999e-05 ***
# Residual                     280    59378 0.33595                     
# Total                        551   176747 1.00000                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## store adonis results in data frame
inff_R2_ID <- inff_a1$R2[3] #extract the R2 value for infant_ID
inff_R2_tp <- inff_a1$R2[4] #extract the R2 value for time_point (factor)
inff_p_ID <- inff_a1$'Pr(>F)'[3] #extract the p value for infant_ID
inff_p_tp <- inff_a1$'Pr(>F)'[4] #extract the p value for time_point (factor)

adonis2.results.inff <- data.frame(dataset=c(rep("infant_faeces",2)), statistic=c("adonis2_correction_DNAmethod_reads_ID", "adonis2_correction_DNAmethod_reads_timepoint"),
                                   column=c("time_point", "infant_ID"), n=c(rep(nrow(inffdf),2)), R2=c(inff_R2_tp, inff_R2_ID), p=c(inff_p_tp, inff_p_ID))
adonis2.results.inff

# save results
write.table(adonis2.results.inff, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/beta_diversity/231212_infant_faeces_adonis_results_corr_isolmethod_reads_by_time_point_and_ID.txt", sep="\t", row.names=F, quote=F)


# ## additionally, include corrections for infant_birth_delivery_mode and infant_ffq_feeding_mode
# 
# ## check if there are NAs in infant_birth_delivery_mode and infant_ffq_feeding_mode columns
# nrow(inffRAclr[is.na(inffRAclr$infant_birth_delivery_mode),]) #9
# nrow(inffRAclr[is.na(inffRAclr$infant_ffq_feeding_mode),]) #120
# 
# inffRAclr_cor_delmode <- inffRAclr[!is.na(inffRAclr$infant_birth_delivery_mode),]
# 
# ## create Aitchison distance matrix and add metadata
# inffait <- vegdist(inffRAclr_cor_delmode[,grep("g__", colnames(inffRAclr_cor_delmode))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
# inffaitm <- cmdscale(inffait, k=3) #create distance matrix for first 3 components
# inffaitd <- data.frame(inffaitm)
# 
# inffdf <- merge(inffRAclr_cor_delmode, inffaitd, by="row.names")
# inffdf <- inffdf[,-1]
# rownames(inffdf) <- inffdf$seq_16S_sample_ID
# 
# inff.beta.cmd <- cmdscale(as.matrix(inffait), k=2, eig=T)
# inffPCoA1 <- inff.beta.cmd$eig[1]/sum(inff.beta.cmd$eig)*100
# inffPCoA2 <- inff.beta.cmd$eig[2]/sum(inff.beta.cmd$eig)*100
# 
# ## run adonis
# inff_a2 <- adonis2(inffait ~ inffdf$DNA_isolation_method + inffdf$seq_16S_n_reads_clean + inffdf$infant_ID + inffdf$time_point + inffdf$infant_birth_delivery_mode, by="margin") #default: 999 permutations
# # Permutation test for adonis under reduced model
# # Marginal effects of terms
# # Permutation: free
# # Number of permutations: 999
# # 
# # adonis2(formula = inffait ~ inffdf$DNA_isolation_method + inffdf$seq_16S_n_reads_clean + inffdf$infant_ID + inffdf$time_point + inffdf$infant_birth_delivery_mode, by = "margin")
# # Df SumOfSqs      R2      F Pr(>F)    
# # inffdf$DNA_isolation_method         1      371 0.00212 1.7414  0.017 *  
# # inffdf$seq_16S_n_reads_clean        1      259 0.00149 1.2189  0.165    
# # inffdf$infant_ID                  260   104555 0.59910 1.8900  0.001 ***
# # inffdf$time_point                   3     4560 0.02613 7.1444  0.001 ***
# # inffdf$infant_birth_delivery_mode   0        0 0.00000    Inf           
# # Residual                          276    58724 0.33649                  
# # Total                             542   174519 1.00000                  
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# 
# inffRAclr_cor_delmode_feeding <- inffRAclr[!is.na(inffRAclr$infant_birth_delivery_mode) & !is.na(inffRAclr$infant_ffq_feeding_mode),]
# 
# ## create Aitchison distance matrix and add metadata
# inffait <- vegdist(inffRAclr_cor_delmode_feeding[,grep("g__", colnames(inffRAclr_cor_delmode_feeding))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
# inffaitm <- cmdscale(inffait, k=3) #create distance matrix for first 3 components
# inffaitd <- data.frame(inffaitm)
# 
# inffdf <- merge(inffRAclr_cor_delmode_feeding, inffaitd, by="row.names")
# inffdf <- inffdf[,-1]
# rownames(inffdf) <- inffdf$seq_16S_sample_ID
# 
# inff.beta.cmd <- cmdscale(as.matrix(inffait), k=2, eig=T)
# inffPCoA1 <- inff.beta.cmd$eig[1]/sum(inff.beta.cmd$eig)*100
# inffPCoA2 <- inff.beta.cmd$eig[2]/sum(inff.beta.cmd$eig)*100
# 
# ## run adonis
# inff_a2 <- adonis2(inffait ~ inffdf$DNA_isolation_method + inffdf$seq_16S_n_reads_clean + inffdf$infant_ID + inffdf$time_point + inffdf$infant_birth_delivery_mode + inffdf$infant_ffq_feeding_mode, by="margin") #default: 999 permutations
# # Permutation test for adonis under reduced model
# # Marginal effects of terms
# # Permutation: free
# # Number of permutations: 999
# # 
# # adonis2(formula = inffait ~ inffdf$DNA_isolation_method + inffdf$seq_16S_n_reads_clean + inffdf$infant_ID + inffdf$time_point + inffdf$infant_birth_delivery_mode + inffdf$infant_ffq_feeding_mode, by = "margin")
# # Df SumOfSqs      R2      F Pr(>F)    
# # inffdf$DNA_isolation_method         1      284 0.00205 1.3423  0.115    
# # inffdf$seq_16S_n_reads_clean        1      210 0.00152 0.9927  0.459    
# # inffdf$infant_ID                  219    86224 0.62212 1.8579  0.001 ***
# # inffdf$time_point                   3     3349 0.02416 5.2677  0.001 ***
# # inffdf$infant_birth_delivery_mode   0        0 0.00000   -Inf           
# # inffdf$infant_ffq_feeding_mode      1      180 0.00130 0.8504  0.675    
# # Residual                          201    42596 0.30733                  
# # Total                             427   138598 1.00000                  
# # ---
# #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


### ===== 16.2. INFANT FAECES PLOT PCOA BY TIME POINT (FACTOR) ===== ###

inff_pcoa_tp <- ggplot(inffdf, aes(x=X1, y=X2, color=time_point))+
  geom_point(alpha=0.5)+
  stat_ellipse(aes(group=time_point))+
  labs(x=paste0("PC1 (", signif(inffPCoA1,5), "%)"), y=paste0("PC2 (", signif(inffPCoA2,5), "%)"), color="Time point",
       title=paste0("Infant faeces, n=", nrow(inffdf), "\n", "Beta diversity by time point (factor)", "\n",
                    "Adonis with correction: R2=", signif(inff_R2_tp,3), " and p=", signif(inff_p_tp,3)))+
  scale_color_manual(values = c("1_month" = "lightgreen",
                                "2_months" = "limegreen",
                                "3_months" = "#009E73",
                                "6_months" = "darkolivegreen"))+
  theme_bw()
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/beta_diversity/231212_infant_faeces_PCoA_PC1_PC2_by_time_point.pdf"), inff_pcoa_tp, dpi=500, width=20, height=15, units="cm")


### ===== 16.3. INFANT FAECES PLOT PCOA BY INDIVIDUAL ID ===== ###

inff_pcoa_ID <- ggplot(inffdf, aes(x=X1, y=X2, color=infant_ID))+
  geom_point(alpha=0.5)+
  labs(x=paste0("PC1 (", signif(inffPCoA1,5), "%)"), y=paste0("PC2 (", signif(inffPCoA2,5), "%)"),
       title=paste0("Infant faeces, n=", nrow(inffdf), "\n", "Beta diversity by infant ID", "\n",
                    "Adonis with correction: R2=", signif(inff_R2_ID,3), " and p=", signif(inff_p_ID,3)))+
  theme_bw()+
  theme(legend.position = "none")
ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/beta_diversity/231212_infant_faeces_PCoA_PC1_PC2_by_individual_ID.pdf"), inff_pcoa_ID, dpi=500, width=20, height=15, units="cm")


##### =========================== 17. PREPARE DATA FRAME FOR COMPARISON OF AITCHISON DISTANCES THAT SHOULD INCLUDE MATERNAL DATA FOR SECOND TWIN INFANT =========================== #####

## For analysis of the distances between maternal faeces - milk:
# Work with data frame 'dfRAclr', which contains non-duplicated maternal data, to ensure that mothers are only included once.
# This data frame is ready to use for analysis.

## For analysis of the distances between milk - infant faeces and maternal faeces and infant faeces:
# Work with data frame 'df', which contains duplicated maternal data, to ensure that for mothers of twins, data for both infants is included.
# -> This data frame needs to be prepared for analysis.


### ===== 17.1. SELECT NEW DATA FRAME ===== ###

## save data frame under new name before adapting it
adf <- df

## check sample type distribution
table(adf$sample_origin_type, useNA="ifany")
# 737 human milk
# 552 infant faeces
# 172 maternal faeces

## drop empty levels
adf <- droplevels(adf)


### ===== 17.2. CALCULATE RELATIVE ABUNDANCES (BASED ON GENERA) ===== ###

## save only columns with absolute bacterial abundances
rownames(adf) <- paste0(adf$sample_origin_type, "__", adf$time_point, "__", adf$seq_16S_sample_ID)
adfbact <- adf[,531:ncol(adf)] #1461x624

## calculate relative bacterial abundances
adfbact_relab <- (adfbact/rowSums(adfbact))
table(rowSums(adfbact_relab), useNA="ifany") #all samples have 1

## merge back with metadata
adfRA <- merge(adf[,1:530], adfbact_relab, by="row.names")
rownames(adfRA) <- adfRA$Row.names
adfRA <- adfRA[,-1]


### ===== 17.3. CLR-TRANSFORM RELATIVE ABUNDANCES ===== ###

adfRAclr <- adfRA[,-1154] #exclude column 'Other'
adfRAclr[,531:ncol(adfRAclr)] <- as.data.frame(do_clr_externalWeighting(adfRAclr[,531:ncol(adfRAclr)], adfRAclr[,531:ncol(adfRAclr)]))

## update colnames to show transformation
colnames(adfRAclr)[531:ncol(adfRAclr)] <- c(paste0(colnames(adfRAclr)[531:ncol(adfRAclr)], "_clr"))


### ===== 17.4. ADD COLUMN COMBINING SAMPLE TYPE AND TIME POINT ===== ###

## Add column combining sample type and time point
adfRAclr$sample_origin_type_time_point <- NA

adfRAclr[adfRAclr$sample_origin_type=="mother_human_milk" & adfRAclr$time_point=="1_month", "sample_origin_type_time_point"] <- "human_milk_1_month"
adfRAclr[adfRAclr$sample_origin_type=="mother_human_milk" & adfRAclr$time_point=="2_months", "sample_origin_type_time_point"] <- "human_milk_2_months"
adfRAclr[adfRAclr$sample_origin_type=="mother_human_milk" & adfRAclr$time_point=="3_months", "sample_origin_type_time_point"] <- "human_milk_3_months"
adfRAclr[adfRAclr$sample_origin_type=="mother_human_milk" & adfRAclr$time_point=="6_months", "sample_origin_type_time_point"] <- "human_milk_6_months"

adfRAclr[adfRAclr$sample_origin_type=="mother_faeces" & adfRAclr$time_point=="3_months", "sample_origin_type_time_point"] <- "maternal_faeces_3_months"

adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="1_month", "sample_origin_type_time_point"] <- "infant_faeces_1_month"
adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="2_months", "sample_origin_type_time_point"] <- "infant_faeces_2_months"
adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="3_months", "sample_origin_type_time_point"] <- "infant_faeces_3_months"
adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="6_months", "sample_origin_type_time_point"] <- "infant_faeces_6_months"

table(adfRAclr$sample_origin_type_time_point, useNA="ifany")


##### =========================== 18. HUMAN MILK - INFANT FAECES AITCHISON DISTANCES BY TIME =========================== #####

## Aim: Extract Aitchison distance for corresponding human milk - infant faecal samples, per time point
##      Investigate if infant faeces becomes more dissimilar to milk with increasing infant age

## Note: Work with data frame 'adfRAclr' to ensure that all twin infants + their mother's data are included!


### ===== 18.1. CREATE DATA FRAME WITH ONLY CORRESPONDING HUMAN MILK AND INFANT FAECES DATA ===== ###

## select human milk and infant faecal samples into separate data frames (per sample type and per time)
dfclrmilkM1 <- adfRAclr[adfRAclr$sample_origin_type=="mother_human_milk" & adfRAclr$time_point=="1_month",]
dfclrinffM1 <- adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="1_month",]

dfclrmilkM2 <- adfRAclr[adfRAclr$sample_origin_type=="mother_human_milk" & adfRAclr$time_point=="2_months",]
dfclrinffM2 <- adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="2_months",]

dfclrmilkM3 <- adfRAclr[adfRAclr$sample_origin_type=="mother_human_milk" & adfRAclr$time_point=="3_months",]
dfclrinffM3 <- adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="3_months",]

dfclrmilkM6 <- adfRAclr[adfRAclr$sample_origin_type=="mother_human_milk" & adfRAclr$time_point=="6_months",]
dfclrinffM6 <- adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="6_months",]


## add a column (identifier) combining family ID, NEXT participation number, infant number and time point
dfclrmilkM1$identifier <- c(paste0(dfclrmilkM1$family_ID, "__", dfclrmilkM1$NEXT_participation_number, "__", dfclrmilkM1$infant_number, "__", dfclrmilkM1$time_point)) #303
dfclrinffM1$identifier <- c(paste0(dfclrinffM1$family_ID, "__", dfclrinffM1$NEXT_participation_number, "__", dfclrinffM1$infant_number, "__", dfclrinffM1$time_point)) #229

dfclrmilkM2$identifier <- c(paste0(dfclrmilkM2$family_ID, "__", dfclrmilkM2$NEXT_participation_number, "__", dfclrmilkM2$infant_number, "__", dfclrmilkM2$time_point)) #52
dfclrinffM2$identifier <- c(paste0(dfclrinffM2$family_ID, "__", dfclrinffM2$NEXT_participation_number, "__", dfclrinffM2$infant_number, "__", dfclrinffM2$time_point)) #38

dfclrmilkM3$identifier <- c(paste0(dfclrmilkM3$family_ID, "__", dfclrmilkM3$NEXT_participation_number, "__", dfclrmilkM3$infant_number, "__", dfclrmilkM3$time_point)) #271
dfclrinffM3$identifier <- c(paste0(dfclrinffM3$family_ID, "__", dfclrinffM3$NEXT_participation_number, "__", dfclrinffM3$infant_number, "__", dfclrinffM3$time_point)) #195

dfclrmilkM6$identifier <- c(paste0(dfclrmilkM6$family_ID, "__", dfclrmilkM6$NEXT_participation_number, "__", dfclrmilkM6$infant_number, "__", dfclrmilkM6$time_point)) #111
dfclrinffM6$identifier <- c(paste0(dfclrinffM6$family_ID, "__", dfclrinffM6$NEXT_participation_number, "__", dfclrinffM6$infant_number, "__", dfclrinffM6$time_point)) #90


## extract identifier from each data frame
dfclrmilkM1fams <- dfclrmilkM1$identifier #303
dfclrinffM1fams <- dfclrinffM1$identifier #229

dfclrmilkM2fams <- dfclrmilkM2$identifier #52
dfclrinffM2fams <- dfclrinffM2$identifier #38

dfclrmilkM3fams <- dfclrmilkM3$identifier #271
dfclrinffM3fams <- dfclrinffM3$identifier #195

dfclrmilkM6fams <- dfclrmilkM6$identifier #111
dfclrinffM6fams <- dfclrinffM6$identifier #90


## for each time point, check which families have both sample types
dfclrmilkinffM1fams <- intersect(dfclrmilkM1fams, dfclrinffM1fams) #210
dfclrmilkinffM2fams <- intersect(dfclrmilkM2fams, dfclrinffM2fams) #34
dfclrmilkinffM3fams <- intersect(dfclrmilkM3fams, dfclrinffM3fams) #190
dfclrmilkinffM6fams <- intersect(dfclrmilkM6fams, dfclrinffM6fams) #75


## keep only data from families which have both sample types
dfclrmilkM1c <- dfclrmilkM1[dfclrmilkM1$identifier %in% dfclrmilkinffM1fams,]
dfclrinffM1c <- dfclrinffM1[dfclrinffM1$identifier %in% dfclrmilkinffM1fams,]

dfclrmilkM2c <- dfclrmilkM2[dfclrmilkM2$identifier %in% dfclrmilkinffM2fams,]
dfclrinffM2c <- dfclrinffM2[dfclrinffM2$identifier %in% dfclrmilkinffM2fams,]

dfclrmilkM3c <- dfclrmilkM3[dfclrmilkM3$identifier %in% dfclrmilkinffM3fams,]
dfclrinffM3c <- dfclrinffM3[dfclrinffM3$identifier %in% dfclrmilkinffM3fams,]

dfclrmilkM6c <- dfclrmilkM6[dfclrmilkM6$identifier %in% dfclrmilkinffM6fams,]
dfclrinffM6c <- dfclrinffM6[dfclrinffM6$identifier %in% dfclrmilkinffM6fams,]


## merge data frames into 1 data frame
dfclrmilkinffc <- as.data.frame(rbind(dfclrmilkM1c,dfclrinffM1c,dfclrmilkM2c,dfclrinffM2c,dfclrmilkM3c,dfclrinffM3c,dfclrmilkM6c,dfclrinffM6c))
dfclrmilkinffc <- droplevels(dfclrmilkinffc)


## add identifiers to rownames
rownames(dfclrmilkinffc) <- c(paste0(dfclrmilkinffc$identifier, "___", dfclrmilkinffc$sample_origin_type, "__", dfclrmilkinffc$seq_16S_sample_ID))


### ===== 18.2. CALCULATE AITCHISON DISTANCES ===== ###

## Calculate Aitchison distance matrix
milkinffait <- vegdist(dfclrmilkinffc[,grep("g__", colnames(dfclrmilkinffc))], method="euclidian")

# save as data frame with rownames and colnames allowing to identify samples
milkinffaitd <- as.data.frame(as.matrix(milkinffait))

# only keep human milk samples in the rows and only keep infant faecal samples in the columns
milkinffaitd2 <- milkinffaitd[grep("mother_human_milk", rownames(milkinffaitd)), grep("infant_faeces", colnames(milkinffaitd))]


### ===== 18.3. FUNCTION TO EXTRACT CORRESPONDING AITCHISON DISTANCES ===== ###

## function to extract Aitchison distances of corresponding samples from the same family
extract.corresponding.distances <- function(inputdata){
  
  identifier <- c()
  aitdist <- c()
  
  for (i in rownames(inputdata)){
    print(paste0("Starting analysis for ", i))
    
    # extract family ID from rowname and save it to vector
    iden <- gsub("___.*", "", i)
    identifier <- c(identifier, iden)
    
    # extract distance from corresponding samples and save it to vector
    aitdist <- c(aitdist, inputdata[i, grep(iden, colnames(inputdata))])
  }
  
  # combine results in data frame
  print(paste0("Combining results in data frame"))
  
  res <- data.frame(identifier=identifier, aitchison_distance=aitdist)
  
  return(res)
}

## extract distances
milk_inff_distances <- extract.corresponding.distances(inputdata=milkinffaitd2)


## separate data by time points
milk_inff_distances_M1 <- milk_inff_distances[grep("1_month", milk_inff_distances$identifier),]
milk_inff_distances_M2 <- milk_inff_distances[grep("2_months", milk_inff_distances$identifier),]
milk_inff_distances_M3 <- milk_inff_distances[grep("3_months", milk_inff_distances$identifier),]
milk_inff_distances_M6 <- milk_inff_distances[grep("6_months", milk_inff_distances$identifier),]

## add time point column
milk_inff_distances_M1$time_point <- "1_month"
milk_inff_distances_M2$time_point <- "2_months"
milk_inff_distances_M3$time_point <- "3_months"
milk_inff_distances_M6$time_point <- "6_months"

## merge data frames back into 1 data frame
milk_inff_distances_all <- as.data.frame(rbind(milk_inff_distances_M1, milk_inff_distances_M2, milk_inff_distances_M3, milk_inff_distances_M6))

## fix data structure
for (i in c(1,3)){milk_inff_distances_all[,i] <- as.factor(as.character(milk_inff_distances_all[,i]))}


### ===== 18.4. PLOT HUMAN MILK - INFANT FAECES AITCHISON DISTANCES PER TIME POINT ===== ###

milkinffdisplot <- ggplot(milk_inff_distances_all, aes(x=time_point, y=aitchison_distance))+
  geom_jitter(alpha=0.5, color="#7DCDC8")+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  labs(x="Month(s) postpartum", y="Aitchison distance", title="Human milk - infant faeces")+
  theme(axis.text.x = element_text(angle=90),
        panel.grid = element_blank())+
  scale_y_continuous(limits = c(0,50), breaks = c(0,25,50))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/milk_and_infant_faeces/240222_Aitchison_distances_betw_corresp_milk_and_infant_faeces_by_time_point.pdf", milkinffdisplot, dpi=500, width=15, height=15, units="cm")


### ===== 18.5. HISTOGRAMS FOR AITCHISON DISTANCES ===== ###

create.histograms <- function(inputdata, timepoint){
  histrichness <- ggplot(inputdata, aes(x=aitchison_distance)) + geom_histogram()
  ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/milk_and_infant_faeces/checks/240222_Aitchison_distances_betw_corresp_milk_and_infant_faeces_histogram_", timepoint, "_n", nrow(inputdata), ".pdf"), histrichness, dpi=500, width=20, height=20, units="cm", useDingbats=F)
}

create.histograms(milk_inff_distances_all[milk_inff_distances_all$time_point=="1_month",], timepoint="M1")
create.histograms(milk_inff_distances_all[milk_inff_distances_all$time_point=="2_months",], timepoint="M2")
create.histograms(milk_inff_distances_all[milk_inff_distances_all$time_point=="3_months",], timepoint="M3")
create.histograms(milk_inff_distances_all[milk_inff_distances_all$time_point=="6_months",], timepoint="M6")

## -> At each time point, the Aitchison distances are normally distributed.


### ===== 18.6. COMPARE HUMAN MILK - INFANT FAECES AITCHISON DISTANCES BY TIME ===== ###

## add numeric time point column to data frame
milk_inff_distances_all$time_point_numeric <- milk_inff_distances_all$time_point
levels(milk_inff_distances_all$time_point_numeric) <- c(1,2,3,6)
milk_inff_distances_all$time_point_numeric <- as.numeric(as.character(milk_inff_distances_all$time_point_numeric))

milkinff_m1 <- lm(aitchison_distance ~ time_point_numeric, data=milk_inff_distances_all)
summary(milkinff_m1)
# Call:
#   lm(formula = aitchison_distance ~ time_point_numeric, data = milk_inff_distances_all)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -14.7295  -2.9775  -0.0357   2.9996  13.8194 
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         26.1615     0.3692  70.868  < 2e-16 ***
# time_point_numeric   1.0121     0.1208   8.381 5.22e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 4.592 on 507 degrees of freedom
# Multiple R-squared:  0.1217,    Adjusted R-squared:  0.1199 
# F-statistic: 70.23 on 1 and 507 DF,  p-value: 5.219e-16


##### =========================== 19. INFANT FAECES - MATERNAL FAECES AITCHISON DISTANCES BY TIME =========================== #####

## Aim: Extract Aitchison distance for corresponding infant (M1, M2, M3, M6) - maternal (M3) faecal samples, for each time point combination (M1-M3, M2-M3, M3-M3, M6-M3)
##      Investigate if infant faeces becomes more similar to maternal faeces within the first 6 months of life

## Note: Work with data frame 'adfRAclr' to ensure that all twin infants + their mother's data are included!


### ===== 19.1. CREATE DATA FRAME WITH ONLY CORRESPONDING HUMAN MILK AND INFANT FAECES DATA ===== ###

## select infant and maternal faecal samples into separate data frames (per sample type and per time)
df2clrinffM1 <- adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="1_month",]
df2clrinffM2 <- adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="2_months",]
df2clrinffM3 <- adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="3_months",]
df2clrinffM6 <- adfRAclr[adfRAclr$sample_origin_type=="infant_faeces" & adfRAclr$time_point=="6_months",]

df2clrmatfM3 <- adfRAclr[adfRAclr$sample_origin_type=="mother_faeces" & adfRAclr$time_point=="3_months",]


## add a column (identifier) combining family ID, NEXT participation number, infant number and time point
df2clrinffM1$identifier <- c(paste0(df2clrinffM1$family_ID, "__", df2clrinffM1$NEXT_participation_number, "__", df2clrinffM1$infant_number)) #229
df2clrinffM2$identifier <- c(paste0(df2clrinffM2$family_ID, "__", df2clrinffM2$NEXT_participation_number, "__", df2clrinffM2$infant_number)) #38
df2clrinffM3$identifier <- c(paste0(df2clrinffM3$family_ID, "__", df2clrinffM3$NEXT_participation_number, "__", df2clrinffM3$infant_number)) #195
df2clrinffM6$identifier <- c(paste0(df2clrinffM6$family_ID, "__", df2clrinffM6$NEXT_participation_number, "__", df2clrinffM6$infant_number)) #90

df2clrmatfM3$identifier <- c(paste0(df2clrmatfM3$family_ID, "__", df2clrmatfM3$NEXT_participation_number, "__", df2clrmatfM3$infant_number)) #172


## extract identifier from each data frame
df2clrinffM1fams <- df2clrinffM1$identifier #229
df2clrinffM2fams <- df2clrinffM2$identifier #38
df2clrinffM3fams <- df2clrinffM3$identifier #195
df2clrinffM6fams <- df2clrinffM6$identifier #90

df2clrmatfM3fams <- df2clrmatfM3$identifier #172


## for each time point combination, check which families have both sample types
df2clrinffM1matfM6fams <- intersect(df2clrinffM1fams, df2clrmatfM3fams) #153
df2clrinffM2matfM6fams <- intersect(df2clrinffM2fams, df2clrmatfM3fams) #2
df2clrinffM3matfM6fams <- intersect(df2clrinffM3fams, df2clrmatfM3fams) #158
df2clrinffM6matfM6fams <- intersect(df2clrinffM6fams, df2clrmatfM3fams) #60


## keep only data from families which have both sample types
df2clrinffM1c <- df2clrinffM1[df2clrinffM1$identifier %in% df2clrinffM1matfM6fams,]
df2clrinffM2c <- df2clrinffM2[df2clrinffM2$identifier %in% df2clrinffM2matfM6fams,]
df2clrinffM3c <- df2clrinffM3[df2clrinffM3$identifier %in% df2clrinffM3matfM6fams,]
df2clrinffM6c <- df2clrinffM6[df2clrinffM6$identifier %in% df2clrinffM6matfM6fams,]

df2clrmatfM3forinffM1c <- df2clrmatfM3[df2clrmatfM3$identifier %in% df2clrinffM1matfM6fams,]
df2clrmatfM3forinffM2c <- df2clrmatfM3[df2clrmatfM3$identifier %in% df2clrinffM2matfM6fams,]
df2clrmatfM3forinffM3c <- df2clrmatfM3[df2clrmatfM3$identifier %in% df2clrinffM3matfM6fams,]
df2clrmatfM3forinffM6c <- df2clrmatfM3[df2clrmatfM3$identifier %in% df2clrinffM6matfM6fams,]


## merge data frames into 1 data frame for each to-be-investigated time point combination
# note that this creates 4 data frames, to comapre (1) inff M1 - matf M3, (2) inff M2 - matf M3, (3) inff M3 - matf M3 and (4) inff M6 - matf M3
dfclrinffM1matfM3c <- as.data.frame(rbind(df2clrinffM1c,df2clrmatfM3forinffM1c))
dfclrinffM1matfM3c <- droplevels(dfclrinffM1matfM3c)

dfclrinffM2matfM3c <- as.data.frame(rbind(df2clrinffM2c,df2clrmatfM3forinffM2c))
dfclrinffM2matfM3c <- droplevels(dfclrinffM2matfM3c)

dfclrinffM3matfM3c <- as.data.frame(rbind(df2clrinffM3c,df2clrmatfM3forinffM3c))
dfclrinffM3matfM3c <- droplevels(dfclrinffM3matfM3c)

dfclrinffM6matfM3c <- as.data.frame(rbind(df2clrinffM6c,df2clrmatfM3forinffM6c))
dfclrinffM6matfM3c <- droplevels(dfclrinffM6matfM3c)


## add identifiers to rownames
rownames(dfclrinffM1matfM3c) <- c(paste0(dfclrinffM1matfM3c$identifier, "___", dfclrinffM1matfM3c$sample_origin_type, "__", dfclrinffM1matfM3c$seq_16S_sample_ID))
rownames(dfclrinffM2matfM3c) <- c(paste0(dfclrinffM2matfM3c$identifier, "___", dfclrinffM2matfM3c$sample_origin_type, "__", dfclrinffM2matfM3c$seq_16S_sample_ID))
rownames(dfclrinffM3matfM3c) <- c(paste0(dfclrinffM3matfM3c$identifier, "___", dfclrinffM3matfM3c$sample_origin_type, "__", dfclrinffM3matfM3c$seq_16S_sample_ID))
rownames(dfclrinffM6matfM3c) <- c(paste0(dfclrinffM6matfM3c$identifier, "___", dfclrinffM6matfM3c$sample_origin_type, "__", dfclrinffM6matfM3c$seq_16S_sample_ID))


### ===== 19.2. CALCULATE AITCHISON DISTANCES ===== ###

### INFANT FAECES M1 - MATERNAL FAECES M3

## Calculate Aitchison distance matrix
inffM1matfM3ait <- vegdist(dfclrinffM1matfM3c[,grep("g__", colnames(dfclrinffM1matfM3c))], method="euclidian")

# save as data frame with rownames and colnames allowing to identify samples
inffM1matfM3aitd <- as.data.frame(as.matrix(inffM1matfM3ait))

# only keep human milk samples in the rows and only keep infant faecal samples in the columns
inffM1matfM3aitd2 <- inffM1matfM3aitd[grep("mother_faeces", rownames(inffM1matfM3aitd)), grep("infant_faeces", colnames(inffM1matfM3aitd))]


### INFANT FAECES M2 - MATERNAL FAECES M3

## Calculate Aitchison distance matrix
inffM2matfM3ait <- vegdist(dfclrinffM2matfM3c[,grep("g__", colnames(dfclrinffM2matfM3c))], method="euclidian")

# save as data frame with rownames and colnames allowing to identify samples
inffM2matfM3aitd <- as.data.frame(as.matrix(inffM2matfM3ait))

# only keep human milk samples in the rows and only keep infant faecal samples in the columns
inffM2matfM3aitd2 <- inffM2matfM3aitd[grep("mother_faeces", rownames(inffM2matfM3aitd)), grep("infant_faeces", colnames(inffM2matfM3aitd))]


### INFANT FAECES M3 - MATERNAL FAECES M3

## Calculate Aitchison distance matrix
inffM3matfM3ait <- vegdist(dfclrinffM3matfM3c[,grep("g__", colnames(dfclrinffM3matfM3c))], method="euclidian")

# save as data frame with rownames and colnames allowing to identify samples
inffM3matfM3aitd <- as.data.frame(as.matrix(inffM3matfM3ait))

# only keep human milk samples in the rows and only keep infant faecal samples in the columns
inffM3matfM3aitd2 <- inffM3matfM3aitd[grep("mother_faeces", rownames(inffM3matfM3aitd)), grep("infant_faeces", colnames(inffM3matfM3aitd))]


### INFANT FAECES M6 - MATERNAL FAECES M3

## Calculate Aitchison distance matrix
inffM6matfM3ait <- vegdist(dfclrinffM6matfM3c[,grep("g__", colnames(dfclrinffM6matfM3c))], method="euclidian")

# save as data frame with rownames and colnames allowing to identify samples
inffM6matfM3aitd <- as.data.frame(as.matrix(inffM6matfM3ait))

# only keep human milk samples in the rows and only keep infant faecal samples in the columns
inffM6matfM3aitd2 <- inffM6matfM3aitd[grep("mother_faeces", rownames(inffM6matfM3aitd)), grep("infant_faeces", colnames(inffM6matfM3aitd))]


### ===== 19.3. FUNCTION TO EXTRACT CORRESPONDING AITCHISON DISTANCES ===== ###

## function to extract Aitchison distances of corresponding samples from the same family
extract.corresponding.distances <- function(inputdata){
  
  identifier <- c()
  aitdist <- c()
  
  for (i in rownames(inputdata)){
    print(paste0("Starting analysis for ", i))
    
    # extract family ID from rowname and save it to vector
    iden <- gsub("___.*", "", i)
    identifier <- c(identifier, iden)
    
    # extract distance from corresponding samples and save it to vector
    aitdist <- c(aitdist, inputdata[i, grep(iden, colnames(inputdata))])
  }
  
  # combine results in data frame
  print(paste0("Combining results in data frame"))
  
  res <- data.frame(identifier=identifier, aitchison_distance=aitdist)
  
  return(res)
}

## extract distances
inffM1_matfM3_distances <- extract.corresponding.distances(inputdata=inffM1matfM3aitd2)
inffM2_matfM3_distances <- extract.corresponding.distances(inputdata=inffM2matfM3aitd2)
inffM3_matfM3_distances <- extract.corresponding.distances(inputdata=inffM3matfM3aitd2)
inffM6_matfM3_distances <- extract.corresponding.distances(inputdata=inffM6matfM3aitd2)

## add time point column
inffM1_matfM3_distances$time_point <- "1_month"
inffM2_matfM3_distances$time_point <- "2_months"
inffM3_matfM3_distances$time_point <- "3_months"
inffM6_matfM3_distances$time_point <- "6_months"

## merge data frames back into 1 data frame
inff_matf_distances_all <- as.data.frame(rbind(inffM1_matfM3_distances, inffM2_matfM3_distances, inffM3_matfM3_distances, inffM6_matfM3_distances))

## fix data structure
for (i in c(1,3)){inff_matf_distances_all[,i] <- as.factor(as.character(inff_matf_distances_all[,i]))}


### ===== 19.4. PLOT INFANT FAECES - MATERNAL FAECES AITCHISON DISTANCES PER TIME POINT ===== ###

inffmatfdisplot <- ggplot(inff_matf_distances_all, aes(x=time_point, y=aitchison_distance))+
  geom_jitter(alpha=0.5, color="#5BB4E5")+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  labs(x="Month(s) postpartum", y="Aitchison distance", title="Infant faeces - Maternal faeces")+
  theme(axis.text.x = element_text(angle=90),
        panel.grid = element_blank())+
  scale_y_continuous(limits = c(0,50), breaks = c(0,25,50))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/maternal_and_infant_faeces/240222_Aitchison_distances_betw_corresp_infant_and_maternal_faeces_by_time_point.pdf", inffmatfdisplot, dpi=500, width=15, height=15, units="cm")


### ===== 19.5. HISTOGRAMS FOR AITCHISON DISTANCES ===== ###

create.histograms <- function(inputdata, timepoint){
  histrichness <- ggplot(inputdata, aes(x=aitchison_distance)) + geom_histogram()
  ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/maternal_and_infant_faeces/checks/240222_Aitchison_distances_betw_corresp_infant_and_maternal_faeces_histogram_", timepoint, "_n", nrow(inputdata), ".pdf"), histrichness, dpi=500, width=20, height=20, units="cm", useDingbats=F)
}

create.histograms(inff_matf_distances_all[inff_matf_distances_all$time_point=="1_month",], timepoint="M1")
create.histograms(inff_matf_distances_all[inff_matf_distances_all$time_point=="2_months",], timepoint="M2")
create.histograms(inff_matf_distances_all[inff_matf_distances_all$time_point=="3_months",], timepoint="M3")
create.histograms(inff_matf_distances_all[inff_matf_distances_all$time_point=="6_months",], timepoint="M6")

## -> At each time point, the Aitchison distances are normally distributed.


### ===== 19.6. COMPARE INFANT FAECES - MATERNAL FAECES AITCHISON DISTANCES BY TIME ===== ###

## add numeric time point column to data frame
inff_matf_distances_all$time_point_numeric <- inff_matf_distances_all$time_point
levels(inff_matf_distances_all$time_point_numeric) <- c(1,2,3,6)
inff_matf_distances_all$time_point_numeric <- as.numeric(as.character(inff_matf_distances_all$time_point_numeric))

inffmatf_m1 <- lm(aitchison_distance ~ time_point_numeric, data=inff_matf_distances_all)
summary(inffmatf_m1)
# Call:
#   lm(formula = aitchison_distance ~ time_point_numeric, data = inff_matf_distances_all)
# 
# Residuals:
#   Min      1Q  Median      3Q     Max 
# -8.4369 -1.7257  0.2655  1.8586  7.0917 
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        40.64062    0.26584 152.875   <2e-16 ***
# time_point_numeric  0.15658    0.08392   1.866   0.0629 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.796 on 371 degrees of freedom
# Multiple R-squared:  0.009296,  Adjusted R-squared:  0.006625 
# F-statistic: 3.481 on 1 and 371 DF,  p-value: 0.06286


##### =========================== 20. MATERNAL FAECES - HUMAN MILK AITCHISON DISTANCES BY TIME =========================== #####

## Aim: Extract Aitchison distance for corresponding milk (M1, M2, M3, M6) - maternal (M3) faecal samples, for each time point combination (M1-M3, M2-M3, M3-M3, M6-M3)
##      Investigate how similar maternal faeces and human milk are within the first 6 months postpartum

## Note: Work with data frame 'dfRAclr' to NOT include any duplicated maternal data!


### ===== 20.1. CREATE DATA FRAME WITH ONLY CORRESPONDING HUMAN MILK AND MATERNAL FAECES DATA ===== ###

## select human milk and maternal faecal samples into separate data frames (per sample type and per time)
df3clrmilkM1 <- dfRAclr[dfRAclr$sample_origin_type=="mother_human_milk" & dfRAclr$time_point=="1_month",]
df3clrmilkM2 <- dfRAclr[dfRAclr$sample_origin_type=="mother_human_milk" & dfRAclr$time_point=="2_months",]
df3clrmilkM3 <- dfRAclr[dfRAclr$sample_origin_type=="mother_human_milk" & dfRAclr$time_point=="3_months",]
df3clrmilkM6 <- dfRAclr[dfRAclr$sample_origin_type=="mother_human_milk" & dfRAclr$time_point=="6_months",]

df3clrmatfM3 <- dfRAclr[dfRAclr$sample_origin_type=="mother_faeces" & dfRAclr$time_point=="3_months",]


## add a column (identifier) combining family ID, NEXT participation number, infant number and time point
df3clrmilkM1$identifier <- c(paste0(df3clrmilkM1$family_ID, "__", df3clrmilkM1$NEXT_participation_number, "__", df3clrmilkM1$infant_number)) #300
df3clrmilkM2$identifier <- c(paste0(df3clrmilkM2$family_ID, "__", df3clrmilkM2$NEXT_participation_number, "__", df3clrmilkM2$infant_number)) #50
df3clrmilkM3$identifier <- c(paste0(df3clrmilkM3$family_ID, "__", df3clrmilkM3$NEXT_participation_number, "__", df3clrmilkM3$infant_number)) #268
df3clrmilkM6$identifier <- c(paste0(df3clrmilkM6$family_ID, "__", df3clrmilkM6$NEXT_participation_number, "__", df3clrmilkM6$infant_number)) #107

df3clrmatfM3$identifier <- c(paste0(df3clrmatfM3$family_ID, "__", df3clrmatfM3$NEXT_participation_number, "__", df3clrmatfM3$infant_number)) #170


## extract identifier from each data frame
df3clrmilkM1fams <- df3clrmilkM1$identifier #300
df3clrmilkM2fams <- df3clrmilkM2$identifier #50
df3clrmilkM3fams <- df3clrmilkM3$identifier #268
df3clrmilkM6fams <- df3clrmilkM6$identifier #107

df3clrmatfM3fams <- df3clrmatfM3$identifier #170


## for each time point combination, check which families have both sample types
df3clrmilkM1matfM6fams <- intersect(df3clrmilkM1fams, df3clrmatfM3fams) #153
df3clrmilkM2matfM6fams <- intersect(df3clrmilkM2fams, df3clrmatfM3fams) #0
df3clrmilkM3matfM6fams <- intersect(df3clrmilkM3fams, df3clrmatfM3fams) #164
df3clrmilkM6matfM6fams <- intersect(df3clrmilkM6fams, df3clrmatfM3fams) #72


## keep only data from families which have both sample types
df3clrmilkM1c <- df3clrmilkM1[df3clrmilkM1$identifier %in% df3clrmilkM1matfM6fams,]
df3clrmilkM2c <- df3clrmilkM2[df3clrmilkM2$identifier %in% df3clrmilkM2matfM6fams,]
df3clrmilkM3c <- df3clrmilkM3[df3clrmilkM3$identifier %in% df3clrmilkM3matfM6fams,]
df3clrmilkM6c <- df3clrmilkM6[df3clrmilkM6$identifier %in% df3clrmilkM6matfM6fams,]

df3clrmatfM3formilkM1c <- df3clrmatfM3[df3clrmatfM3$identifier %in% df3clrmilkM1matfM6fams,]
df3clrmatfM3formilkM2c <- df3clrmatfM3[df3clrmatfM3$identifier %in% df3clrmilkM2matfM6fams,]
df3clrmatfM3formilkM3c <- df3clrmatfM3[df3clrmatfM3$identifier %in% df3clrmilkM3matfM6fams,]
df3clrmatfM3formilkM6c <- df3clrmatfM3[df3clrmatfM3$identifier %in% df3clrmilkM6matfM6fams,]


## merge data frames into 1 data frame for each to-be-investigated time point combination
# note that this creates 3 data frames, to compare (1) milk M1 - matf M3, (2) milk M3 - matf M3 and (3) milk M6 - matf M3
# there is no overlap in maternal faecal M3 and milk M2 samples, so milk M2 - matf M3 cannot be investigated
dfclrmilkM1matfM3c <- as.data.frame(rbind(df3clrmilkM1c,df3clrmatfM3formilkM1c))
dfclrmilkM1matfM3c <- droplevels(dfclrmilkM1matfM3c)

# dfclrmilkM2matfM3c <- as.data.frame(rbind(df3clrmilkM2c,df3clrmatfM3formilkM2c))
# dfclrmilkM2matfM3c <- droplevels(dfclrmilkM2matfM3c)

dfclrmilkM3matfM3c <- as.data.frame(rbind(df3clrmilkM3c,df3clrmatfM3formilkM3c))
dfclrmilkM3matfM3c <- droplevels(dfclrmilkM3matfM3c)

dfclrmilkM6matfM3c <- as.data.frame(rbind(df3clrmilkM6c,df3clrmatfM3formilkM6c))
dfclrmilkM6matfM3c <- droplevels(dfclrmilkM6matfM3c)


## add identifiers to rownames
rownames(dfclrmilkM1matfM3c) <- c(paste0(dfclrmilkM1matfM3c$identifier, "___", dfclrmilkM1matfM3c$sample_origin_type, "__", dfclrmilkM1matfM3c$seq_16S_sample_ID))
# rownames(dfclrmilkM2matfM3c) <- c(paste0(dfclrmilkM2matfM3c$identifier, "___", dfclrmilkM2matfM3c$sample_origin_type, "__", dfclrmilkM2matfM3c$seq_16S_sample_ID))
rownames(dfclrmilkM3matfM3c) <- c(paste0(dfclrmilkM3matfM3c$identifier, "___", dfclrmilkM3matfM3c$sample_origin_type, "__", dfclrmilkM3matfM3c$seq_16S_sample_ID))
rownames(dfclrmilkM6matfM3c) <- c(paste0(dfclrmilkM6matfM3c$identifier, "___", dfclrmilkM6matfM3c$sample_origin_type, "__", dfclrmilkM6matfM3c$seq_16S_sample_ID))


### ===== 20.2. CALCULATE AITCHISON DISTANCES ===== ###

### HUMAN MILK M1 - MATERNAL FAECES M3

## Calculate Aitchison distance matrix
milkM1matfM3ait <- vegdist(dfclrmilkM1matfM3c[,grep("g__", colnames(dfclrmilkM1matfM3c))], method="euclidian")

# save as data frame with rownames and colnames allowing to identify samples
milkM1matfM3aitd <- as.data.frame(as.matrix(milkM1matfM3ait))

# only keep human milk samples in the rows and only keep infant faecal samples in the columns
milkM1matfM3aitd2 <- milkM1matfM3aitd[grep("mother_faeces", rownames(milkM1matfM3aitd)), grep("mother_human_milk", colnames(milkM1matfM3aitd))]


# ### HUMAN MILK M2 - MATERNAL FAECES M3
# 
# ## Calculate Aitchison distance matrix
# milkM2matfM3ait <- vegdist(dfclrmilkM2matfM3c[,grep("g__", colnames(dfclrmilkM2matfM3c))], method="euclidian")
# 
# # save as data frame with rownames and colnames allowing to identify samples
# milkM2matfM3aitd <- as.data.frame(as.matrix(milkM2matfM3ait))
# 
# # only keep human milk samples in the rows and only keep infant faecal samples in the columns
# milkM2matfM3aitd2 <- milkM2matfM3aitd[grep("mother_faeces", rownames(milkM2matfM3aitd)), grep("mother_human_milk", colnames(milkM2matfM3aitd))]


### HUMAN MILK M3 - MATERNAL FAECES M3

## Calculate Aitchison distance matrix
milkM3matfM3ait <- vegdist(dfclrmilkM3matfM3c[,grep("g__", colnames(dfclrmilkM3matfM3c))], method="euclidian")

# save as data frame with rownames and colnames allowing to identify samples
milkM3matfM3aitd <- as.data.frame(as.matrix(milkM3matfM3ait))

# only keep human milk samples in the rows and only keep infant faecal samples in the columns
milkM3matfM3aitd2 <- milkM3matfM3aitd[grep("mother_faeces", rownames(milkM3matfM3aitd)), grep("mother_human_milk", colnames(milkM3matfM3aitd))]


### HUMAN MILK M6 - MATERNAL FAECES M3

## Calculate Aitchison distance matrix
milkM6matfM3ait <- vegdist(dfclrmilkM6matfM3c[,grep("g__", colnames(dfclrmilkM6matfM3c))], method="euclidian")

# save as data frame with rownames and colnames allowing to identify samples
milkM6matfM3aitd <- as.data.frame(as.matrix(milkM6matfM3ait))

# only keep human milk samples in the rows and only keep infant faecal samples in the columns
milkM6matfM3aitd2 <- milkM6matfM3aitd[grep("mother_faeces", rownames(milkM6matfM3aitd)), grep("mother_human_milk", colnames(milkM6matfM3aitd))]


### ===== 20.3. FUNCTION TO EXTRACT CORRESPONDING AITCHISON DISTANCES ===== ###

## function to extract Aitchison distances of corresponding samples from the same family
extract.corresponding.distances <- function(inputdata){
  
  identifier <- c()
  aitdist <- c()
  
  for (i in rownames(inputdata)){
    print(paste0("Starting analysis for ", i))
    
    # extract family ID from rowname and save it to vector
    iden <- gsub("___.*", "", i)
    identifier <- c(identifier, iden)
    
    # extract distance from corresponding samples and save it to vector
    aitdist <- c(aitdist, inputdata[i, grep(iden, colnames(inputdata))])
  }
  
  # combine results in data frame
  print(paste0("Combining results in data frame"))
  
  res <- data.frame(identifier=identifier, aitchison_distance=aitdist)
  
  return(res)
}

## extract distances
milkM1_matfM3_distances <- extract.corresponding.distances(inputdata=milkM1matfM3aitd2)
# milkM2_matfM3_distances <- extract.corresponding.distances(inputdata=milkM2matfM3aitd2)
milkM3_matfM3_distances <- extract.corresponding.distances(inputdata=milkM3matfM3aitd2)
milkM6_matfM3_distances <- extract.corresponding.distances(inputdata=milkM6matfM3aitd2)

## add time point column
milkM1_matfM3_distances$time_point <- "1_month"
# milkM2_matfM3_distances$time_point <- "2_months"
milkM3_matfM3_distances$time_point <- "3_months"
milkM6_matfM3_distances$time_point <- "6_months"

## merge data frames back into 1 data frame
milk_matf_distances_all <- as.data.frame(rbind(milkM1_matfM3_distances, milkM3_matfM3_distances, milkM6_matfM3_distances))

## fix data structure
for (i in c(1,3)){milk_matf_distances_all[,i] <- as.factor(as.character(milk_matf_distances_all[,i]))}


### ===== 20.4. PLOT MATERNAL FAECES - HUMAN MILK AITCHISON DISTANCES PER TIME POINT ===== ###

milkmatfdisplot <- ggplot(milk_matf_distances_all, aes(x=time_point, y=aitchison_distance))+
  geom_jitter(alpha=0.5, color="#5BB4E5")+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  labs(x="Month(s) postpartum", y="Aitchison distance", title="Maternal faeces - Human milk")+
  theme(axis.text.x = element_text(angle=90),
        panel.grid = element_blank())+
  scale_y_continuous(limits = c(0,60), breaks = c(0,20,40,60))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/maternal_faeces_and_milk/240222_Aitchison_distances_betw_corresp_maternal_faeces_and_milk_by_time_point.pdf", milkmatfdisplot, dpi=500, width=15, height=15, units="cm")


### ===== 20.5. HISTOGRAMS FOR AITCHISON DISTANCES ===== ###

create.histograms <- function(inputdata, timepoint){
  histrichness <- ggplot(inputdata, aes(x=aitchison_distance)) + geom_histogram()
  ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/maternal_faeces_and_milk/checks/240222_Aitchison_distances_betw_corresp_maternal_faeces_and_milk_histogram_", timepoint, "_n", nrow(inputdata), ".pdf"), histrichness, dpi=500, width=20, height=20, units="cm", useDingbats=F)
}

create.histograms(milk_matf_distances_all[milk_matf_distances_all$time_point=="1_month",], timepoint="M1")
# create.histograms(milk_matf_distances_all[milk_matf_distances_all$time_point=="2_months",], timepoint="M2")
create.histograms(milk_matf_distances_all[milk_matf_distances_all$time_point=="3_months",], timepoint="M3")
create.histograms(milk_matf_distances_all[milk_matf_distances_all$time_point=="6_months",], timepoint="M6")

## -> At each time point, the Aitchison distances are normally distributed.


### ===== 20.6. COMPARE MATERNAL FAECES - HUMAN MILK AITCHISON DISTANCES BY TIME ===== ###

## add numeric time point column to data frame
milk_matf_distances_all$time_point_numeric <- milk_matf_distances_all$time_point
levels(milk_matf_distances_all$time_point_numeric) <- c(1,3,6)
milk_matf_distances_all$time_point_numeric <- as.numeric(as.character(milk_matf_distances_all$time_point_numeric))

milkmatf_m1 <- lm(aitchison_distance ~ time_point_numeric, data=milk_matf_distances_all)
summary(milkmatf_m1)
# Call:
#   lm(formula = aitchison_distance ~ time_point_numeric, data = milk_matf_distances_all)
# 
# Residuals:
#   Min     1Q Median     3Q    Max 
# -8.282 -2.035 -0.196  1.883 11.054 
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        43.67139    0.28801  151.63   <2e-16 ***
# time_point_numeric  0.14247    0.08743    1.63    0.104    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 3.078 on 387 degrees of freedom
# Multiple R-squared:  0.006815,  Adjusted R-squared:  0.004248 
# F-statistic: 2.655 on 1 and 387 DF,  p-value: 0.104


##### =========================== 21. CREATE COMBINED PLOT OF AITCHISON DISTANCES BETWEEN CORRESPONDING SAMPLES BY SAMPLE TYPE COMPARISON AND TIME (FIGURE 4B) =========================== #####

### ===== 21.1. MERGE AND PREPARE DATA FRAMES ===== ###

## add column showing comparison to data frames
milk_inff_distances_all$comparison <- "human_milk_infant_faeces"
inff_matf_distances_all$comparison <- "infant_faeces_maternal_faeces"
milk_matf_distances_all$comparison <- "maternal_faeces_human_milk"

## combine data frames
distdf <- as.data.frame(rbind(milk_inff_distances_all, inff_matf_distances_all, milk_matf_distances_all))
distdf$comparison <- as.factor(as.character(distdf$comparison))

## resort levels in comparison column
distdf$comparison <- factor(distdf$comparison, levels = c("maternal_faeces_human_milk", "human_milk_infant_faeces", "infant_faeces_maternal_faeces"))


### ===== 21.2. COMBINED PLOT ===== ###

combdisplot <- ggplot(distdf, aes(x=time_point, y=aitchison_distance, color=comparison))+
  geom_jitter(alpha=0.5)+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  labs(x="Month(s) postpartum", y="Aitchison distance")+
  scale_x_discrete(labels = c("1","2","3","6"))+
  scale_y_continuous(limits = c(0,60), breaks = c(0,20,40,60))+
  scale_color_manual(values = c("maternal_faeces_human_milk" = c("turquoise"),
                                "human_milk_infant_faeces" = "#CC79A7", 
                                "infant_faeces_maternal_faeces" = "#5BB4E5"))+
  facet_grid(.~comparison)+
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.background = element_rect(color = "white", fill="white"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/all_sample_types/240222_Aitchison_distances_betw_corresp_samples_by_comparison_and_time_point_v2.pdf", combdisplot, dpi=500, width=20, height=12, units="cm")


### ===== 21.3. COMBINE STATISTICS RESULTS AND ADJUST FOR MULTIPLE TESTING ===== ###

## save statistics results from models in data frames
milkinff_sum1 <- as.data.frame(summary(milkinff_m1)$coef)
milkinff_p <- milkinff_sum1[rownames(milkinff_sum1)=="time_point_numeric", "Pr(>|t|)"]

inffmatf_sum1 <- as.data.frame(summary(inffmatf_m1)$coef)
inffmatf_p <- inffmatf_sum1[rownames(inffmatf_sum1)=="time_point_numeric", "Pr(>|t|)"]

milkmatf_sum1 <- as.data.frame(summary(milkmatf_m1)$coef)
milkmatf_p <- milkmatf_sum1[rownames(milkmatf_sum1)=="time_point_numeric", "Pr(>|t|)"]


##  create data frame with results from distance comparison
distresdf <- data.frame(comparison=as.factor(c("human_milk_infant_faeces", "infant_faeces_maternal_faeces", "maternal_faeces_human_milk")),
                        n=c(nrow(milk_inff_distances_all), nrow(inff_matf_distances_all), nrow(milk_matf_distances_all)),
                        p=c(milkinff_p, inffmatf_p, milkmatf_p))

## correct for multiple testing
distresdf$FDR_BH <- p.adjust(distresdf$p, method="BH")

## save results data frame
write.table(distresdf, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/all_sample_types/240222_lm_results_Aitchison_distances_betw_corresp_samples_by_time_point_v2.txt", sep="\t", row.names=F, quote=F)


### ===== 21.4. COMPARE PAIR-WISE DISTANCES BETWEEN CORRESPONDING SAMPLES FROM ALL TIME POINTS ===== ###

## Aim: extract distances between pooled samples from all sample types to show that human milk and infant faeces are closer (smaller distances) than infant - maternal faeces and human milk - maternal faeces


## run 3 pair-wise comparisons
comp1_matf_milk <- lm(aitchison_distance ~ comparison, data=distdf[distdf$comparison=="maternal_faeces_human_milk" | distdf$comparison=="human_milk_infant_faeces",])
summary(comp1_matf_milk)

comp2_matf_inff <- lm(aitchison_distance ~ comparison, data=distdf[distdf$comparison=="maternal_faeces_human_milk" | distdf$comparison=="infant_faeces_maternal_faeces",])
summary(comp2_matf_inff)

comp3_milk_inff <- lm(aitchison_distance ~ comparison, data=distdf[distdf$comparison=="human_milk_infant_faeces" | distdf$comparison=="infant_faeces_maternal_faeces",])
summary(comp3_milk_inff)


## save statistics results from models in data frames
comp1_matf_milk_sum <- as.data.frame(summary(comp1_matf_milk)$coef)
comp1_matf_milk_p <- comp1_matf_milk_sum[rownames(comp1_matf_milk_sum)=="comparisonhuman_milk_infant_faeces", "Pr(>|t|)"]

comp2_matf_inff_sum <- as.data.frame(summary(comp2_matf_inff)$coef)
comp2_matf_inff_p <- comp2_matf_inff_sum[rownames(comp2_matf_inff_sum)=="comparisoninfant_faeces_maternal_faeces", "Pr(>|t|)"]

comp3_milk_inff_sum <- as.data.frame(summary(comp3_milk_inff)$coef)
comp3_milk_inff_p <- comp3_milk_inff_sum[rownames(comp3_milk_inff_sum)=="comparisoninfant_faeces_maternal_faeces", "Pr(>|t|)"]


##  create data frame with results from distance comparison
distresdfcomp <- data.frame(comparison=as.factor(c("comparison_matf_milk_to_milk_inff", "comparison_matf_milk_to_inff_matf", "comparison_milk_inff_to_inff_matf")),
                            n=c(nrow(distdf[distdf$comparison=="maternal_faeces_human_milk" | distdf$comparison=="human_milk_infant_faeces",]),
                                nrow(distdf[distdf$comparison=="maternal_faeces_human_milk" | distdf$comparison=="infant_faeces_maternal_faeces",]),
                                nrow(distdf[distdf$comparison=="human_milk_infant_faeces" | distdf$comparison=="infant_faeces_maternal_faeces",])),
                            p=c(comp1_matf_milk_p, comp2_matf_inff_p, comp3_milk_inff_p))

## correct for multiple testing
distresdfcomp$FDR_BH <- p.adjust(distresdfcomp$p, method="BH")

## save results data frame
write.table(distresdfcomp, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/5_microbiota_by_sample_type_and_time/beta_diversity/all_sample_types/240222_lm_results_Aitchison_distances_betw_corresp_samples_by_comparison_and_time_point_v2.txt", sep="\t", row.names=F, quote=F)








