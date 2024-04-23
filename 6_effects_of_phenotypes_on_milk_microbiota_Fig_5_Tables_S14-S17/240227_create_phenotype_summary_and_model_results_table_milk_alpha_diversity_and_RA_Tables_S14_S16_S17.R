######################################################################################################################################################
### CREATE PHENOTYPES SUMMARY STATISTICS + MODEL RESULTS TABLE FOR MILK MICROBIOTA DATA (ALPHA DIV AND REL ABUNDANCES) FOR MILK COMPOSITION PAPER  ###
######################################################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessions type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 0. DATA IMPORT

# 1. SELECT DATA OF INTEREST
#    1.1 MILK MICROBIOTA DATA AND MATERNAL PHENOTYPE DATA
#    1.2 MILK MICROBIOTA DATA AND INFANT PHENOTYPE DATA

# 2. ENSURE CORRECT DATA STRUCTURE
#    2.1 MATERNAL PHENOTYPES
#    2.2 INFANT PHENOTYPES

# 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
#    3.1 MATERNAL PHENOTYPES
#    3.2 INFANT PHENOTYPES

# 4. CHANGE SELECTED CATEGORICAL PHENOTYPE DATA TO NUMERIC FOR USE IN ASSOCIATION STUDIES
#    4.1 MATERNAL PHENOTYPES
#    4.2 INFANT PHENOTYPES

# 5. CREATE TABLE SHOWING PHENOTYPE DATA TYPE (FACTOR OR NUMERIC/INTEGER) AND USE
#    5.1 FUNCTION TO CREATE OVERVIEW TABLES
#    5.2 MATERNAL PHENOTYPES
#    5.3 INFANT PHENOTYPES

# 6. CREATE BARPLOTS AND CREATE HISTOGRAMS TO SHOW DATA DISTRIBUTION TO GUIDE DECISION ON PHENOTYPE DATA TRANSFORMATION
#    6.1 DATA DISTRIBUTION OF MATERNAL PHENOTYPES
#    6.2 DATA DISTRIBUTION OF INFANT PHENOTYPES

# 7. CREATE SUMMARY STATISTICS FOR CATEGORICAL PHENOTYPES
#    7.1 FUNCTION TO CREATE SUMMARY STATISTICS
#    7.2 SUMMARY STATISTICS FOR CATEGORICAL MATERNAL PHENOTYPES
#    7.3 SUMMARY STATISTICS FOR CATEGORICAL INFANT PHENOTYPES

# 8. CREATE SUMMARY STATISTICS FOR NUMERIC/INTEGER PHENOTYPES
#    8.1 FUNCTION TO CREATE SUMMARY STATISTICS
#    8.2 SUMMARY STATISTICS FOR MATERNAL PHENOTYPES
#    8.3 SUMMARY STATISTICS FOR INFANT PHENOTYPES

# 9. COMBINE SUMMARY STATISTICS TABLES FOR ALL PHENOTYPES
#    9.1 PREPARE SUMMARY STATISTICS TABLE WITH CATEGORICAL PHENOTYPES
#    9.2 PREPARE SUMMARY STATISTICS TABLE WITH NUMERIC PHENOTYPES
#    9.3 COMBINE SUMMARY STATISTICS TABLES WITH CATEGORICAL AND NUMERIC PHENOTYPES

# 10. INVERSE-RANK TRANSFORM NUMERIC PHENOTYPES AND MILK ALPHA DIVERSITY DATA
#    10.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
#    10.2 INVERSE-RANK TRANSFORM NUMERIC PHENOTYPE DATA
#    10.3 INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY DATA

# 11. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED NUMERIC PHENOTYPES
#    11.1 FUNCTION TO MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED NUMERIC PHENOTYPES
#    11.2 MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED NUMERIC PHENOTYPES

# 12. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED MILK ALPHA DIVERSITY DATA
#    12.1 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED MILK ALPHA DIVERSITY DATA FROM MATERNAL PHENOTYPE FILE
#    12.2 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED MILK ALPHA DIVERSITY DATA FROM INFANT PHENOTYPE FILE

# 13. CALCULATE MILK RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE
#    13.1 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA FOR DATA FRAME WITH MATERNAL PHENOTYPES
#    13.2 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA FOR DATA FRAME WITH INFANT PHENOTYPES

# 14. CLR-TRANSFORM MILK RELATIVE ABUNDANCES
#    14.1 FUNCTION FOR CLR-TRANSFORMATION
#    14.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES IN DATA FRAME WITH MATERNAL PHENOTYPES
#    14.3 CLR-TRANSFORM MILK RELATIVE ABUNDANCES IN DATA FRAME WITH INFANT PHENOTYPES

# 15. PREPARE MATERNAL AND INFANT DATA SETS FOR ASSOCIATION OF PHENOTYPES WITH MILK MICROBIOTA DATA
#    15.1 OVERVIEW OF CURRENT DATA FRAMES PREPARED FOR ASSOCIATIONS OF MILK MICROBIOTA DATA AND PHENOTYPES
#    15.2 PREPARE DATA FRAMES WITH MATERNAL PHENOTYPES
#    15.3 PREPARE DATA FRAMES WITH INFANT PHENOTYPES

# 16. MODEL FOR ASSOCIATION OF MULTIPLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA WITH CORRECTION FOR DNA ISOLATION BATCH, SEQUENCING READS, TIME AND ID
#    16.1 PREPARE DATA FRAMES WITH MATERNAL STATIC PHENOTYPES
#    16.2 PREPARE DATA FRAMES WITH INFANT STATIC PHENOTYPES
#    16.3 FUNCTION FOR ASSOCIATION OF STATIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR DNA ISOLATION BATCH, SEQUENCING READS, TIME AND ID
#    16.4 RUN ASSOCIATION OF MATERNAL STATIC PHENOTYPES WITH MILK ALPHA DIVERSITY MEASURES AND CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES
#    16.5 RUN ASSOCIATION OF INFANT STATIC PHENOTYPES WITH MILK ALPHA DIVERSITY MEASURES AND CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES

# 17. MODEL FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR SEQUENCING READS
#    17.1 PREPARE DATA FRAMES WITH MATERNAL 1-TIME-POINT PHENOTYPES
#    17.2 PREPARE DATA FRAMES WITH INFANT 1-TIME-POINT PHENOTYPES
#    17.3 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR DNA ISOLATION BATCH AND SEQUENCING READS
#    17.4 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR SEQUENCING READS
#    17.5 RUN ASSOCIATION OF MATERNAL SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA
#    17.6 RUN ASSOCIATION OF INFANT SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA

# 18. MODEL FOR ASSOCIATION OF MULTIPLE-TIME-POINT DYNAMIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR DNA ISOLATION BATCH, SEQUENCING READS, TIME AND ID
#    18.1 PREPARE DATA FRAMES WITH MATERNAL DYNAMIC PHENOTYPES
#    18.2 PREPARE DATA FRAMES WITH INFANT DYNAMIC PHENOTYPES
#    18.3 FUNCTION FOR ASSOCIATION OF DYNAMIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR DNA ISOLATION BATCH, SEQUENCING READS, TIME AND ID
#    18.4 FUNCTION FOR ASSOCIATION OF DYNAMIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR SEQUENCING READS, TIME AND ID
#    18.5 RUN ASSOCIATION OF MATERNAL DYNAMIC PHENOTYPES WITH WITH MILK ALPHA DIVERSITY MEASURES AND CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES
#    18.6 RUN ASSOCIATION OF INFANT DYNAMIC PHENOTYPES WITH WITH MILK ALPHA DIVERSITY MEASURES AND CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES

# 19. COMBINE RESULTS AND CORRECT FOR MULTIPLE TESTING
#    19.1 AVAILABLE RESULTS DATA FRAMES FOR MERGING
#    19.2 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM STATIC PHENOTYPE ASSOCIATIONS
#    19.3 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM 1-TIME-POINT PHENOTYPE ASSOCIATIONS
#    19.4 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM DYNAMIC PHENOTYPE ASSOCIATIONS
#    19.5 MERGE ALL RESULTS DATA FRAMES
#    19.6 EXCLUDE PHENOTYPES
#    19.7 ADD PHENOTYPE GROUPS FOR PLOTTING
#    19.8 ADD PRETTY LABELS FOR PHENOTYPES FOR PLOTTING
#    19.9 ENSURE CORRECT N FOR PHENOTYPES WITH >2 CATEGORIES
#    19.10 SPLIT BY OUTPUT TYPE (ALPHA DIVERSITY / RELATIVE ABUNDANCES)
#    19.11 CORRECT FOR MULTIPLE TESTING AND SAVE RESULTS TABLES

# 20. PREPARE SUPPLEMENTARY TABLES (PHENOTYPE SUMMARY STATISTICS TABLE AND MODEL RESULTS TABLES)
#    20.1 PHENOTYPE SUMMARY STATISTICS TABLE (TABLE S14)
#    20.2 MILK ALPHA DIVERSITY - MODEL RESULTS TABLE (TABLE S16)
#    20.3 MILK RELATIVE BACTERIAL ABUNDANCES - MODEL RESULTS TABLE (TABLE S17)


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/")

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

## resort columns
data <- data[,c(2:14,1,15:ncol(data))]


##### =========================== 1. SELECT DATA OF INTEREST =========================== #####

### ===== 1.1 MILK MICROBIOTA DATA AND MATERNAL PHENOTYPE DATA ===== ###

## select (1) milk samples and (2) every mother only 1x (i.e. remove duplication of data generated for twin babies) and select IDs, maternal phenotype columns and bacterial data
mmilk_tmp <- data[data$sample_type=="human_milk",]
mmilk <- mmilk_tmp[grep(";1", mmilk_tmp$mother_sample_ID), c(1:426,472:1154)] #725x1109


### ===== 1.2 MILK MICROBIOTA DATA AND INFANT PHENOTYPE DATA ===== ###

## select milk samples and select IDs, infant phenotype columns (which can potentially affect maternal milk microbiota) and bacterial data
imilk <- data[data$sample_type=="human_milk", c(1:15,427:435,442:445,472:1154)] #737x711


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

### ===== 2.1 MATERNAL PHENOTYPES ===== ###

## Maternal phenotypes: check that data structure is correct
str(mmilk[,1:50])
str(mmilk[,51:100])
str(mmilk[,101:150])
str(mmilk[,151:200])
str(mmilk[,201:250])
str(mmilk[,251:300])
str(mmilk[,301:350])
str(mmilk[,351:400])
str(mmilk[,401:450])
str(mmilk[,c(451:486,1109)])
# -> data structure is correct


### ===== 2.2 INFANT PHENOTYPES ===== ###

## Infant phenotypes: check that data structure is correct
str(imilk[,1:50])
str(imilk[,c(50:88,711)])
# -> data structure is correct


##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

### ===== 3.1 MATERNAL PHENOTYPES ===== ###

## Maternal phenotypes: for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
mmilk$mother_milk_collection_season <- factor(mmilk$mother_milk_collection_season, levels = c("Spring", "Summer", "Autumn", "Winter"))
mmilk$mother_milk_collection_month <- factor(mmilk$mother_milk_collection_month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",  "Oct",  "Nov", "Dec"))
mmilk$mother_milk_collection_notes <- factor(mmilk$mother_milk_collection_notes, levels = c("No_records_of_incorrect_sample_handling", levels(mmilk$mother_milk_collection_notes)[c(1,2,4)]))
mmilk$mother_milk_collection_breasts_sampled <- factor(mmilk$mother_milk_collection_breasts_sampled, levels=c("one_breast", "both_breasts"))
mmilk$mother_genetics_blood_group_genotype <- factor(mmilk$mother_genetics_blood_group_genotype, levels = c("O/O", "A/O", "A/A", "B/O", "B/B", "A/B"))
mmilk$mother_blood_group <- factor(mmilk$mother_blood_group, levels = c("O", "A", "B", "AB"))
mmilk$mother_exp_living_situation <- factor(mmilk$mother_exp_living_situation, levels = c("big_house", "small_house", "flat", "farm", "other"))
mmilk$mother_birth_gestational_age_categories <- factor(mmilk$mother_birth_gestational_age_categories, levels = c("term", "preterm"))
mmilk$mother_breastpump_brand <- factor(mmilk$mother_breastpump_brand, levels = c("Philips", "Medela", "Ardo", "Other"))

# also change it for milk groups
mmilk$mother_milk_HMO_milk_group <- factor(mmilk$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


### ===== 3.2 INFANT PHENOTYPES ===== ###

## Infant phenotypes: for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
imilk$infant_birth_delivery_mode <- factor(imilk$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
imilk$infant_birth_delivery_mode_detailed <- factor(imilk$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                        "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
imilk$infant_ffq_feeding_type_delivery <- factor(imilk$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))

# also change it for milk groups
imilk$mother_milk_HMO_milk_group <- factor(imilk$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 4. CHANGE SELECTED CATEGORICAL PHENOTYPE DATA TO NUMERIC FOR USE IN ASSOCIATION STUDIES  =========================== #####

## Some categorical phenotypes will be used as numeric phenotypes in models. These phenotypes are here changed to numeric before making the summary statistics tables so that the tables show the data as used for associations.
## I will later add info to the supplementary table on which factors/categories the numbers correspond to.


### ===== 4.1 MATERNAL PHENOTYPES ===== ###

## vector with the names of the columns that will be changed to numeric
matnumcatphenos <- c("mother_milk_collection_time","mother_education_level","mother_net_income","mother_health_dental")

## check data before change
# head(mmilk[,colnames(mmilk) %in% matnumcatphenos])
# tail(mmilk[,colnames(mmilk) %in% matnumcatphenos])
# table(mmilk$mother_milk_collection_time, useNA="ifany")
# table(mmilk$mother_education_level, useNA="ifany")
# table(mmilk$mother_net_income, useNA="ifany")
# table(mmilk$mother_health_dental, useNA="ifany")

## change data entries for numcatphenos to numeric
for (i in matnumcatphenos){mmilk[,i] <- as.numeric(as.character(gsub("__.*", "", mmilk[,i])))}

## show conversion from categorical to numeric data structure in column name
a <- c()
for (i in matnumcatphenos){a <- c(a, grep(i, colnames(mmilk)))}
length(a) #5
length(matnumcatphenos) #4 -> one phenotype too much was pulled from the data frame
setdiff(colnames(mmilk)[a], matnumcatphenos) # -> remove "mother_health_dental_bleeding_gums_during_preg" from vector a
colnames(mmilk)[a]
a <- a[-5]
setdiff(colnames(mmilk)[a], matnumcatphenos) 
length(intersect(colnames(mmilk)[a], matnumcatphenos))

for (i in a){colnames(mmilk)[i] <- paste0(colnames(mmilk)[i], "_numeric")}

## check data after change
# head(mmilk[,a])
# tail(mmilk[,a])
# table(mmilk$mother_milk_collection_time_numeric, useNA="ifany")
# table(mmilk$mother_education_level_numeric, useNA="ifany")
# table(mmilk$mother_net_income_numeric, useNA="ifany")
# table(mmilk$mother_health_dental_numeric, useNA="ifany")


### ===== 4.2 INFANT PHENOTYPES ===== ###

## No infant phenotypes were changed to numeric.


##### =========================== 5. CREATE TABLE SHOWING PHENOTYPE DATA TYPE (FACTOR OR NUMERIC/INTEGER) AND USE  =========================== #####

### ===== 5.1 FUNCTION TO CREATE OVERVIEW TABLES ===== ###

create.overview.table <- function(inputdata){
  # create empty vectors
  print("Creating empty vectors")
  column_name <- c()
  data_type <- c()
  participant <- c()
  
  # add data to vectors
  print("Adding data to vectors")
  for (i in 1:ncol(inputdata)){
    print(paste0("Adding data from column ", i, " ", colnames(inputdata)[i]))
    column_name <- c(column_name, colnames(inputdata)[i])
    data_type <- c(data_type, class(inputdata[,i]))
    participant <- c(participant, paste0(gsub("_.*", "", colnames(inputdata)[i])))
  }
  
  # combine in 1 data frame
  print("Combining data in data frame")
  my_df <- data.frame(column_name=column_name,
                      data_type=data_type,
                      participant=participant)
  
  my_df[my_df$participant=="time","participant"] <- NA
  
  print("Saving data frame")
  return(my_df)
}


### ===== 5.2 MATERNAL PHENOTYPES ===== ###

## Maternal phenotypes: create overview
overview_mmilk <- create.overview.table(mmilk[,c(8:10,16:426,428:458)])

## add info in data_use column to indicate whether a used phenotype is static or dynamic or 1-time-point
overview_mmilk$model_use <- c(rep("dynamic",length(1:2)),           #time_point and time_point_numeric
                              rep("static",length(3)),              #type_pregnancy
                              rep("dynamic",length(4:9)),           #milk collection -> exclude selected phenotypes from models
                              rep("static",length(10:384)),         #maternal genetics (exclude genotypes for FUT2 and FUT3?), maternal anthropometrics, food pref, ffq (exclude items from 1st analysis), stress, education, exp living and pets, mat health, preg/birth comp
                              rep("dynamic",length(385:386)),       #subclinical mastitis and Na/K ratio
                              rep("static",length(387:398)),        #maternal GIT diary, ROME, pre-preg med use
                              rep("1_time_point",length(399:404)),  #preg and birth med use
                              rep("dynamic",length(405:414)),       #maternal postpartum med use, breastpumping
                              rep("static",length(415:417)),        #Le/Se milk group
                              rep("dynamic",length(418:445)))       #milk HMOs

## add info on which phenotypes to exclude as x from association models
overview_mmilk$model_exclusion <- c(rep("exclude",length(1:2)),    #time_point and time_point_numeric (exclude as it will be corrected for)
                                   rep("include",length(3:8)),     #type_pregnancy, milk collection
                                   rep("exclude",length(9:14)),    #milk collection notes, maternal FUT2 and FUT3 genotypes
                                   rep("include",length(15:445)))  #maternal anthropometrics, food pref, food groups, food items, nutrients, stress, education, exp living and pets, mat health, preg/birth comp, subclinical mastitis and Na/K ratio, maternal GIT diary, ROME, preg and birth med use, maternal postpartum med use, breastpumping, Le/Se milk group, milk HMOs

## add info on which microbiota time points phenotypes were/will be linked to
overview_mmilk$model_time_point <- c(rep(NA,length(1:2)),                      #time_point and time_point_numeric (exclude as it will be corrected for)
                                    rep("1_to_6_month(s)",length(3:8)),        #type_pregnancy, milk collection
                                    rep(NA,length(9:14)),                      #milk collection notes, maternal FUT2 and FUT3 genotypes
                                    rep("1_to_6_month(s)",length(15:384)),     #maternal anthropometrics, food pref, food groups, food items, nutrients, stress, education, exp living and pets, mat health, preg/birth comp
                                    rep("1_to_3_month(s)",length(385:386)),    #subclinical mastitis and Na/K ratio,
                                    rep("1_to_6_month(s)",length(387:398)),    #maternal GIT diary, ROME, prior preg med use
                                    rep("1_month",length(399:404)),            #preg and birth med use
                                    rep("1_and_3_month(s)",length(405:412)),   #maternal postpartum med use
                                    rep("1_to_3_month(s)",length(413:414)),    #breastpumping (no data for M6 was included)
                                    rep("1_to_6_month(s)",length(415:445)))    #Le/Se milk group, milk HMOs

# test <- mmilk[,c(8:10,16:426,428:458)]
# for (i in 1:ncol(test)){
#   print(paste0(i, "__", colnames(test)[i]))
#   print(table(test[!is.na(test[,i]), "time_point"]))
# }

## separate phenotypes by data type
factor_mmilk <- overview_mmilk[overview_mmilk$data_type=="factor",] #70 categorical phenotypes
numeric_mmilk <- overview_mmilk[overview_mmilk$data_type=="numeric" | overview_mmilk$data_type=="integer",] #375 numeric phenotypes


### ===== 5.3 INFANT PHENOTYPES ===== ###

## Infant phenotypes: create overview
overview_imilk <- create.overview.table(imilk[,c(8:9,16:28)])

## add info in data_use column to indicate whether a used phenotype is static or dynamic or dynamic or 1-time-point
overview_imilk$model_use <- c(rep("dynamic",length(1:2)),    #time_point and time_point_numeric
                             rep("static",length(3:11)),     #infant genetics, infant sex, infant birth place and mode, birthweight
                             rep("1_time_point",length(12:13)), #infant feeding birth
                             rep("dynamic",length(14)),         #infant feeding mode
                             rep("1_time_point",length(15)))    #infant food introduction

## add info on which phenotypes to exclude as x from association models
overview_imilk$model_exclusion <- c(rep("exclude",length(1:3)),     #time_point and time_point_numeric (exclude as it will be corrected for), infant genetics (FUT2)
                                   rep("include",length(4:15)))    #infant Se status, infant sex, infant birth place and mode, birthweight, infant feeding

## add info on which microbiota time points phenotypes were/will be linked to
overview_imilk$model_time_point <- c(rep(NA,length(1:3)),                    #time_point and time_point_numeric, infant genetics (FUT2)
                                    rep("1_to_6_month(s)",length(4:11)),     #infant Se status, infant sex, infant birth place and mode, birthweight
                                    rep("1_month",length(12:13)),         #infant feeding birth
                                    rep("1_to_6_month(s)",length(14)),       #infant feeding mode
                                    rep("6_months",length(15)))            #infant food introduction

# test <- imilk[,c(8:9,16:28)]
# for (i in 1:ncol(test)){
#   print(paste0(i, "__", colnames(test)[i]))
#   print(table(test[!is.na(test[,i]), "time_point"]))
# }

## separate phenotypes by data type
factor_imilk <- overview_imilk[overview_imilk$data_type=="factor",] #12 categorical phenotypes
numeric_imilk <- overview_imilk[overview_imilk$data_type=="numeric" | overview_imilk$data_type=="integer",] #3 numeric phenotypes


##### =========================== 6. CREATE BARPLOTS AND CREATE HISTOGRAMS TO SHOW DATA DISTRIBUTION TO GUIDE DECISION ON PHENOTYPE DATA TRANSFORMATION  =========================== #####

library(ggplot2)

### ===== 6.1 DATA DISTRIBUTION OF MATERNAL PHENOTYPES ===== ###

## create bar plots for all categorical phenotypes
pdf(paste0("phenotype_data_distribution/240130_milk_microbiota_categorical_maternal_phenotypes_data_distribution_barplots_1_n70.pdf"))
for (i in factor_mmilk$column_name){
  print(paste0("Creating barplots for ", i))
  barplots <- ggplot(mmilk[!is.na(mmilk[,i]),], aes(x=mmilk[!is.na(mmilk[,i]),i]))+
    geom_bar()+
    labs(i)+
    ggtitle(paste0(i))+
    theme(axis.text.x = element_text(angle = 90))
  print(barplots)
}
dev.off()

## create histograms for all numeric phenotypes
pdf(paste0("phenotype_data_distribution/240130_milk_microbiota_numeric_maternal_phenotypes_data_distribution_histograms_1_n375.pdf"))
for (i in numeric_mmilk$column_name){
  print(paste0("Creating histogram for ", i))
  histo <- ggplot(mmilk[!is.na(mmilk[,i]),], aes(x=mmilk[!is.na(mmilk[,i]),i]))+
    geom_histogram(bins=50)+
    labs(i)+
    ggtitle(paste0(i))
  print(histo)
}
dev.off()


### ===== 6.2 DATA DISTRIBUTION OF INFANT PHENOTYPES ===== ###

## create bar plots for all categorical phenotypes
pdf(paste0("phenotype_data_distribution/240130_milk_microbiota_categorical_infant_phenotypes_data_distribution_barplots_1_n12.pdf"))
for (i in factor_imilk$column_name){
  print(paste0("Creating barplots for ", i))
  barplots <- ggplot(imilk[!is.na(imilk[,i]),], aes(x=imilk[!is.na(imilk[,i]),i]))+
    geom_bar()+
    labs(i)+
    ggtitle(paste0(i))+
    theme(axis.text.x = element_text(angle = 90))
  print(barplots)
}
dev.off()

## create histograms for all numeric phenotypes
pdf(paste0("phenotype_data_distribution/240130_milk_microbiota_numeric_infant_phenotypes_data_distribution_histograms_1_n3.pdf"))
for (i in numeric_imilk$column_name){
  print(paste0("Creating histogram for ", i))
  histo <- ggplot(imilk[!is.na(imilk[,i]),], aes(x=imilk[!is.na(imilk[,i]),i]))+
    geom_histogram(bins=50)+
    labs(i)+
    ggtitle(paste0(i))
  print(histo)
}
dev.off()


##### =========================== 7. CREATE SUMMARY STATISTICS FOR CATEGORICAL PHENOTYPES  =========================== #####

### ===== 7.1 FUNCTION TO CREATE SUMMARY STATISTICS ===== ###

create.summary.factors <- function(inputdata){
  # create empty vectors
  print("Creating empty vectors")
  column_name <- c()
  
  n_total <- c()
  n_answered <- c()
  perc_answered <- c()
  n_missing <- c()
  perc_missing <- c()
  
  all_levels <- c()
  n_per_level <- c()
  
  # add data to vectors
  print("Adding data to vectors")
  for (i in 1:ncol(inputdata)){
    print(paste0("Adding data from column ", i, " ", colnames(inputdata)[i]))
    
    column_name <- c(column_name, colnames(inputdata)[i])
    
    n_total <- c(n_total, nrow(inputdata))
    n_answered    <- c(n_answered,    nrow(inputdata[!is.na(inputdata[,i]),]))
    perc_answered <- c(perc_answered, nrow(inputdata[!is.na(inputdata[,i]),])/nrow(inputdata))
    n_missing     <- c(n_missing,     nrow(inputdata[is.na(inputdata[,i]),]))
    perc_missing  <- c(perc_missing,  nrow(inputdata[is.na(inputdata[,i]),])/nrow(inputdata))
    
    all_levels <- c(all_levels, paste(collapse=";",levels(as.factor(as.character(inputdata[,i])))))
    n_per_level <- c(n_per_level, paste(collapse=";",table(as.factor(as.character(inputdata[,i])))))
  }
  
  # combine in 1 data frame
  print("Combining data in data frame")
  my_df <- data.frame(column_name=column_name,
                      all_levels=all_levels,
                      n_per_level=n_per_level,
                      n_total=n_total,
                      n_answered=n_answered,
                      perc_answered=perc_answered,
                      n_missing=n_missing,
                      perc_missing=perc_missing)
  
  print("Saving data frame")
  return(my_df)
}


### ===== 7.2 SUMMARY STATISTICS FOR CATEGORICAL MATERNAL PHENOTYPES ===== ###

## create summary
factor_mmilk_sum <- create.summary.factors(mmilk[,colnames(mmilk) %in% factor_mmilk$column_name])

## polish perc
factor_mmilk_sum$perc_answered <- signif(factor_mmilk_sum$perc_answered, 2)
factor_mmilk_sum$perc_missing <- signif(factor_mmilk_sum$perc_missing, 2)

## change colnames
colnames(factor_mmilk_sum)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== 7.3 SUMMARY STATISTICS FOR CATEGORICAL INFANT PHENOTYPES ===== ###

## create summary
factor_imilk_sum <- create.summary.factors(imilk[,colnames(imilk) %in% factor_imilk$column_name])

## polish perc
factor_imilk_sum$perc_answered <- signif(factor_imilk_sum$perc_answered, 2)
factor_imilk_sum$perc_missing <- signif(factor_imilk_sum$perc_missing, 2)

## change colnames
colnames(factor_imilk_sum)[1:3] <- c("Phenotype", "Categories", "n_per_category")


##### =========================== 8. CREATE SUMMARY STATISTICS FOR NUMERIC/INTEGER PHENOTYPES  =========================== #####

### ===== 8.1 FUNCTION TO CREATE SUMMARY STATISTICS ===== ###

create.summary.numeric <- function(inputdata){
  
  # calculate summary statistics
  print(paste0("Calculating summary statistics"))

  SumStatsTmp <- NULL
  for (i in 1:ncol(inputdata)){
    summary.data <- as.vector(summary(inputdata[,i]))
    if (length(summary.data) == 6){
      summary.data <- c(summary.data, 0)
    }
    SumStatsTmp <- as.data.frame(rbind(SumStatsTmp, summary.data))
  }

  # fix row- and colnames
  print(paste0("Setting rownames and colnames in summary statistics table"))
  rownames(SumStatsTmp) <- 1:nrow(SumStatsTmp)
  colnames(SumStatsTmp) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.", "NA's")
  
  # add colname as column
  print(paste0("Adding phenotype in column"))
  SumStatsTmp$column_name <- colnames(inputdata)
  
  # add n_total
  print(paste0("Adding n_total"))
  SumStatsTmp$n_total <- nrow(inputdata)
  
  # add perc_missing
  print(paste0("Adding perc_missing"))
  colnames(SumStatsTmp)[7] <- "n_missing"
  SumStatsTmp$perc_missing <- c(SumStatsTmp$n_missing/SumStatsTmp$n_total)
  
  # add n_answered
  print(paste0("Adding n_answered"))
  SumStatsTmp$n_answered <- nrow(inputdata)-SumStatsTmp$n_missing
  
  # add perc_answered
  print(paste0("Adding perc_answered"))
  SumStatsTmp$perc_answered <- c(SumStatsTmp$n_answered/SumStatsTmp$n_total)
  
  # add SD
  print(paste0("Adding SD"))
  sd <- c()
  for (i in 1:ncol(inputdata)){
    sd <- c(sd, sd(inputdata[,i], na.rm=T))
  }
  SumStatsTmp$SD <- sd
  
  # fix colnames
  print(paste0("Prettifying colnames"))
  colnames(SumStatsTmp)[c(1:6,8)] <- c("Minimum", "1st_quartile", "Median", "Mean", "3rd_quartile", "Maximum", "Phenotype")
  
  # resort columns
  print(paste0("Resorting colums"))
  SumStats <- SumStatsTmp[,c(8,4,13,3,2,5,1,6,9,11,12,7,10)]
  
  print(paste0("Returning summary statistics table"))
  return(SumStats)
}


### ===== 8.2 SUMMARY STATISTICS FOR MATERNAL PHENOTYPES ===== ###

## Calculate summary statistics
numeric_mmilk_sum <- create.summary.numeric(mmilk[,colnames(mmilk) %in% numeric_mmilk$column_name])


## Add a 'Note' column to indicate levels of categorical phenotypes that were converted to numeric for associations
numeric_mmilk_sum$Note <- NA

# matnumcatphenos

## Add note for mother_milk_collection_time
data$mother_milk_collection_time <- factor(data$mother_milk_collection_time, levels = c("1__20h", "2__21h", "3__22h", "4__23h", "5__00h",
                                                                                      "6__01h", "7__02h", "8__03h", "9__04h", "10__05h",
                                                                                      "11__06h", "12__07h", "13__08h", "14__09h", "15__10h",
                                                                                      "16__11h", "17__12h", "18__13h", "19__14h", "20__15h",
                                                                                      "21__16h", "22__17h", "23__18h", "24__19h")) #first correctly sort the levels in the data data frame, from which we pull the original levels
numeric_mmilk_sum[numeric_mmilk_sum$Phenotype=="mother_milk_collection_time_numeric", "Note"] <- c(paste0("This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: ",
                                                                                                   paste(collapse=";", levels(data$mother_milk_collection_time)))) # add note and levels

## Add note for mother_education_level
numeric_mmilk_sum[numeric_mmilk_sum$Phenotype=="mother_education_level_numeric", "Note"] <- c(paste0("This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: ",
                                                                                              paste(collapse=";", levels(data$mother_education_level))))

## Add note for mother_net_income
numeric_mmilk_sum[numeric_mmilk_sum$Phenotype=="mother_net_income_numeric", "Note"] <- c(paste0("This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: ",
                                                                                         paste(collapse=";", levels(data$mother_net_income))))

## Add note for mother_health_dental
numeric_mmilk_sum[numeric_mmilk_sum$Phenotype=="mother_health_dental_numeric", "Note"] <- c(paste0("This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: ",
                                                                                            paste(collapse=";", levels(data$mother_health_dental))))

## Check notes
# numeric_mmilk_sum[!is.na(numeric_mmilk_sum$Note),]


### ===== 8.3 SUMMARY STATISTICS FOR INFANT PHENOTYPES ===== ###

## Calculate summary statistics
numeric_imilk_sum <- create.summary.numeric(imilk[,colnames(imilk) %in% numeric_imilk$column_name])


## Add a 'Note' column to indicate levels of categorical phenotypes that were converted to numeric for associations
numeric_imilk_sum$Note <- NA


##### =========================== 9. COMBINE SUMMARY STATISTICS TABLES FOR ALL PHENOTYPES  =========================== #####

### ===== 9.1 PREPARE SUMMARY STATISTICS TABLE WITH CATEGORICAL PHENOTYPES ===== ###

## combine maternal and infant summary statistics tables, excluding the time_point columns #80x8
factor_sum <- as.data.frame(rbind(factor_mmilk_sum[factor_mmilk_sum$Phenotype!="time_point",],
                                  factor_imilk_sum[factor_imilk_sum$Phenotype!="time_point",]))


## add info on data type, model and linked microbiota time points from factor overview tables
# combine maternal and infant factor overview tables
factor_all <- as.data.frame(rbind(factor_mmilk, factor_imilk))
colnames(factor_all)[1] <- "Phenotype"

# add info
library(dplyr)
factor_sum2 <- left_join(factor_sum, factor_all[,c(1:2,4:6)], by="Phenotype")


### ===== 9.2 PREPARE SUMMARY STATISTICS TABLE WITH NUMERIC PHENOTYPES ===== ###

## combine maternal and infant summary statistics tables, excluding the time_point_numeric columns #378x6
numeric_sum <- as.data.frame(rbind(numeric_mmilk_sum[numeric_mmilk_sum$Phenotype!="time_point_numeric",],
                                   numeric_imilk_sum[numeric_imilk_sum$Phenotype!="time_point_numeric",]))

## add info on data type, model and linked microbiota time points from numeric overview tables
# combine maternal and infant numeric overview tables
numeric_all <- as.data.frame(rbind(numeric_mmilk, numeric_imilk))
colnames(numeric_all)[1] <- "Phenotype"

# add info
# library(dplyr)
numeric_sum2 <- left_join(numeric_sum, numeric_all[,c(1:2,4:6)], by="Phenotype")


### ===== 9.3 COMBINE SUMMARY STATISTICS TABLES WITH CATEGORICAL AND NUMERIC PHENOTYPES ===== ###

## combine tables #456x20
library(plyr)
all_sum_tmp <- rbind.fill(factor_sum2, numeric_sum2)

## resort columns
all_sum <- all_sum_tmp[,c(1,9,4:8,2:3,13:19,11,10,12,20)] #456x20

## save phenotype summary statistics table
# write.table(all_sum, file="phenotype_summary_statistics/240301_milk_microbiota_all_maternal_and_infant_phenotypes_summary_statistics_n456.txt", row.names=F, col.names=T, sep="\t", quote=F)


##### =========================== 10. INVERSE-RANK TRANSFORM NUMERIC PHENOTYPES AND MILK ALPHA DIVERSITY DATA =========================== #####

### ===== 10.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 10.2 INVERSE-RANK TRANSFORM NUMERIC PHENOTYPE DATA ===== ###

## Maternal phenotypes
mmilk_invr <- mmilk
mmilk_invr[,numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"]] <- as.data.frame(apply(mmilk_invr[,numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"]], 2, invrank))

colnames(mmilk_invr)[which(colnames(mmilk_invr) %in% numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"])] <- paste0(colnames(mmilk_invr)[which(colnames(mmilk_invr) %in% numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"])], "_invr")


## Infant phenotypes
imilk_invr <- imilk
imilk_invr[,numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"]] <- as.data.frame(apply(imilk_invr[,numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"]], 2, invrank))

colnames(imilk_invr)[which(colnames(imilk_invr) %in% numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"])] <- paste0(colnames(imilk_invr)[which(colnames(imilk_invr) %in% numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"])], "_invr")


### ===== 10.3 INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY DATA ===== ###

## Maternal data frame
mmilk_invr[,grep("alpha_div",colnames(mmilk_invr))] <- as.data.frame(apply(mmilk_invr[,grep("alpha_div",colnames(mmilk_invr))], 2, invrank))
colnames(mmilk_invr)[grep("alpha_div", colnames(mmilk_invr))] <- paste0(colnames(mmilk_invr)[grep("alpha_div", colnames(mmilk_invr))], "_invr")

## Infant data frame
imilk_invr[,grep("alpha_div",colnames(imilk_invr))] <- as.data.frame(apply(imilk_invr[,grep("alpha_div", colnames(imilk_invr))], 2, invrank))
colnames(imilk_invr)[grep("alpha_div", colnames(imilk_invr))] <- paste0(colnames(imilk_invr)[grep("alpha_div", colnames(imilk_invr))], "_invr")


##### =========================== 11. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED NUMERIC PHENOTYPES =========================== #####

# library(ggplot2)

### ===== 11.1 FUNCTION TO MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED NUMERIC PHENOTYPES ===== ###

## make histograms for inverse-rank transformed numeric phenotype data
check.distribution.invr.pheno <- function(inputdata, phenotype_columns, outputname){
  pdf(paste0("phenotype_data_distribution/", outputname, ".pdf"))
  
  for (i in c(phenotype_columns)){
    print(paste0("Creating histogram for ", i))
    milkhist <- ggplot(inputdata,
                      aes(x=as.numeric(inputdata[,i]))) + geom_histogram() + labs(x=i)
    print(milkhist)
  }
  
  dev.off()
}


### ===== 11.2 MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED NUMERIC PHENOTYPES ===== ###

## Maternal phenotypes
mmilk_invr_phenotype_columns <- paste0(numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"], "_invr")
check.distribution.invr.pheno(mmilk_invr,
                              phenotype_columns = mmilk_invr_phenotype_columns,
                              outputname = "240130_milk_microbiota_numeric_maternal_phenotypes_data_distribution_histograms_2_after_invr_transformation_n374")

# Infant phenotypes
imilk_invr_phenotype_columns <- paste0(numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"], "_invr")
check.distribution.invr.pheno(imilk_invr,
                              phenotype_columns = imilk_invr_phenotype_columns,
                              outputname = "240130_milk_microbiota_numeric_infant_phenotypes_data_distribution_histograms_2_after_invr_transformation_n2")


##### =========================== 12. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED MILK ALPHA DIVERSITY DATA =========================== #####

### ===== 12.1 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED MILK ALPHA DIVERSITY DATA FROM MATERNAL PHENOTYPE FILE ===== ###

## 1. Create histograms for inverse-rank transformed alpha diversity measures from all time points combined
pdf(paste0("microbiota_data_distribution/240130_milk_alpha_diversity_measures_for_maternal_phenotypes_histograms.pdf"))
for (i in c(483:485)){
  milkhist <- ggplot(mmilk_invr,
                     aes(x=as.numeric(mmilk_invr[,i]))) + geom_histogram() + labs(x=colnames(mmilk_invr)[i])
  print(milkhist)
}
dev.off()


### ===== 12.2 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED MILK ALPHA DIVERSITY DATA FROM INFANT PHENOTYPE FILE ===== ###

## 1. Create histograms for inverse-rank transformed alpha diversity measures from all time points combined
pdf(paste0("microbiota_data_distribution/240130_milk_alpha_diversity_measures_for_infant_phenotypes_histograms.pdf"))
for (i in c(85:87)){
  milkhist <- ggplot(imilk_invr,
                     aes(x=as.numeric(imilk_invr[,i]))) + geom_histogram() + labs(x=colnames(imilk_invr)[i])
  print(milkhist)
}
dev.off()


##### =========================== 13. CALCULATE MILK RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE =========================== #####

### ===== 13.1 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA FOR DATA FRAME WITH MATERNAL PHENOTYPES ===== ###

## save only columns with absolute bacterial abundances
rownames(mmilk_invr) <- mmilk_invr$seq_16S_sample_ID
mmilkbact <- mmilk_invr[,486:ncol(mmilk_invr)] #725x624

## calculate relative bacterial abundances
mmilkbact_relab <- (mmilkbact/rowSums(mmilkbact))
table(rowSums(mmilkbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
mmilkbact_relab_tmp <- mmilkbact_relab[,-624]

## select only bacterial genera with ≥10% prevalence
mmilkbact_relab_tmp2 <- mmilkbact_relab_tmp[,colSums(mmilkbact_relab_tmp>0)>=73]
length(mmilkbact_relab_tmp2) #36
# colnames(mmilkbact_relab_tmp2)

## merge back with metadata
mmilk_invr_RA <- merge(mmilk_invr[,1:485], mmilkbact_relab_tmp2, by="row.names")
rownames(mmilk_invr_RA) <- mmilk_invr_RA$Row.names
mmilk_invr_RA <- mmilk_invr_RA[,-1]


### ===== 13.2 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA FOR DATA FRAME WITH INFANT PHENOTYPES ===== ###

## save only columns with absolute bacterial abundances
rownames(imilk_invr) <- imilk_invr$seq_16S_sample_ID
imilkbact <- imilk_invr[,88:ncol(imilk_invr)] #737x624

## calculate relative bacterial abundances
imilkbact_relab <- (imilkbact/rowSums(imilkbact))
table(rowSums(imilkbact_relab), useNA="ifany") #all samples have 1

## select only bacterial genera that were included for association with the maternal phenotypes
imilkbact_relab_tmp <- imilkbact_relab[,colnames(imilkbact_relab) %in% colnames(mmilkbact_relab_tmp2)]
length(imilkbact_relab_tmp) #36
# colnames(imilkbact_relab_tmp)

## merge back with metadata
imilk_invr_RA <- merge(imilk_invr[,1:87], imilkbact_relab_tmp, by="row.names")
rownames(imilk_invr_RA) <- imilk_invr_RA$Row.names
imilk_invr_RA <- imilk_invr_RA[,-1]


##### =========================== 14. CLR-TRANSFORM MILK RELATIVE ABUNDANCES =========================== #####

### ===== 14.1 FUNCTION FOR CLR-TRANSFORMATION ===== ###

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


### ===== 14.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES IN DATA FRAME WITH MATERNAL PHENOTYPES ===== ###

mmilk_invr_RAclr <- mmilk_invr_RA
mmilk_invr_RAclr[,486:ncol(mmilk_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(mmilk_invr_RAclr[,486:ncol(mmilk_invr_RAclr)], mmilk_invr_RAclr[,486:ncol(mmilk_invr_RAclr)]))

colnames(mmilk_invr_RAclr)[486:ncol(mmilk_invr_RAclr)] <- c(paste0(colnames(mmilk_invr_RAclr)[486:ncol(mmilk_invr_RAclr)], "_clr"))


### ===== 14.3 CLR-TRANSFORM MILK RELATIVE ABUNDANCES IN DATA FRAME WITH INFANT PHENOTYPES ===== ###

imilk_invr_RAclr <- imilk_invr_RA
imilk_invr_RAclr[,88:ncol(imilk_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(imilk_invr_RAclr[,88:ncol(imilk_invr_RAclr)], imilk_invr_RAclr[,88:ncol(imilk_invr_RAclr)]))

colnames(imilk_invr_RAclr)[88:ncol(imilk_invr_RAclr)] <- c(paste0(colnames(imilk_invr_RAclr)[88:ncol(imilk_invr_RAclr)], "_clr"))


##### =========================== 15. PREPARE MATERNAL AND INFANT DATA SETS FOR ASSOCIATION OF PHENOTYPES WITH MILK MICROBIOTA DATA =========================== #####

### ===== 15.1 OVERVIEW OF CURRENT DATA FRAMES PREPARED FOR ASSOCIATIONS OF MILK MICROBIOTA DATA AND PHENOTYPES ===== ###

## Now the data frames are ready to be used for the association models:
# Data                                          Name              Size
# Maternal phenotypes, milk microbiota data:    mmilk_invr_RAclr: 725x521
# Infant phenotypes,   milk microbiota data:    imilk_invr_RAclr: 737x123


### ===== 15.2 PREPARE DATA FRAMES WITH MATERNAL PHENOTYPES ===== ###

## Create data frame from mmilk_invr_RAclr for use in association models for maternal static phenotypes
rownames(mmilk_invr_RAclr) <- mmilk_invr_RAclr$mother_sample_ID

overview_mmilk[overview_mmilk$model_exclusion=="include" & (overview_mmilk$data_type=="numeric" | overview_mmilk$data_type=="integer"), "column_name"] <-
  paste0(overview_mmilk[overview_mmilk$model_exclusion=="include" & (overview_mmilk$data_type=="numeric" | overview_mmilk$data_type=="integer"), "column_name"], "_invr")


### ===== 15.3 PREPARE DATA FRAMES WITH INFANT PHENOTYPES ===== ###

## Create data frame from imilk_invr_RAclr for use in association models for infant static phenotypes
rownames(imilk_invr_RAclr) <- imilk_invr_RAclr$infant_sample_ID

overview_imilk[overview_imilk$model_exclusion=="include" & (overview_imilk$data_type=="numeric" | overview_imilk$data_type=="integer"), "column_name"] <-
  paste0(overview_imilk[overview_imilk$model_exclusion=="include" & (overview_imilk$data_type=="numeric" | overview_imilk$data_type=="integer"), "column_name"], "_invr")


##### =========================== 16. MODEL FOR ASSOCIATION OF MULTIPLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA WITH CORRECTION FOR DNA ISOLATION BATCH, SEQUENCING READS, TIME AND ID =========================== #####

### ===== 16.1 PREPARE DATA FRAMES WITH MATERNAL STATIC PHENOTYPES ===== ###

# selected basic info
basic_mmilk_invr_RAclr <- mmilk_invr_RAclr[,c(1:9,11:15,466,481,483:485,486:521)] #IDs, DNA isolation batch, clean sequencing reads, invr-transformed milk alpha diversity measures and clr-transformed RA

# select maternal static phenotypes that should be included
v_static_mmilk_invr_RAclr <- overview_mmilk[overview_mmilk$participant=="mother" & overview_mmilk$model_use=="static" & overview_mmilk$model_exclusion=="include","column_name"]
pheno_mmilk_invr_RAclr_static <- mmilk_invr_RAclr[,colnames(mmilk_invr_RAclr) %in% v_static_mmilk_invr_RAclr]

# combine in 1 data frame that is to be used for models
mmilk_invr_RAclr_static <- merge(basic_mmilk_invr_RAclr, pheno_mmilk_invr_RAclr_static, by="row.names")
mmilk_invr_RAclr_static <- mmilk_invr_RAclr_static[,-1] #remove Row.names column
#725x441, of which 386 are phenotypes
#str(mmilk_invr_RAclr_static)

## Maternal data set (mmilk_invr_RAclr_static)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56:441 = maternal static phenotypes


### ===== 16.2 PREPARE DATA FRAMES WITH INFANT STATIC PHENOTYPES ===== ###

# selected basic info
basic_imilk_invr_RAclr <- imilk_invr_RAclr[,c(1:9,11:15,68,83,85:87,88:123)] #IDs, DNA isolation batch, clean sequencing reads, invr-transformed milk alpha diversity measures and clr-transformed RA

# select infant static phenotypes that should be included
v_static_imilk_invr_RAclr <- overview_imilk[overview_imilk$participant=="infant" & overview_imilk$model_use=="static" & overview_imilk$model_exclusion=="include","column_name"]
pheno_imilk_invr_RAclr_static <- imilk_invr_RAclr[,colnames(imilk_invr_RAclr) %in% v_static_imilk_invr_RAclr]

# combine in 1 data frame that is to be used for models
imilk_invr_RAclr_static <- merge(basic_imilk_invr_RAclr, pheno_imilk_invr_RAclr_static, by="row.names")
imilk_invr_RAclr_static <- imilk_invr_RAclr_static[,-1] #remove Row.names column
#737x63, of which 8 are phenotypes
#str(imilk_invr_RAclr_static)

## Infant data set (imilk_invr_RAclr_static)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56:63  = infant static phenotypes


### ===== 16.3 FUNCTION FOR ASSOCIATION OF STATIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR DNA ISOLATION BATCH, SEQUENCING READS, TIME AND ID ===== ###

library(mmrm)

##### FUNCTIONS (version 4.7 from 23 January 2024, from Alex) ##########

.run_mmrm_single = function(bacteria,phenotype,time,time_cat,NEXT_ID,covariates) {
  if(class(bacteria) !="numeric") stop("bacteria should be numeric vector!")
  if(class(time)!="numeric") stop("time should be 'double'!")
  if(class(time_cat)!='factor') stop("time_cat should be 'factor'")
  if(class(NEXT_ID)!='factor') stop("NEXT_ID should be 'factor'")
  
  data.fit = try(data.frame(bac = bacteria,trait = phenotype,Time = time,Time_cat = time_cat,NEXT_ID = NEXT_ID,covariates))
  if(class(data.fit)[1]!="data.frame") stop("check that your input data is syncronized, i.e. all inputs have the same samples in the same order")
  data.fit = data.fit[complete.cases(data.fit),]
  
  if(length(table(as.character(data.fit$Time_cat))) > 3) {
    formula.null.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + Time3 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + Time3 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + Time3 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + Time3 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
    
    data.fit2 = data.frame(data.fit, poly(data.fit$Time,3))
    colnames(data.fit2)[(ncol(data.fit2)-2): ncol(data.fit2)] = c("Time1","Time2","Time3")
  } else if (length(table(as.character(data.fit$Time_cat))) > 2) {
    formula.null.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + Time2 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + Time2 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
    
    data.fit2 = data.frame(data.fit, poly(data.fit$Time,2))
    colnames(data.fit2)[(ncol(data.fit2)-1): ncol(data.fit2)] = c("Time1","Time2")
  } else if (length(table(as.character(data.fit$Time_cat))) > 1){
    formula.null.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + trait + us(Time_cat|NEXT_ID)")
    formula.run.us = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + trait + trait:Time1 + us(Time_cat|NEXT_ID)")
    formula.null.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                            "+ Time1 + trait + cs(Time_cat|NEXT_ID)")
    formula.run.cs = paste("bac ~ ", paste(colnames(covariates),collapse = "+"),
                           "+ Time1 + trait + trait:Time1 + cs(Time_cat|NEXT_ID)")
    
    data.fit2 = data.frame(data.fit, poly(data.fit$Time,1))
    colnames(data.fit2)[ncol(data.fit2)] = c("Time1")
  } else {
    stop ("Data has only one timepoint! Swtich to lm() to run the models")
  }
  #run ML models with US matrix
  nullModel.reml = try(mmrm(as.formula(formula.null.us), data= data.fit2,reml = T))
  runModel.reml = try(mmrm(as.formula(formula.run.us), data= data.fit2,reml = T))
  
  if(class(runModel.reml)[1] != "try-error") summary.runModel.reml <-try(summary(runModel.reml))
  if(class(nullModel.reml)[1] != "try-error") summary.nullModel.reml <-try(summary(nullModel.reml))
  if(!exists("summary.runModel.reml")) {
    summary.runModel.reml = NA
    class(summary.runModel.reml) = "try-error"
  }
  if(!exists("summary.nullModel.reml")) {
    summary.nullModel.reml = NA
    class(summary.nullModel.reml) = "try-error"
  }
  
  if (class(nullModel.reml)[1] == "try-error"|
      class(runModel.reml)[1] == "try-error"|
      class(summary.nullModel.reml)[1] == "try-error"|
      class(summary.runModel.reml)[1] == "try-error") {
    runModel.reml = try(mmrm(as.formula(formula.run.cs), data= data.fit2,reml = T))
    nullModel.reml = try(mmrm(as.formula(formula.null.cs), data= data.fit2,reml = T))
    if(class(runModel.reml)[1] != "try-error") summary.runModel.reml <-try(summary(runModel.reml))
    if(class(nullModel.reml)[1] != "try-error") summary.nullModel.reml <-try(summary(nullModel.reml))
    model.type = "CompoundSymmetry"
  } else {model.type = "Unstructured"}
  
  if(class(nullModel.reml)[1]=="try-error"|
     class(runModel.reml)[1]=="try-error") {
    trait.results = data.frame(row.names = NULL,
                               type = "Failure",
                               Covar.type = model.type,
                               bac = NA,
                               trait = NA,
                               N = nrow(data.fit2),
                               levels = NA,
                               "Estimate" = NA,
                               "Std. Error" = NA,
                               "df" = NA,
                               "t value"= NA,
                               "Pr(>|t|)" = NA 
    )
    timeInt.results = data.frame(row.names = NULL,
                                 type = "Failure",
                                 Covar.type = model.type,
                                 bac = NA,
                                 trait = NA,
                                 N = nrow(data.fit2),
                                 levels = NA,
                                 "Estimate" = NA,
                                 "Std. Error" = NA,
                                 "df" = NA,
                                 "t value"= NA,
                                 "Pr(>|t|)" = NA
    )
  } else {
    trait.results = data.frame(row.names = NULL,
                               type = "Success",
                               Covar.type = model.type,
                               bac = NA,
                               trait = NA,
                               N = nrow(data.fit2),
                               levels = sub("trait",
                                            "",
                                            rownames(summary.nullModel.reml$coef)[
                                              grep("trait",rownames(summary.nullModel.reml$coef))]),
                               summary.nullModel.reml$coef[
                                 grep("trait",rownames(summary.nullModel.reml$coef)),,drop =F])
    timeInt.results = data.frame(row.names = NULL,
                                 type = "Success",
                                 Covar.type = model.type,
                                 bac = NA,
                                 trait = NA,
                                 N = nrow(data.fit2),
                                 levels = sub("trait",
                                              "",
                                              rownames(summary.runModel.reml$coef)[
                                                grep(":",rownames(summary.runModel.reml$coef))]),
                                 summary.runModel.reml$coef[grep(":",rownames(summary.runModel.reml$coef)),,drop = F])
  }
  
  list(trait = trait.results,
       time.trait = timeInt.results)
}

run_mmrm_analysis = function( bacteria, phenotypes, time, time_cat, NEXT_ID,covariates){
  output_names = c("trait","time.trait")
  output = sapply(output_names, function(x) NULL)
  for (i in 1:ncol(bacteria)) {
    print (i)
    for(j in 1:ncol(phenotypes)){
      print(j)
      singleAssoc = .run_mmrm_single(bacteria[,i],phenotypes[,j],time = time,time_cat = time_cat, NEXT_ID = NEXT_ID,covariates = covariates)
      singleAssoc$trait[,"bac"] = colnames(bacteria)[i]
      singleAssoc$trait[,"trait"] = colnames(phenotypes)[j]
      singleAssoc$time.trait[,"bac"] = colnames(bacteria)[i]
      singleAssoc$time.trait[,"trait"] = colnames(phenotypes)[j]
      output$trait = rbind(output$trait,singleAssoc$trait)
      output$time.trait = rbind(output$time.trait,singleAssoc$time.trait)
    }
    
  }
  colnames(output$trait)[8] = "SE"
  colnames(output$time.trait)[8] = "SE"
  colnames(output$time.trait)[11] = "P"
  colnames(output$trait)[11] = "P"
  
  
  output$trait$FDR = p.adjust(output$trait$P)
  output
  
}


### ===== 16.4 RUN ASSOCIATION OF MATERNAL STATIC PHENOTYPES WITH MILK ALPHA DIVERSITY MEASURES AND CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES ===== ###

## Maternal data set (mmilk_invr_RAclr_static)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56:441 = maternal static phenotypes

mmilk_invr_RAclr_static_results <- run_mmrm_analysis(bacteria   = mmilk_invr_RAclr_static[,c(17:19,20:55)],  #invr-transformed milk alpha diversity measures and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
                                                     phenotypes = mmilk_invr_RAclr_static[,c(56:441)],       #maternal static phenotypes 56:438
                                                     time       = mmilk_invr_RAclr_static[,9],  #time_point_numeric
                                                     time_cat   = mmilk_invr_RAclr_static[,8],  #time_point (factor)
                                                     NEXT_ID    = mmilk_invr_RAclr_static[,6],  #mother_ID (do not use family_ID here - doesn't work)
                                                     covariates = mmilk_invr_RAclr_static[,15:16,drop = F]) #correct for DNA_isolation_batch and seq_16S_n_reads_clean; try removing drop=F if it doesn't work
mmilk_invr_RAclr_static_results_trait <- mmilk_invr_RAclr_static_results$trait #15483x12
mmilk_invr_RAclr_static_results_time_trait <- mmilk_invr_RAclr_static_results$time.trait #15483x11
# write.table(mmilk_invr_RAclr_static_results_trait, file="association_model_results/240130_milk_microbiota_maternal_static_phenotypes_model_results_by_level.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(mmilk_invr_RAclr_static_results_time_trait, file="association_model_results/240130_milk_microbiota_maternal_static_phenotypes_model_results_by_level_and_time.txt", row.names=F, col.names=T, sep="\t", quote=F)
# mmilk_invr_RAclr_static_results_trait <- read.table(file="association_model_results/240130_milk_microbiota_maternal_static_phenotypes_model_results_by_level.txt", header=T, sep="\t", stringsAsFactors=T)
# mmilk_invr_RAclr_static_results_time_trait <- read.table(file="association_model_results/240130_milk_microbiota_maternal_static_phenotypes_model_results_by_level_and_time.txt", header=T, sep="\t", stringsAsFactors=T)

mmilk_invr_RAclr_static_results_trait$dataset <- "mother_milk_microbiota_static"


### ===== 16.5 RUN ASSOCIATION OF INFANT STATIC PHENOTYPES WITH MILK ALPHA DIVERSITY MEASURES AND CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES ===== ###

## Infant data set (imilk_invr_RAclr_static)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56:63  = infant static phenotypes

imilk_invr_RAclr_static_results <- run_mmrm_analysis(bacteria   = imilk_invr_RAclr_static[,c(17:19,20:55)], #invr-transformed milk alpha diversity measures and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
                                                     phenotypes = imilk_invr_RAclr_static[,c(56:63)], #infant static phenotypes 56:63
                                                     time       = imilk_invr_RAclr_static[,9],  #time_point_numeric
                                                     time_cat   = imilk_invr_RAclr_static[,8],  #time_point (factor)
                                                     NEXT_ID    = imilk_invr_RAclr_static[,7],  #infant_ID (do not use family_ID here - doesn't work)
                                                     covariates = imilk_invr_RAclr_static[,15:16,drop = F]) #correct for DNA_isolation_batch and seq_16S_n_reads_clean; try removing drop=F if it doesnt work
imilk_invr_RAclr_static_results_trait <- imilk_invr_RAclr_static_results$trait #585x12
imilk_invr_RAclr_static_results_time_trait <- imilk_invr_RAclr_static_results$time.trait #585x11 (no FDR column compared to trait table)
# write.table(imilk_invr_RAclr_static_results_trait, file="association_model_results/240130_milk_microbiota_infant_static_phenotypes_model_results_by_level.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(imilk_invr_RAclr_static_results_time_trait, file="association_model_results/240130_milk_microbiota_infant_static_phenotypes_model_results_by_level_and_time.txt", row.names=F, col.names=T, sep="\t", quote=F)
# imilk_invr_RAclr_static_results_trait <- read.table(file="association_model_results/240130_milk_microbiota_infant_static_phenotypes_model_results_by_level.txt", header=T, sep="\t", stringsAsFactors=T)
# imilk_invr_RAclr_static_results_time_trait <- read.table(file="association_model_results/240130_milk_microbiota_infant_static_phenotypes_model_results_by_level_and_time.txt", header=T, sep="\t", stringsAsFactors=T)

imilk_invr_RAclr_static_results_trait$dataset <- "infant_milk_microbiota_static"


##### =========================== 17. MODEL FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR SEQUENCING READS =========================== #####

### ===== 17.1 PREPARE DATA FRAMES WITH MATERNAL 1-TIME-POINT PHENOTYPES ===== ###

## Create data frame from mmilk_invr_RAclr for use in association models for infant static phenotypes
rownames(mmilk_invr_RAclr) <- mmilk_invr_RAclr$mother_sample_ID

# selected basic info
basic_mmilk_invr_RAclr <- mmilk_invr_RAclr[,c(1:9,11:15,466,481,483:485,486:521)] #IDs, DNA isolation batch, clean sequencing reads, invr-transformed milk alpha diversity measures and clr-transformed RA

# select maternal 1-time-point phenotypes that should be included -> 6 phenotypes
v_1tp_mmilk_invr_RAclr <- overview_mmilk[overview_mmilk$participant=="mother" & overview_mmilk$model_use=="1_time_point" & overview_mmilk$model_exclusion=="include","column_name"]
pheno_mmilk_invr_RAclr_1tp <- mmilk_invr_RAclr[,colnames(mmilk_invr_RAclr) %in% v_1tp_mmilk_invr_RAclr]

# combine in 1 data frame that is to be used for models
mmilk_invr_RAclr_1tp <- merge(basic_mmilk_invr_RAclr, pheno_mmilk_invr_RAclr_1tp, by="row.names")
mmilk_invr_RAclr_1tp <- mmilk_invr_RAclr_1tp[,-1] #remove Row.names column
#725x61, of which 6 are phenotypes
#str(mmilk_invr_RAclr_1tp)

## Maternal data set (mmilk_invr_RAclr_1tp)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56:61  = maternal 1-time-point phenotypes


### ===== 17.2 PREPARE DATA FRAMES WITH INFANT 1-TIME-POINT PHENOTYPES ===== ###

## Create data frame from imilk_invr_RAclr for use in association models for infant static phenotypes
rownames(imilk_invr_RAclr) <- imilk_invr_RAclr$infant_sample_ID

# selected basic info
basic_imilk_invr_RAclr <- imilk_invr_RAclr[,c(1:9,11:15,68,83,85:87,88:123)] #IDs, DNA isolation batch, clean sequencing reads, invr-transformed milk alpha diversity measures and clr-transformed RA

# select maternal 1-time-point phenotypes that should be included -> 6 phenotypes
v_1tp_imilk_invr_RAclr <- overview_imilk[overview_imilk$participant=="infant" & overview_imilk$model_use=="1_time_point" & overview_imilk$model_exclusion=="include","column_name"]
pheno_imilk_invr_RAclr_1tp <- imilk_invr_RAclr[,colnames(imilk_invr_RAclr) %in% v_1tp_imilk_invr_RAclr]

# combine in 1 data frame that is to be used for models
imilk_invr_RAclr_1tp <- merge(basic_imilk_invr_RAclr, pheno_imilk_invr_RAclr_1tp, by="row.names")
imilk_invr_RAclr_1tp <- imilk_invr_RAclr_1tp[,-1] #remove Row.names column
# #737x58, of which 3 are phenotypes
# #str(imilk_invr_RAclr_1tp)

## Infant data set (imilk_invr_RAclr_1tp)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56:58  = infant 1-time-point phenotype(s)


### ===== 17.3 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR DNA ISOLATION BATCH AND SEQUENCING READS ===== ###

run.lmer.cor.batch.reads.1tp <- function(datasetname, inputdata, xcolumn, ycolumn){
  p_models <- c()
  
  levels <- c()
  # n_levels <- c()
  e <- c()
  p <- c()
  
  d <- c()
  x <- c()
  y <- c()
  n <- c()
  stat <- c()
  
  for (i in xcolumn){
    
    for (j in ycolumn){
      print(paste0("i: ",i,"j: ",j))
      
      #print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #mixed model with correction for milk_group as fixed effects and WITH the phenotype of interest (i)
      m1 <- lm(my_df[,j] ~ DNA_isolation_batch + seq_16S_n_reads_clean + my_df[,i], data=my_df)
      sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
      sum1$rows <- rownames(sum1)
      
      levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, rep(paste0("lm_1_time_point_correction_DNAbatch_reads"), length(sum1[grep("my_df",sum1$rows),"rows"])))
      d <- c(d, rep(datasetname, length(sum1[grep("my_df",sum1$rows),"rows"])))
      x <- c(x, rep(colnames(inputdata)[i], length(sum1[grep("my_df",sum1$rows),"rows"])))
      y <- c(y, rep(colnames(inputdata)[j], length(sum1[grep("my_df",sum1$rows),"rows"])))
      n <- c(n, rep(nrow(my_df), length(sum1[grep("my_df",sum1$rows),"rows"])))
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, estimate=e, p=p)
  
  #remove "my_df[, i]" from levels
  res$levels <- substring(res$levels, 11)
  
  return(res)
}


### ===== 17.4 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR SEQUENCING READS ===== ###

run.lmer.cor.reads.1tp <- function(datasetname, inputdata, xcolumn, ycolumn){
  p_models <- c()
  
  levels <- c()
  e <- c()
  p <- c()
  
  d <- c()
  x <- c()
  y <- c()
  n <- c()
  stat <- c()
  
  for (i in xcolumn){
    
    for (j in ycolumn){
      print(paste0("i: ",i,"j: ",j))
      
      #print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #mixed model with correction for milk_group as fixed effects and WITH the phenotype of interest (i)
      m1 <- lm(my_df[,j] ~ seq_16S_n_reads_clean + my_df[,i], data=my_df)
      sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
      sum1$rows <- rownames(sum1)
      
      levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, rep(paste0("lm_1_time_point_correction_reads"), length(sum1[grep("my_df",sum1$rows),"rows"])))
      d <- c(d, rep(datasetname, length(sum1[grep("my_df",sum1$rows),"rows"])))
      x <- c(x, rep(colnames(inputdata)[i], length(sum1[grep("my_df",sum1$rows),"rows"])))
      y <- c(y, rep(colnames(inputdata)[j], length(sum1[grep("my_df",sum1$rows),"rows"])))
      n <- c(n, rep(nrow(my_df), length(sum1[grep("my_df",sum1$rows),"rows"])))
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, estimate=e, p=p)
  
  #remove "my_df[, i]" from levels
  res$levels <- substring(res$levels, 11)
  
  return(res)
}


### ===== 17.5 RUN ASSOCIATION OF MATERNAL SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA ===== ###

mmilk_1tp_results <- run.lmer.cor.batch.reads.1tp(datasetname = "mother_milk_microbiota_1_time_point",
                                                  inputdata   = mmilk_invr_RAclr_1tp,
                                                  xcolumn     = c(56:61),   #single time point maternal phenotypes
                                                  ycolumn     = c(17:55))   #invr-transformed milk alpha diversity measures and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
mmilk_1tp_results$Note <- NA
# write.table(mmilk_1tp_results, file="association_model_results/240130_milk_microbiota_maternal_1_time_point_phenotypes_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# mmilk_1tp_results <- read.table(file="association_model_results/240130_milk_microbiota_maternal_1_time_point_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 17.6 RUN ASSOCIATION OF INFANT SINGLE-TIME-POINT STATIC PHENOTYPES WITH MILK MICROBIOTA DATA ===== ###

## first run infant phenotypes that were present in both DNA isolation batches #fix function to allow 3 levels for phenotype 56!
imilk_1tp_results_1 <- run.lmer.cor.batch.reads.1tp(datasetname = "infant_milk_microbiota_1_time_point",
                                                  inputdata   = imilk_invr_RAclr_1tp,
                                                  xcolumn     = c(56:57),   #single time point maternal phenotypes
                                                  ycolumn     = c(17:55))   #invr-transformed milk alpha diversity measures and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
imilk_1tp_results_1$Note <- NA

## ! Note: this infant phenotype was only present for samples from 1 batch, so the correction for batch effect was removed here.
imilk_1tp_results_2 <- run.lmer.cor.reads.1tp(datasetname = "infant_milk_microbiota_1_time_point",
                                            inputdata   = imilk_invr_RAclr_1tp,
                                            xcolumn     = c(58),      #single time point infant phenotypes
                                            ycolumn     = c(17:55))   #invr-transformed milk alpha diversity measures and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
imilk_1tp_results_2$Note <- "No_correction_for_DNA_isolation_batch_applied_as_all_individuals_who_had_data_for_this_phenotype_were_from_the_same_DNA_isolation_batch"

## combine results for infant phenotypes
imilk_1tp_results <- as.data.frame(rbind(imilk_1tp_results_1, imilk_1tp_results_2))

# write.table(imilk_1tp_results, file="association_model_results/240130_milk_microbiota_infant_1_time_point_phenotypes_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# imilk_1tp_results <- read.table(file="association_model_results/240130_milk_microbiota_infant_1_time_point_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


##### =========================== 18. MODEL FOR ASSOCIATION OF MULTIPLE-TIME-POINT DYNAMIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR DNA ISOLATION BATCH, SEQUENCING READS, TIME AND ID =========================== #####

### ===== 18.1 PREPARE DATA FRAMES WITH MATERNAL DYNAMIC PHENOTYPES ===== ###

# selected basic info
basic_mmilk_invr_RAclr <- mmilk_invr_RAclr[,c(1:9,11:15,466,481,483:485,486:521)] #IDs, DNA isolation batch, clean sequencing reads, invr-transformed milk alpha diversity measures and clr-transformed RA

# select maternal dynamic phenotypes that should be included
v_dynamic_mmilk_invr_RAclr <- overview_mmilk[overview_mmilk$participant=="mother" & overview_mmilk$model_use=="dynamic" & overview_mmilk$model_exclusion=="include","column_name"]
pheno_mmilk_invr_RAclr_dynamic <- mmilk_invr_RAclr[,colnames(mmilk_invr_RAclr) %in% v_dynamic_mmilk_invr_RAclr]

# combine in 1 data frame that is to be used for models
mmilk_invr_RAclr_dynamic <- merge(basic_mmilk_invr_RAclr, pheno_mmilk_invr_RAclr_dynamic, by="row.names")
mmilk_invr_RAclr_dynamic <- mmilk_invr_RAclr_dynamic[,-1] #remove Row.names column
#725x100, of which 45 are phenotypes
#str(mmilk_invr_RAclr_dynamic)

## Maternal data set (mmilk_invr_RAclr_dynamic)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56:100  = maternal dynamic phenotypes


### ===== 18.2 PREPARE DATA FRAMES WITH INFANT DYNAMIC PHENOTYPES ===== ###

# select infant dynamic phenotypes that should be included: only 1 phenpotype: infant_ffq_feeding_mode
v_dynamic_imilk_invr_RAclr <- overview_imilk[overview_imilk$participant=="infant" & overview_imilk$model_use=="dynamic" & overview_imilk$model_exclusion=="include","column_name"]

# selected basic info and the 1 dynamic infant phenotype (col 27)
imilk_invr_RAclr_dynamic <- imilk_invr_RAclr[,c(1:9,11:15,68,83,85:87,88:123,27)] #IDs, DNA isolation batch, clean sequencing reads, invr-transformed milk alpha diversity measures and clr-transformed RA
#737x56, of which 1 is a phenotype
#str(imilk_invr_RAclr_dynamic)

## Infant data set (imilk_invr_RAclr_dynamic)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56     = infant dynamic phenotype(s)


### ===== 18.3 FUNCTION FOR ASSOCIATION OF DYNAMIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR DNA ISOLATION BATCH, SEQUENCING READS, TIME AND ID ===== ###

run.lmer.cor.DNAbatch.reads.time.ID <- function(datasetname, inputdata, xcolumn, ycolumn){
  
  levels <- c()
  e <- c()
  p <- c()
  
  d <- c()
  x <- c()
  y <- c()
  n <- c()
  stat <- c()
  
  for (i in xcolumn){
    
    for (j in ycolumn){
      print(paste0("i: ",i,"j: ",j))
      
      #print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #run mixed model with correction for milk_group and time_point as fixed effects and WITH the phenotype of interest (i)
      m1 <- lmerTest::lmer(my_df[,j] ~ DNA_isolation_batch + seq_16S_n_reads_clean + time_point_numeric + (1|mother_ID) + my_df[,i], REML=F, data=my_df)
      sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
      sum1$rows <- rownames(sum1)
      # print(sum1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, rep(paste0("lmer_correction_DNAbatch_reads_time_ID"), length(sum1[grep("my_df",sum1$rows),"rows"])))
      d <- c(d, rep(datasetname, length(sum1[grep("my_df",sum1$rows),"rows"])))
      x <- c(x, rep(colnames(inputdata)[i], length(sum1[grep("my_df",sum1$rows),"rows"])))
      y <- c(y, rep(colnames(inputdata)[j], length(sum1[grep("my_df",sum1$rows),"rows"])))
      n <- c(n, rep(nrow(my_df), length(sum1[grep("my_df",sum1$rows),"rows"])))
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, estimate=e, p=p)
  
  #remove "my_df[, i]" from levels
  res$levels <- substring(res$levels, 11)
  
  return(res)
}


### ===== 18.4 FUNCTION FOR ASSOCIATION OF DYNAMIC PHENOTYPES WITH MILK MICROBIOTA DATA WITH CORRECTION FOR SEQUENCING READS, TIME AND ID ===== ###

run.lmer.cor.reads.time.ID <- function(datasetname, inputdata, xcolumn, ycolumn){
  
  levels <- c()
  e <- c()
  p <- c()
  
  d <- c()
  x <- c()
  y <- c()
  n <- c()
  stat <- c()
  
  for (i in xcolumn){
    
    for (j in ycolumn){
      print(paste0("i: ",i,"j: ",j))
      
      #print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #run mixed model with correction for milk_group and time_point as fixed effects and WITH the phenotype of interest (i)
      m1 <- lmerTest::lmer(my_df[,j] ~ seq_16S_n_reads_clean + time_point_numeric + (1|mother_ID) + my_df[,i], REML=F, data=my_df)
      sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
      sum1$rows <- rownames(sum1)
      # print(sum1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for covariates and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, rep(paste0("lmer_correction_reads_time_ID"), length(sum1[grep("my_df",sum1$rows),"rows"])))
      d <- c(d, rep(datasetname, length(sum1[grep("my_df",sum1$rows),"rows"])))
      x <- c(x, rep(colnames(inputdata)[i], length(sum1[grep("my_df",sum1$rows),"rows"])))
      y <- c(y, rep(colnames(inputdata)[j], length(sum1[grep("my_df",sum1$rows),"rows"])))
      n <- c(n, rep(nrow(my_df), length(sum1[grep("my_df",sum1$rows),"rows"])))
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, estimate=e, p=p)
  
  #remove "my_df[, i]" from levels
  res$levels <- substring(res$levels, 11)
  
  return(res)
}


### ===== 18.5 RUN ASSOCIATION OF MATERNAL DYNAMIC PHENOTYPES WITH WITH MILK ALPHA DIVERSITY MEASURES AND CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES ===== ###

## Maternal data set (mmilk_invr_RAclr_dynamic)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56:100  = maternal dynamic phenotypes

## run association for maternal dynamic phenotypes that have data in both DNA isolation batches
mmilk_dynamic_results_1 <- run.lmer.cor.DNAbatch.reads.time.ID(datasetname = "mother_milk_microbiota_dynamic",
                                                               inputdata   = mmilk_invr_RAclr_dynamic,
                                                               xcolumn     = c(56:70,73:100),        #dynamic maternal phenotypes (with correction for DNA isolation batch)
                                                               ycolumn     = c(17:19,20:55))   #invr-transformed milk alpha diversity measures and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
mmilk_dynamic_results_1$Note <- NA

## run association for maternal dynamic phenotypes that have data in only 1 DNA isolation batch
mmilk_dynamic_results_2 <- run.lmer.cor.reads.time.ID(datasetname = "mother_milk_microbiota_dynamic",
                                                      inputdata   = mmilk_invr_RAclr_dynamic,
                                                      xcolumn     = c(71:72),         #dynamic maternal phenotypes (WITHOUT correction for DNA isolation batch)
                                                      ycolumn     = c(17:19,20:55))   #invr-transformed milk alpha diversity measures and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
mmilk_dynamic_results_2$Note <- "No_correction_for_DNA_isolation_batch_applied_as_all_individuals_who_had_data_for_this_phenotype_were_from_the_same_DNA_isolation_batch"

## combine association results
mmilk_dynamic_results <- as.data.frame(rbind(mmilk_dynamic_results_1, mmilk_dynamic_results_2))
# write.table(mmilk_dynamic_results, file="association_model_results/240130_milk_microbiota_maternal_dynamic_phenotypes_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# mmilk_dynamic_results <- read.table(file="association_model_results/240130_milk_microbiota_maternal_dynamic_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 18.6 RUN ASSOCIATION OF INFANT DYNAMIC PHENOTYPES WITH WITH MILK ALPHA DIVERSITY MEASURES AND CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES ===== ###

## Infant data set (imilk_invr_RAclr_dynamic)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15:16  = DNA_isolation_batch and seq_16S_n_reads_clean
## 17:19  = invr-transformed milk alpha diversity measures
## 20:55  = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
## 56     = infant dynamic phenotypes

imilk_dynamic_results <- run.lmer.cor.DNAbatch.reads.time.ID(datasetname = "infant_milk_microbiota_dynamic",
                                                             inputdata   = imilk_invr_RAclr_dynamic,
                                                             xcolumn     = c(56),           #dynamic infant phenotypes
                                                             ycolumn     = c(17:19,20:55))  #invr-transformed milk alpha diversity measures and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
imilk_dynamic_results$Note <- NA
# write.table(imilk_dynamic_results, file="association_model_results/240130_milk_microbiota_infant_dynamic_phenotypes_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# imilk_dynamic_results <- read.table(file="association_model_results/240130_milk_microbiota_infant_dynamic_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


##### =========================== 19. COMBINE RESULTS AND CORRECT FOR MULTIPLE TESTING =========================== #####

### ===== 19.1 AVAILABLE RESULTS DATA FRAMES FOR MERGING ===== ###

## Results data frames                      Size
# mmilk_invr_RAclr_static_results_trait 15483 x 13
# imilk_invr_RAclr_static_results_trait   585 x 13
# mmilk_1tp_results                       234 x  9
# imilk_1tp_results                       156 x  9
# mmilk_dynamic_results                  2301 x  9
# imilk_dynamic_results                    39 x  9

## if necessary, import these data frames:

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/")

mmilk_invr_RAclr_static_results_trait <- read.table(file="association_model_results/240130_milk_microbiota_maternal_static_phenotypes_model_results_by_level.txt", header=T, sep="\t", stringsAsFactors=T)
mmilk_invr_RAclr_static_results_trait$dataset <- "mother_milk_microbiota_static"
# mmilk_invr_RAclr_static_results_time_trait <- read.table(file="association_model_results/240130_milk_microbiota_maternal_static_phenotypes_model_results_by_level_and_time.txt", header=T, sep="\t", stringsAsFactors=T)

imilk_invr_RAclr_static_results_trait <- read.table(file="association_model_results/240130_milk_microbiota_infant_static_phenotypes_model_results_by_level.txt", header=T, sep="\t", stringsAsFactors=T)
imilk_invr_RAclr_static_results_trait$dataset <- "infant_milk_microbiota_static"
# imilk_invr_RAclr_static_results_time_trait <- read.table(file="association_model_results/240130_milk_microbiota_infant_static_phenotypes_model_results_by_level_and_time.txt", header=T, sep="\t", stringsAsFactors=T)

mmilk_1tp_results <- read.table(file="association_model_results/240130_milk_microbiota_maternal_1_time_point_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)
imilk_1tp_results <- read.table(file="association_model_results/240130_milk_microbiota_infant_1_time_point_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)

mmilk_dynamic_results <- read.table(file="association_model_results/240130_milk_microbiota_maternal_dynamic_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)
imilk_dynamic_results <- read.table(file="association_model_results/240130_milk_microbiota_infant_dynamic_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 19.2 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM STATIC PHENOTYPE ASSOCIATIONS ===== ###

results_static <- as.data.frame(rbind(mmilk_invr_RAclr_static_results_trait, imilk_invr_RAclr_static_results_trait))

## select relevant rows and columns
# exclude results from phenotype that could not always run (model could not be built)
# table(results_static[is.na(results_static$P),"trait"])
results_static_1 <- results_static[results_static$trait!="mother_weight_gain_preg_kg_invr",] #this excludes 39 rows
results_static_2 <- results_static_1[!is.na(results_static_1$P),] #this excludes 117 rows (39 each for mother_ffq_item_spreads_yeast_extract_grams_per_day_invr, mother_ffq_item_sugar_in_tea_grams_per_day_invr and mother_ffq_item_wine_fort_ml_per_day_invr)

# only select relevant columns
results_static_3 <- results_static_2[,c(2:7,11,13)]

## fix column names
colnames(results_static_3) <- c("statistic", "y", "x", "n_total", "levels", "estimate", "p", "dataset")

## add missing columns
results_static_3$Note <- NA

## resort columns
sel_results_static <- results_static_3[,c(8,1,3,2,4:7,9)] #15912x9
sel_results_static$Note <- as.factor(sel_results_static$Note)


### ===== 19.3 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM 1-TIME-POINT PHENOTYPE ASSOCIATIONS ===== ###

results_1tp <- as.data.frame(rbind(mmilk_1tp_results, imilk_1tp_results))
sel_results_1tp <- results_1tp #390x9


### ===== 19.4 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM DYNAMIC PHENOTYPE ASSOCIATIONS ===== ###

results_dynamic <- as.data.frame(rbind(mmilk_dynamic_results, imilk_dynamic_results))
sel_results_dynamic <- results_dynamic #2340x9


### ===== 19.5 MERGE ALL RESULTS DATA FRAMES ===== ###

## merge all results
# library(plyr)
model_results <- as.data.frame(rbind.fill(sel_results_static, sel_results_1tp, sel_results_dynamic)) #18642x9


### ===== 19.6 EXCLUDE PHENOTYPES ===== ###

## exclude all food preferences and all individual food items
model_results_short_1 <- model_results[-grep("mother_food_pref", model_results$x),] #exclude 85 food preferences x39 (this excludes 3315 rows)
model_results_short_2 <- model_results_short_1[-grep("mother_ffq_item", model_results_short_1$x),] #exclude 160 individual food items x39 (this excludes 6240 rows)

## exclude selected detailed phenotypes (n=14) and mother_milk_collection_volume_pumped_invr as it is highly correlated with time point postpartum
exc_phenos <- c("mother_genetics_blood_group_genotype",
                "mother_ffq_fish_seafood_consumption_yes_no", "mother_ffq_meat_consumption_yes_no",
                "mother_exp_pets_cat", "mother_exp_pets_dog",
                "mother_birth_gestational_age_categories",
                "mother_milk_HMO_Le", "mother_milk_HMO_Se",
                "infant_genetics_FUT2_secretor", "infant_Se_status_same_as_mother",
                "infant_birth_delivery_place_detailed", "infant_birth_delivery_mode_detailed",
                "mother_milk_collection_month",
                "mother_dis_subclinical_mastitis_Na_K_ratio_invr","mother_milk_collection_volume_pumped_invr")

model_results_short_3 <- model_results_short_2[!(model_results_short_2$x %in% exc_phenos),] #exclude similar or detailed phenotypes (this excludes 1365 rows) and phenotypes correlated with time postpartum (this excludes 39 rows)

model_results_short <- model_results_short_3 ## 7683x9


### ===== 19.7 ADD PHENOTYPE GROUPS FOR PLOTTING ===== ###

## show all included unique phentypes
# length(unique(model_results_short$x)) #185
# unique(model_results_short$x)

## add phenotype groups to cluster results easily
model_results_short$phenotype_group <- NA

model_results_short[model_results_short$x=="mother_type_pregnancy", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"
model_results_short[model_results_short$x=="mother_blood_group", "phenotype_group"] <- "Maternal_genetics"
model_results_short[model_results_short$x=="mother_age_birth_invr", "phenotype_group"] <- "Maternal_and_infant_anthropometrics"                                                 
model_results_short[model_results_short$x=="mother_previous_preg", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"                                                  
model_results_short[model_results_short$x=="mother_gravida_invr", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"                                                   
model_results_short[model_results_short$x=="mother_para_invr", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"                                                      
model_results_short[model_results_short$x=="mother_height_cm_invr", "phenotype_group"] <- "Maternal_and_infant_anthropometrics"                                                
model_results_short[model_results_short$x=="mother_weight_prior_preg_kg_invr", "phenotype_group"] <- "Maternal_and_infant_anthropometrics"            
model_results_short[model_results_short$x=="mother_weight_early_preg_kg_invr", "phenotype_group"] <- "Maternal_and_infant_anthropometrics"                                      
model_results_short[model_results_short$x=="mother_BMI_early_preg_invr", "phenotype_group"] <- "Maternal_and_infant_anthropometrics"
model_results_short[model_results_short$x=="mother_weight_just_before_delivery_kg_invr", "phenotype_group"] <- "Maternal_and_infant_anthropometrics"
model_results_short[model_results_short$x=="mother_food_alcohol_during_first_trimester", "phenotype_group"] <- "Maternal_diet"
model_results_short[model_results_short$x=="mother_food_avoid_dairy_during_preg", "phenotype_group"] <- "Maternal_diet"                                   
model_results_short[model_results_short$x=="mother_food_avoid_gluten_during_preg", "phenotype_group"] <- "Maternal_diet"
model_results_short[grep("mother_ffq_group", model_results_short$x), "phenotype_group"] <- "Maternal_diet"
model_results_short[model_results_short$x=="mother_ffq_energy_kcal_invr", "phenotype_group"] <- "Maternal_diet"
model_results_short[grep("mother_ffq_nutrients", model_results_short$x), "phenotype_group"] <- "Maternal_diet"
model_results_short[model_results_short$x=="mother_ffq_adapted_Mediterranean_diet_score_invr", "phenotype_group"] <- "Maternal_diet"
model_results_short[model_results_short$x=="mother_stress_mean_pregnancy_stress_invr", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="mother_stress_mean_postpartum_stress_invr", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="mother_education_level_numeric_invr", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="mother_net_income_numeric_invr", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[grep("mother_exp", model_results_short$x), "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="mother_health_prior_preg_eating_disorder", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_reprod_age_menarche_invr", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_reprod_prior_preg_regular_menstr_pattern", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_reprod_prior_preg_contraceptive_use", "phenotype_group"] <- "Maternal_medication_use"
model_results_short[model_results_short$x=="mother_health_reprod_prior_preg_contraceptive_use_hormonal", "phenotype_group"] <- "Maternal_medication_use"
model_results_short[model_results_short$x=="mother_health_reprod_prior_preg_endometriosis", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_reprod_prior_preg_PCOS", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_reprod_prior_preg_ovarian_cysts", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_reprod_prior_preg_STD", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_reprod_preg_IVF", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"
model_results_short[model_results_short$x=="mother_health_dental_numeric_invr", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_dental_bleeding_gums_during_preg", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_health_preg_hyperemesis_gravidarum_severe_nausea_early_preg", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"
model_results_short[grep("mother_preg_comp", model_results_short$x), "phenotype_group"] <- "Pregnancy_and_birth_characteristics"
model_results_short[model_results_short$x=="mother_birth_comp_placenta_birth_comp", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"
model_results_short[model_results_short$x=="mother_birth_gestational_age_weeks_invr", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"
model_results_short[model_results_short$x=="mother_birth_water_spontaneously_broken", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"
model_results_short[grep("mother_dis_diary", model_results_short$x), "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_dis_ROME3_functional_bloating", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_dis_ROME3_IBS", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[model_results_short$x=="mother_med_prior_preg_folic_acid_use", "phenotype_group"] <- "Maternal_medication_use"
model_results_short[model_results_short$x=="mother_milk_HMO_milk_group", "phenotype_group"] <- "Maternal_genetics"
model_results_short[model_results_short$x=="infant_sex", "phenotype_group"] <- "Maternal_and_infant_anthropometrics"
model_results_short[model_results_short$x=="infant_birth_delivery_place", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"                                           
model_results_short[model_results_short$x=="infant_birth_delivery_mode", "phenotype_group"] <- "Pregnancy_and_birth_characteristics"                                            
model_results_short[model_results_short$x=="infant_birthweight_g_invr", "phenotype_group"] <- "Maternal_and_infant_anthropometrics"
model_results_short[grep("mother_med", model_results_short$x), "phenotype_group"] <- "Maternal_medication_use"
model_results_short[model_results_short$x=="infant_ffq_feeding_type_delivery", "phenotype_group"] <- "Infant_feeding"
model_results_short[model_results_short$x=="infant_ffq_breastfeeding_behaviour_after_birth_1bad_5excellent_invr", "phenotype_group"] <- "Infant_feeding"
model_results_short[model_results_short$x=="infant_ffq_solid_food_introduced", "phenotype_group"] <- "Infant_feeding"
model_results_short[grep("mother_milk_collection", model_results_short$x), "phenotype_group"] <- "Human_milk_collection"
model_results_short[model_results_short$x=="mother_dis_subclinical_mastitis", "phenotype_group"] <- "Maternal_health_and_diseases"
model_results_short[grep("ugml", model_results_short$x), "phenotype_group"] <- "Human_milk_HMO_concentrations"                
model_results_short[model_results_short$x=="mother_breastpump_brand", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="mother_breastpump_times_per_day_invr", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="infant_ffq_feeding_mode", "phenotype_group"] <- "Infant_feeding"

# table(model_results_short$phenotype_group, useNA="ifany")
# length(unique(model_results_short$x)) #185
# length(unique(model_results_short$phenotype_group)) #10

phenos_info <- model_results_short[,c(3,10)]
phenos_info <- phenos_info[!duplicated(phenos_info),]
table(phenos_info$phenotype_group, useNA="ifany")
# Human_milk_collection              Human_milk_HMO_concentrations 
# 3                                  28 
# Infant_feeding                     Maternal_and_infant_anthropometrics 
# 4                                  8 
# Maternal_diet                      Maternal_genetics 
# 72                                 2 
# Maternal_health_and_diseases       Maternal_lifestyle_and_exposures 
# 21                                 11 
# Maternal_medication_use            Pregnancy_and_birth_characteristics 
# 17                                 19

unique(model_results_short[model_results_short$phenotype_group=="Human_milk_collection", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Human_milk_HMO_concentrations", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Infant_feeding", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_and_infant_anthropometrics", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_diet", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_genetics", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_health_and_diseases", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_lifestyle_and_exposures", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_medication_use", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Pregnancy_and_birth_characteristics", "x"])

write.table(phenos_info, file="association_model_results/240229_phenotypes_associated_with_milk_microbiota_real_n185.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 19.8 ADD PRETTY LABELS FOR PHENOTYPES FOR PLOTTING ===== ###

## ensure results for traits with >2 categories are correctly displayed (combine x with level)
model_results_short$trait_name_in_plot <- NA
model_results_short[model_results_short$levels!="", "trait_name_in_plot"] <- c(paste0(model_results_short[model_results_short$levels!="", "x"],
                                                                                         "_", model_results_short[model_results_short$levels!="", "levels"]))
model_results_short[is.na(model_results_short$trait_name_in_plot), "trait_name_in_plot"] <- as.character(as.factor(model_results_short[is.na(model_results_short$trait_name_in_plot), "x"]))
# length(unique(model_results_short$trait_name_in_plot)) #197

## create column with polished trait names for plotting
model_results_short$trait_name_in_plot <- gsub("_invr", "", model_results_short$trait_name_in_plot)
model_results_short$trait_name_in_plot <- as.factor(as.character(model_results_short$trait_name_in_plot))

levels(model_results_short$trait_name_in_plot) <- c("infant_birth_delivery_mode_C-section" = "delivery_mode_(C-section_vs_vaginal_birth)",
                                                            "infant_birth_delivery_place_hospital" = "delivery_place_(hospital_vs_home)",
                                                            "infant_birthweight_g" = "infant_birth_weight_(grams)",
                                                            "infant_ffq_breastfeeding_behaviour_after_birth_1bad_5excellent" = "excellent_infant_breastfeeding_behaviour_after_birth",
                                                            "infant_ffq_feeding_mode_mixed_feeding" = "feeding_mode_(mixed_feeding_vs_human_milk_feeding)",
                                                            "infant_ffq_feeding_type_delivery_breast_milk_via_bottle"  = "feeding_mode_after_birth_(breast_milk_via_bottle_vs_breastfeeding)",              
                                                            "infant_ffq_feeding_type_delivery_mixed_feeding" = "feeding_mode_after_birth_(mixed_feeding_vs_breastfeeding)",
                                                            "infant_ffq_solid_food_introduced_yes" = "solid_food_introduced",
                                                            "infant_sex_male" = "infant_sex_(male_vs_female)",
                                                            "mother_age_birth" = "maternal_age_at_birth_(years)",
                                                            "mother_birth_comp_placenta_birth_comp_yes" = "placenta_birth_complication",
                                                            "mother_birth_gestational_age_weeks" = "gestational_age_(weeks)",
                                                            "mother_birth_water_spontaneously_broken_yes" = "water_spontaneously_broken",
                                                            "mother_blood_group_A" = "maternal_blood_group_(A_vs_O)",
                                                            "mother_blood_group_AB" = "maternal_blood_group_(AB_vs_O)",
                                                            "mother_blood_group_B" = "maternal_blood_group_(B_vs_O)",
                                                            "mother_BMI_early_preg" = "maternal_BMI_during_early_pregnancy_(kg/m2)",
                                                            "mother_breastpump_brand_Ardo" = "breast_pump_brand_(Ardo_vs_Philips)",
                                                            "mother_breastpump_brand_Medela" = "breast_pump_brand_(Medela_vs_Philips)",
                                                            "mother_breastpump_brand_Other" = "breast_pump_brand_('other'_vs_Philips)",
                                                            "mother_breastpump_times_per_day" = "breast_pump_use_(times/day)",
                                                            "mother_dis_diary_abdominal_pain_mean" = "abdominal_pain",
                                                            "mother_dis_diary_bloated_mean" = "bloated",
                                                            "mother_dis_diary_burping_mean" = "burping",
                                                            
                                                            "mother_dis_diary_constipation_mean" = "constipation",
                                                            "mother_dis_diary_diarrhea_mean" = "diarrhea",
                                                            "mother_dis_diary_flatulence_mean" = "flatulence",
                                                            "mother_dis_diary_freq_mean" = "stool_frequency_(stools/day)",
                                                            "mother_dis_diary_nausea_mean" = "nausea",
                                                            "mother_dis_diary_not_well_mean" = "uneasy_abdominal_feeling",
                                                            "mother_dis_ROME3_functional_bloating_yes" = "functional_bloating_(ROME3)",
                                                            "mother_dis_ROME3_IBS_yes" = "irritable_bowel_syndrome_(ROME3)",
                                                            "mother_dis_subclinical_mastitis_yes" = "subclinical_mastitis",
                                                            "mother_education_level_numeric" = "higher_maternal_education",
                                                            "mother_exp_living_situation_farm" = "living_situation_(farm_vs_big_house)",
                                                            "mother_exp_living_situation_flat" = "living_situation_(flat_vs_big_house)",
                                                            "mother_exp_living_situation_other" = "living_situation_('other'_vs_big_house)",
                                                            "mother_exp_living_situation_small_house" = "living_situation_(small_house_vs_big_house)",
                                                            "mother_exp_pets_any_yes" = "pets",
                                                            "mother_exp_preg_any_smoking_during_preg_yes" = "maternal_smoking_during_pregnancy",
                                                            "mother_exp_preg_secondhand_smoke_during_first_trimester_yes" = "secondhand_smoking_exposure_during_first_trimester",
                                                            "mother_exp_prior_preg_smoked_one_whole_year_yes" = "maternal_smoking_for_at_least_1_year_before_pregnancy",
                                                            "mother_ffq_adapted_Mediterranean_diet_score" = "Mediterranean_diet",
                                                            
                                                            "mother_ffq_energy_kcal" = "total_energy_intake_(kcal)",
                                                            "mother_ffq_group_alcohol_ml_per_day" = "alcohol_(ml/day)",
                                                            "mother_ffq_group_cereals_grams_per_day" = "cereals_grams/day)",
                                                            "mother_ffq_group_cheese_grams_per_day" = "cheese_(grams/day",
                                                            "mother_ffq_group_coffee_ml_per_day" = "coffee_(grams/day)",
                                                            "mother_ffq_group_dairy_grams_per_day" = "dairy_(grams/day)",
                                                            "mother_ffq_group_eggs_grams_per_day" = "eggs_(grams/day)",
                                                            "mother_ffq_group_fats_and_oils_grams_per_day" = "fats_and_oils_(grams/day)",
                                                            "mother_ffq_group_fish_seafood_grams_per_day" = "fish_and_seafood_(grams/day)",
                                                            "mother_ffq_group_fruits_grams_per_day" = "fruits_(grams/day)",
                                                            "mother_ffq_group_legumes_grams_per_day" = "legumes_(grams/day)",
                                                            "mother_ffq_group_meat_grams_per_day" = "meat_(grams/day)",
                                                            "mother_ffq_group_meat_substitute_grams_per_day" = "meat_substitutes_(grams/day)",
                                                            "mother_ffq_group_nonrefined_high_fiber_breads_grams_per_day" = "not_refined_high_fiber_breads_(grams/day)",
                                                            "mother_ffq_group_nuts_and_seeds_grams_per_day" = "nuts_and_seeds_(grams/day)",
                                                            "mother_ffq_group_pasta_grams_per_day" = "pasta_(grams/day)",
                                                            "mother_ffq_group_pastry_grams_per_day" = "pastry_(grams/day)",
                                                            "mother_ffq_group_potatoes_grams_per_day" = "potatoes_(grams/day)",
                                                            "mother_ffq_group_prepared_meal_grams_per_day" = "prepared_meals_(grams/day)",
                                                            
                                                            "mother_ffq_group_refined_white_breads_grams_per_day" = "refined_white_breads_(grams/day)",
                                                            "mother_ffq_group_rice_grams_per_day" = "rice_(grams/day)",
                                                            "mother_ffq_group_sandwich_fillers_grams_per_day" = "sandwich_fillers_(grams/day)",
                                                            "mother_ffq_group_sauces_condiments_grams_per_day" = "sauces_and_condiments_(grams/day)",
                                                            "mother_ffq_group_savoury_snacks_grams_per_day" = "savoury_snacks_(grams/day)",
                                                            "mother_ffq_group_sugar_sweets_grams_per_day" = "sugar_and_sweets_(grams/day)",
                                                            "mother_ffq_group_sweet_drinks_and_no_or_low_alc_drinks_ml_per_day" = "sweet_and_non-alcoholic_and_low-alcohol_drinks_(ml/day)",
                                                            "mother_ffq_group_tea_ml_per_day" = "tea_(ml/day)",
                                                            "mother_ffq_group_vegetables_grams_per_day" = "vegetables_(grams/day)",
                                                            "mother_ffq_group_water_ml_per_day" = "water_(ml/day)",
                                                            "mother_ffq_nutrients_ALA" = "alpha-linolenic_acid",
                                                            "mother_ffq_nutrients_alcohol_total" = "alcohol_(total)",
                                                            "mother_ffq_nutrients_beta.carotene" = "beta-carotene",
                                                            "mother_ffq_nutrients_calcium" = "calcium",
                                                            "mother_ffq_nutrients_carbohydrates_total" = "carbohydrates_(total)",
                                                            "mother_ffq_nutrients_cholesterol" = "cholesterol",
                                                            "mother_ffq_nutrients_DHA" = "docosahexaenoic_acid",
                                                            "mother_ffq_nutrients_dietary_fiber_total" = "dietary_fiber_(total)",
                                                            "mother_ffq_nutrients_dietary_folate_equivalents" = "dietary_folate_equivalents",
                                                            
                                                            "mother_ffq_nutrients_EPA" = "eicosapentaenoic_acid",
                                                            "mother_ffq_nutrients_fat_total" = "fat_(total)",
                                                            "mother_ffq_nutrients_fatty_acids_total_monounsaturated" = "monounsaturated_fatty_acids_(total)",
                                                            "mother_ffq_nutrients_fatty_acids_total_polyunsaturated" = "polyunsaturated_fatty_acids_(total)",
                                                            "mother_ffq_nutrients_fatty_acids_total_saturated" = "saturated_fatty_acids_(total)",
                                                            "mother_ffq_nutrients_fatty_acids_total_trans" = "trans-fatty_acids_(total)",
                                                            "mother_ffq_nutrients_folic_acid" = "folic_acid",
                                                            "mother_ffq_nutrients_iron_haem" = "haem_iron",
                                                            "mother_ffq_nutrients_iron_non_haem" = "non-haem_iron",
                                                            "mother_ffq_nutrients_iron_total" = "iron_(total)",
                                                            "mother_ffq_nutrients_linoleic_acid" = "linoleic_acid",
                                                            "mother_ffq_nutrients_lycopene" = "lycopene",
                                                            "mother_ffq_nutrients_magnesium" = "magnesium",
                                                            "mother_ffq_nutrients_mono_and_disaccharides_total" = "mono-_and_disaccharides_(total)",
                                                            "mother_ffq_nutrients_polysaccharides_total" = "polysaccharides_(total)",
                                                            "mother_ffq_nutrients_potassium" = "potassium",
                                                            "mother_ffq_nutrients_protein_animal" = "animal_protein",
                                                            "mother_ffq_nutrients_protein_total" = "protein_(total)",
                                                            "mother_ffq_nutrients_protein_vegetable" = "vegetable_protein",
                                                            
                                                            "mother_ffq_nutrients_retinol" = "retinol",
                                                            "mother_ffq_nutrients_retinol_activity_equivalents" = "retinol_activity_equivalents",
                                                            "mother_ffq_nutrients_vitamin_B1" = "vitamin_B1",
                                                            "mother_ffq_nutrients_vitamin_B12" = "vitamin_B12",
                                                            "mother_ffq_nutrients_vitamin_B2" = "vitamin_B2",
                                                            "mother_ffq_nutrients_vitamin_B6" = "vitamin_B6",
                                                            "mother_ffq_nutrients_vitamin_C" = "vitamin_C",
                                                            "mother_ffq_nutrients_vitamin_D" = "vitamin_D",
                                                            "mother_ffq_nutrients_vitamin_E" = "vitamin_E",
                                                            "mother_ffq_nutrients_water_total" = "water_(total)",
                                                            "mother_ffq_nutrients_zinc" = "zinc",
                                                            "mother_food_alcohol_during_first_trimester_yes" = "alcohol_consumption_during_first_trimester",
                                                            "mother_food_avoid_dairy_during_preg_yes" = "avoiding_dairy_during_pregnancy",
                                                            "mother_food_avoid_gluten_during_preg_yes" = "avoiding_gluten_during_pregnancy",
                                                            "mother_gravida" = "gravida",
                                                            "mother_health_dental_bleeding_gums_during_preg_yes" = "bleeding_gums_during_pregnancy",
                                                            "mother_health_dental_numeric" = "poor_dental_health",
                                                            "mother_health_preg_hyperemesis_gravidarum_severe_nausea_early_preg_yes" = "self-reported_hyperemesis_gravidarum_or_severe_nausea_during_early_pregnancy",
                                                            "mother_health_prior_preg_eating_disorder_yes" = "eating_disorder_before_pregnancy",
                                                            
                                                            "mother_health_reprod_age_menarche" = "age_at_menarche_(years)",
                                                            "mother_health_reprod_preg_IVF_yes" = "in_vitro_fertilization",
                                                            "mother_health_reprod_prior_preg_contraceptive_use_hormonal_yes" = "hormonal_contraceptive_use_before_pregnancy",
                                                            "mother_health_reprod_prior_preg_contraceptive_use_yes" = "contraceptive_use_before_pregnancy",
                                                            "mother_health_reprod_prior_preg_endometriosis_yes" = "endometriosis_before_pregnancy",
                                                            "mother_health_reprod_prior_preg_ovarian_cysts_yes" = "ovarian_cysts_before_pregnancy",
                                                            "mother_health_reprod_prior_preg_PCOS_yes" = "polycystic_ovary_syndrome_before_pregnancy",
                                                            "mother_health_reprod_prior_preg_regular_menstr_pattern_yes" = "regular_menstrual_pattern_before_pregnancy",
                                                            "mother_health_reprod_prior_preg_STD_yes" = "sexually_transmitted_disease_before_pregnancy",
                                                            "mother_height_cm" = "maternal_height_(cm)",
                                                            "mother_med_after_birth_oxytocin_yes" = "oxytocin_after_birth",
                                                            "mother_med_analgesics_N02_yes" = "analgesics_(N02)",
                                                            "mother_med_anti_hemorrhoids_C05A_yes" = "antihemorrhoids_(C05A)",
                                                            "mother_med_anti_thrombotic_agents_B01A_yes" = "antithrombotic_agents_(B01A)",
                                                            "mother_med_antibiotics_J01_yes" = "antibiotics_(J01)",
                                                            "mother_med_antifungals_topical_D01A_yes" = "topical_antifungals_(D01A)",
                                                            "mother_med_beta_lactam_antibacterials_penicillins_J01C_yes" = "beta-lactam_antibacterials_penicillins_(J01)",
                                                            "mother_med_birth_epidural_yes" = "epidural_during_birth",
                                                            "mother_med_birth_IV_contraction_inducers_yes" = "intravenous_contraction_inducers_during_birth",
                                                            
                                                            "mother_med_birth_IV_pump_remifentanil_yes" = "intravenous_remifentanil_pump_during_birth",
                                                            "mother_med_birth_oxytocin_yes" = "oxytocin_during_birth",
                                                            "mother_med_birth_spinal_yes" = "spinal_block_during_birth",
                                                            "mother_med_iron_preparations_B03A_yes" = "iron_preparations_(B03A)",
                                                            "mother_med_osmotically_acting_laxatives_A06AD_yes" = "osmotically_acting_laxatives_(A06AD)",
                                                            "mother_med_prior_preg_folic_acid_use_yes" = "folic_acid_use_before_pregnancy",
                                                            "mother_milk_collection_breasts_sampled_both_breasts" = "sampled_from_both_breasts",
                                                            "mother_milk_collection_season_Autumn" = "collection_season_(autumn_vs_spring)",
                                                            "mother_milk_collection_season_Summer" = "collection_season_(summer_vs_spring)",
                                                            "mother_milk_collection_season_Winter" = "collection_season_(winter_vs_spring)",
                                                            "mother_milk_collection_time_numeric" = "time_of_the_day",
                                                            # "mother_milk_collection_volume_pumped" = "pumped_milk_volume_(ml)",
                                                            "mother_milk_HMO_2FL_ugml" = "2'FL_(µg/ml)",
                                                            "mother_milk_HMO_3F3SL_ugml" = "3'F3'SL_(µg/ml)",
                                                            "mother_milk_HMO_3FL_ugml" = "3'FL_(µg/ml)",
                                                            "mother_milk_HMO_3GL_ugml" = "3'GL_(µg/ml)",
                                                            "mother_milk_HMO_3SL_ugml" = "3'SL_(µg/ml)",
                                                            "mother_milk_HMO_6GL_ugml" = "6'GL_(µg/ml)",
                                                            "mother_milk_HMO_6SL_ugml" = "6'SL_(µg/ml)",
                                                            
                                                            "mother_milk_HMO_A_tetra_ugml" = "A-tetra_(µg/ml)",
                                                            "mother_milk_HMO_DFLNHa_ugml" = "DFLNHα_(µg/ml)",
                                                            "mother_milk_HMO_DSLNT_ugml" = "DSLNT_(µg/ml)",
                                                            "mother_milk_HMO_Fuc_ugml" = "fucosylatted_HMOs_(µg/ml)",
                                                            "mother_milk_HMO_LDFT_ugml" = "LDFT_(µg/ml)",
                                                            "mother_milk_HMO_LNDH_I_ugml" = "LNDH_I_(µg/ml)",
                                                            "mother_milk_HMO_LNFP_I_ugml" = "LNFP_I_(µg/ml)",
                                                            "mother_milk_HMO_LNFP_II_ugml" = "LNFP_II_(µg/ml)",
                                                            "mother_milk_HMO_LNFP_III_ugml" = "LNFP_III_(µg/ml)",
                                                            "mother_milk_HMO_LNFP_V_ugml" = "LNFP_V_(µg/ml)",
                                                            "mother_milk_HMO_LNH_ugml" = "LNH_(µg/ml)",
                                                            "mother_milk_HMO_LNnDFH_ugml" = "LNnDFH_(µg/ml)",
                                                            "mother_milk_HMO_LNnFP_V_ugml" = "LNnFP_V_(µg/ml)",
                                                            "mother_milk_HMO_LNnT_ugml" = "LNnT_(µg/ml)",
                                                            "mother_milk_HMO_LNT_ugml" = "LNT_(µg/ml)",
                                                            "mother_milk_HMO_LSTb_ugml" = "LSTb_(µg/ml)",
                                                            "mother_milk_HMO_LSTc_ugml" = "LSTc_(µg/ml)",
                                                            "mother_milk_HMO_MFLNH_III_ugml" = "MFLNH_III_(µg/ml)",
                                                            "mother_milk_HMO_milk_group_Le-Se-" = "maternal_milk_group_(Le-Se-_vs_Le+Se+)",
                                                            
                                                            "mother_milk_HMO_milk_group_Le-Se+" = "maternal_milk_group_(Le-Se+_vs_Le+Se+)",
                                                            "mother_milk_HMO_milk_group_Le+Se-" = "maternal_milk_group_(Le+Se-_vs_Le+Se+)",
                                                            "mother_milk_HMO_Neut_ugml" = "neutral_HMOs_(µg/ml)",
                                                            "mother_milk_HMO_Sia_ugml" = "sialylated_HMOs_(µg/ml)",
                                                            "mother_milk_HMO_Total_ugml" = "total_HMOs_(µg/ml)",
                                                            "mother_net_income_numeric" = "higher_maternal_net_income",
                                                            "mother_para" = "para",
                                                            "mother_preg_comp_breech_yes" = "breech",
                                                            "mother_preg_comp_cholestatis_yes" = "cholestatis",
                                                            "mother_preg_comp_gestational_diabetes_yes" = "gestational_diabetes",
                                                            "mother_preg_comp_hyperemesis_gravidarum_yes" = "hyperemesis_gravidarum",
                                                            "mother_preg_comp_hypertension_yes" = "hypertension",
                                                            "mother_preg_comp_past_C_section_yes" = "past_C-section",
                                                            "mother_preg_comp_preeclampsia_yes" = "preeclampsia",
                                                            "mother_preg_comp_prolonged_rupture_membranes_yes" = "prolonged_rupture_of_membranes",
                                                            "mother_previous_preg_yes" = "previous_pregnancy",
                                                            "mother_stress_mean_postpartum_stress" = "higher_maternal_postpartum_stress",
                                                            "mother_stress_mean_pregnancy_stress" = "higher_maternal_stress_during_pregnancy",
                                                            "mother_type_pregnancy_twin_pregnancy" = "twin_pregnancy",
                                                            "mother_weight_early_preg_kg" = "maternal_weight_during_early_pregnancy_(kg)",
                                                            "mother_weight_just_before_delivery_kg" = "maternal_weight_just_before_birth_(kg)",
                                                            "mother_weight_prior_preg_kg" = "maternal_weight_before_pregnancy_(kg)")

## make first letter capitalized for each trait
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
model_results_short$trait_name_in_plot <- as.factor(as.character(firstup(as.character(as.factor(model_results_short$trait_name_in_plot)))))


### ===== 19.9 ENSURE CORRECT N FOR PHENOTYPES WITH >2 CATEGORIES ===== ###

for (i in c(1,9,10)){model_results_short[,i] <- as.factor(as.character(model_results_short[,i]))}
levels(model_results_short$levels)[1] <- NA
model_results_short$levels <- droplevels(model_results_short$levels)
levels(model_results_short$levels)


### MATERNAL BLOOD GROUP ###

all_sum[all_sum$Phenotype=="mother_blood_group",c(1,4,8,9)]
# Phenotype                 n_answered  Categories  n_per_category
# mother_blood_group        620         A;AB;B;O    280;20;44;276

model_results_short[model_results_short$x=="mother_blood_group" & model_results_short$levels=="A", "n_total"] <- 280+276
model_results_short[model_results_short$x=="mother_blood_group" & model_results_short$levels=="AB", "n_total"] <- 20+276
model_results_short[model_results_short$x=="mother_blood_group" & model_results_short$levels=="B", "n_total"] <- 44+276


### BREAST PUMP BRAND

all_sum[all_sum$Phenotype=="mother_breastpump_brand",c(1,4,8,9)]
# Phenotype                       n_answered   Categories                     n_per_category
# mother_breastpump_brand         49           Ardo;Medela;Other;Philips      4;17;6;22

model_results_short[model_results_short$x=="mother_breastpump_brand" & model_results_short$levels=="Ardo", "n_total"] <- 4+22
model_results_short[model_results_short$x=="mother_breastpump_brand" & model_results_short$levels=="Medela", "n_total"] <- 17+22
model_results_short[model_results_short$x=="mother_breastpump_brand" & model_results_short$levels=="Other", "n_total"] <- 6+22


### MILK COLLECTION SEASON

all_sum[all_sum$Phenotype=="mother_milk_collection_season",c(1,4,8,9)]
# Phenotype                            n_answered   Categories                   n_per_category
# mother_milk_collection_season        724          Autumn;Spring;Summer;Winter  163;212;173;176

model_results_short[model_results_short$x=="mother_milk_collection_season" & model_results_short$levels=="Autumn", "n_total"] <- 163+212
model_results_short[model_results_short$x=="mother_milk_collection_season" & model_results_short$levels=="Summer", "n_total"] <- 173+212
model_results_short[model_results_short$x=="mother_milk_collection_season" & model_results_short$levels=="Winter", "n_total"] <- 176+212


### INFANT FEEDING TYPE AFTER DELIVERY

all_sum[all_sum$Phenotype=="infant_ffq_feeding_type_delivery",c(1,4,8,9)]
# Phenotype                               n_answered  Categories                                               n_per_category
# infant_ffq_feeding_type_delivery        216         breast_milk_via_bottle;breastfeeding;mixed_feeding       2;204;10

model_results_short[model_results_short$x=="infant_ffq_feeding_type_delivery" & model_results_short$levels=="breast_milk_via_bottle", "n_total"] <- 2+204
model_results_short[model_results_short$x=="infant_ffq_feeding_type_delivery" & model_results_short$levels=="mixed_feeding", "n_total"] <- 10+204


### LIVING SITUATION

all_sum[all_sum$Phenotype=="mother_exp_living_situation",c(1,4,8,9)]
# Phenotype                          n_answered   Categories                             n_per_category
# mother_exp_living_situation        605          big_house;farm;flat;other;small_house  268;24;15;57;241

model_results_short[model_results_short$x=="mother_exp_living_situation" & model_results_short$levels=="farm", "n_total"] <- 268+24
model_results_short[model_results_short$x=="mother_exp_living_situation" & model_results_short$levels=="flat", "n_total"] <- 268+15
model_results_short[model_results_short$x=="mother_exp_living_situation" & model_results_short$levels=="other", "n_total"] <- 268+57
model_results_short[model_results_short$x=="mother_exp_living_situation" & model_results_short$levels=="small_house", "n_total"] <- 268+241


### MILK GROUP

all_sum[all_sum$Phenotype=="mother_milk_HMO_milk_group",c(1,4,8,9)]
# Phenotype                         n_answered     Categories                   n_per_category
# mother_milk_HMO_milk_group        725            Le-Se-;Le-Se+;Le+Se-;Le+Se+  13;41;177;494

model_results_short[model_results_short$x=="mother_milk_HMO_milk_group" & model_results_short$levels=="Le-Se-", "n_total"] <- 13+494
model_results_short[model_results_short$x=="mother_milk_HMO_milk_group" & model_results_short$levels=="Le-Se+", "n_total"] <- 41+494
model_results_short[model_results_short$x=="mother_milk_HMO_milk_group" & model_results_short$levels=="Le+Se-", "n_total"] <- 177+494


### ===== 19.10 SPLIT BY OUTPUT TYPE (ALPHA DIVERSITY / RELATIVE ABUNDANCES) ===== ###

## separate results for milk alpha diversity and for milk relative bacterial abundances
model_results_alpha <- model_results_short[grep("alpha_div", model_results_short$y),] #591x11
model_results_RA <- model_results_short[grep("g__", model_results_short$y),] #7092x11


### ===== 19.11 CORRECT FOR MULTIPLE TESTING AND SAVE RESULTS TABLES ===== ###

### MILK ALPHA DIVERSITY

## exclude results for Simpson diversity
model_results_alpha <- model_results_alpha[model_results_alpha$y!="alpha_div_simpson_invr",]

## calculate FDR
model_results_alpha$FDR <- p.adjust(model_results_alpha$p, method="BH")

## reorder data frame by FDR results
model_results_alpha <- model_results_alpha[order(model_results_alpha$FDR, model_results_alpha$p),]

## check results
nrow(model_results_alpha[model_results_alpha$FDR<0.05,]) # ->  0 FDR significant associations
nrow(model_results_alpha[model_results_alpha$p<0.05,])   # -> 30 nominally significant associations
length(unique(model_results_alpha$x)) #185

## save results
# write.table(model_results_alpha, file="association_model_results/240229_model_results_milk_alpha_diversity_by_phenotypes_n185.txt", row.names=F, col.names=T, sep="\t", quote=F)
# saveRDS(model_results_alpha, file="association_model_results/240229_model_results_milk_alpha_diversity_by_phenotypes_n185.rds")
# model_results_alpha <- readRDS(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/association_model_results/240229_model_results_milk_alpha_diversity_by_phenotypes_n185.rds")


### MILK RELATIVE BACTERIAL ABUNDANCES

## calculate FDR
model_results_RA$FDR <- p.adjust(model_results_RA$p, method="BH")

## reorder data frame by FDR results
model_results_RA <- model_results_RA[order(model_results_RA$FDR, model_results_RA$p),]

## check results
nrow(model_results_RA[model_results_RA$FDR<0.05,]) # ->   8 FDR significant associations
nrow(model_results_RA[model_results_RA$p<0.05,])   # -> 613 nominally significant associations
length(unique(model_results_RA$x)) #185

## save results
# write.table(model_results_RA, file="association_model_results/240229_model_results_milk_relative_abundances_by_phenotypes_n185.txt", row.names=F, col.names=T, sep="\t", quote=F)
# saveRDS(model_results_RA, file="association_model_results/240229_model_results_milk_relative_abundances_by_phenotypes_n185.rds")


##### =========================== 20. PREPARE SUPPLEMENTARY TABLES (PHENOTYPE SUMMARY STATISTICS TABLE AND MODEL RESULTS TABLES)  =========================== #####

### ===== 20.1 PHENOTYPE SUMMARY STATISTICS TABLE (TABLE S14) ===== ###

## desired rows in phenotype summary statistics supplementary table
## -> only include rows for phenotypes that were included in the analysis / model results table

## select only relevant rows from all_sum table
# length(unique(model_results_alpha$x)) #185
# length(unique(model_results_RA$x)) #185
# length(intersect(unique(model_results_alpha$x), unique(model_results_RA$x))) #185
p <- unique(model_results_alpha$x) #vector with all 185 included phenotypes
p <- gsub("_invr", "", p) #ensure the phenotypes have the same colnames as in the all_sum table (remove '_invr' extension for transformed phenotypes)
psum <- all_sum[all_sum$Phenotype %in% p,] #185x20


## current columns in phenotype summary statistics table
# colnames(all_sum)
# [1] "Phenotype"        "data_type"        "n_total"          "n_answered"      
# [5] "perc_answered"    "n_missing"        "perc_missing"     "Categories"      
# [9] "n_per_category"   "Mean"             "SD"               "Median"          
# [13] "1st_quartile"     "3rd_quartile"     "Minimum"          "Maximum"         
# [17] "model_exclusion"  "model_use"        "model_time_point" "Note" 

## desired columns in phenotype summary statistics supplementary table
# [1] "Phenotype"        "data_type"        "n_total"          "n_answered"      
# [5] "perc_answered"    "n_missing"        "perc_missing"     "Categories"      
# [9] "n_per_category"   "Mean"             "SD"               "Median"          
# [13] "1st_quartile"     "3rd_quartile"     "Minimum"          "Maximum"         
# [17]   "model_use"        "model_time_point" "Note" 
## -> keep all from all_sum except "model_exclusion"
## -> add 'phenotype' group column from the model results tables

## create an intermediate table that links the 'Phenotype' to the 'phenotype_group'
pgroups <- model_results_alpha[,c(3,10)]
pgroups <- pgroups[!duplicated(pgroups),] #185x2, removed duplicated rows

## add pgroups to psum table
pgroups$x <- gsub("_invr", "", pgroups$x)
colnames(pgroups)[1] <- "Phenotype"
psum2 <- left_join(psum, pgroups, by="Phenotype")

## remove "model_exclusion" and "model_use" columns from psum2
psum3 <- psum2[,-c(17:18)]

## resort columns in psum3 ##note: this also excludes model_time_point!
psum4 <- psum3[,c(19,1:16,18)] #185x18

## correct entries for stool diary
psum4$Phenotype <- gsub("mother_dis_diary_freq_mean", "mother_dis_diary_stool_freq_mean", psum4$Phenotype)
psum4$Phenotype <- gsub("mother_dis_diary_not_well_mean", "mother_dis_diary_uneasy_abdominal_feeling_mean", psum4$Phenotype)

## save polished phenotype summary statistics table
write.table(psum4, file="phenotype_summary_statistics/240301_milk_microbiota_included_maternal_and_infant_phenotypes_summary_statistics_n185.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 20.2 MILK ALPHA DIVERSITY - MODEL RESULTS TABLE (TABLE S16) ===== ###

resalpha_tmp <- model_results_alpha

# add 'phenotype_type column
resalpha_tmp$phenotype_type <- as.factor(as.character(gsub(".*microbiota_", "", resalpha_tmp$dataset)))

# fix entries in 'model' ('statistic') column
# resalpha_tmp$statistic <- as.factor(as.character(resalpha_tmp$statistic))
levels(resalpha_tmp$statistic) <- c("mmrm_compound_symmetry_correction_for_DNAbatch_reads_time_ID",
                                    "mmrm_unstructured_correction_for_DNAbatch_reads_time_ID",
                                    "lm_correction_for_DNAbatch_reads",
                                    "lm_correction_for_reads",
                                    "lmer_correction_for_DNAbatch_reads_time_ID",
                                    "lmer_correction_for_reads_time_ID")

# remove extensions for transformations from x and y columns (only mention i table description later)
resalpha_tmp$x <- gsub("_invr", "", resalpha_tmp$x)
resalpha_tmp$y <- gsub("_invr", "", resalpha_tmp$y)

# add columns from all_sum (data_type, model_time_point)
all_sum_tmp <- all_sum[all_sum$Phenotype %in% resalpha_tmp$x, c(1,2,19)] #185x3
colnames(all_sum_tmp)[1] <- "x"
resalpha_tmp2 <- left_join(resalpha_tmp, all_sum_tmp, by="x")

# # ensure NAs show for numeric/integer phenotypes in investigated category column
# resalpha_tmp2$levels <- as.factor(as.character(resalpha_tmp2$levels))
# levels(resalpha_tmp2$levels)[1] <- NA

# add column showing the reference category for categorical phenotypes
resalpha_tmp2$reference_category <- resalpha_tmp2$levels
levels(resalpha_tmp2$reference_category) <- c("O","O","O",
                                            "big_house","big_house","Le+Se+",
                                            "Le+Se+","Le+Se+","big_house",
                                            "big_house","single_pregnancy","no",
                                            "vaginal_birth","home","female",
                                            "breastfeeding","breastfeeding","Philips",
                                            "Spring","one_breast","Philips",
                                            "Philips","Spring","Spring")

# remove columns 'dataset' and resort columns in resalpha_tmp2
resalpha_tmp3 <- resalpha_tmp2[,c(10,3,14,16,6,11,4,5,7,8,12,15,13,2,9)]

# fix colnames in resalpha_tmp3
colnames(resalpha_tmp3) <- c("Phenotype_group", "Phenotype", "Data_type", "Reference_category", "Tested_category", "Phenotype_label_in_plots",
                             "Alpha_diversity_measure",
                             "n", "Estimate", "p", "FDR",
                             "Outcome_time_points_linked_to_phenotype", "Phenotype_type", "Model", "Note")

# fix entries for alpha diversity measures
resalpha_tmp3$Alpha_diversity_measure <- gsub("alpha_div_", "", resalpha_tmp3$Alpha_diversity_measure)

# ## fix maternal education label
# resalpha_tmp3$Phenotype_label_in_plots <- as.factor(as.character(gsub("Higer_maternal_education", "Higher_maternal_education", resalpha_tmp3$Phenotype_label_in_plots)))


resalpha <- resalpha_tmp3

## save polished phenotype summary statistics table
write.table(resalpha, file="association_model_results/240302_polished_model_results_milk_alpha_diversity_by_phenotypes_n185.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 20.3 MILK RELATIVE BACTERIAL ABUNDANCES - MODEL RESULTS TABLE (TABLE S17) ===== ###

resRA_tmp <- model_results_RA

# add 'phenotype_type column
resRA_tmp$phenotype_type <- as.factor(as.character(gsub(".*microbiota_", "", resRA_tmp$dataset)))

# fix entries in 'model' ('statistic') column
# resRA_tmp$statistic <- as.factor(as.character(resRA_tmp$statistic))
levels(resRA_tmp$statistic) <- c("mmrm_compound_symmetry_correction_for_DNAbatch_reads_time_ID",
                                    "mmrm_unstructured_correction_for_DNAbatch_reads_time_ID",
                                    "lm_correction_for_DNAbatch_reads",
                                    "lm_correction_for_reads",
                                    "lmer_correction_for_DNAbatch_reads_time_ID",
                                    "lmer_correction_for_reads_time_ID")


# remove extensions for transformations from x and y columns (only mention i table description later)
resRA_tmp$x <- gsub("_invr", "", resRA_tmp$x)
resRA_tmp$y <- gsub("g__", "", resRA_tmp$y)
resRA_tmp$y <- gsub("_clr", "", resRA_tmp$y)

# add columns from all_sum (data_type, model_time_point)
# all_sum_tmp <- all_sum[all_sum$Phenotype %in% resRA_tmp$x, c(1,2,19)] #185x3
# colnames(all_sum_tmp)[1] <- "x"
resRA_tmp2 <- left_join(resRA_tmp, all_sum_tmp, by="x")

# # ensure NAs show for numeric/integer phenotypes in investigated category column
# resRA_tmp2$levels <- as.factor(as.character(resRA_tmp2$levels))
# levels(resRA_tmp2$levels)[1] <- NA

# add column showing the reference category for categorical phenotypes
resRA_tmp2$reference_category <- resRA_tmp2$levels
levels(resRA_tmp2$reference_category) <- c("O","O","O",
                                              "big_house","big_house","Le+Se+",
                                              "Le+Se+","Le+Se+","big_house",
                                              "big_house","single_pregnancy","no",
                                              "vaginal_birth","home","female",
                                              "breastfeeding","breastfeeding","Philips",
                                              "Spring","one_breast","Philips",
                                              "Philips","Spring","Spring")

# remove columns 'dataset' and resort columns in resRA_tmp2
resRA_tmp3 <- resRA_tmp2[,c(10,3,14,16,6,11,4,5,7,8,12,15,13,2,9)]

# fix colnames in resRA_tmp3
colnames(resRA_tmp3) <- c("Phenotype_group", "Phenotype", "Data_type", "Reference_category", "Tested_category", "Phenotype_label_in_plots",
                             "Relative_bacterial_abundance",
                             "n", "Estimate", "p", "FDR",
                             "Outcome_time_points_linked_to_phenotype", "Phenotype_type", "Model", "Note")

# ## fix maternal education label
# resRA_tmp3$Phenotype_label_in_plots <- as.factor(as.character(gsub("Higer_maternal_education", "Higher_maternal_education", resRA_tmp3$Phenotype_label_in_plots)))

resRA <- resRA_tmp3

## save polished phenotype summary statistics table
write.table(resRA, file="association_model_results/240301_polished_model_results_milk_relative_abundances_by_phenotypes_n185.txt", row.names=F, col.names=T, sep="\t", quote=F)








