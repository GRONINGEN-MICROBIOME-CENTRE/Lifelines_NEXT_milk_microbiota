################################################################################################################
### CREATE PHENOTYPES SUMMARY STATISTICS + MODEL RESULTS TABLE FOR REAL HMO DATA FOR MILK COMPOSITION PAPER  ###
################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 0. DATA IMPORT

# 1. SELECT DATA OF INTEREST
#    1.1 MATERNAL PHENOTYPE DATA
#    1.2 INFANT PHENOTYPE DATA

# 2. ENSURE CORRECT DATA STRUCTURE
#    2.1 MATERNAL PHENOTYPES
#    2.2 INFANT PHENOTYPES

# 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
#    3.1 MATERNAL PHENOTYPES
#    3.2 INFANT PHENOTYPES

# 4. CHANGE SELECTED CATEGORICAL PHENOTYPE DATA TO NUMERIC FOR USE IN ASSOCIATION STUDIES
#    4.1 MATERNAL PHENOTYPES
#    4.2 INFANT PHENOTYPES

# 5. CREATE TABLE SHOWING PHENOTYPE DATA TYPE (FACTOR OR NUMERIC/INTEGER)
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

# 10. INVERSE-RANK TRANSFORM NUMERIC PHENOTYPES AND HMO DATA
#    10.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
#    10.2 INVERSE-RANK TRANSFORM NUMERIC PHENOTYPE DATA
#    10.3 INVERSE-RANK TRANSFORM HMO DATA

# 11. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED NUMERIC PHENOTYPES
#    11.1 FUNCTION TO MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED NUMERIC PHENOTYPES
#    11.2 MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED NUMERIC PHENOTYPES

# 12. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED HMO DATA
#    12.1 FUNCTION TO MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED HMO CONCENTRATIONS PER MILK GROUP
#    12.2 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED HMO DATA FROM MATERNAL PHENOTYPE FILE
#    12.3 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED HMO DATA FROM INFANT PHENOTYPE FILE

# 13. PREPARE MATERNAL AND INFANT DATA SETS FOR ASSOCIATION OF PHENOTYPES WITH HMOs
#    13.1 OVERVIEW OF CURRENT DATA FRAMES PREPARED FOR ASSOCIATIONS OF HMOs AND PHENOTYPE
#    13.2 PREPARE DATA FRAMES WITH MATERNAL PHENOTYPES
#    13.3 PREPARE DATA FRAMES WITH INFANT PHENOTYPES

# 14. MODEL FOR ASSOCIATION OF MULTIPLE-TIME-POINT STATIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID
#    14.1 PREPARE DATA FRAMES WITH MATERNAL STATIC PHENOTYPES
#    14.2 PREPARE DATA FRAMES WITH INFANT STATIC PHENOTYPES
#    14.3 FUNCTION FOR ASSOCIATION OF STATIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID
#    14.4 RUN ASSOCIATION OF MATERNAL STATIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS
#    14.5 RUN ASSOCIATION OF INFANT STATIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS

# 15. MODEL FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP
#    15.1 PREPARE DATA FRAMES WITH MATERNAL 1-TIME-POINT PHENOTYPES
#    15.2 PREPARE DATA FRAMES WITH INFANT 1-TIME-POINT PHENOTYPES
#    15.3 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP
#    15.4 RUN ASSOCIATION OF MATERNAL SINGLE-TIME-POINT STATIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS
#    15.5 RUN ASSOCIATION OF INFANT SINGLE-TIME-POINT STATIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS

# 16. MODEL FOR ASSOCIATION OF MULTIPLE-TIME-POINT DYNAMIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID
#    16.1 PREPARE DATA FRAMES WITH MATERNAL DYNAMIC PHENOTYPES
#    16.2 PREPARE DATA FRAMES WITH INFANT DYNAMIC PHENOTYPES
#    16.3 FUNCTION FOR ASSOCIATION OF DYNAMIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID
#    16.4 RUN ASSOCIATION OF MATERNAL DYNAMIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS
#    16.5 RUN ASSOCIATION OF INFANT DYNAMIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS

# 17. COMBINE RESULTS AND CORRECT FOR MULTIPLE TESTING
#    17.1 PREPARE DATA FRAMES FOR MERGING
#    17.2 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM STATIC PHENOTYPE ASSOCIATIONS
#    17.3 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM 1-TIME-POINT PHENOTYPE ASSOCIATIONS
#    17.4 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM DYNAMIC PHENOTYPE ASSOCIATIONS
#    17.5 MERGE ALL RESULTS DATA FRAMES AND SPLIT BY OUTPUT TYPE (ALPHA DIVERSITY / RELATIVE ABUNDANCES)
#    17.6 EXCLUDE PHENOTYPES
#    17.7 ADD PHENOTYPE GROUPS FOR PLOTTING
#    17.8 ADD PRETTY LABELS FOR PHENOTYPES FOR PLOTTING
#    17.9 ENSURE CORRECT N FOR PHENOTYPES WITH >2 CATEGORIES
#    17.10 CORRECT FOR MULTIPLE TESTING AND SAVE RESULTS TABLES

# 18. PREPARE SUPPLEMENTARY TABLES (PHENOTYPE SUMMARY STATISTICS TABLE AND MODEL RESULTS TABLES)
#    18.1 PHENOTYPE SUMMARY STATISTICS TABLE (TABLE S7) 
#    18.2 MILK HMOs - MODEL RESULTS TABLE (TABLE S8)


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/")

## import file
hmo <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/phenotypes/240108_milk_composition_paper_sel_phenotypes_hmo_data_n1563.txt", header=T, sep="\t", stringsAsFactors=T)
#1563x512

# Notes:
# This file contains the following:
# - Phenotypes linked to measured HMO data.
#   !! Note that maternal phenotype and HMO data was duplicated for twins. This is indicated by ";1" (single baby/first baby) and ";2" (twin baby/duplicated data) in the mother_ID and mother_sample_ID.
#      -> For association with maternal phenotypes (i.e. non-duplicated data), grep only ";1"
#      -> For association with infant phenotypes, use both ";1" and ";2"
# - Real measured levels of 24 single HMOs and 4 grouped HMOs in µg/ml and HMO-based Le and Se status and milk groups
# - A quality control sample was included in all HMO batches and if there was less than <15% variation, the quality of the run was considered good. No need to correct for a batch effect from HMO measurements.


##### =========================== 1. SELECT DATA OF INTEREST =========================== #####

### ===== 1.1 MATERNAL PHENOTYPE DATA ===== ###

## select every mother only 1x (i.e. remove duplication of data generated for twin babies) and select IDs, maternal phenotype columns and HMOs
mhmo <- hmo[grep(";1", hmo$mother_sample_ID),c(1:434,481:512)] #1542x466


### ===== 1.2 INFANT PHENOTYPE DATA ===== ###

## create data frames to use for infants: select IDs, infant phenotype columns (which can potentially affect maternal milk HMO concentrations) and HMOs
ihmo <- hmo[,c(1:15,435:443,450:453,481:512)] #1563x60

# [435] "infant_genetics_FUT2_G_to_A_rs601338"                                             
# [436] "infant_genetics_FUT2_secretor"                                                    
# [437] "infant_Se_status_same_as_mother"  
# [438] "infant_sex"                                                                       
# [439] "infant_birth_delivery_place"                                                      
# [440] "infant_birth_delivery_place_detailed"                                             
# [441] "infant_birth_delivery_mode"                                                       
# [442] "infant_birth_delivery_mode_detailed"                                              
# [443] "infant_birthweight_g" 
# 
# [450] "infant_ffq_feeding_type_delivery"                                                 
# [451] "infant_ffq_breastfeeding_behaviour_after_birth_1bad_5excellent"
# [452] "infant_ffq_feeding_mode" 
# [453] "infant_ffq_solid_food_introduced"


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

### ===== 2.1 MATERNAL PHENOTYPES ===== ###

## Maternal phenotypes: check that data structure is correct
str(mhmo[,1:50])
str(mhmo[,51:100])
str(mhmo[,101:150])
str(mhmo[,151:200])
str(mhmo[,201:250])
str(mhmo[,251:300])
str(mhmo[,301:350])
str(mhmo[,351:400])
str(mhmo[,400:ncol(mhmo)])
# -> data structure is correct


### ===== 2.2 INFANT PHENOTYPES ===== ###

## Infant phenotypes: check that data structure is correct
str(ihmo[,1:50])
str(ihmo[,50:ncol(ihmo)])
# -> data structure is correct


##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

### ===== 3.1 MATERNAL PHENOTYPES ===== ###

## Maternal phenotypes: for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
mhmo$mother_milk_collection_season <- factor(mhmo$mother_milk_collection_season, levels = c("Spring", "Summer", "Autumn", "Winter"))
mhmo$mother_milk_collection_month <- factor(mhmo$mother_milk_collection_month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",  "Oct",  "Nov", "Dec"))
mhmo$mother_milk_collection_notes <- factor(mhmo$mother_milk_collection_notes, levels = c("No_records_of_incorrect_sample_handling", levels(mhmo$mother_milk_collection_notes)[c(3,5,7,6,2,1)]))
mhmo$mother_milk_collection_breasts_sampled <- factor(mhmo$mother_milk_collection_breasts_sampled, levels=c("one_breast", "both_breasts"))
mhmo$mother_genetics_blood_group_genotype <- factor(mhmo$mother_genetics_blood_group_genotype, levels = c("O/O", "A/O", "A/A", "B/O", "B/B", "A/B"))
mhmo$mother_blood_group <- factor(mhmo$mother_blood_group, levels = c("O", "A", "B", "AB"))
mhmo$mother_exp_living_situation <- factor(mhmo$mother_exp_living_situation, levels = c("big_house", "small_house", "flat", "farm", "other"))
mhmo$mother_birth_gestational_age_categories <- factor(mhmo$mother_birth_gestational_age_categories, levels = c("term", "preterm"))
mhmo$mother_breastpump_brand <- factor(mhmo$mother_breastpump_brand, levels = c("Philips", "Medela", "Ardo", "Other"))

# also change it for milk groups
mhmo$mother_milk_HMO_milk_group <- factor(mhmo$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


### ===== 3.2 INFANT PHENOTYPES ===== ###

## Infant phenotypes: for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
ihmo$infant_birth_delivery_mode <- factor(ihmo$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
ihmo$infant_birth_delivery_mode_detailed <- factor(ihmo$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                        "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
ihmo$infant_ffq_feeding_type_delivery <- factor(ihmo$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))

# also change it for milk groups
ihmo$mother_milk_HMO_milk_group <- factor(ihmo$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 4. CHANGE SELECTED CATEGORICAL PHENOTYPE DATA TO NUMERIC FOR USE IN ASSOCIATION STUDIES  =========================== #####

## Some categorical phenotypes will be used as numeric phenotypes in models. These phenotypes are here changed to numeric before making the summary statistics tables so that the tables show the data as used for associations.
## I will later add info to the supplementary table on which factors/categories the numbers correspond to.


### ===== 4.1 MATERNAL PHENOTYPES ===== ###

## vector with the names of the columns that will be changed to numeric
matnumcatphenos <- c("mother_milk_collection_time","mother_education_level","mother_net_income","mother_health_dental")

## check data before change
# head(mhmo[,colnames(mhmo) %in% matnumcatphenos])
# tail(mhmo[,colnames(mhmo) %in% matnumcatphenos])
# table(mhmo$mother_milk_collection_time, useNA="ifany")
# table(mhmo$mother_education_level, useNA="ifany")
# table(mhmo$mother_net_income, useNA="ifany")
# table(mhmo$mother_health_dental, useNA="ifany")

## change data entries for numcatphenos to numeric
for (i in matnumcatphenos){mhmo[,i] <- as.numeric(as.character(gsub("__.*", "", mhmo[,i])))}

## show conversion from categorical to numeric data structure in column name
a <- c()
for (i in matnumcatphenos){a <- c(a, grep(i, colnames(mhmo)))}
length(a) #5
length(matnumcatphenos) #4 -> one phenotype too much was pulled from the data frame
setdiff(colnames(mhmo)[a], matnumcatphenos) # -> remove "mother_health_dental_bleeding_gums_during_preg" from vector a
colnames(mhmo)[a]
a <- a[-5]
setdiff(colnames(mhmo)[a], matnumcatphenos) 
length(intersect(colnames(mhmo)[a], matnumcatphenos))

for (i in a){colnames(mhmo)[i] <- paste0(colnames(mhmo)[i], "_numeric")}

## check data after change
# head(mhmo[,a])
# tail(mhmo[,a])
# table(mhmo$mother_milk_collection_time_numeric, useNA="ifany")
# table(mhmo$mother_education_level_numeric, useNA="ifany")
# table(mhmo$mother_net_income_numeric, useNA="ifany")
# table(mhmo$mother_health_dental_numeric, useNA="ifany")


### ===== 4.2 INFANT PHENOTYPES ===== ###

## No infant phenotypes were changed to numeric.


##### =========================== 5. CREATE TABLE SHOWING PHENOTYPE DATA TYPE (FACTOR OR NUMERIC/INTEGER)  =========================== #####

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
overview_mhmo <- create.overview.table(mhmo[,c(8:10,16:434)])

## add info in data_use column to indicate whether a used phenotype is static or dynamic or 1-time-point
overview_mhmo$model_use <- c(rep("dynamic",length(1:2)),     #time_point and time_point_numeric
                             rep("static",length(3)),        #type_pregnancy
                             rep("dynamic",length(4:9)),     #milk collection -> exclude selected phenotypes from models
                             rep("static",length(10:386)),   #maternal genetics (exclude genotypes for FUT2 and FUT3?), maternal anthropometrics, food pref, ffq (exclude items from 1st analysis), stress, education, exp living and pets, mat health, preg/birth comp
                             rep("dynamic",length(387:388)), #subclinical mastitis and Na/K ratio
                             rep("static",length(389:410)),  #maternal GIT diary, ROME, preg and birth med use
                             rep("dynamic",length(411:422))) #maternal postpartum med use, breastpumping

## add info on which phenotypes to exclude as x from association models
overview_mhmo$model_exclusion <- c(rep("exclude",length(1:2)),     #time_point and time_point_numeric (exclude as it will be corrected for)
                                   rep("include",length(3:8)),     #type_pregnancy, milk collection
                                   rep("exclude",length(9:16)),    #milk collection notes, maternal FUT2 and FUT3 genotypes
                                   rep("include",length(17:422)))  #maternal anthropometrics, food pref, food groups, food items, nutrients, stress, education, exp living and pets, mat health, preg/birth comp, subclinical mastitis and Na/K ratio, maternal GIT diary, ROME, preg and birth med use, maternal postpartum med use, breastpumping

## add info on which HMO time points phenotypes were/will be linked to
overview_mhmo$model_time_point <- c(rep(NA,length(1:2)),                       #time_point and time_point_numeric (exclude as it will be corrected for)
                                    rep("0.5_to_6_months",length(3:8)),        #type_pregnancy, milk collection
                                    rep(NA,length(9:16)),                      #milk collection notes, maternal FUT2 and FUT3 genotypes
                                    rep("0.5_to_6_months",length(17:386)),     #maternal anthropometrics, food pref, food groups, food items, nutrients, stress, education, exp living and pets, mat health, preg/birth comp
                                    rep("1_to_3_month(s)",length(387:388)),    #subclinical mastitis and Na/K ratio,
                                    rep("0.5_to_6_months",length(389:401)),    #maternal GIT diary, ROME, prior preg med use
                                    rep("0.5_and_1_month(s)",length(402:410)), #preg and birth med use
                                    rep("1_and_3_month(s)",length(411:419)),   #maternal postpartum med use
                                    rep("0.5_to_6_months",length(420:422)))    #breastpumping

# test <- mhmo[,c(8:10,16:434)]
# for (i in 1:ncol(test)){
#   print(paste0(i, "__", colnames(test)[i]))
#   print(table(test[!is.na(test[,i]), "time_point"]))
# }

## separate phenotypes by data type
factor_mhmo <- overview_mhmo[overview_mhmo$data_type=="factor",] #75 categorical phenotypes
numeric_mhmo <- overview_mhmo[overview_mhmo$data_type=="numeric" | overview_mhmo$data_type=="integer",] #347 numeric phenotypes


### ===== 5.3 INFANT PHENOTYPES ===== ###

## Infant phenotypes: create overview
overview_ihmo <- create.overview.table(ihmo[,c(8:9,16:28)])

## add info in data_use column to indicate whether a used phenotype is static or dynamic or dynamic or 1-time-point
overview_ihmo$model_use <- c(rep("dynamic",length(1:2)),     #time_point and time_point_numeric
                             rep("static",length(3:13)),     #infant genetics, infant sex, infant birth place and mode, birthweight, infant feeding birth
                             rep("dynamic",length(14)),      #infant feeding mode
                             rep("1_time_point",length(15))) #infant food introduction

## add info on which phenotypes to exclude as x from association models
overview_ihmo$model_exclusion <- c(rep("exclude",length(1:3)),     #time_point and time_point_numeric (exclude as it will be corrected for), infant genetics (FUT2)
                                   rep("include",length(4:15)))    #infant Se status, infant sex, infant birth place and mode, birthweight, infant feeding

## add info on which HMO time points phenotypes were/will be linked to
overview_ihmo$model_time_point <- c(rep(NA,length(1:3)),                     #time_point and time_point_numeric, infant genetics (FUT2)
                                    rep("0.5_to_6_months",length(4:11)),     #infant Se status, infant sex, infant birth place and mode, birthweight
                                    rep("0.5_and_1_month(s)",length(12:13)), #infant feeding birth
                                    rep("0.5_to_6_months",length(14)),       #infant feeding mode
                                    rep("6_month(s)",length(15)))            #infant food introduction

# test <- ihmo[,c(8:9,16:28)]
# for (i in 1:ncol(test)){
#   print(paste0(i, "__", colnames(test)[i]))
#   print(table(test[!is.na(test[,i]), "time_point"]))
# }

## separate phenotypes by data type
factor_ihmo <- overview_ihmo[overview_ihmo$data_type=="factor",] #12 categorical phenotypes
numeric_ihmo <- overview_ihmo[overview_ihmo$data_type=="numeric" | overview_ihmo$data_type=="integer",] #3 numeric phenotypes


##### =========================== 6. CREATE BARPLOTS AND CREATE HISTOGRAMS TO SHOW DATA DISTRIBUTION TO GUIDE DECISION ON PHENOTYPE DATA TRANSFORMATION  =========================== #####

library(ggplot2)

### ===== 6.1 DATA DISTRIBUTION OF MATERNAL PHENOTYPES ===== ###

## create bar plots for all categorical phenotypes
pdf(paste0("phenotype_data_distribution/240122_HMO_real_categorical_maternal_phenotypes_data_distribution_barplots_2_n75.pdf"))
for (i in factor_mhmo$column_name){
  print(paste0("Creating barplots for ", i))
  barplots <- ggplot(mhmo[!is.na(mhmo[,i]),], aes(x=mhmo[!is.na(mhmo[,i]),i]))+
    geom_bar()+
    labs(i)+
    ggtitle(paste0(i))+
    theme(axis.text.x = element_text(angle = 90))
  print(barplots)
}
dev.off()

## create histograms for all numeric phenotypes
pdf(paste0("phenotype_data_distribution/240122_HMO_real_numeric_maternal_phenotypes_data_distribution_histograms_2_n347.pdf"))
for (i in numeric_mhmo$column_name){
  print(paste0("Creating histogram for ", i))
  histo <- ggplot(mhmo[!is.na(mhmo[,i]),], aes(x=mhmo[!is.na(mhmo[,i]),i]))+
    geom_histogram(bins=50)+
    labs(i)+
    ggtitle(paste0(i))
  print(histo)
}
dev.off()


### ===== 6.2 DATA DISTRIBUTION OF INFANT PHENOTYPES ===== ###

## create bar plots for all categorical phenotypes
pdf(paste0("phenotype_data_distribution/240122_HMO_real_categorical_infant_phenotypes_data_distribution_barplots_2_n12.pdf"))
for (i in factor_ihmo$column_name){
  print(paste0("Creating barplots for ", i))
  barplots <- ggplot(ihmo[!is.na(ihmo[,i]),], aes(x=ihmo[!is.na(ihmo[,i]),i]))+
    geom_bar()+
    labs(i)+
    ggtitle(paste0(i))+
    theme(axis.text.x = element_text(angle = 90))
  print(barplots)
}
dev.off()

## create histograms for all numeric phenotypes
pdf(paste0("phenotype_data_distribution/240122_HMO_real_numeric_infant_phenotypes_data_distribution_histograms_2_n3.pdf"))
for (i in numeric_ihmo$column_name){
  print(paste0("Creating histogram for ", i))
  histo <- ggplot(ihmo[!is.na(ihmo[,i]),], aes(x=ihmo[!is.na(ihmo[,i]),i]))+
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
factor_mhmo_sum <- create.summary.factors(mhmo[,colnames(mhmo) %in% factor_mhmo$column_name])

## polish perc
factor_mhmo_sum$perc_answered <- signif(factor_mhmo_sum$perc_answered, 2)
factor_mhmo_sum$perc_missing <- signif(factor_mhmo_sum$perc_missing, 2)

## change colnames
colnames(factor_mhmo_sum)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== 7.3 SUMMARY STATISTICS FOR CATEGORICAL INFANT PHENOTYPES ===== ###

## create summary
factor_ihmo_sum <- create.summary.factors(ihmo[,colnames(ihmo) %in% factor_ihmo$column_name])

## polish perc
factor_ihmo_sum$perc_answered <- signif(factor_ihmo_sum$perc_answered, 2)
factor_ihmo_sum$perc_missing <- signif(factor_ihmo_sum$perc_missing, 2)

## change colnames
colnames(factor_ihmo_sum)[1:3] <- c("Phenotype", "Categories", "n_per_category")


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
numeric_mhmo_sum <- create.summary.numeric(mhmo[,colnames(mhmo) %in% numeric_mhmo$column_name])


## Add a 'Note' column to indicate levels of categorical phenotypes that were converted to numeric for associations
numeric_mhmo_sum$Note <- NA

# matnumcatphenos

## Add note for mother_milk_collection_time
hmo$mother_milk_collection_time <- factor(hmo$mother_milk_collection_time, levels = c("1__20h", "2__21h", "3__22h", "4__23h", "5__00h",
                                                                                      "6__01h", "7__02h", "8__03h", "9__04h", "10__05h",
                                                                                      "11__06h", "12__07h", "13__08h", "14__09h", "15__10h",
                                                                                      "16__11h", "17__12h", "18__13h", "19__14h", "20__15h",
                                                                                      "21__16h", "22__17h", "23__18h", "24__19h")) # first correctly sort the levels in the hmo data frame, from which we pull the original levels
numeric_mhmo_sum[numeric_mhmo_sum$Phenotype=="mother_milk_collection_time_numeric", "Note"] <- c(paste0("This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: ",
                                                                                                        paste(collapse=";", levels(hmo$mother_milk_collection_time)))) # add note and levels

## Add note for mother_education_level
numeric_mhmo_sum[numeric_mhmo_sum$Phenotype=="mother_education_level_numeric", "Note"] <- c(paste0("This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: ",
                                                                                                   paste(collapse=";", levels(hmo$mother_education_level))))

## Add note for mother_net_income
numeric_mhmo_sum[numeric_mhmo_sum$Phenotype=="mother_net_income_numeric", "Note"] <- c(paste0("This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: ",
                                                                                              paste(collapse=";", levels(hmo$mother_net_income))))

## Add note for mother_health_dental
numeric_mhmo_sum[numeric_mhmo_sum$Phenotype=="mother_health_dental_numeric", "Note"] <- c(paste0("This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: ",
                                                                                                 paste(collapse=";", levels(hmo$mother_health_dental))))

## Check notes
# numeric_mhmo_sum[!is.na(numeric_mhmo_sum$Note),]


### ===== 8.3 SUMMARY STATISTICS FOR INFANT PHENOTYPES ===== ###

## Calculate summary statistics
numeric_ihmo_sum <- create.summary.numeric(ihmo[,colnames(ihmo) %in% numeric_ihmo$column_name])


## Add a 'Note' column to indicate levels of categorical phenotypes that were converted to numeric for associations
numeric_ihmo_sum$Note <- NA


##### =========================== 9. COMBINE SUMMARY STATISTICS TABLES FOR ALL PHENOTYPES  =========================== #####

### ===== 9.1 PREPARE SUMMARY STATISTICS TABLE WITH CATEGORICAL PHENOTYPES ===== ###

## combine maternal and infant summary statistics tables, excluding the time_point columns #85x9
factor_sum <- as.data.frame(rbind(factor_mhmo_sum[factor_mhmo_sum$Phenotype!="time_point",],
                                  factor_ihmo_sum[factor_ihmo_sum$Phenotype!="time_point",]))


## add info on data type, model and linked HMO time points from factor overview tables
# combine maternal and infant factor overview tables
factor_all <- as.data.frame(rbind(factor_mhmo, factor_ihmo))
colnames(factor_all)[1] <- "Phenotype"

# add info
library(dplyr)
factor_sum2 <- left_join(factor_sum, factor_all[,c(1:2,4:6)], by="Phenotype")


### ===== 9.2 PREPARE SUMMARY STATISTICS TABLE WITH NUMERIC PHENOTYPES ===== ###

## combine maternal and infant summary statistics tables, excluding the time_point_numeric columns #348x15
numeric_sum <- as.data.frame(rbind(numeric_mhmo_sum[numeric_mhmo_sum$Phenotype!="time_point_numeric",],
                                   numeric_ihmo_sum[numeric_ihmo_sum$Phenotype!="time_point_numeric",]))

## add info on data type, model and linked HMO time points from numeric overview tables
# combine maternal and infant numeric overview tables
numeric_all <- as.data.frame(rbind(numeric_mhmo, numeric_ihmo))
colnames(numeric_all)[1] <- "Phenotype"

# add info
# library(dplyr)
numeric_sum2 <- left_join(numeric_sum, numeric_all[,c(1:2,4:6)], by="Phenotype")


### ===== 9.3 COMBINE SUMMARY STATISTICS TABLES WITH CATEGORICAL AND NUMERIC PHENOTYPES ===== ###

## combine tables
library(plyr)
all_sum_tmp <- rbind.fill(factor_sum2, numeric_sum2)

## resort columns
all_sum <- all_sum_tmp[,c(1,9,4:8,2:3,13:19,11,10,12,20)] #433x20

## save phenotype summary statistics table
# write.table(all_sum, file="phenotype_summary_statistics/240122_HMOs_real_all_maternal_and_infant_phenotypes_summary_statistics_n433.txt", row.names=F, col.names=T, sep="\t", quote=F)


##### =========================== 10. INVERSE-RANK TRANSFORM NUMERIC PHENOTYPES AND HMO DATA =========================== #####

### ===== 10.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 10.2 INVERSE-RANK TRANSFORM NUMERIC PHENOTYPE DATA ===== ###

## Maternal phenotypes
mhmo_invr <- mhmo
mhmo_invr[,numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"]] <- as.data.frame(apply(mhmo_invr[,numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"]], 2, invrank))

colnames(mhmo_invr)[which(colnames(mhmo_invr) %in% numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"])] <- paste0(colnames(mhmo_invr)[which(colnames(mhmo_invr) %in% numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"])], "_invr")


## Infant phenotypes
ihmo_invr <- ihmo
ihmo_invr[,numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"]] <- as.data.frame(apply(ihmo_invr[,numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"]], 2, invrank))

colnames(ihmo_invr)[which(colnames(ihmo_invr) %in% numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"])] <- paste0(colnames(ihmo_invr)[which(colnames(ihmo_invr) %in% numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"])], "_invr")


### ===== 10.3 INVERSE-RANK TRANSFORM HMO DATA ===== ###

## Maternal data frame
mhmo_invr[,grep("_ugml",colnames(mhmo_invr))] <- as.data.frame(apply(mhmo_invr[,grep("ugml", colnames(mhmo_invr))], 2, invrank))
colnames(mhmo_invr)[grep("_ugml", colnames(mhmo_invr))] <- paste0(colnames(mhmo_invr)[grep("_ugml", colnames(mhmo_invr))], "_invr")

## Infant data frame
ihmo_invr[,grep("_ugml",colnames(ihmo_invr))] <- as.data.frame(apply(ihmo_invr[,grep("ugml", colnames(ihmo_invr))], 2, invrank))
colnames(ihmo_invr)[grep("_ugml", colnames(ihmo_invr))] <- paste0(colnames(ihmo_invr)[grep("_ugml", colnames(ihmo_invr))], "_invr")


##### =========================== 11. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED NUMERIC PHENOTYPES =========================== #####

# library(ggplot2)

### ===== 11.1 FUNCTION TO MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED NUMERIC PHENOTYPES ===== ###

## make histograms for inverse-rank transformed numeric phenotype data
check.distribution.invr.pheno <- function(inputdata, phenotype_columns, outputname){
  pdf(paste0("phenotype_data_distribution/", outputname, ".pdf"))
  
  for (i in c(phenotype_columns)){
    print(paste0("Creating histogram for ", i))
    hmohist <- ggplot(inputdata,
                      aes(x=as.numeric(inputdata[,i]))) + geom_histogram() + labs(x=i)
    print(hmohist)
  }
  
  dev.off()
}


### ===== 11.2 MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED NUMERIC PHENOTYPES ===== ###

## Maternal phenotypes
mhmo_invr_phenotype_columns <- paste0(numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"], "_invr")
check.distribution.invr.pheno(mhmo_invr,
                              phenotype_columns = mhmo_invr_phenotype_columns,
                              outputname = "240122_HMO_real_numeric_maternal_phenotypes_data_distribution_histograms_3_after_invr_transformation_n346")

# Infant phenotypes
ihmo_invr_phenotype_columns <- paste0(numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"], "_invr")
check.distribution.invr.pheno(ihmo_invr,
                              phenotype_columns = ihmo_invr_phenotype_columns,
                              outputname="240122_HMO_real_numeric_infant_phenotypes_data_distribution_histograms_3_after_invr_transformation_n2")


##### =========================== 12. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED HMO DATA =========================== #####

### ===== 12.1 FUNCTION TO MAKE HISTOGRAMS FOR INVERSE-RANK TRANSFORMED HMO CONCENTRATIONS PER MILK GROUP ===== ###

check.distribution.invr.hmo <- function(inputdata, milk_group, time_point, hmo_columns, outputname){
  print(paste0("Creating histograms for ", milk_group, " at ", time_point, " time point(s)"))
  
  pdf(paste0("HMO_data_distribution/", outputname, "_", milk_group, "_", time_point, ".pdf"))
  
  for (i in c(hmo_columns)){
    hmohist <- ggplot(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group,],
                      aes(x=as.numeric(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group, i]))) + geom_histogram() + labs(x=colnames(inputdata)[i])
    print(hmohist)
  }
  
  dev.off()
}


### ===== 12.2 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED HMO DATA FROM MATERNAL PHENOTYPE FILE ===== ###

## 1. Create histograms for inverse-rank transformed HMO concentrations from all milk groups combined and all time points combined
pdf(paste0("HMO_data_distribution/240122_HMO_invr_real_levels_for_maternal_phenotypes_histograms_all_milk_groups.pdf"))
for (i in c(439:466)){
  hmohist <- ggplot(mhmo_invr,
                    aes(x=as.numeric(mhmo_invr[,i]))) + geom_histogram() + labs(x=colnames(mhmo_invr)[i])
  print(hmohist)
}
dev.off()

## 2. Create histograms for inverse-rank transformed HMO concentrations per milk group, with all time points combined
check.distribution.invr.hmo(mhmo_invr, milk_group="Le-Se-", time_point="all", hmo_columns=c(439:466), outputname="240122_HMO_invr_real_levels_for_maternal_phenotypes_histograms")
check.distribution.invr.hmo(mhmo_invr, milk_group="Le-Se+", time_point="all", hmo_columns=c(439:466), outputname="240122_HMO_invr_real_levels_for_maternal_phenotypes_histograms")
check.distribution.invr.hmo(mhmo_invr, milk_group="Le+Se-", time_point="all", hmo_columns=c(439:466), outputname="240122_HMO_invr_real_levels_for_maternal_phenotypes_histograms")
check.distribution.invr.hmo(mhmo_invr, milk_group="Le+Se+", time_point="all", hmo_columns=c(439:466), outputname="240122_HMO_invr_real_levels_for_maternal_phenotypes_histograms")


### ===== 12.3 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED HMO DATA FROM INFANT PHENOTYPE FILE ===== ###

## 1. Create histograms for inverse-rank transformed HMO concentrations from all milk groups combined and all time points combined
pdf(paste0("HMO_data_distribution/240122_HMO_invr_real_levels_for_infant_phenotypes_histograms_all_milk_groups.pdf"))
for (i in c(33:60)){
  hmohist <- ggplot(ihmo_invr,
                    aes(x=as.numeric(ihmo_invr[,i]))) + geom_histogram() + labs(x=colnames(ihmo_invr)[i])
  print(hmohist)
}
dev.off()

## 2. Create histograms for inverse-rank transformed HMO concentrations per milk group, with all time points combined
check.distribution.invr.hmo(ihmo_invr, milk_group="Le-Se-", time_point="all", hmo_columns=c(33:60), outputname="240122_HMO_invr_real_levels_for_infant_phenotypes_histograms")
check.distribution.invr.hmo(ihmo_invr, milk_group="Le-Se+", time_point="all", hmo_columns=c(33:60), outputname="240122_HMO_invr_real_levels_for_infant_phenotypes_histograms")
check.distribution.invr.hmo(ihmo_invr, milk_group="Le+Se-", time_point="all", hmo_columns=c(33:60), outputname="240122_HMO_invr_real_levels_for_infant_phenotypes_histograms")
check.distribution.invr.hmo(ihmo_invr, milk_group="Le+Se+", time_point="all", hmo_columns=c(33:60), outputname="240122_HMO_invr_real_levels_for_infant_phenotypes_histograms")


##### =========================== 13. PREPARE MATERNAL AND INFANT DATA SETS FOR ASSOCIATION OF PHENOTYPES WITH HMOs =========================== #####

### ===== 13.1 OVERVIEW OF CURRENT DATA FRAMES PREPARED FOR ASSOCIATIONS OF HMOs AND PHENOTYPE ===== ###

## Now the data frames are ready to be used for the association models:
# Data                                   Name       Size
# Maternal phenotypes, real HMO data:    mhmo_invr: 1542x466
# Infant phenotypes,   real HMO data:    ihmo_invr: 1563x60


### ===== 13.2 PREPARE DATA FRAMES WITH MATERNAL PHENOTYPES ===== ###

## Create data frame from mhmo_invr for use in association models for maternal static phenotypes
rownames(mhmo_invr) <- mhmo_invr$mother_sample_ID

overview_mhmo[overview_mhmo$model_exclusion=="include" & (overview_mhmo$data_type=="numeric" | overview_mhmo$data_type=="integer"), "column_name"] <-
  paste0(overview_mhmo[overview_mhmo$model_exclusion=="include" & (overview_mhmo$data_type=="numeric" | overview_mhmo$data_type=="integer"), "column_name"], "_invr")


### ===== 13.3 PREPARE DATA FRAMES WITH INFANT PHENOTYPES ===== ###

## Create data frame from ihmo_invr for use in association models for infant static phenotypes
rownames(ihmo_invr) <- ihmo_invr$infant_sample_ID

overview_ihmo[overview_ihmo$model_exclusion=="include" & (overview_ihmo$data_type=="numeric" | overview_ihmo$data_type=="integer"), "column_name"] <-
  paste0(overview_ihmo[overview_ihmo$model_exclusion=="include" & (overview_ihmo$data_type=="numeric" | overview_ihmo$data_type=="integer"), "column_name"], "_invr")


##### =========================== 14. MODEL FOR ASSOCIATION OF MULTIPLE-TIME-POINT STATIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID =========================== #####

### ===== 14.1 PREPARE DATA FRAMES WITH MATERNAL STATIC PHENOTYPES ===== ###

# selected basic info
basic_mhmo_invr <- mhmo_invr[,c(1:9,11:15,436,439:466)] #IDs, milk group and invr-transformed HMO concentrations

# select maternal static phenotypes that should be included
v_static_mhmo_invr <- overview_mhmo[overview_mhmo$participant=="mother" & overview_mhmo$model_use=="static" & overview_mhmo$model_exclusion=="include","column_name"]
pheno_mhmo_invr_static <- mhmo_invr[,colnames(mhmo_invr) %in% v_static_mhmo_invr]

# combine in 1 data frame that is to be used for models
mhmo_invr_static <- merge(basic_mhmo_invr, pheno_mhmo_invr_static, by="row.names")
mhmo_invr_static <- mhmo_invr_static[,-1] #remove Row.names column
#1542x436, of which 393 are phenotypes
#str(mhmo_invr_static)

## Maternal data set (mhmo_invr_static)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = milk HMO conc.
## 44:436 = maternal static phenotypes


### ===== 14.2 PREPARE DATA FRAMES WITH INFANT STATIC PHENOTYPES ===== ###

# selected basic info
basic_ihmo_invr <- ihmo_invr[,c(1:9,11:15,30,33:60)] #IDs, milk group and invr-transformed HMO concentrations

# select infant static phenotypes that should be included
v_static_ihmo_invr <- overview_ihmo[overview_ihmo$participant=="infant" & overview_ihmo$model_use=="static" & overview_ihmo$model_exclusion=="include","column_name"]
pheno_ihmo_invr_static <- ihmo_invr[,colnames(ihmo_invr) %in% v_static_ihmo_invr]

# combine in 1 data frame that is to be used for models
ihmo_invr_static <- merge(basic_ihmo_invr, pheno_ihmo_invr_static, by="row.names")
ihmo_invr_static <- ihmo_invr_static[,-1] #remove Row.names column
#1563x53, of which 10 are phenotypes
#str(ihmo_invr_static)

## Infant data set (ihmo_invr_static)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = milk HMO conc.
## 44:53 = infant static phenotypes


### ===== 14.3 FUNCTION FOR ASSOCIATION OF STATIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID ===== ###

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


### ===== 14.4 RUN ASSOCIATION OF MATERNAL STATIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS ===== ###

## Maternal data set (mhmo_invr_static)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = milk HMO conc.
## 44:436 = maternal static phenotypes (maternal food groups: 147:174)

mhmo_static_results <- run_mmrm_analysis(bacteria   = mhmo_invr_static[,c(16:43)], #inverse-rank transformed HMO concentrations, all: 16:43
                                         phenotypes = mhmo_invr_static[,c(44:436)], #maternal static phenotypes, all 44:436
                                         time       = mhmo_invr_static[,9],  #time_point_numeric
                                         time_cat   = mhmo_invr_static[,8],  #time_point (factor)
                                         NEXT_ID    = mhmo_invr_static[,6],  #mother_ID (do not use family_ID here - doesn't work)
                                         covariates = mhmo_invr_static[,15,drop = F]) #correct for mother_milk_HMO_milk_group; try removing drop=F if it doesn't work
mhmo_static_results_trait <- mhmo_static_results$trait #11256x12
mhmo_static_results_time.trait <- mhmo_static_results$time.trait #11256x11
write.table(mhmo_static_results_trait, file="association_model_results/240123_HMOs_invr_real_maternal_static_phenotypes_model_results_by_level.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(mhmo_static_results_time.trait, file="association_model_results/240123_HMOs_invr_real_maternal_static_phenotypes_model_results_by_level_and_time.txt", row.names=F, col.names=T, sep="\t", quote=F)
# mhmo_static_results_trait <- read.table("association_model_results/240123_HMOs_invr_real_maternal_static_phenotypes_model_results_by_level.txt", header=T, sep="\t", stringsAsFactors=T)
# mhmo_static_results_time.trait <- read.table("association_model_results/240123_HMOs_invr_real_maternal_static_phenotypes_model_results_by_level_and_time.txt", header=T, sep="\t", stringsAsFactors=T)

mhmo_static_results_trait$dataset <- "mother_real_HMOs_static"


### ===== 14.5 RUN ASSOCIATION OF INFANT STATIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS ===== ###

## Infant data set (ihmo_invr_static)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43 = milk HMO conc.
## 44:53 = infant static phenotypes

ihmo_static_results <- run_mmrm_analysis(bacteria   = ihmo_invr_static[,c(16:43)], #inverse-rank transformed HMO concentrations
                                         phenotypes = ihmo_invr_static[,c(44:53)], #infant static phenotypes 44:53
                                         time       = ihmo_invr_static[,9],  #time_point_numeric
                                         time_cat   = ihmo_invr_static[,8],  #time_point (factor)
                                         NEXT_ID    = ihmo_invr_static[,7],  #infant_ID (do not use family_ID here - doesn't work)
                                         covariates = ihmo_invr_static[,15,drop = F]) #correct for mother_milk_HMO_milk_group; try removing drop=F if it doesnt work
ihmo_static_results_trait <- ihmo_static_results$trait #504x12
ihmo_static_results_trait.time <- ihmo_static_results$time.trait #504x11
write.table(ihmo_static_results_trait, file="association_model_results/240123_HMOs_invr_real_infant_static_phenotypes_model_results_by_level.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(ihmo_static_results_trait.time, file="association_model_results/240123_HMOs_invr_real_infant_static_phenotypes_model_results_by_level_and_time.txt", row.names=F, col.names=T, sep="\t", quote=F)
# ihmo_static_results_trait <- read.table("association_model_results/240123_HMOs_invr_real_infant_static_phenotypes_model_results_by_level.txt", header=T, sep="\t", stringsAsFactors=T)
# ihmo_static_results_trait.time <- read.table("association_model_results/240123_HMOs_invr_real_infant_static_phenotypes_model_results_by_level_and_time.txt", header=T, sep="\t", stringsAsFactors=T)

ihmo_static_results_trait$dataset <- "infant_real_HMOs_static"


##### =========================== 15. MODEL FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP =========================== #####

### ===== 15.1 PREPARE DATA FRAMES WITH MATERNAL 1-TIME-POINT PHENOTYPES ===== ###

## No maternal phenotypes are only to be associated with milk HMO levels at 1 time point.


### ===== 15.2 PREPARE DATA FRAMES WITH INFANT 1-TIME-POINT PHENOTYPES ===== ###

## Create data frame from ihmo_invr for use in association models for infant static phenotypes
rownames(ihmo_invr) <- ihmo_invr$infant_sample_ID

# select infant 1-time-point phenotypes that should be included -> onlly 1 column: infant_ffq_solid_food_introduced is only linked to milk HMOs at 6 months
v_1tp_ihmo_invr <- overview_ihmo[overview_ihmo$participant=="infant" & overview_ihmo$model_use=="1_time_point" & overview_ihmo$model_exclusion=="include","column_name"]

# selected basic info and 1-time-point phenotype
ihmo_invr_1tp <- ihmo_invr[,c(1:9,11:15,30,33:60,28)] #IDs, milk group and invr-transformed HMO concentrations, 28=infant_ffq_solid_food_introduced
#1563x44, of which 1 is a phenotype
#str(ihmo_invr_1tp)

## Infant data set (ihmo_invr_1tp)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = invr-transformed milk HMOs
## 44     = infant 1-time-point phenotype(s)


### ===== 15.3 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT STATIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP ===== ###

run.lmer.cor.milkgroup.1tp <- function(datasetname, inputdata, xcolumn, ycolumn){
  
  levels <- c()
  n_levels <- c()
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
      
      #second mixed model with correction for milk_group as fixed effects and WITH the phenotype of interest (i)
      m1 <- lm(my_df[,j] ~ mother_milk_HMO_milk_group + my_df[,i], data=my_df)
      sum1 <- summary(m1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      levels <- c(levels, paste(collapse=";",levels(my_df[,i])))
      n_levels <- c(n_levels, paste(collapse=";",table(my_df[,i])))
      e <- c(e, c(sum1$coefficients[grep("my_df",rownames(sum1$coefficients)), "Estimate"])) #first rows are for milk groups, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      p <- c(p, c(sum1$coefficients[grep("my_df",rownames(sum1$coefficients)), "Pr(>|t|)"])) #first rows are for milk groups, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, paste0("lm_1_time_point_correction_milkgroup"))
      d <- c(d, datasetname)
      x <- c(x, colnames(inputdata)[i])
      y <- c(y, colnames(inputdata)[j])
      n <- c(n, nrow(my_df))
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, n_levels=n_levels, estimate=e, p=p)
  
  return(res)
}


### ===== 15.4 RUN ASSOCIATION OF MATERNAL SINGLE-TIME-POINT STATIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS ===== ###

## No maternal phenotypes are only to be associated with milk HMO levels at 1 time point.


### ===== 15.5 RUN ASSOCIATION OF INFANT SINGLE-TIME-POINT STATIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS ===== ###

## Infant data set (ihmo_invr_1tp)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = invr-transformed milk HMOs
## 44     = infant 1-time-point phenotype(s)

ihmo_1tp_results <- run.lmer.cor.milkgroup.1tp(datasetname = "infant_real_HMOs_1_time_point",
                                               inputdata   = ihmo_invr_1tp,
                                               xcolumn     = c(44),      #single time point infant phenotypes (categorical, only 2 levels)
                                               ycolumn     = c(16:43))   #HMO concentrations
ihmo_1tp_results$Note <- NA
# write.table(ihmo_1tp_results, file="association_model_results/240123_HMOs_invr_real_infant_1_time_point_phenotypes_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# ihmo_1tp_results <- read.table(file="association_model_results/240123_HMOs_invr_real_infant_1_time_point_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


##### =========================== 16. MODEL FOR ASSOCIATION OF MULTIPLE-TIME-POINT DYNAMIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID =========================== #####

### ===== 16.1 PREPARE DATA FRAMES WITH MATERNAL DYNAMIC PHENOTYPES ===== ###

# selected basic info
basic_mhmo_invr <- mhmo_invr[,c(1:9,11:15,436,439:466)] #IDs, milk group and invr-transformed HMO concentrations

# select maternal dynamic phenotypes that should be included
v_dynamic_mhmo_invr <- overview_mhmo[overview_mhmo$participant=="mother" & overview_mhmo$model_use=="dynamic" & overview_mhmo$model_exclusion=="include","column_name"]
pheno_mhmo_invr_dynamic <- mhmo_invr[,colnames(mhmo_invr) %in% v_dynamic_mhmo_invr]

# combine in 1 data frame that is to be used for models
mhmo_invr_dynamic <- merge(basic_mhmo_invr, pheno_mhmo_invr_dynamic, by="row.names")
mhmo_invr_dynamic <- mhmo_invr_dynamic[,-1] #remove Row.names column
#1542x62, of which 19 are phenotypes
#str(mhmo_invr_dynamic)

## Maternal data set (mhmo_invr_dynamic)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = milk HMO conc.
## 44:62 = maternal dynamic phenotypes


### ===== 16.2 PREPARE DATA FRAMES WITH INFANT DYNAMIC PHENOTYPES ===== ###

# select infant dynamic phenotypes that should be included: only 1 phenpotype: infant_ffq_feeding_mode
v_dynamic_ihmo_invr <- overview_ihmo[overview_ihmo$participant=="infant" & overview_ihmo$model_use=="dynamic" & overview_ihmo$model_exclusion=="include","column_name"]

# selected basic info and the 1 dynamic infant phenotype (col 27)
ihmo_invr_dynamic <- ihmo_invr[,c(1:9,11:15,30,33:60,27)] #IDs, milk group and invr-transformed HMO concentrations
#1563x44, of which 1 is a phenotype
#str(ihmo_invr_dynamic)

## Infant data set (ihmo_invr_dynamic)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = milk HMO conc.
## 44     = infant dynamic phenotypes


### ===== 16.3 FUNCTION FOR ASSOCIATION OF DYNAMIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID ===== ###

run.lmer.cor.milkgroup.time.ID <- function(datasetname, inputdata, xcolumn, ycolumn){
  
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
      m1 <- lmerTest::lmer(my_df[,j] ~ mother_milk_HMO_milk_group + time_point_numeric + (1|mother_ID) + my_df[,i], REML=F, data=my_df)
      sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
      sum1$rows <- rownames(sum1)
      # print(sum1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for milk groups and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for milk groups and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for milk groups and time points, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, rep(paste0("lmer_correction_milkgroup_time_ID"), length(sum1[grep("my_df",sum1$rows),"rows"])))
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


### ===== 16.4 RUN ASSOCIATION OF MATERNAL DYNAMIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS ===== ###

## Maternal data set (mhmo_invr_dynamic)
## Columns:
## 1:14 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = milk HMO conc.
## 44:62 = maternal dynamic phenotypes

mhmo_dynamic_results <- run.lmer.cor.milkgroup.time.ID(datasetname = "mother_real_HMOs_dynamic",
                                                       inputdata   = mhmo_invr_dynamic,
                                                       xcolumn     = c(44:62),   #dynamic maternal phenotypes (3 have >2 categories: mother_milk_collection_season, mother_milk_collection_month, mother_breastpump_brand)
                                                       ycolumn     = c(16:43))   #HMO concentrations
mhmo_dynamic_results$Note <- NA
# write.table(mhmo_dynamic_results, file="association_model_results/240123_HMOs_invr_real_maternal_dynamic_phenotypes_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# mhmo_dynamic_results <- read.table(file="association_model_results/240123_HMOs_invr_real_maternal_dynamic_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 16.5 RUN ASSOCIATION OF INFANT DYNAMIC PHENOTYPES WITH REAL MILK HMO CONCENTRATIONS ===== ###

## Infant data set (ihmo_invr_dynamic)
## Columns:
## 1:14 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 15   = milk_group
## 16:43  = milk HMO conc.
## 44     = infant dynamic phenotypes

ihmo_dynamic_results <- run.lmer.cor.milkgroup.time.ID(datasetname = "infant_real_HMOs_dynamic",
                                                       inputdata   = ihmo_invr_dynamic,
                                                       xcolumn     = c(44),     #dynamic infant phenotypes (categorical, only 2 levels)
                                                       ycolumn     = c(16:43))  #HMO concentrations
ihmo_dynamic_results$Note <- NA
# write.table(ihmo_dynamic_results, file="association_model_results/240123_HMOs_invr_real_infant_dynamic_phenotypes_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# ihmo_dynamic_results <- read.table(file="association_model_results/240123_HMOs_invr_real_infant_dynamic_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


##### =========================== 17. COMBINE RESULTS AND CORRECT FOR MULTIPLE TESTING =========================== #####

### ===== 17.1 PREPARE DATA FRAMES FOR MERGING ===== ###

## Results data frames      Size
# mhmo_static_results_trait 11256 x 13
# ihmo_static_results_trait   504 x 13
# 
# #mhmo_1tp_results
# ihmo_1tp_results             28 x 10
# 
# mhmo_dynamic_results        924 x  9
# ihmo_dynamic_results         28 x  9

## if necessary, import these data frames:

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/")

mhmo_static_results_trait <- read.table("association_model_results/240123_HMOs_invr_real_maternal_static_phenotypes_model_results_by_level.txt", header=T, sep="\t", stringsAsFactors=T)
mhmo_static_results_trait$dataset <- "mother_real_HMOs_static"
# mhmo_static_results_time.trait <- read.table("association_model_results/240123_HMOs_invr_real_maternal_static_phenotypes_model_results_by_level_and_time.txt", header=T, sep="\t", stringsAsFactors=T)

ihmo_static_results_trait <- read.table("association_model_results/240123_HMOs_invr_real_infant_static_phenotypes_model_results_by_level.txt", header=T, sep="\t", stringsAsFactors=T)
ihmo_static_results_trait$dataset <- "infant_real_HMOs_static"
# ihmo_static_results_trait.time <- read.table("association_model_results/240123_HMOs_invr_real_infant_static_phenotypes_model_results_by_level_and_time.txt", header=T, sep="\t", stringsAsFactors=T)

ihmo_1tp_results <- read.table(file="association_model_results/240123_HMOs_invr_real_infant_1_time_point_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)

mhmo_dynamic_results <- read.table(file="association_model_results/240123_HMOs_invr_real_maternal_dynamic_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)
ihmo_dynamic_results <- read.table(file="association_model_results/240123_HMOs_invr_real_infant_dynamic_phenotypes_model_results.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 17.2 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM STATIC PHENOTYPE ASSOCIATIONS ===== ###

results_static <- as.data.frame(rbind(mhmo_static_results_trait, ihmo_static_results_trait)) #11760 x 13

## select relevant rows and columns
# exclude results from phenotype that could not always run (model could not be built)
# results_static[is.na(results_static$P),]
results_static_1 <- results_static[results_static$trait!="mother_ffq_item_sugar_in_tea_grams_per_day_invr",] #this excludes 28 rows, 11732 x 13
results_static_2 <- results_static_1[results_static_1$trait!="mother_weight_gain_preg_kg_invr",] #this excludes 28 rows, 11704 x 13
# results_static_2[is.na(results_static_2$P),]

# only select relevant columns
results_static_3 <- results_static_2[,c(2:7,11,13)]

## fix column names
colnames(results_static_3) <- c("statistic", "y", "x", "n_total", "levels", "estimate", "p", "dataset")

## add missing columns
results_static_3$Note <- NA

## resort columns
sel_results_static <- results_static_3[,c(8,1,3,2,4:7,9)] #11704x9


### ===== 17.3 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM 1-TIME-POINT PHENOTYPE ASSOCIATIONS ===== ###

# results_1tp <- as.data.frame(rbind(mhmo_1tp_results, ihmo_1tp_results))
results_1tp <- ihmo_1tp_results

## ensure only the tested level shows in the levels column and not the reference
results_1tp$levels <- gsub("no;", "", results_1tp$levels)

## select relevant columns
sel_results_1tp <- results_1tp[,c(1:6,8:10)] #28x9


### ===== 17.4 MERGE AND PREPARE DATA FRAMES WITH RESULTS FROM DYNAMIC PHENOTYPE ASSOCIATIONS ===== ###

results_dynamic <- as.data.frame(rbind(mhmo_dynamic_results, ihmo_dynamic_results))
sel_results_dynamic <- results_dynamic #952x9


### ===== 17.5 MERGE ALL RESULTS DATA FRAMES AND SPLIT BY OUTPUT TYPE (ALPHA DIVERSITY / RELATIVE ABUNDANCES) ===== ###

## merge all results
model_results <- as.data.frame(rbind.fill(sel_results_static, sel_results_1tp, sel_results_dynamic)) #12684x9


### ===== 17.6 EXCLUDE PHENOTYPES ===== ###

## exclude all food preferences and all individual food items
model_results_short_1 <- model_results[-grep("mother_food_pref", model_results$x),] #exclude 85 food preferences x28 (this excludes 2380 rows)
model_results_short_2 <- model_results_short_1[-grep("mother_ffq_item", model_results_short_1$x),] #exclude 162 individual food items x28 (this excludes 4536 rows)

## exclude selected detailed phenotypes (n=12) and exclude mother_milk_collection_volume_pumped_invr as it is highly correlated with time postpartum
exc_phenos <- c("mother_genetics_blood_group_genotype",
                "mother_ffq_fish_seafood_consumption_yes_no", "mother_ffq_meat_consumption_yes_no",
                "mother_exp_pets_cat", "mother_exp_pets_dog",
                "mother_birth_gestational_age_categories",
                "infant_genetics_FUT2_secretor", "infant_Se_status_same_as_mother",
                "infant_birth_delivery_place_detailed", "infant_birth_delivery_mode_detailed",
                "mother_milk_collection_month",
                "mother_dis_subclinical_mastitis_Na_K_ratio_invr","mother_milk_collection_volume_pumped_invr")

model_results_short_3 <- model_results_short_2[!(model_results_short_2$x %in% exc_phenos),] #exclude similar or detailed phenotypes (this excludes 924 rows) and phenotypes correlated with time point postpartum (this excludes another 28 rows)

model_results_short <- model_results_short_3 ## 4816x9


### ===== 17.7 ADD PHENOTYPE GROUPS FOR PLOTTING ===== ###

## show all included unique phentypes
# length(unique(model_results_short$x)) #163
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
model_results_short[model_results_short$x=="mother_breastpump_yes_no", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="mother_breastpump_brand", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="mother_breastpump_times_per_day_invr", "phenotype_group"] <- "Maternal_lifestyle_and_exposures"
model_results_short[model_results_short$x=="infant_ffq_feeding_mode", "phenotype_group"] <- "Infant_feeding"


# table(model_results_short$phenotype_group, useNA="ifany")
# length(unique(model_results_short$x)) #162
# length(unique(model_results_short$phenotype_group)) #9

phenos_info <- model_results_short[,c(3,10)]
phenos_info <- phenos_info[!duplicated(phenos_info),]
table(phenos_info$phenotype_group, useNA="ifany")
# Human_milk_collection               Infant_feeding 
# 3                                   4 
# Maternal_and_infant_anthropometrics Maternal_diet 
# 8                                   72 
# Maternal_genetics                   Maternal_health_and_diseases 
# 1                                   21 
# Maternal_lifestyle_and_exposures    Maternal_medication_use 
# 12                                  22 
# Pregnancy_and_birth_characteristics 
# 19

unique(model_results_short[model_results_short$phenotype_group=="Human_milk_collection", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Infant_feeding", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_and_infant_anthropometrics", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_diet", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_genetics", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_health_and_diseases", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_lifestyle_and_exposures", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Maternal_medication_use", "x"])
unique(model_results_short[model_results_short$phenotype_group=="Pregnancy_and_birth_characteristics", "x"])

write.table(phenos_info, file="association_model_results/240229_phenotypes_associated_with_milk_HMOs_real_n162.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 17.8 ADD PRETTY LABELS FOR PHENOTYPES FOR PLOTTING ===== ###

## ensure results for traits with >2 categories are correctly displayed (combine x with level)
model_results_short$trait_name_in_plot <- NA
model_results_short[model_results_short$levels!="", "trait_name_in_plot"] <- c(paste0(model_results_short[model_results_short$levels!="", "x"],
                                                                                      "_", model_results_short[model_results_short$levels!="", "levels"]))
model_results_short[is.na(model_results_short$trait_name_in_plot), "trait_name_in_plot"] <- as.character(as.factor(model_results_short[is.na(model_results_short$trait_name_in_plot), "x"]))
# length(unique(model_results_short$trait_name_in_plot)) #172

## create column with polished trait names for plotting
model_results_short$trait_name_in_plot <- gsub("_invr", "", model_results_short$trait_name_in_plot)
model_results_short$trait_name_in_plot <- as.factor(as.character(model_results_short$trait_name_in_plot))

levels(model_results_short$trait_name_in_plot) <- c("infant_birth_delivery_mode_C-section" = "delivery_mode_(C-section_vs_vaginal_birth)",
                                                    "infant_birth_delivery_place_hospital" = "delivery_place_(hospital_vs_home)",
                                                    "infant_birthweight_g" = "infant_birth_weight_(grams)",
                                                    "infant_ffq_breastfeeding_behaviour_after_birth_1bad_5excellent" = "excellent_infant_breastfeeding_behaviour_after_birth",
                                                    "infant_ffq_feeding_mode_mixed_feeding" = "feeding_mode_(mixed_feeding_vs_human_milk_feeding)",
                                                    "infant_ffq_feeding_type_delivery_breast_milk_via_bottle" = "feeding_mode_after_birth_(breast_milk_via_bottle_vs_breastfeeding)",             
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
                                                    "mother_breastpump_yes_no_yes" = "breast_pump_use",
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
                                                    "mother_med_after_birth_iron_preparations_B03A_yes" = "iron_preparations_after_birth_(B03A)",
                                                    "mother_med_after_birth_oxytocin_yes" = "oxytocin_after_birth",
                                                    "mother_med_after_birth_paracetamol_N02BE01_yes" = "paracetamol_after_birth_(N02BE01)",
                                                    "mother_med_after_birth_uterotonics_G02A_yes" = "uterotonics_after_birth_(G02A)",
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
                                                    "mother_med_gynaecological_antiinfectives_G01A_yes" = "gynaecological_antiinfectives_(G01A)",
                                                    "mother_med_iron_preparations_B03A_yes" = "iron_preparations_(B03A)",
                                                    "mother_med_osmotically_acting_laxatives_A06AD_yes" = "osmotically_acting_laxatives_(A06AD)",
                                                    "mother_med_preg_folic_acid_use_during_first_trimester_yes" = "folic_acid_use_during_first_trimester",
                                                    "mother_med_prior_preg_folic_acid_use_yes" = "folic_acid_use_before_pregnancy",
                                                    "mother_milk_collection_breasts_sampled_both_breasts" = "sampled_from_both_breasts",
                                                    "mother_milk_collection_season_Autumn" = "collection_season_(autumn_vs_spring)",
                                                    "mother_milk_collection_season_Summer" = "collection_season_(summer_vs_spring)",
                                                    "mother_milk_collection_season_Winter" = "collection_season_(winter_vs_spring)",
                                                    "mother_milk_collection_time_numeric" = "time_of_the_day",
                                                    # "mother_milk_collection_volume_pumped" = "pumped_milk_volume_(ml)",
                                                  
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


### ===== 17.9 ENSURE CORRECT N FOR PHENOTYPES WITH >2 CATEGORIES ===== ###

for (i in c(1,6,10)){model_results_short[,i] <- as.factor(as.character(model_results_short[,i]))}
levels(model_results_short$levels)[1] <- NA

### MATERNAL BLOOD GROUP ###

all_sum[all_sum$Phenotype=="mother_blood_group",c(1,4,8,9)]
# Phenotype                n_answered Categories n_per_category
# mother_blood_group       1204       A;AB;B;O   538;37;90;539
  
model_results_short[model_results_short$x=="mother_blood_group" & model_results_short$levels=="A", "n_total"] <- 538+539
model_results_short[model_results_short$x=="mother_blood_group" & model_results_short$levels=="AB", "n_total"] <- 37+539
model_results_short[model_results_short$x=="mother_blood_group" & model_results_short$levels=="B", "n_total"] <- 90+539


### BREAST PUMP BRAND

all_sum[all_sum$Phenotype=="mother_breastpump_brand",c(1,4,8,9)]
# Phenotype                      n_answered      Categories                  n_per_category
# mother_breastpump_brand        253             Ardo;Medela;Other;Philips   20;60;25;148

model_results_short[model_results_short$x=="mother_breastpump_brand" & model_results_short$levels=="Ardo", "n_total"] <- 20+148
model_results_short[model_results_short$x=="mother_breastpump_brand" & model_results_short$levels=="Medela", "n_total"] <- 60+148
model_results_short[model_results_short$x=="mother_breastpump_brand" & model_results_short$levels=="Other", "n_total"] <- 25+148


### MILK COLLECTION SEASON

all_sum[all_sum$Phenotype=="mother_milk_collection_season",c(1,4,8,9)]
# Phenotype                           n_answered     Categories                  n_per_category
# mother_milk_collection_season       1540           Autumn;Spring;Summer;Winter 393;422;354;371

model_results_short[model_results_short$x=="mother_milk_collection_season" & model_results_short$levels=="Autumn", "n_total"] <- 393+422
model_results_short[model_results_short$x=="mother_milk_collection_season" & model_results_short$levels=="Summer", "n_total"] <- 354+422
model_results_short[model_results_short$x=="mother_milk_collection_season" & model_results_short$levels=="Winter", "n_total"] <- 371+422


### INFANT FEEDING TYPE AFTER DELIVERY

all_sum[all_sum$Phenotype=="infant_ffq_feeding_type_delivery",c(1,4,8,9)]
# Phenotype                               n_answered          Categories                                               n_per_category
# infant_ffq_feeding_type_delivery        430                 breast_milk_via_bottle;breastfeeding;mixed_feeding       3;404;23

model_results_short[model_results_short$x=="infant_ffq_feeding_type_delivery" & model_results_short$levels=="breast_milk_via_bottle", "n_total"] <- 3+404
model_results_short[model_results_short$x=="infant_ffq_feeding_type_delivery" & model_results_short$levels=="mixed_feeding", "n_total"] <- 23+404


### LIVING SITUATION

all_sum[all_sum$Phenotype=="mother_exp_living_situation",c(1,4,8,9)]
# Phenotype                         n_answered     Categories                              n_per_category
# mother_exp_living_situation       1167           big_house;farm;flat;other;small_house   534;49;45;109;430

model_results_short[model_results_short$x=="mother_exp_living_situation" & model_results_short$levels=="farm", "n_total"] <- 534+49
model_results_short[model_results_short$x=="mother_exp_living_situation" & model_results_short$levels=="flat", "n_total"] <- 534+45
model_results_short[model_results_short$x=="mother_exp_living_situation" & model_results_short$levels=="other", "n_total"] <- 534+109
model_results_short[model_results_short$x=="mother_exp_living_situation" & model_results_short$levels=="small_house", "n_total"] <- 534+430


### ===== 17.10 CORRECT FOR MULTIPLE TESTING AND SAVE RESULTS TABLES ===== ###

## calculate FDR
model_results_short$FDR <- p.adjust(model_results_short$p, method="BH")

## reorder data frame by FDR results
model_results_short <- model_results_short[order(model_results_short$FDR, model_results_short$p),] 

## check results
nrow(model_results_short[model_results_short$FDR<0.05,]) # ->  12 FDR significant associations
nrow(model_results_short[model_results_short$p<0.05,])   # -> 419 nominally significant associations
length(unique(model_results_short$x)) #162

## save results #4844x12
# write.table(model_results_short, file="association_model_results/240229_model_results_milk_HMOs_real_by_phenotypes_n162.txt", row.names=F, col.names=T, sep="\t", quote=F)
# saveRDS(model_results_short, file="association_model_results/240229_model_results_milk_HMOs_real_by_phenotypes_n162.rds")
# model_results_short <- readRDS(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results/240229_model_results_milk_HMOs_real_by_phenotypes_n162.rds")


##### =========================== 18. PREPARE SUPPLEMENTARY TABLES (PHENOTYPE SUMMARY STATISTICS TABLE AND MODEL RESULTS TABLES)  =========================== #####

### ===== 18.1 PHENOTYPE SUMMARY STATISTICS TABLE (TABLE S7) ===== ###

## desired rows in phenotype summary statistics supplementary table
## -> only include rows for phenotypes that were included in the analysis / model results table

## select only relevant rows from all_sum table
# length(unique(model_results_short$x)) #162
p <- unique(model_results_short$x) #vector with all 162 included phenotypes
p <- gsub("_invr", "", p) #ensure the phenotypes have the same colnames as in the all_sum table (remove '_invr' extension for transformed phenotypes)
psum <- all_sum[all_sum$Phenotype %in% p,] #162x20


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
pgroups <- model_results_short[,c(3,10)]
pgroups <- pgroups[!duplicated(pgroups),] #162x2, removed duplicated rows

## add pgroups to psum table
pgroups$x <- gsub("_invr", "", pgroups$x)
colnames(pgroups)[1] <- "Phenotype"
psum2 <- left_join(psum, pgroups, by="Phenotype")

## remove "model_exclusion" and "model_use" columns from psum2
psum3 <- psum2[,-c(17:19)]

## resort columns in psum3
psum4 <- psum3[,c(18,1:17)] #162x18

## correct entries for stool diary
psum4$Phenotype <- gsub("mother_dis_diary_freq_mean", "mother_dis_diary_stool_freq_mean", psum4$Phenotype)
psum4$Phenotype <- gsub("mother_dis_diary_not_well_mean", "mother_dis_diary_uneasy_abdominal_feeling_mean", psum4$Phenotype)

## save polished phenotype summary statistics table
write.table(psum4, file="phenotype_summary_statistics/240229_milk_HMOs_real_included_maternal_and_infant_phenotypes_summary_statistics_n162.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 18.2 MILK HMOs - MODEL RESULTS TABLE (TABLE S8)  ===== ###

resHMO_tmp <- model_results_short

# add 'phenotype_type' column
resHMO_tmp$phenotype_type <- resHMO_tmp$dataset
resHMO_tmp$phenotype_type <- gsub("infant_real_HMOs_", "", resHMO_tmp$phenotype_type)
resHMO_tmp$phenotype_type <- gsub("mother_real_HMOs_", "", resHMO_tmp$phenotype_type)

# fix entries in 'model' ('statistic') column
resHMO_tmp$statistic <- as.factor(as.character(resHMO_tmp$statistic))
levels(resHMO_tmp$statistic) <- c("mmrm_compound_symmetry_correction_for_milkgroup_time_ID",
                                  "lm_correction_for_milkgroup",
                                  "lmer_correction_for_milkgroup_time_ID",
                                  "mmrm_unstructured_correction_for_milkgroup_time_ID")

# remove extensions for transformations from x and y columns (only mention i table description later)
resHMO_tmp$x <- gsub("_invr", "", resHMO_tmp$x)
resHMO_tmp$y <- gsub("_invr", "", resHMO_tmp$y)

# add columns from all_sum (data_type, model_time_point)
all_sum_tmp <- all_sum[all_sum$Phenotype %in% resHMO_tmp$x, c(1,2,19)] #162x3
colnames(all_sum_tmp)[1] <- "x"
resHMO_tmp2 <- left_join(resHMO_tmp, all_sum_tmp, by="x")

# # ensure NAs show for numeric/integer phenotypes in investigated category column
# resHMO_tmp2$levels <- as.factor(as.character(resHMO_tmp2$levels))
# levels(resHMO_tmp2$levels)[1] <- NA

# add column showing the reference category for categorical phenotypes
resHMO_tmp2$reference_category <- resHMO_tmp2$levels
levels(resHMO_tmp2$reference_category) <- c("O","O","Philips",
                                            "Spring","O","one_breast",
                                            "breastfeeding","vaginal_birth","big_house",
                                            "big_house","home","female",
                                            "Philips","breastfeeding","big_house",
                                            "Philips","big_house","Spring",
                                            "single_pregnancy","Spring","no")

# remove columns 'dataset' and resort columns in resHMO_tmp2
resHMO_tmp3 <- resHMO_tmp2[,c(10,3,14,16,6,11,4,5,7,8,12,15,13,2,9)]

# fix colnames in resHMO_tmp3
colnames(resHMO_tmp3) <- c("Phenotype_group", "Phenotype", "Data_type", "Reference_category", "Tested_category", "Phenotype_label_in_plots",
                           "HMO",
                           "n", "Estimate", "p", "FDR",
                           "Outcome_time_points_linked_to_phenotype", "Phenotype_type", "Model", "Note")

# fix entries for HMO diversity measures
resHMO_tmp3$HMO <- gsub("mother_milk_HMO_", "", resHMO_tmp3$HMO)
resHMO_tmp3$HMO <- gsub("_ugml", "", resHMO_tmp3$HMO)

resHMO <- resHMO_tmp3

## save polished results table
write.table(resHMO, file="association_model_results/240229_polished_model_results_milk_HMOs_real_by_phenotypes_n162.txt", row.names=F, col.names=T, sep="\t", quote=F)







