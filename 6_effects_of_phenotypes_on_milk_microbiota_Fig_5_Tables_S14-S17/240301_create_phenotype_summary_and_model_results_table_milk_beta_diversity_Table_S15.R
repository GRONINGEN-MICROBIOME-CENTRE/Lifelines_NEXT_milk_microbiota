########################################################################################################################################
### CREATE PHENOTYPES SUMMARY STATISTICS + MODEL RESULTS TABLE FOR MILK MICROBIOTA DATA (BETA DIVERSITY) FOR MILK COMPOSITION PAPER  ###
########################################################################################################################################

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

# 6. INVERSE-RANK TRANSFORM NUMERIC PHENOTYPES AND MILK ALPHA DIVERSITY DATA
#    6.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
#    6.2 INVERSE-RANK TRANSFORM DATA IN DATA FRAME WITH MATERNAL PHENOTYPES
#    6.3 INVERSE-RANK TRANSFORM DATA IN DATA FRAME WITH INFANT PHENOTYPES

# 7. CALCULATE MILK RELATIVE ABUNDANCES
#    7.1 CALCULATE MILK RELATIVE ABUNDANCES FOR DATA FRAME WITH MATERNAL PHENOTYPES
#    7.2 CALCULATE MILK RELATIVE ABUNDANCES FOR DATA FRAME WITH INFANT PHENOTYPES

# 8. CLR-TRANSFORM MILK RELATIVE ABUNDANCES
#    8.1 FUNCTION FOR CLR-TRANSFORMATION
#    8.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES IN DATA FRAME WITH MATERNAL PHENOTYPES
#    8.3 CLR-TRANSFORM MILK RELATIVE ABUNDANCES IN DATA FRAME WITH INFANT PHENOTYPES

# 9. ASSOCIATE PHENOTYPES WITH MILK BETA DIVERSITY
#    9.1 FUNCTION FOR ADONIS WITH CORRECTION FOR DNA ISOLATION BATCH AND READS
#    9.2 ASSOCIATE MATERNAL PHENOTYPES WITH MILK BETA DIVERSITY WITH CORRECTION FOR DNA ISOLATION BATCH AND READS
#    9.3 ASSOCIATE INFANT PHENOTYPES WITH MILK BETA DIVERSITY WITH CORRECTION FOR DNA ISOLATION BATCH AND READS
#    9.4 FUNCTION FOR ADONIS WITH CORRECTION FOR READS
#    9.5 ASSOCIATE MATERNAL PHENOTYPES WITH MILK BETA DIVERSITY WITH CORRECTION FOR READS
#    9.6 ASSOCIATE INFANT PHENOTYPES WITH MILK BETA DIVERSITY WITH CORRECTION FOR READS

# 10. COMBINE RESULTS AND CORRECT FOR MULTIPLE TESTING
#    10.1 EXCLUDE PHENOTYPES
#    10.2 ADD TIME POINTS TO RESULTS FROM ADONIS
#    10.3 ADD INFO ON CORRECTION FACTORS TO RESULTS FROM ADONIS
#    10.4 COMBINE ADONIS RESULTS PER TIME POINT AND CORRECT FOR MULTIPLE TESTING
#    10.5 POLISH RESULTS TABLES AND SAVE THEM (TABLES S15A-D)
#    10.6 CHECK RESULTS
#    [10.7 COMBINE ALL RESULTS FROM ADONIS AND CORRECT FOR MULTIPLE TESTING]


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

## add info on which phenotypes to exclude as x from association models
overview_mmilk$model_exclusion <- c(rep("exclude",length(1:2)),    #time_point and time_point_numeric (exclude as it will be corrected for)
                                   rep("include",length(3:4)),     #type_pregnancy, milk collection
                                   rep("exclude",length(5)),       #milk collection month
                                   rep("include",length(6:8)),     #milk collection
                                   rep("exclude",length(9:15)),    #milk collection notes, maternal FUT2 and FUT3 genotypes, maternal blood group genotype
                                   rep("include",length(16:24)),   #maternal anthropometrics
                                   rep("exclude",length(25)),      #mother_weight_gain_preg_kg (small n)
                                   rep("include",length(26:29)),   #maternal anthropometrics
                                   rep("exclude",length(30:116)),  #food pref, meat and fish yes/no
                                   rep("include",length(117:144)), #food groups
                                   rep("exclude",length(145:307)), #food items
                                   rep("include",length(308:354)), #nutrients, aMED, stress, education, exp living and pets
                                   rep("exclude",length(355:356)), #cats and dogs
                                   rep("include",length(357:381)), #smoking, mat health, preg/birth comp
                                   rep("exclude",length(382)),     #mother_birth_gestational_age_categories
                                   rep("include",length(383:385)), #preg/birth comp, subclinical mastitis                               
                                   rep("exclude",length(386)),     #Na/K ratio
                                   rep("include",length(387:415)), #maternal GIT diary, ROME, preg and birth med use, maternal postpartum med use, breastpumping, milk group
                                   rep("exclude",length(416:417)), #Le/Se status
                                   rep("include",length(418:445))) #milk HMOs

## exclude all rows from phenotypes that should not be included in the models
overview_mmilk2 <- overview_mmilk[overview_mmilk$model_exclusion=="include",]
rownames(overview_mmilk2) <- 1:nrow(overview_mmilk2)

# rownames(mmilk) <- mmilk$seq_16S_sample_ID
# test <- mmilk[,colnames(mmilk) %in% overview_mmilk2$column_name]
# test$seq_16S_sample_ID <- rownames(test)
# test2 <- merge(test, mmilk[,c(8,459,466)], by="seq_16S_sample_ID")
# for (i in 1:ncol(test2)){
#   print(paste0(i, "__", colnames(test2)[i]))
#   print(table(test2[!is.na(test2[,i]), "DNA_isolation_batch"]))
# }

# only in batch 2:
# 150__mother_breastpump_brand
# 151__mother_breastpump_times_per_day

## add info in data_use column to indicate whether a used phenotype is present (and should thus be corrected for) both DNA isolation batches
overview_mmilk2$model_use <- c(rep("adonis_with_correction_for_DNAbatch_reads",length(1:147)),
                              rep("adonis_with_correction_for_reads",length(148:149)),
                              rep("adonis_with_correction_for_DNAbatch_reads",length(150:178)))

# for (i in 1:ncol(test2)){
#   print(paste0(i, "__", colnames(test2)[i]))
#   print(table(test2[!is.na(test2[,i]), "time_point"]))
# }

## add info on which microbiota time points phenotypes were/will be linked to
overview_mmilk2$model_time_point <- c(rep("1_to_6_month(s)",length(1:120)),
                                     rep("1_to_3_month(s)",length(121)),        #subclinical mastitis
                                     rep("1_to_6_month(s)",length(122:133)),
                                     rep("1_month",length(134:139)),            #preg and birth med use
                                     rep("1_and_3_month(s)",length(140:147)),   #maternal postpartum med use
                                     rep("1_to_3_month(s)",length(148:149)),    #breastpumping (no data for M6 was included)
                                     rep("1_to_6_month(s)",length(150:178))) 

## separate phenotypes by data type
factor_mmilk <- overview_mmilk2[overview_mmilk2$data_type=="factor",] #54 categorical phenotypes
numeric_mmilk <- overview_mmilk2[overview_mmilk2$data_type=="numeric" | overview_mmilk2$data_type=="integer",] #124 numeric phenotypes


### ===== 5.3 INFANT PHENOTYPES ===== ###

## Infant phenotypes: create overview
overview_imilk <- create.overview.table(imilk[,c(19,20,22,24:28)])

## add info on which phenotypes to exclude as x from association models
overview_imilk$model_exclusion <- c(rep("include",length(1:8)))

# rownames(imilk) <- imilk$seq_16S_sample_ID
# itest <- imilk[,colnames(imilk) %in% overview_imilk$column_name]
# itest$seq_16S_sample_ID <- rownames(itest)
# itest2 <- merge(itest, imilk[,c(8,61,68)], by="seq_16S_sample_ID")
# for (i in 1:ncol(itest2)){
#   print(paste0(i, "__", colnames(itest2)[i]))
#   print(table(itest2[!is.na(itest2[,i]), "DNA_isolation_batch"]))
# }

# only in batch 2:
# 9__infant_ffq_solid_food_introduced

## add info in data_use column to indicate whether a used phenotype is present (and should thus be corrected for) both DNA isolation batches
overview_imilk$model_use <- c(rep("adonis_with_correction_for_DNAbatch_reads",length(1:7)),
                             rep("adonis_with_correction_for_reads",length(8))) #infant food introduction

# for (i in 1:ncol(itest2)){
#   print(paste0(i, "__", colnames(itest2)[i]))
#   print(table(itest2[!is.na(itest2[,i]), "time_point"]))
# }

## add info on which microbiota time points phenotypes were/will be linked to
overview_imilk$model_time_point <- c(rep("1_to_6_month(s)",length(1:4)), 
                                    rep("1_month",length(5:6)),         #infant feeding birth
                                    rep("1_to_6_month(s)",length(7)),   #infant feeding mode
                                    rep("6_months",length(8)))          #infant food introduction


## separate phenotypes by data type
factor_imilk <- overview_imilk[overview_imilk$data_type=="factor",] #6 categorical phenotypes
numeric_imilk <- overview_imilk[overview_imilk$data_type=="numeric" | overview_imilk$data_type=="integer",] #2 numeric phenotypes


##### =========================== 6. INVERSE-RANK TRANSFORM NUMERIC PHENOTYPES AND MILK ALPHA DIVERSITY DATA =========================== #####

### ===== 6.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 6.2 INVERSE-RANK TRANSFORM DATA IN DATA FRAME WITH MATERNAL PHENOTYPES ===== ###

## Maternal phenotypes
mmilk_invr <- mmilk
mmilk_invr[,numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"]] <- as.data.frame(apply(mmilk_invr[,numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"]], 2, invrank))

colnames(mmilk_invr)[which(colnames(mmilk_invr) %in% numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"])] <- paste0(colnames(mmilk_invr)[which(colnames(mmilk_invr) %in% numeric_mmilk[numeric_mmilk$model_exclusion=="include","column_name"])], "_invr")

# also adapt the column names in the overview_mmilk2 table
overview_mmilk2[overview_mmilk2$data_type=="numeric" | overview_mmilk2$data_type=="integer", "column_name"] <- paste0(overview_mmilk2[overview_mmilk2$data_type=="numeric" | overview_mmilk2$data_type=="integer", "column_name"], "_invr")


### ===== 6.3 INVERSE-RANK TRANSFORM DATA IN DATA FRAME WITH INFANT PHENOTYPES ===== ###

## Infant phenotypes
imilk_invr <- imilk
imilk_invr[,numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"]] <- as.data.frame(apply(imilk_invr[,numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"]], 2, invrank))

colnames(imilk_invr)[which(colnames(imilk_invr) %in% numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"])] <- paste0(colnames(imilk_invr)[which(colnames(imilk_invr) %in% numeric_imilk[numeric_imilk$model_exclusion=="include","column_name"])], "_invr")

# also adapt the column names in the overview_mmilk2 table
overview_imilk[overview_imilk$data_type=="numeric" | overview_imilk$data_type=="integer", "column_name"] <- paste0(overview_imilk[overview_imilk$data_type=="numeric" | overview_imilk$data_type=="integer", "column_name"], "_invr")


##### =========================== 7. CALCULATE MILK RELATIVE ABUNDANCES =========================== #####

### ===== 7.1 CALCULATE MILK RELATIVE ABUNDANCES FOR DATA FRAME WITH MATERNAL PHENOTYPES ===== ###

## save only columns with absolute bacterial abundances
rownames(mmilk_invr) <- mmilk_invr$seq_16S_sample_ID
mmilkbact <- mmilk_invr[,486:ncol(mmilk_invr)] #725x624

## calculate relative bacterial abundances
mmilkbact_relab <- (mmilkbact/rowSums(mmilkbact))
table(rowSums(mmilkbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
mmilkbact_relab_tmp <- mmilkbact_relab[,-624]

## exclude columns for bacteria that were not identified in any sample
mmilkbact_relab_tmp2 <- mmilkbact_relab_tmp[,colSums(mmilkbact_relab_tmp)>0]
dim(mmilkbact_relab_tmp2)

## merge back with metadata
mmilk_invr_RA <- merge(mmilk_invr[,1:485], mmilkbact_relab_tmp2, by="row.names")
rownames(mmilk_invr_RA) <- mmilk_invr_RA$Row.names
mmilk_invr_RA <- mmilk_invr_RA[,-1]


### ===== 7.2 CALCULATE MILK RELATIVE ABUNDANCES FOR DATA FRAME WITH INFANT PHENOTYPES ===== ###

## save only columns with absolute bacterial abundances
rownames(imilk_invr) <- imilk_invr$seq_16S_sample_ID
imilkbact <- imilk_invr[,88:ncol(imilk_invr)] #737x624

## calculate relative bacterial abundances
imilkbact_relab <- (imilkbact/rowSums(imilkbact))
table(rowSums(imilkbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
imilkbact_relab_tmp <- imilkbact_relab[,-624]

## exclude columns for bacteria that were not identified in any sample
imilkbact_relab_tmp2 <- imilkbact_relab_tmp[,colSums(imilkbact_relab_tmp)>0]
dim(imilkbact_relab_tmp2)

## merge back with metadata
imilk_invr_RA <- merge(imilk_invr[,1:87], imilkbact_relab_tmp2, by="row.names")
rownames(imilk_invr_RA) <- imilk_invr_RA$Row.names
imilk_invr_RA <- imilk_invr_RA[,-1]


##### =========================== 8. CLR-TRANSFORM MILK RELATIVE ABUNDANCES =========================== #####

### ===== 8.1 FUNCTION FOR CLR-TRANSFORMATION ===== ###

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


### ===== 8.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES IN DATA FRAME WITH MATERNAL PHENOTYPES ===== ###

mmilk_invr_RAclr <- mmilk_invr_RA
mmilk_invr_RAclr[,486:ncol(mmilk_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(mmilk_invr_RAclr[,486:ncol(mmilk_invr_RAclr)], mmilk_invr_RAclr[,486:ncol(mmilk_invr_RAclr)]))

colnames(mmilk_invr_RAclr)[486:ncol(mmilk_invr_RAclr)] <- c(paste0(colnames(mmilk_invr_RAclr)[486:ncol(mmilk_invr_RAclr)], "_clr"))


### ===== 8.3 CLR-TRANSFORM MILK RELATIVE ABUNDANCES IN DATA FRAME WITH INFANT PHENOTYPES ===== ###

imilk_invr_RAclr <- imilk_invr_RA
imilk_invr_RAclr[,88:ncol(imilk_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(imilk_invr_RAclr[,88:ncol(imilk_invr_RAclr)], imilk_invr_RAclr[,88:ncol(imilk_invr_RAclr)]))

colnames(imilk_invr_RAclr)[88:ncol(imilk_invr_RAclr)] <- c(paste0(colnames(imilk_invr_RAclr)[88:ncol(imilk_invr_RAclr)], "_clr"))


##### =========================== 9. ASSOCIATE PHENOTYPES WITH MILK BETA DIVERSITY =========================== #####

library(vegan)


### ===== 9.1 FUNCTION FOR ADONIS WITH CORRECTION FOR DNA ISOLATION BATCH AND READS ===== ###

## function to run Adonis with correction for DNA isolation batch and clean read counts
run.adonis2.cor.DNAbatch.reads <- function(mydata, my_nontransformed_data, columns, nontransformed_columns){
  column <- c()
  R2 <- c()
  p <- c()
  n <- c()
  all_levels <- c()
  n_per_level <- c()
  my_mean <- c()
  my_sd <- c()
  
  note <- c()
  
  for (i in columns){
    print(paste0("Starting analysis for ", colnames(mydata)[i]))
    
    d <- mydata[!is.na(mydata[,i]),]

    ## AITCHISON DISTANCE MATRIX + METADATA
    print("Creating Aitchison distance matrix")
    ait <- vegdist(d[,grep("g__", colnames(d))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
    aitm <- cmdscale(ait, k=3) #create distance matrix for first 3 components
    aitd <- data.frame(aitm)
    
    print("Creating data frame with Aitchison distances and metadata")
    new_df <- merge(d, aitd, by="row.names")
    new_df <- new_df[,-1]
    rownames(new_df) <- new_df$seq_16S_sample_ID
    
    column <- c(column, colnames(new_df)[i])
    n <- c(n, nrow(new_df))
    
    
    ## ADONIS
    print("Running adonis")
    a <- table(new_df[,i])/nrow(new_df)
    
    if (any(a>0.90)){ #set to NA if >90% of the data is the same (this includes if there is only 1 factor because that is 100% the same)
      R2 <- c(R2, NA)
      p <- c(p, NA)
      note <- c(note, paste0("Adonis_could_not_be_run_as_there_was_too_little_variation_in_the_data"))
    }else{ ##only run adonis if there is sufficient variation in the answers (at least 3% need to have different answer), this also means that >=2 levels are present for the variable of interest
      a1 <- adonis2(ait ~ new_df$DNA_isolation_batch + new_df$seq_16S_n_reads_clean + new_df[,i], by="margin")
      R2 <- c(R2, a1$R2[3]) #we pick the 3rd value as the first 2 are for the factors we correct for
      p <- c(p, a1$'Pr(>F)'[3]) #we pick the 3rd value as the first 2 are for the factors we correct for
      note <- c(note, NA)
    }
    }
  
  
  for (j in nontransformed_columns){
    d_nontf <- my_nontransformed_data[!is.na(my_nontransformed_data[,j]),]
    
    if(class(d_nontf[,j]) == "factor" | class(d_nontf[,j]) == "character"){
      print("Retrieving summary statistics for categorical phenotype")
      all_levels <- c(all_levels, paste(collapse=";",levels(as.factor(as.character(d_nontf[,j])))))
      n_per_level <- c(n_per_level, paste(collapse=";",table(as.factor(as.character(d_nontf[,j])))))
      my_mean <- c(my_mean, NA)
      my_sd <- c(my_sd, NA)
    }else{
      print("Retrieving summary statistics for numeric/integer phenotype")
      all_levels <- c(all_levels, NA)
      n_per_level <- c(n_per_level, NA)
      my_mean <- c(my_mean, mean(d_nontf[,j]))
      my_sd <- c(my_sd, sd(d_nontf[,j]))
    }
  }

  
  print("Preparing results data frame")
  b <- data.frame(column=column, n=n, categories=all_levels, n_per_category=n_per_level, mean=my_mean, SD=my_sd, R2=R2, pvalue=p, note=note)
  return(b)
}


### ===== 9.2 ASSOCIATE MATERNAL PHENOTYPES WITH MILK BETA DIVERSITY WITH CORRECTION FOR DNA ISOLATION BATCH AND READS ===== ###

## M1 - phenotypes that are present in both DNA isolation batches
m_phenos_M1 <- overview_mmilk2[overview_mmilk2$model_use=="adonis_with_correction_for_DNAbatch_reads" & overview_mmilk2$model_exclusion=="include", "column_name"]
m_phenos_M1_nontransformed <- gsub("_invr", "", m_phenos_M1)

m_adonis2_results_M1_1 <- run.adonis2.cor.DNAbatch.reads(mydata = mmilk_invr_RAclr[mmilk_invr_RAclr$time_point=="1_month",],
                                                         my_nontransformed_data = mmilk[mmilk$time_point=="1_month",],
                                                         columns = which(colnames(mmilk_invr_RAclr) %in% colnames(mmilk_invr_RAclr)[colnames(mmilk_invr_RAclr) %in% m_phenos_M1]),
                                                         nontransformed_columns = which(colnames(mmilk) %in% colnames(mmilk)[colnames(mmilk) %in% m_phenos_M1_nontransformed]))


## M2 - phenotypes that are present in both DNA isolation batches
m_phenos_M2 <- overview_mmilk2[overview_mmilk2$model_use=="adonis_with_correction_for_DNAbatch_reads" & overview_mmilk2$model_exclusion=="include" & 
                                 (overview_mmilk2$model_time_point=="1_to_6_month(s)" | overview_mmilk2$model_time_point=="1_to_3_month(s)"), #include model_time_point for selection!
                               "column_name"] 
m_phenos_M2_nontransformed <- gsub("_invr", "", m_phenos_M2)

m_adonis2_results_M2_1 <- run.adonis2.cor.DNAbatch.reads(mydata = mmilk_invr_RAclr[mmilk_invr_RAclr$time_point=="2_months",],
                                                         my_nontransformed_data = mmilk[mmilk$time_point=="2_months",],
                                                         columns = which(colnames(mmilk_invr_RAclr) %in% colnames(mmilk_invr_RAclr)[colnames(mmilk_invr_RAclr) %in% m_phenos_M2]),
                                                         nontransformed_columns = which(colnames(mmilk) %in% colnames(mmilk)[colnames(mmilk) %in% m_phenos_M2_nontransformed]))


## M3 - phenotypes that are present in both DNA isolation batches
m_phenos_M3 <- overview_mmilk2[overview_mmilk2$model_use=="adonis_with_correction_for_DNAbatch_reads" & overview_mmilk2$model_exclusion=="include" & 
                                 (overview_mmilk2$model_time_point=="1_to_6_month(s)" | overview_mmilk2$model_time_point=="1_to_3_month(s)" | overview_mmilk2$model_time_point=="1_and_3_month(s)"), #include model_time_point for selection!
                               "column_name"]
m_phenos_M3_nontransformed <- gsub("_invr", "", m_phenos_M3)

m_adonis2_results_M3_1 <- run.adonis2.cor.DNAbatch.reads(mydata = mmilk_invr_RAclr[mmilk_invr_RAclr$time_point=="3_months",],
                                                         my_nontransformed_data = mmilk[mmilk$time_point=="3_months",],
                                                         columns = which(colnames(mmilk_invr_RAclr) %in% colnames(mmilk_invr_RAclr)[colnames(mmilk_invr_RAclr) %in% m_phenos_M3]),
                                                         nontransformed_columns = which(colnames(mmilk) %in% colnames(mmilk)[colnames(mmilk) %in% m_phenos_M3_nontransformed]))


## M6 - phenotypes that are present in both DNA isolation batches
m_phenos_M6 <- overview_mmilk2[overview_mmilk2$model_use=="adonis_with_correction_for_DNAbatch_reads" & overview_mmilk2$model_exclusion=="include" & 
                                 overview_mmilk2$model_time_point=="1_to_6_month(s)", #include model_time_point for selection!
                               "column_name"]
m_phenos_M6_nontransformed <- gsub("_invr", "", m_phenos_M6)

m_adonis2_results_M6_1 <- run.adonis2.cor.DNAbatch.reads(mydata = mmilk_invr_RAclr[mmilk_invr_RAclr$time_point=="6_months",],
                                                         my_nontransformed_data = mmilk[mmilk$time_point=="6_months",],
                                                         columns = which(colnames(mmilk_invr_RAclr) %in% colnames(mmilk_invr_RAclr)[colnames(mmilk_invr_RAclr) %in% m_phenos_M6]),
                                                         nontransformed_columns = which(colnames(mmilk) %in% colnames(mmilk)[colnames(mmilk) %in% m_phenos_M6_nontransformed]))


### ===== 9.3 ASSOCIATE INFANT PHENOTYPES WITH MILK BETA DIVERSITY WITH CORRECTION FOR DNA ISOLATION BATCH AND READS ===== ###

## M1 - phenotypes that are present in both DNA isolation batches
i_phenos_M1 <- overview_imilk[overview_imilk$model_use=="adonis_with_correction_for_DNAbatch_reads" & overview_imilk$model_exclusion=="include"& 
                                (overview_imilk$model_time_point=="1_to_6_month(s)" | overview_imilk$model_time_point=="1_month"), #include model_time_point for selection!
                              "column_name"]
i_phenos_M1_nontransformed <- gsub("_invr", "", i_phenos_M1)

i_adonis2_results_M1_1 <- run.adonis2.cor.DNAbatch.reads(mydata = imilk_invr_RAclr[imilk_invr_RAclr$time_point=="1_month",],
                                                         my_nontransformed_data = imilk[imilk$time_point=="1_month",],
                                                         columns = which(colnames(imilk_invr_RAclr) %in% colnames(imilk_invr_RAclr)[colnames(imilk_invr_RAclr) %in% i_phenos_M1]),
                                                         nontransformed_columns = which(colnames(imilk) %in% colnames(imilk)[colnames(imilk) %in% i_phenos_M1_nontransformed]))


## M2 - phenotypes that are present in both DNA isolation batches
i_phenos_M2 <- overview_imilk[overview_imilk$model_use=="adonis_with_correction_for_DNAbatch_reads" & overview_imilk$model_exclusion=="include" & 
                                 overview_imilk$model_time_point=="1_to_6_month(s)", #include model_time_point for selection!
                               "column_name"] 
i_phenos_M2_nontransformed <- gsub("_invr", "", i_phenos_M2)

i_adonis2_results_M2_1 <- run.adonis2.cor.DNAbatch.reads(mydata = imilk_invr_RAclr[imilk_invr_RAclr$time_point=="2_months",],
                                                         my_nontransformed_data = imilk[imilk$time_point=="2_months",],
                                                         columns = which(colnames(imilk_invr_RAclr) %in% colnames(imilk_invr_RAclr)[colnames(imilk_invr_RAclr) %in% i_phenos_M2]),
                                                         nontransformed_columns = which(colnames(imilk) %in% colnames(imilk)[colnames(imilk) %in% i_phenos_M2_nontransformed]))


## M3 - phenotypes that are present in both DNA isolation batches
i_phenos_M3 <- overview_imilk[overview_imilk$model_use=="adonis_with_correction_for_DNAbatch_reads" & overview_imilk$model_exclusion=="include" & 
                                 overview_imilk$model_time_point=="1_to_6_month(s)", #include model_time_point for selection!
                               "column_name"]
i_phenos_M3_nontransformed <- gsub("_invr", "", i_phenos_M3)

i_adonis2_results_M3_1 <- run.adonis2.cor.DNAbatch.reads(mydata = imilk_invr_RAclr[imilk_invr_RAclr$time_point=="3_months",],
                                                         my_nontransformed_data = imilk[imilk$time_point=="3_months",],
                                                         columns = which(colnames(imilk_invr_RAclr) %in% colnames(imilk_invr_RAclr)[colnames(imilk_invr_RAclr) %in% i_phenos_M3]),
                                                         nontransformed_columns = which(colnames(imilk) %in% colnames(imilk)[colnames(imilk) %in% i_phenos_M3_nontransformed]))


## M6 - phenotypes that are present in both DNA isolation batches
i_phenos_M6 <- overview_imilk[overview_imilk$model_use=="adonis_with_correction_for_DNAbatch_reads" & overview_imilk$model_exclusion=="include" & 
                                 (overview_imilk$model_time_point=="1_to_6_month(s)" | overview_imilk$model_time_point=="6_months"), #include model_time_point for selection!
                               "column_name"]
i_phenos_M6_nontransformed <- gsub("_invr", "", i_phenos_M6)

i_adonis2_results_M6_1 <- run.adonis2.cor.DNAbatch.reads(mydata = imilk_invr_RAclr[imilk_invr_RAclr$time_point=="6_months",],
                                                         my_nontransformed_data = imilk[imilk$time_point=="6_months",],
                                                         columns = which(colnames(imilk_invr_RAclr) %in% colnames(imilk_invr_RAclr)[colnames(imilk_invr_RAclr) %in% i_phenos_M6]),
                                                         nontransformed_columns = which(colnames(imilk) %in% colnames(imilk)[colnames(imilk) %in% i_phenos_M6_nontransformed]))


### ===== 9.4 FUNCTION FOR ADONIS WITH CORRECTION FOR READS ===== ###

## function to run Adonis with correction for DNA isolation batch and clean read counte
run.adonis2.cor.reads <- function(mydata, my_nontransformed_data, columns, nontransformed_columns){
  column <- c()
  R2 <- c()
  p <- c()
  n <- c()
  all_levels <- c()
  n_per_level <- c()
  my_mean <- c()
  my_sd <- c()
  
  note <- c()
  
  for (i in columns){
    print(paste0("Starting analysis for ", colnames(mydata)[i]))
    
    d <- mydata[!is.na(mydata[,i]),]
    print(paste0("The data frame for this phenotype has these dimensions: ", dim(d)))
    
    ## AITCHISON DISTANCE MATRIX + METADATA
    print("Creating Aitchison distance matrix")
    ait <- vegdist(d[,grep("g__", colnames(d))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
    aitm <- cmdscale(ait, k=3) #create distance matrix for first 3 components
    aitd <- data.frame(aitm)
    
    print("Creating data frame with Aitchison distances and metadata")
    new_df <- merge(d, aitd, by="row.names")
    new_df <- new_df[,-1]
    rownames(new_df) <- new_df$seq_16S_sample_ID
    
    column <- c(column, colnames(new_df)[i])
    n <- c(n, nrow(new_df))
    
    ## ADONIS
    print("Running adonis")
    a <- table(new_df[,i])/nrow(new_df)
    
    if (any(a>0.90)){ #set to NA if >90% of the data is the same (this includes if there is only 1 factor because that is 100% the same)
      R2 <- c(R2, NA)
      p <- c(p, NA)
      note <- c(note, paste0("Adonis_could_not_be_run_as_there_was_too_little_variation_in_the_data"))
    }else{ ##only run adonis if there is sufficient variation in the answers (at least 3% need to have different answer), this also means that >=2 levels are present for the variable of interest
      a1 <- adonis2(ait ~ new_df$seq_16S_n_reads_clean + new_df[,i], by="margin")
      print(a1)
      R2 <- c(R2, a1$R2[2]) #we pick the 3rd value as the first 2 are for the factors we correct for
      p <- c(p, a1$'Pr(>F)'[2]) #we pick the 3rd value as the first 2 are for the factors we correct for
      note <- c(note, NA)
    }
  }
  
  
  for (j in nontransformed_columns){
    d_nontf <- my_nontransformed_data[!is.na(my_nontransformed_data[,j]),]
    
    if(class(d_nontf[,j]) == "factor" | class(d_nontf[,j]) == "character"){
      print("Retrieving summary statistics for categorical phenotype")
      all_levels <- c(all_levels, paste(collapse=";",levels(as.factor(as.character(d_nontf[,j])))))
      n_per_level <- c(n_per_level, paste(collapse=";",table(as.factor(as.character(d_nontf[,j])))))
      my_mean <- c(my_mean, NA)
      my_sd <- c(my_sd, NA)
    }else{
      print("Retrieving summary statistics for numeric/integer phenotype")
      all_levels <- c(all_levels, NA)
      n_per_level <- c(n_per_level, NA)
      my_mean <- c(my_mean, mean(d_nontf[,j]))
      my_sd <- c(my_sd, sd(d_nontf[,j]))
    }
  }
  
  
  print("Preparing results data frame")
  b <- data.frame(column=column, n=n, categories=all_levels, n_per_category=n_per_level, mean=my_mean, SD=my_sd, R2=R2, pvalue=p, note=note)
  return(b)
}


### ===== 9.5 ASSOCIATE MATERNAL PHENOTYPES WITH MILK BETA DIVERSITY WITH CORRECTION FOR READS ===== ###

## M1 - phenotypes that are present in only 1 DNA isolation batch
m_phenos_correads_M1 <- overview_mmilk2[overview_mmilk2$model_use=="adonis_with_correction_for_reads" & overview_mmilk2$model_exclusion=="include", "column_name"]
m_phenos_correads_M1_nontransformed <- gsub("_invr", "", m_phenos_correads_M1)

m_adonis2_results_M1_2 <- run.adonis2.cor.reads(mydata = mmilk_invr_RAclr[mmilk_invr_RAclr$time_point=="1_month",],
                                                my_nontransformed_data = mmilk[mmilk$time_point=="1_month",],
                                                columns = which(colnames(mmilk_invr_RAclr) %in% colnames(mmilk_invr_RAclr)[colnames(mmilk_invr_RAclr) %in% m_phenos_correads_M1]),
                                                nontransformed_columns = which(colnames(mmilk) %in% colnames(mmilk)[colnames(mmilk) %in% m_phenos_correads_M1_nontransformed]))


# ## M2 - phenotypes that are present in only 1 DNA isolation batch
# m_phenos_correads_M2 <- overview_mmilk2[overview_mmilk2$model_use=="adonis_with_correction_for_reads" & overview_mmilk2$model_exclusion=="include" & 
#                                  (overview_mmilk2$model_time_point=="1_to_6_month(s)" | overview_mmilk2$model_time_point=="1_to_3_month(s)"), #include model_time_point for selection!
#                                "column_name"] 
# m_phenos_correads_M2_nontransformed <- gsub("_invr", "", m_phenos_correads_M2)
# 
# m_adonis2_results_M2_2 <- run.adonis2.cor.reads(mydata = mmilk_invr_RAclr[mmilk_invr_RAclr$time_point=="2_months",],
#                                                 my_nontransformed_data = mmilk[mmilk$time_point=="2_months",],
#                                                 columns = which(colnames(mmilk_invr_RAclr) %in% colnames(mmilk_invr_RAclr)[colnames(mmilk_invr_RAclr) %in% m_phenos_correads_M2]),
#                                                 nontransformed_columns = which(colnames(mmilk) %in% colnames(mmilk)[colnames(mmilk) %in% m_phenos_correads_M2_nontransformed]))
## -> No distance matrix can be created for these phenotypes as only 2 individuals have data.


## M3 - phenotypes that are present in only 1 DNA isolation batch
m_phenos_correads_M3 <- overview_mmilk2[overview_mmilk2$model_use=="adonis_with_correction_for_reads" & overview_mmilk2$model_exclusion=="include" & 
                                 (overview_mmilk2$model_time_point=="1_to_6_month(s)" | overview_mmilk2$model_time_point=="1_to_3_month(s)" | overview_mmilk2$model_time_point=="1_and_3_month(s)"), #include model_time_point for selection!
                               "column_name"]
m_phenos_correads_M3_nontransformed <- gsub("_invr", "", m_phenos_correads_M3)

m_adonis2_results_M3_2 <- run.adonis2.cor.reads(mydata = mmilk_invr_RAclr[mmilk_invr_RAclr$time_point=="3_months",],
                                                my_nontransformed_data = mmilk[mmilk$time_point=="3_months",],
                                                columns = which(colnames(mmilk_invr_RAclr) %in% colnames(mmilk_invr_RAclr)[colnames(mmilk_invr_RAclr) %in% m_phenos_correads_M3]),
                                                nontransformed_columns = which(colnames(mmilk) %in% colnames(mmilk)[colnames(mmilk) %in% m_phenos_correads_M3_nontransformed]))


# ## M6 - phenotypes that are present in only 1 DNA isolation batch
# m_phenos_correads_M6 <- overview_mmilk2[overview_mmilk2$model_use=="adonis_with_correction_for_reads" & overview_mmilk2$model_exclusion=="include" & 
#                                  overview_mmilk2$model_time_point=="1_to_6_month(s)", #include model_time_point for selection!
#                                "column_name"]
# m_phenos_correads_M6_nontransformed <- gsub("_invr", "", m_phenos_correads_M6)
#
# m_adonis2_results_M6_2 <- run.adonis2.cor.reads(mydata = mmilk_invr_RAclr[mmilk_invr_RAclr$time_point=="6_months",],
#                                                 my_nontransformed_data = mmilk[mmilk$time_point=="6_months",],
#                                                 columns = which(colnames(mmilk_invr_RAclr) %in% colnames(mmilk_invr_RAclr)[colnames(mmilk_invr_RAclr) %in% m_phenos_correads_M6]),
#                                                 nontransformed_columns = which(colnames(mmilk) %in% colnames(mmilk)[colnames(mmilk) %in% m_phenos_correads_M6_nontransformed]))
## -> No phenotypes for M6.


### ===== 9.6 ASSOCIATE INFANT PHENOTYPES WITH MILK BETA DIVERSITY WITH CORRECTION FOR READS ===== ###

# ## M1 - phenotypes that are present in only 1 DNA isolation batch
# i_phenos_correads_M1 <- overview_imilk[overview_imilk$model_use=="adonis_with_correction_for_reads" & overview_imilk$model_exclusion=="include"& 
#                                 (overview_imilk$model_time_point=="1_to_6_month(s)" | overview_imilk$model_time_point=="1_month"), #include model_time_point for selection!
#                               "column_name"]
# 
# i_adonis2_results_M1_2 <- run.adonis2.cor.reads(mydata = imilk_invr_RAclr[imilk_invr_RAclr$time_point=="1_month",],
#                                                          columns = which(colnames(imilk_invr_RAclr) %in% colnames(imilk_invr_RAclr)[colnames(imilk_invr_RAclr) %in% i_phenos_correads_M1]))
# ## -> No phenotypes for M1.


# ## M2 - phenotypes that are present in only 1 DNA isolation batch
# i_phenos_correads_M2 <- overview_imilk[overview_imilk$model_use=="adonis_with_correction_for_reads" & overview_imilk$model_exclusion=="include" & 
#                                 overview_imilk$model_time_point=="1_to_6_month(s)", #include model_time_point for selection!
#                               "column_name"] 
# 
# i_adonis2_results_M2_2 <- run.adonis2.cor.reads(mydata = imilk_invr_RAclr[imilk_invr_RAclr$time_point=="2_months",],
#                                                          columns = which(colnames(imilk_invr_RAclr) %in% colnames(imilk_invr_RAclr)[colnames(imilk_invr_RAclr) %in% i_phenos_correads_M2]))
# # ## -> No phenotypes for M2.


# ## M3 - phenotypes that are present in only 1 DNA isolation batch
# i_phenos_correads_M3 <- overview_imilk[overview_imilk$model_use=="adonis_with_correction_for_reads" & overview_imilk$model_exclusion=="include" & 
#                                 overview_imilk$model_time_point=="1_to_6_month(s)", #include model_time_point for selection!
#                               "column_name"]
# 
# i_adonis2_results_M3_2 <- run.adonis2.cor.reads(mydata = imilk_invr_RAclr[imilk_invr_RAclr$time_point=="3_months",],
#                                                          columns = which(colnames(imilk_invr_RAclr) %in% colnames(imilk_invr_RAclr)[colnames(imilk_invr_RAclr) %in% i_phenos_correads_M3]))
# # ## -> No phenotypes for M3.


## M6 - phenotypes that are present in only 1 DNA isolation batch
i_phenos_correads_M6 <- overview_imilk[overview_imilk$model_use=="adonis_with_correction_for_reads" & overview_imilk$model_exclusion=="include" & 
                                (overview_imilk$model_time_point=="1_to_6_month(s)" | overview_imilk$model_time_point=="6_months"), #include model_time_point for selection!
                              "column_name"]
i_phenos_correads_M6_nontransformed <- gsub("_invr", "", i_phenos_correads_M6)

i_adonis2_results_M6_2 <- run.adonis2.cor.reads(mydata = imilk_invr_RAclr[imilk_invr_RAclr$time_point=="6_months",],
                                                my_nontransformed_data = imilk[imilk$time_point=="3_months",],
                                                columns = which(colnames(imilk_invr_RAclr) %in% colnames(imilk_invr_RAclr)[colnames(imilk_invr_RAclr) %in% i_phenos_correads_M6]),
                                                nontransformed_columns = which(colnames(imilk) %in% colnames(imilk)[colnames(imilk) %in% i_phenos_correads_M6_nontransformed]))


##### =========================== 10. COMBINE RESULTS AND CORRECT FOR MULTIPLE TESTING =========================== #####

### ===== 10.1 EXCLUDE PHENOTYPES ===== ###

m_adonis2_results_M1_1 <- m_adonis2_results_M1_1[m_adonis2_results_M1_1$column!="mother_milk_collection_volume_pumped_invr",]
m_adonis2_results_M2_1 <- m_adonis2_results_M2_1[m_adonis2_results_M2_1$column!="mother_milk_collection_volume_pumped_invr",]
m_adonis2_results_M3_1 <- m_adonis2_results_M3_1[m_adonis2_results_M3_1$column!="mother_milk_collection_volume_pumped_invr",]
m_adonis2_results_M6_1 <- m_adonis2_results_M6_1[m_adonis2_results_M6_1$column!="mother_milk_collection_volume_pumped_invr",]


### ===== 10.2 ADD TIME POINTS TO RESULTS FROM ADONIS ===== ###

m_adonis2_results_M1_1$time_point <- "1_month"
m_adonis2_results_M2_1$time_point <- "2_months"
m_adonis2_results_M3_1$time_point <- "3_months"
m_adonis2_results_M6_1$time_point <- "6_months"

i_adonis2_results_M1_1$time_point <- "1_month"
i_adonis2_results_M2_1$time_point <- "2_months"
i_adonis2_results_M3_1$time_point <- "3_months"
i_adonis2_results_M6_1$time_point <- "6_months"

m_adonis2_results_M1_2$time_point <- "1_month"
m_adonis2_results_M3_2$time_point <- "3_months"

i_adonis2_results_M6_2$time_point <- "6_months"


### ===== 10.3 ADD INFO ON CORRECTION FACTORS TO RESULTS FROM ADONIS ===== ###

m_adonis2_results_M1_1$model <- "adonis_correction_for_DNA_isolation_batch_and_reads"
m_adonis2_results_M2_1$model <- "adonis_correction_for_DNA_isolation_batch_and_reads"
m_adonis2_results_M3_1$model <- "adonis_correction_for_DNA_isolation_batch_and_reads"
m_adonis2_results_M6_1$model <- "adonis_correction_for_DNA_isolation_batch_and_reads"

i_adonis2_results_M1_1$model <- "adonis_correction_for_DNA_isolation_batch_and_reads"
i_adonis2_results_M2_1$model <- "adonis_correction_for_DNA_isolation_batch_and_reads"
i_adonis2_results_M3_1$model <- "adonis_correction_for_DNA_isolation_batch_and_reads"
i_adonis2_results_M6_1$model <- "adonis_correction_for_DNA_isolation_batch_and_reads"

m_adonis2_results_M1_2$model <- "adonis_correction_for_reads"
m_adonis2_results_M3_2$model <- "adonis_correction_for_reads"

i_adonis2_results_M6_2$model <- "adonis_correction_for_reads"


### ===== 10.4 COMBINE ADONIS RESULTS PER TIME POINT AND CORRECT FOR MULTIPLE TESTING ===== ###

## combine results per time point
adonis_results_M1 <- as.data.frame(rbind(m_adonis2_results_M1_1, i_adonis2_results_M1_1, m_adonis2_results_M1_2)) #184x11
adonis_results_M2 <- as.data.frame(rbind(m_adonis2_results_M2_1, i_adonis2_results_M2_1))                         #166x11
adonis_results_M3 <- as.data.frame(rbind(m_adonis2_results_M3_1, i_adonis2_results_M3_1, m_adonis2_results_M3_2)) #176x11
adonis_results_M6 <- as.data.frame(rbind(m_adonis2_results_M6_1, i_adonis2_results_M6_1, i_adonis2_results_M6_2)) #166x11

## correct for multiple testing
adonis_results_M1$FDR <- p.adjust(adonis_results_M1$p, method="BH")
adonis_results_M2$FDR <- p.adjust(adonis_results_M2$p, method="BH")
adonis_results_M3$FDR <- p.adjust(adonis_results_M3$p, method="BH")
adonis_results_M6$FDR <- p.adjust(adonis_results_M6$p, method="BH")

## sort results by FDR and p
adonis_results_M1 <- adonis_results_M1[order(adonis_results_M1$FDR, adonis_results_M1$p),]
adonis_results_M2 <- adonis_results_M2[order(adonis_results_M2$FDR, adonis_results_M2$p),]
adonis_results_M3 <- adonis_results_M3[order(adonis_results_M3$FDR, adonis_results_M3$p),]
adonis_results_M6 <- adonis_results_M6[order(adonis_results_M6$FDR, adonis_results_M6$p),]


### ===== 10.5 POLISH RESULTS TABLES AND SAVE THEM (TABLES S15A-D) ===== ###

## resort columns
adonis_results_M1 <- adonis_results_M1[,c(10,1:8,12,11,9)]
adonis_results_M2 <- adonis_results_M2[,c(10,1:8,12,11,9)]
adonis_results_M3 <- adonis_results_M3[,c(10,1:8,12,11,9)]
adonis_results_M6 <- adonis_results_M6[,c(10,1:8,12,11,9)]

## save results per time point
write.table(adonis_results_M1, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/beta_diversity/240301_milk_beta_diversity_by_maternal_and_infant_phenotypes_adonis_results_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis_results_M2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/beta_diversity/240301_milk_beta_diversity_by_maternal_and_infant_phenotypes_adonis_results_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis_results_M3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/beta_diversity/240301_milk_beta_diversity_by_maternal_and_infant_phenotypes_adonis_results_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis_results_M6, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/beta_diversity/240301_milk_beta_diversity_by_maternal_and_infant_phenotypes_adonis_results_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 10.6 CHECK RESULTS ===== ###

## check results
#M1
nrow(adonis_results_M1[adonis_results_M1$FDR<0.05 & !is.na(adonis_results_M1$FDR),]) # 0 FDR significant associations
nrow(adonis_results_M1[adonis_results_M1$p<0.05 & !is.na(adonis_results_M1$p),])     #21 nominally significant associations

#M2
nrow(adonis_results_M2[adonis_results_M2$FDR<0.05 & !is.na(adonis_results_M2$FDR),]) # 0 FDR significant associations
nrow(adonis_results_M2[adonis_results_M2$p<0.05 & !is.na(adonis_results_M2$p),])     #13 nominally significant associations

#M3
nrow(adonis_results_M3[adonis_results_M3$FDR<0.05 & !is.na(adonis_results_M3$FDR),]) # 0 FDR significant associations
nrow(adonis_results_M3[adonis_results_M3$p<0.05 & !is.na(adonis_results_M3$p),])     # 6 nominally significant associations

#M6
nrow(adonis_results_M6[adonis_results_M6$FDR<0.05 & !is.na(adonis_results_M6$FDR),]) # 0 FDR significant associations
nrow(adonis_results_M6[adonis_results_M6$p<0.05 & !is.na(adonis_results_M6$p),])     # 7 nominally significant associations


## check if some phenotypes have nominally significant effects at M1 and M3 (the time points with the largest n)
nom_sign_M1_M3 <- intersect(adonis_results_M1[adonis_results_M1$p<0.05 & !is.na(adonis_results_M1$p),"column"], adonis_results_M3[adonis_results_M3$p<0.05 & !is.na(adonis_results_M3$p),"column"])
nom_sign_M1_M3
# "infant_ffq_feeding_mode"

## check if some phenotypes have nominally significant effects at M1, M3 and M6 (the time points with the largest n)
nom_sign_M1_M3_M6 <- intersect(nom_sign_M1_M3, adonis_results_M6[adonis_results_M6$p<0.05 & !is.na(adonis_results_M6$p),"column"])
nom_sign_M1_M3_M6
# none

## check if some phenotypes have nominally significant effects at M1 and M3 (the time points with the largest n)
nom_sign_M1_M2_M3 <- intersect(nom_sign_M1_M3, adonis_results_M2[adonis_results_M2$p<0.05 & !is.na(adonis_results_M2$p),"column"])
nom_sign_M1_M2_M3
# "infant_ffq_feeding_mode"

## check results for infant_ffq_feeding_mode at each time point
adonis_results_M1[adonis_results_M1$column=="infant_ffq_feeding_mode",]
adonis_results_M2[adonis_results_M2$column=="infant_ffq_feeding_mode",]
adonis_results_M3[adonis_results_M3$column=="infant_ffq_feeding_mode",]
adonis_results_M6[adonis_results_M6$column=="infant_ffq_feeding_mode",]


# ### ===== 10.7 COMBINE ALL RESULTS FROM ADONIS AND CORRECT FOR MULTIPLE TESTING ===== ###
# 
# ## combine data frames
# adonis_results <- as.data.frame(rbind(adonis_results_M1, adonis_results_M2, adonis_results_M3, adonis_results_M6))
# 
# ## correct for multiple testing
# adonis_results$FDR <- p.adjust(adonis_results$p, method="BH")
# 
# ## sort results by FDR and p
# adonis_results <- adonis_results[order(adonis_results$FDR, adonis_results$p),]
# 
# ## check results
# nrow(adonis_results[adonis_results$FDR<0.05 & !is.na(adonis_results$FDR),]) # 0 FDR significant associations
# nrow(adonis_results[adonis_results$p<0.05 & !is.na(adonis_results$p),])     #54 nominally significant associations
# 
# ## save results
# write.table(adonis_results, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/beta_diversity/240206_milk_beta_diversity_by_maternal_and_infant_phenotypes_adonis_results_all.txt", row.names=F, col.names=T, sep="\t", quote=F)

