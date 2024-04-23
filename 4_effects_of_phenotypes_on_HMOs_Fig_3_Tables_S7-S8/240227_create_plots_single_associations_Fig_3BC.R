################################################################################################################
### CHECK CORRELATION OF DYNAMIC PHENOTYPES WITH TIME POINT POSTPARTUM AND CREATE PLOTS FOR FIGURE 3B AND 3C ###
################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 0. DATA IMPORT
# 
# 1. SELECT DATA OF INTEREST
# 1.1 MATERNAL PHENOTYPE DATA
# 1.2 INFANT PHENOTYPE DATA
# 
# 2. ENSURE CORRECT DATA STRUCTURE
# 2.1 MATERNAL PHENOTYPES
# 2.2 INFANT PHENOTYPES
# 
# 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
# 3.1 MATERNAL PHENOTYPES
# 3.2 INFANT PHENOTYPES
# 
# 4. CHANGE SELECTED CATEGORICAL PHENOTYPE DATA TO NUMERIC FOR USE IN ASSOCIATION STUDIES
# 4.1 MATERNAL PHENOTYPES
# 4.2 INFANT PHENOTYPES
# 
# 5. CREATE TABLE SHOWING PHENOTYPE DATA TYPE (FACTOR OR NUMERIC/INTEGER)
# 5.1 FUNCTION TO CREATE OVERVIEW TABLES
# 5.2 MATERNAL PHENOTYPES
# 5.3 INFANT PHENOTYPES
# 
# 6. INVERSE-RANK TRANSFORM NUMERIC PHENOTYPES AND HMO DATA
# 6.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 6.2 INVERSE-RANK TRANSFORM NUMERIC PHENOTYPE DATA
# 6.3 INVERSE-RANK TRANSFORM HMO DATA
# 
# 7. CHECK CORRELATION OF DYNAMIC PHENOTYPES WITH TIME POINT POSTPARTUM
# 7.1 FUNCTION TO CHECK CORRELATION
# 7.2 CHECK CORRELATION OF MATERNAL DYNAMIC PHENOTYPES WITH TIME POSTPARTUM
# 7.3 CHECK CORRELATION OF INFANT DYNAMIC PHENOTYPES WITH TIME POSTPARTUM
# 
# 8. IMPORT MODEL RESULTS
# 
# 9. PLOT FDR-SIGNIFICANT ASSOCIATIONS
# 9.1 FUNCTION TO CREATE SCATTERPLOTS / BOXPLOTS
# 9.2 CREATE PLOTS FOR FDR-SIGNIFICANT ASSOCIATIONS OF MATERNAL PHENOTYPES WITH MILK HMOs
# 9.3 CREATE PLOTS FOR FDR-SIGNIFICANT ASSOCIATIONS OF INFANT PHENOTYPES WITH MILK HMOs
# 
# 10. PLOT A-TETRA BY SECRETOR STATUS AND BLOOD GROUP (FIGURE 3B)
# 
# 11. TEST EFFECT OF BLACK/GREEN/HERBAL TEA INTAKE ON 3'SL
# 11.1 FUNCTION FOR ASSOCIATION OF DYNAMIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID
# 11.2 RUN ASSOCIATION OF INDIVIDUAL TEA ITEMS WITH 3'SL CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID
# 
# 12. PLOT 3'SL BY TEA INTAKE
# 12.1 CREATE PLOTS FOR GROUPED TEA INTAKE
# 12.2 CREATE PLOTS FOR BLACK/GREEN/HERBAL TEA INTAKE (FIGURE 3C)
# 
# 13 COMBINE PLOTS FOR A-TETRA AND 3'SL (FIGURE 3B AND 3C)






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


##### =========================== 6. INVERSE-RANK TRANSFORM NUMERIC PHENOTYPES AND HMO DATA =========================== #####

### ===== 6.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 6.2 INVERSE-RANK TRANSFORM NUMERIC PHENOTYPE DATA ===== ###

## Maternal phenotypes
mhmo_invr <- mhmo
mhmo_invr[,numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"]] <- as.data.frame(apply(mhmo_invr[,numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"]], 2, invrank))

colnames(mhmo_invr)[which(colnames(mhmo_invr) %in% numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"])] <- paste0(colnames(mhmo_invr)[which(colnames(mhmo_invr) %in% numeric_mhmo[numeric_mhmo$model_exclusion=="include","column_name"])], "_invr")


## Infant phenotypes
ihmo_invr <- ihmo
ihmo_invr[,numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"]] <- as.data.frame(apply(ihmo_invr[,numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"]], 2, invrank))

colnames(ihmo_invr)[which(colnames(ihmo_invr) %in% numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"])] <- paste0(colnames(ihmo_invr)[which(colnames(ihmo_invr) %in% numeric_ihmo[numeric_ihmo$model_exclusion=="include","column_name"])], "_invr")


### ===== 6.3 INVERSE-RANK TRANSFORM HMO DATA ===== ###

## Maternal data frame
mhmo_invr[,grep("_ugml",colnames(mhmo_invr))] <- as.data.frame(apply(mhmo_invr[,grep("ugml", colnames(mhmo_invr))], 2, invrank))
colnames(mhmo_invr)[grep("_ugml", colnames(mhmo_invr))] <- paste0(colnames(mhmo_invr)[grep("_ugml", colnames(mhmo_invr))], "_invr")

## Infant data frame
ihmo_invr[,grep("_ugml",colnames(ihmo_invr))] <- as.data.frame(apply(ihmo_invr[,grep("ugml", colnames(ihmo_invr))], 2, invrank))
colnames(ihmo_invr)[grep("_ugml", colnames(ihmo_invr))] <- paste0(colnames(ihmo_invr)[grep("_ugml", colnames(ihmo_invr))], "_invr")


##### =========================== 7. CHECK CORRELATION OF DYNAMIC PHENOTYPES WITH TIME POINT POSTPARTUM =========================== #####

library(ggplot2)
library(ggpubr)


### ===== 7.1 FUNCTION TO CHECK CORRELATION ===== ###

check.correlation <- function(inputdata, columns){
  
  pheno <- c()
  r <- c()
  p <- c()
  
  for (i in columns){
    print(paste0("Run correlation test for ", i, "_", colnames(inputdata)[i]))
    
    pheno <- c(pheno, colnames(inputdata)[i])
    
    cor_res <- cor.test(inputdata$time_point_numeric, inputdata[,i], method="pearson")
    r <- c(r, cor_res$estimate)
    p <- c(p, cor_res$p.value)
    
    if (p[length(p)]<0.05){
      print(paste0("There was a significant association between time_point_numeric and ", colnames(inputdata)[i]))
      print("Generating correlation plot")
      
      corplot1 <- ggscatter(inputdata[!is.na(inputdata[,i]),], x="time_point_numeric", y=paste0(colnames(inputdata)[i]), cor.method = "pearson", cor.coef = T, add = "reg.line", conf.int = T)
      ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240227_corr_plot_", colnames(inputdata)[i], "by_time_point_numeric.pdf"), corplot1, dpi=500, width=20, height=20, unit="cm")
    }else{
      print(paste0("There was no significant association between time_point_numeric and ", colnames(inputdata)[i]))
    }
  }
  
  print("Combine results in data frame")
  res <- data.frame(phenotype=pheno, r=r, p=p)
  
  print("Correct for multiple testing")
  res$FDR <- p.adjust(res$p, method="BH")
  
  res <- res[order(res$FDR, res$p),]
  
  return(res)
  
}


### ===== 7.2 CHECK CORRELATION OF MATERNAL DYNAMIC PHENOTYPES WITH TIME POSTPARTUM ===== ###

numeric_mhmo_2 <- numeric_mhmo
numeric_mhmo_2[numeric_mhmo_2$model_exclusion=="include" & numeric_mhmo_2$model_use=="dynamic" & (numeric_mhmo_2$data_type=="numeric" | numeric_mhmo_2$data_type=="integer"), "column_name"] <- 
  c(paste0(numeric_mhmo_2[numeric_mhmo_2$model_exclusion=="include" & numeric_mhmo_2$model_use=="dynamic" & (numeric_mhmo_2$data_type=="numeric" | numeric_mhmo_2$data_type=="integer"), "column_name"], "_invr"))

matphenotypes <- colnames(mhmo_invr[,colnames(mhmo_invr) %in% numeric_mhmo_2[numeric_mhmo_2$model_exclusion=="include" & numeric_mhmo_2$model_use=="dynamic" & (numeric_mhmo_2$data_type=="numeric" | numeric_mhmo_2$data_type=="integer"), "column_name"]])

timecorrmatphenos <- check.correlation(mhmo_invr, columns=which(names(mhmo_invr) %in% matphenotypes))
write.table(timecorrmatphenos, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240227_correlation_of_dynamic_phenotypes_with_time_point.txt", row.names=F, col.names=T, sep="\t", quote=F)
timecorrmatphenos
# !! mother_milk_collection_volume_pumped_invr is significantly increasing with time point postpartum (problematic!) !!
# mother_dis_subclinical_mastitis_Na_K_ratio_invr is significantly decreasing with time point postpartum (not relevant - phenotype not used)
# mother_milk_collection_time_numeric_invr and mother_breastpump_times_per_day_invr are not significantly associated with time point postpartum

numeric_mhmo[numeric_mhmo$model_use=="dynamic" & numeric_mhmo$model_exclusion=="include",]
# -> There were no dynamic categorical maternal phenotypes.


### ===== 7.3 CHECK CORRELATION OF INFANT DYNAMIC PHENOTYPES WITH TIME POSTPARTUM ===== ###

numeric_ihmo_2 <- numeric_ihmo
numeric_ihmo_2[numeric_ihmo_2$model_exclusion=="include" & numeric_ihmo_2$model_use=="dynamic" & (numeric_ihmo_2$data_type=="numeric" | numeric_ihmo_2$data_type=="integer"), "column_name"] <- 
  c(paste0(numeric_ihmo_2[numeric_ihmo_2$model_exclusion=="include" & numeric_ihmo_2$model_use=="dynamic" & (numeric_ihmo_2$data_type=="numeric" | numeric_ihmo_2$data_type=="integer"), "column_name"], "_invr"))

infphenotypes <- colnames(ihmo_invr[,colnames(ihmo_invr) %in% numeric_ihmo_2[numeric_ihmo_2$model_exclusion=="include" & numeric_ihmo_2$model_use=="dynamic" & (numeric_ihmo_2$data_type=="numeric" | numeric_ihmo_2$data_type=="integer"), "column_name"]])
infphenotypes
# -> There were no dynamic numeric infant phenotypes.

# timecorrinfphenos <- check.correlation(ihmo_invr, columns=which(names(ihmo_invr) %in% infphenotypes))
# timecorrinfphenos


numeric_ihmo[numeric_ihmo$model_use=="dynamic" & numeric_ihmo$model_exclusion=="include",]
# -> There were no dynamic categorical infant phenotypes.


##### =========================== 8. IMPORT MODEL RESULTS =========================== #####

## import model results
model_results_short <- readRDS(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results/240227_model_results_milk_HMOs_real_by_phenotypes_n162.rds")

## save in new data frame and add a column showing the association (x and y)
model_results <- model_results_short
model_results$association <- c(paste0(model_results$x, "__", model_results$y))

## create vectors with FDR-significant associations, split by mother/infant
maternal_associations <- unique(model_results[model_results$FDR<0.05, "association"])[c(1:11)]
# infant_associations <- unique(model_results[model_results$FDR<0.05, "association"])[c()]


##### =========================== 9. PLOT FDR-SIGNIFICANT ASSOCIATIONS =========================== #####

### ===== 9.1 FUNCTION TO CREATE SCATTERPLOTS / BOXPLOTS ===== ###

create.plots <- function(inputdata, associations){
  
  for (i in associations){
    
    print(paste0("Starting plotting for ", i))
    
    my_x <- model_results[model_results$association==i,"x"]
    my_y <- model_results[model_results$association==i,"y"]
    
    my_class_x <- class(inputdata[,colnames(inputdata)==my_x])
    
    my_df <- inputdata[!is.na(inputdata[,colnames(inputdata)==my_x,]) & !is.na(inputdata[,colnames(inputdata)==my_y,]),]
    print(dim(my_df))
    
    my_x_name <- colnames(my_df)[colnames(my_df)==my_x]
    my_y_name <- colnames(my_df)[colnames(my_df)==my_y]
    print(my_x_name)
    print(my_y_name)
    
    if(my_class_x == "numeric" | my_class_x == "integer"){
      
      print("Generating scatterplot")
      
      scatterplot1 <- ggplot(my_df, aes(x=my_df[,colnames(my_df)==my_x], y=my_df[,colnames(my_df)==my_y]))+
        geom_point(alpha=0.5)+
        theme_bw()+
        geom_smooth(method='lm')+
        labs(x=paste0(my_x_name), y=paste0(my_y_name))+
        facet_grid(.~time_point_numeric)
      
      # ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240227_scatterplot_", my_y, "_by_", my_x, ".pdf"), scatterplot1, dpi=500, width=20, height=20, unit="cm")
      ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/plots_per_timepoint/240227_scatterplot_", my_y, "_by_", my_x, ".pdf"), scatterplot1, dpi=500, width=20, height=20, unit="cm")
      
    }else{
      
      print("Generating boxplot")
      
      boxplot1 <- ggplot(my_df, aes(x=my_df[,colnames(my_df)==my_x], y=my_df[,colnames(my_df)==my_y], color=my_df[,colnames(my_df)==my_x]))+
        geom_jitter(alpha=0.5)+
        geom_boxplot(alpha=0)+
        theme_bw()+
        theme(legend.position = "none")+
        labs(x=paste0(my_x_name), y=paste0(my_y_name))+
        facet_grid(.~time_point_numeric)
      
      # ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240227_boxplot_", my_y_name, "_by_", my_x_name, ".pdf"), boxplot1, dpi=500, width=20, height=20, unit="cm")
      ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/plots_per_timepoint/240227_boxplot_", my_y_name, "_by_", my_x_name, ".pdf"), boxplot1, dpi=500, width=20, height=20, unit="cm")
      
    }
  }
}


### ===== 9.2 CREATE PLOTS FOR FDR-SIGNIFICANT ASSOCIATIONS OF MATERNAL PHENOTYPES WITH MILK HMOs ===== ###

create.plots(mhmo_invr, maternal_associations)


### ===== 9.3 CREATE PLOTS FOR FDR-SIGNIFICANT ASSOCIATIONS OF INFANT PHENOTYPES WITH MILK HMOs ===== ###

# create.plots(ihmo_invr, infant_associations)
# There are no significant associations with infant phenotypes.


##### =========================== 10. PLOT A-TETRA BY SECRETOR STATUS AND BLOOD GROUP (FIGURE 3B) =========================== #####

## plot with invr-transformed A-tetra levels
atetra_plot1 <- ggplot(mhmo_invr[!is.na(mhmo_invr$mother_blood_group),], aes(x=mother_blood_group, y=mother_milk_HMO_A_tetra_ugml_invr, color=mother_blood_group))+
  geom_jitter(alpha=0.5, size=1)+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  facet_grid(.~mother_milk_HMO_Se)+
  scale_color_manual(values = c("O" = "lightcyan3", "A" = "deeppink4", "B" = "gold", "AB" = "chocolate3"))+
  labs(x="Maternal blood group", y="Inverse-rank transformed\nA-tetra abundance [µg/ml]")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_boxplot_mother_milk_HMO_A_tetra_ugml_invr_by_Se_and_bloodgroup.pdf", atetra_plot1, dpi=500, width=10, height=7, unit="cm")


## plot with non-transformed A-tetra levels
atetra_plot2 <- ggplot(mhmo[!is.na(mhmo$mother_blood_group),], aes(x=mother_blood_group, y=mother_milk_HMO_A_tetra_ugml, color=mother_blood_group))+
  geom_jitter(alpha=0.5, size=1)+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  facet_grid(.~mother_milk_HMO_Se)+
  scale_y_continuous(limits = c(-1,700), breaks = c(0,350,700), labels = c(0,350,700))+
  scale_color_manual(values = c("O" = "lightcyan3", "A" = "deeppink4", "B" = "gold", "AB" = "chocolate3"))+
  labs(x="Maternal blood group", y="A-tetra abundance [µg/ml]")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_boxplot_mother_milk_HMO_A_tetra_ugml_by_Se_and_bloodgroup.pdf", atetra_plot2, dpi=500, width=10, height=7, unit="cm")


##### =========================== 11. TEST EFFECT OF BLACK/GREEN/HERBAL TEA INTAKE ON 3'SL =========================== #####

### ===== 11.1 FUNCTION FOR ASSOCIATION OF DYNAMIC PHENOTYPES WITH HMO CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID ===== ###

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
      print(sum1)
      
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


### ===== 11.2 RUN ASSOCIATION OF INDIVIDUAL TEA ITEMS WITH 3'SL CONCENTRATIONS WITH CORRECTION FOR MILK GROUP, TIME AND ID ===== ###

tea_items_results <- run.lmer.cor.milkgroup.time.ID(datasetname = "mother_real_HMOs_dynamic",
                                                    inputdata   = mhmo_invr,
                                                    xcolumn     = c(300:302),   #invr-transformed individual tea items
                                                    ycolumn     = c(445))       #invr-transformed 3'SL concentrations


##### =========================== 12. PLOT 3'SL BY TEA INTAKE =========================== #####

### ===== 12.1 CREATE PLOTS FOR GROUPED TEA INTAKE ===== ###

## plot with invr-transformed 3'SL levels
tea_3SL_plot1 <- ggplot(mhmo_invr[!is.na(mhmo_invr$mother_ffq_group_tea_ml_per_day_invr),], aes(x=mother_ffq_group_tea_ml_per_day_invr, y=mother_milk_HMO_3SL_ugml_invr))+
  geom_point(alpha=0.5, size=1, color="lightblue3")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_smooth(method='lm', color="black", fill="lightblue4")+
  labs(x="Inverse-rank transformed\nmaternal tea intake (ml/day)", y="Inverse-rank transformed\n3'SL abundance [µg/ml]")+
  scale_x_continuous(limits = c(-2.2,3))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_scatterplot_mother_milk_HMO_3SL_ugml_invr_by_mother_ffq_group_tea_ml_per_day_invr_pretty.pdf", tea_3SL_plot1, dpi=500, width=7, height=7, unit="cm")


## plot with non-transformed 3'SL levels
tea_3SL_plot2 <- ggplot(mhmo[!is.na(mhmo$mother_ffq_group_tea_ml_per_day),], aes(x=mother_ffq_group_tea_ml_per_day, y=mother_milk_HMO_3SL_ugml))+
  geom_point(alpha=0.5, size=1, color="lightblue3")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_smooth(method='lm', color="black", fill="lightblue4")+
  labs(x="Maternal tea intake (ml/day)", y="3'SL abundance [µg/ml]")+
  scale_y_continuous(limits = c(0,500))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_scatterplot_mother_milk_HMO_3SL_ugml_by_mother_ffq_group_tea_ml_per_day_pretty.pdf", tea_3SL_plot2, dpi=500, width=7, height=7, unit="cm")


### ===== 12.2 CREATE PLOTS FOR BLACK/GREEN/HERBAL TEA INTAKE (FIGURE 3C) ===== ###

black_tea_3SL_plot2 <- ggplot(mhmo[!is.na(mhmo$mother_ffq_item_tea_black_ml_per_day),], aes(x=(mother_ffq_item_tea_black_ml_per_day/1000), y=mother_milk_HMO_3SL_ugml))+
  geom_point(alpha=0.5, size=1, color="sienna4")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_smooth(method='lm', color="black", fill="sienna4")+
  labs(x="Maternal black tea intake (l/day)", y="3'SL abundance [µg/ml]")+
  scale_x_continuous(limits = c(0,1.5))+
  scale_y_continuous(limits = c(0,500))+
  annotate(geom="text", x=0.75, y=475, label="p=0.16")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_scatterplot_mother_milk_HMO_3SL_ugml_by_mother_ffq_item_tea_black_ml_per_day_pretty.pdf", black_tea_3SL_plot2, dpi=500, width=7, height=7, unit="cm")

green_tea_3SL_plot2 <- ggplot(mhmo[!is.na(mhmo$mother_ffq_item_tea_green_ml_per_day),], aes(x=(mother_ffq_item_tea_green_ml_per_day/1000), y=mother_milk_HMO_3SL_ugml))+
  geom_point(alpha=0.5, size=1, color="darkolivegreen3")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_smooth(method='lm', color="black", fill="darkolivegreen3")+
  labs(x="Maternal green tea intake (l/day)", y="3'SL abundance [µg/ml]")+
  scale_x_continuous(limits = c(0,1.5))+
  scale_y_continuous(limits = c(0,500))+
  annotate(geom="text", x=0.75, y=475, label="p=0.01")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_scatterplot_mother_milk_HMO_3SL_ugml_by_mother_ffq_item_tea_green_ml_per_day_pretty.pdf", green_tea_3SL_plot2, dpi=500, width=7, height=7, unit="cm")

herbs_tea_3SL_plot2 <- ggplot(mhmo[!is.na(mhmo$mother_ffq_item_tea_herbs_ml_per_day),], aes(x=(mother_ffq_item_tea_herbs_ml_per_day/1000), y=mother_milk_HMO_3SL_ugml))+
  geom_point(alpha=0.5, size=1, color="goldenrod3")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_smooth(method='lm', color="black", fill="goldenrod3")+
  labs(x="Maternal herbal tea intake (l/day)", y="3'SL abundance [µg/ml]")+
  scale_x_continuous(limits = c(0,1.5))+
  scale_y_continuous(limits = c(0,500))+
  annotate(geom="text", x=0.75, y=475, label="p=0.03")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_scatterplot_mother_milk_HMO_3SL_ugml_by_mother_ffq_item_tea_herbs_ml_per_day_pretty.pdf", herbs_tea_3SL_plot2, dpi=500, width=7, height=7, unit="cm")


## combine plot for 3 different tea items
library(cowplot)
tea_items_3SL_plots <- plot_grid(black_tea_3SL_plot2, green_tea_3SL_plot2, herbs_tea_3SL_plot2, nrow=1)
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_scatterplot_mother_milk_HMO_3SL_ugml_by_black_green_and_herbal_tea_pretty.pdf", tea_items_3SL_plots, dpi=500, width=21, height=7, unit="cm")


##### =========================== 13 COMBINE PLOTS FOR A-TETRA AND 3'SL (FIGURE 3B AND 3C) =========================== #####

## plot with non-transformed A-tetra levels, without facet grid box
atetra_plot3 <- ggplot(mhmo[!is.na(mhmo$mother_blood_group),], aes(x=mother_blood_group, y=mother_milk_HMO_A_tetra_ugml, color=mother_blood_group))+
  geom_jitter(alpha=0.5, size=1)+
  geom_boxplot(alpha=0, color="black")+
  theme_bw()+
  facet_grid(.~mother_milk_HMO_Se)+
  theme(legend.position = "none",
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text = element_blank())+
  scale_y_continuous(limits = c(-1,700), breaks = c(0,350,700), labels = c(0,350,700))+
  scale_color_manual(values = c("O" = "lightcyan3", "A" = "deeppink4", "B" = "gold", "AB" = "chocolate3"))+
  labs(x="Maternal blood group", y="A-tetra abundance [µg/ml]")

all_plots <- plot_grid(atetra_plot3, tea_items_3SL_plots, rel_widths = c(1,2), nrow=1)
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results_single_plots/240229_A_tetra_and_3SL_plots_pretty.pdf", all_plots, dpi=500, width=35, height=7, unit="cm")















