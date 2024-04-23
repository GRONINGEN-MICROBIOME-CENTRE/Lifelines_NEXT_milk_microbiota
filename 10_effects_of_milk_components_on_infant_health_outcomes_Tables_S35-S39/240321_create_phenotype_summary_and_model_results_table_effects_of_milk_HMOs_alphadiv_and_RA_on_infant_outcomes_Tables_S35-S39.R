####################################################################################################################################################################################
### CREATE PHENOTYPES SUMMARY STATISTICS + MODEL RESULTS TABLES FOR EFFECT OF MILK COMPOSITION (HMOs, MICROBIOTA) ON INFANT OUTCOMES FOR MILK COMPOSITION PAPER (TABLES S35-S39) ###
####################################################################################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessions type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 0. DATA IMPORT
# 0.1 FILE WITH PHENOTYPES + MILK HMOs
# 0.2 FILE WITH PHENOTYPES + MILK MICROBIOTA DATA
# 
# 1. SELECT DATA OF INTEREST
# 1.1 INFANT OUTCOME PHENOTYPES AND MILK HMOs
# 1.2 INFANT OUTCOME PHENOTYPES AND MILK MICROBIOTA DATA
# 
# 2. ENSURE CORRECT DATA STRUCTURE
# 2.1 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK HMOs
# 2.2 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK MICROBIOTA DATA
# 
# 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
# 3.1 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK HMOs
# 3.2 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK MICROBIOTA DATA
# 
# 4. CHANGE SELECTED CATEGORICAL PHENOTYPE DATA TO NUMERIC FOR USE IN ASSOCIATION STUDIES
# 4.1 CHECK WHICH PHENOTYPES TO CONVERT TO NUMERIC
# 4.2 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK HMOs
# 4.3 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK MICROBIOTA DATA
# 
# 5. CREATE HISTOGRAMS TO SHOW DATA DISTRIBUTION TO GUIDE DECISION ON PHENOTYPE DATA TRANSFORMATION
# 5.1 INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs
# 5.2 INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA
# 
# 6. CREATE SUMMARY STATISTICS FOR CATEGORICAL PHENOTYPES
# 6.1 FUNCTION TO CREATE SUMMARY STATISTICS
# 6.2 SUMMARY STATISTICS FOR CATEGORICAL INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs
# 6.3 SUMMARY STATISTICS FOR CATEGORICAL INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA
# 
# 7. CREATE SUMMARY STATISTICS FOR NUMERIC/INTEGER PHENOTYPES
# 7.1 FUNCTION TO CREATE SUMMARY STATISTICS\
# 7.2 SUMMARY STATISTICS FOR INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs
# 7.3 SUMMARY STATISTICS FOR INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA
# [7.4 SUMMARY STATISTICS FOR MILK HMOs]
# 7.5 SUMMARY STATISTICS FOR MILK ALPHA DIVERSITY
# 
# 8. COMBINE SUMMARY STATISTICS TABLES FOR ALL INFANT OUTCOME PHENOTYPES
# 8.1 COMBINE SUMMARY STATISTICS TABLES WITH CATEGORICAL AND NUMERIC PHENOTYPES FOR INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs (TABLE S35A-D)
# 8.2 COMBINE SUMMARY STATISTICS TABLES WITH CATEGORICAL AND NUMERIC PHENOTYPES FOR INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA (TABLE S37A-C)
# 
# 9. CALCULATE MILK RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE
# 9.1 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE
# 9.2 SUMMARY STATISTICS FOR MILK RELATIVE BACTERIAL ABUNDANCES
# 
# 10. INVERSE-RANK TRANSFORM MILK HMOs AND INFANT FAECAL ALPHA DIVERSITY DATA
# 10.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 10.2. INVERSE-RANK TRANSFORM MILK HMOs AND INFANT OUTCOME PHENOTYPES
# 10.3. INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY MEASURES AND INFANT OUTCOME PHENOTYPES
# 
# 11. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED DATA
# 11.1 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs
# 11.2 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA
# 11.3 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED MILK HMOs
# 11.4 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED MILK ALPHA DIVERSITY MEASURES
# 
# 12. CLR-TRANSFORM MILK RELATIVE ABUNDANCES 
# 12.1 FUNCTION FOR CLR-TRANSFORMATION
# 12.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES
# 
# 13. HISTOGRAMS FOR CLR-TRANSFORMED DATA
# 13.1 HISTOGRAMS FOR CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES
# 
# 14. MODEL FUNCTION (LMER/LM)
# 
# 15. DETAILS ABOUT ihmo_num_invr DATA FRAME FOR ASSOCIATION OF BASIC PHENOTYPES AND MILK HMOs WITH INFANT OUTCOMES
# 
# 16. CHECK EFFECTS OF BASIC PHENOTYPES ON INFANT OUTCOMES
# 16.1 RUN MODELS FOR NUMERIC/CATEGORICAL INFANT OUTCOMES ~ NUMERIC BASIC PHENOTYPE
# 16.2 RUN MODELS FOR NUMERIC INFANT OUTCOMES ~ CATEGORICAL BASIC PHENOTYPE
# 16.3 RUN FISHER'S TESTS FOR CATEGORICAL BASIC PHENOTYPE ~ CATEGORICAL INFANT OUTCOMES
# 16.4 COMBINE RESULTS PER TIME POINT AND CORRECT FOR MULTIPLE TESTING
# 16.5 CHECK WHICH RESULTS ARE CONSISTENT ACROSS TIME POINTS
# ## -> When investigating the effect of milk components on infant growth, correct for infant sex and gestational age in the models!
# 
# 17. ASSOCIATIONS OF MILK HMO CONCENTRATIONS WITH INFANT OUTCOME PHENOTYPES
# 17.1 RUN MODELS FOR NUMERIC/CATEGORICAL INFANT OUTCOMES ~ NUMERIC MILK HMO CONCENTRATIONS
# 17.2 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ NUMERIC MILK HMO CONCENTRATIONS
# [17.3 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ NUMERIC MILK HMO CONCENTRATIONS - STRATIFIED BY SEX]
# 
# 18. ASSOCIATIONS OF MILK MICROBIOTA (ALPHA DIVERSITY, RELATIVE ABUNDANCES) WITH INFANT OUTCOME PHENOTYPES
# 18.1 RUN MODELS FOR NUMERIC/CATEGORICAL INFANT OUTCOMES ~ NUMERIC MILK MICROBIOTA DATA
# 18.2 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ NUMERIC MILK MICROBIOTA DATA
# [18.3 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ NUMERIC MILK MICROBIOTA DATA - STRATIFIED BY SEX]
# 
# 19. ASSOCIATIONS OF MILK GROUPS WITH INFANT OUTCOME PHENOTYPES
# 19.1 RUN MODELS FOR NUMERIC INFANT OUTCOMES ~ MILK GROUPS
# 19.2 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ MILK GROUPS
# [19.3 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ MILK GROUPS - STRATIFIED BY INFANT SEX]
# 19.4 RUN FISHER'S TEST FOR CATEGORICAL INFANT OUTCOMES ~ MILK GROUPS
# 
# 20. COMBINE RESULTS FROM ASSOCIATIONS OF MILK GROUPS AND MILK HMOs WITH INFANT OUTCOME PHENOTYPES
# 20.1 COMBINE RESULTS PER TIME POINT
# [20.2 CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS]
# 20.3 MORE FOCUSED ANALYSIS: ONLY USE SELECTED TIME POINTS AND PHENOTYPES AND RECALCULATE FDR (TABLE S36A-D)
# 20.4 MORE FOCUSED ANALYSIS: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS
# 
# 21. COMBINE RESULTS FROM ASSOCIATIONS OF MILK MILK MICROBIOTA DATA WITH INFANT OUTCOME PHENOTYPES
# 21.1 COMBINE RESULTS PER TIME POINT
# [21.2 MILK ALPHA DIVERSITY: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS]
# [21.3 MILK RELATIVE ABUNDANCES: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS]
# 21.4 MILK ALPHA DIVERSITY: MORE FOCUSED ANALYSIS: ONLY USE SELECTED TIME POINTS AND PHENOTYPES AND RECALCULATE FDR (TABLE S38A-C)
# 21.5 MILK ALPHA DIVERSITY: MORE FOCUSED ANALYSIS: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS
# 21.6 MILK RELATIVE ABUNDANCES: MORE FOCUSED ANALYSIS: ONLY USE SELECTED TIME POINTS AND PHENOTYPES AND RECALCULATE FDR (TABLE S39A-C)
# 21.7 MILK RELATIVE ABUNDANCES: MORE FOCUSED ANALYSIS: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes")


### ===== 0.1 FILE WITH PHENOTYPES + MILK HMOs ===== ###

## import file with phenotypes + milk HMO data
hmo <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/phenotypes/240108_milk_composition_paper_sel_phenotypes_hmo_data_n1563.txt", header=T, sep="\t", stringsAsFactors=T)
#1563x512

# Notes:
# This file contains the following:
# - Phenotypes linked to milk, maternal faecal and infant faecal microbiota data.
#   The file also includes real measured levels of 24 single HMOs and 4 grouped HMOs in µg/ml and HMO-based Le and Se status and milk groups as phenotypes.
#   !! Note that maternal phenotype and HMO data was duplicated for twins. This is indicated by ";1" (single baby/first baby) and ";2" (twin baby/duplicated data) in the mother_ID and mother_sample_ID.
#      -> For association with maternal phenotypes (i.e. non-duplicated data), grep only ";1"
#      -> For association with infant phenotypes, use both ";1" and ";2"


### ===== 0.2 FILE WITH PHENOTYPES + MILK MICROBIOTA DATA ===== ###

## import file with phenotypes + microbiota data on genus level
bact <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/231129_milk_composition_paper_sel_phenotypes_sel_microbiota_data_genera_n1515.rds")
#1515x1154

# resort columns
bact <- bact[,c(2:14,1,15:ncol(bact))]


##### =========================== 1. SELECT DATA OF INTEREST =========================== #####

## Note: We keep mothers with duplicated data for twins as we are looking at infant outcomes here!


### ===== 1.1 INFANT OUTCOME PHENOTYPES AND MILK HMOs ===== ###

ihmo <- hmo[,c(1:15,                    #IDs etc.
               25,28,482:484,           #maternal Le/Se status (genetic- and HMO-based)
               435:437,                 #infant Se status (genetic)
               396:397,438:443,450:453, #maternal/infant phenotypes
               485:512,                 #milk HMOs
               444:449,454:480)]        #infant outcomes
#1563x96, 33 of which are infant outcomes


### ===== 1.2 INFANT OUTCOME PHENOTYPES AND MILK MICROBIOTA DATA ===== ###

ibact <- bact[bact$sample_origin_type=="mother_human_milk"
              ,c(1:15,504,507,            #IDs, incl. seq ID and sample type
                 511,526,                 #DNA isolation batch, seq reads
                 25:26,473:475,           #maternal Le/Se status (genetic- and HMO-based)
                 427:429,                 #infant Se status (genetic)
                 394:395,430:434,442:445, #maternal/infant phenotypes (gestational age, infant sex, delivery place and mode, infant feeding)
                 476:503,                 #milk HMOs
                 528:530,                 #alpha diversity
                 531:1154,                #relative bacterial abundances
                 436:441,446:471)]        #infant outcomes
#737x725, 32 of which are infant outcomes


## check which infant outcome is missing for milk microbiota association
setdiff(colnames(ihmo)[64:96], colnames(ibact)[694:725])
#"infant_dis_ROME_colic" is only available for HMOs, too small n for association with milk microbiota


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

### ===== 2.1 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK HMOs ===== ###

## check that data structure is correct
str(ihmo[,1:50])
str(ihmo[,c(51:96)])
# -> data structure is correct


### ===== 2.2 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK MICROBIOTA DATA ===== ###

## check that data structure is correct
str(ibact[,1:50])
str(ibact[,c(51:70)])
str(ibact[,c(692:725)])
# -> data structure is correct


##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

### ===== 3.1 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK HMOs ===== ###

## for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
ihmo$mother_birth_gestational_age_categories <- factor(ihmo$mother_birth_gestational_age_categories, levels = c("term", "preterm"))

ihmo$infant_birth_delivery_mode <- factor(ihmo$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
ihmo$infant_birth_delivery_mode_detailed <- factor(ihmo$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                        "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
ihmo$infant_ffq_feeding_type_delivery <- factor(ihmo$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))

# also change it for milk groups
ihmo$mother_milk_HMO_milk_group <- factor(ihmo$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


### ===== 3.2 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK MICROBIOTA DATA ===== ###

## for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
ibact$mother_birth_gestational_age_categories <- factor(ibact$mother_birth_gestational_age_categories, levels = c("term", "preterm"))

ibact$infant_birth_delivery_mode <- factor(ibact$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
ibact$infant_birth_delivery_mode_detailed <- factor(ibact$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                          "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
ibact$infant_ffq_feeding_type_delivery <- factor(ibact$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))

# also change it for milk groups
ibact$mother_milk_HMO_milk_group <- factor(ibact$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 4. CHANGE SELECTED CATEGORICAL PHENOTYPE DATA TO NUMERIC FOR USE IN ASSOCIATION STUDIES  =========================== #####

## Some categorical phenotypes will be used as numeric phenotypes in models. These phenotypes are here changed to numeric before making the summary statistics tables so that the tables show the data as used for associations.
## I will later add info to the supplementary table on which factors/categories the numbers correspond to.

### ===== 4.1 CHECK WHICH PHENOTYPES TO CONVERT TO NUMERIC ===== ###

## check infant outcomes
for (i in 64:96){
  print(colnames(ihmo)[i])
  print(table(ihmo[,i], useNA="ifany"))
}

numcat_phenos <- c("infant_BITSS", "infant_igsq_2_freq_difficult_stool_past_week", "infant_igsq_3_freq_spitting_milk_normal_day_past_week", 
                   "infant_igsq_5_freq_discomfort_spitting_milk_past_week", "infant_igsq_6_freq_pain_spitting_milk_past_week", "infant_igsq_7_crying_duration_per_day_past_week", 
                   "infant_igsq_8_freq_unstoppable_crying_past_week", "infant_igsq_9_freq_crying_after_feeding_past_week", "infant_igsq_11_freq_unstoppable_restlessness_baby_past_week", 
                   "infant_igsq_12_freq_gas_per_day_past_week", "infant_igsq_13_freq_discomfort_gas_past_week", "infant_dis_ROME_e13_stool_freq_last_month", 
                   "infant_dis_ROME_e14_stool_structure", "infant_dis_ROME_e25_mucus_stool_last_week", "infant_sleep_pattern_problem")
# -> 15 phenotypes will be changed to numeric


### ===== 4.2 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK HMOs ===== ###

## create new data frame
ihmo_num <- ihmo

## change data entries for all numcat_phenos to numeric
for (i in numcat_phenos){ihmo_num[,i] <- as.numeric(as.character(gsub("__.*", "", ihmo_num[,i])))}


### ===== 4.3 DATA FRAME WITH INFANT OUTCOME PHENOTYPES AND MILK MICROBIOTA DATA ===== ###

## create new data frame
ibact_num <- ibact

## change data entries for all numcat_phenos to numeric
for (i in numcat_phenos){ibact_num[,i] <- as.numeric(as.character(gsub("__.*", "", ibact_num[,i])))}


##### =========================== 5. CREATE HISTOGRAMS TO SHOW DATA DISTRIBUTION TO GUIDE DECISION ON PHENOTYPE DATA TRANSFORMATION  =========================== #####

library(ggplot2)


### ===== 5.1 INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs ===== ###

## create histograms for all numeric infant outcome phenotypes
pdf(paste0("phenotype_data_distribution/240222_milk_HMOs_real_infant_outcome_phenotypes_data_distribution_histograms_1_n29.pdf"))
for (i in c(64:69,71:84,87:89,91:96)){
  print(paste0("Creating histogram for ", colnames(ihmo_num)[i]))
  histo <- ggplot(ihmo_num[!is.na(ihmo_num[,i]),], aes(x=ihmo_num[!is.na(ihmo_num[,i]),i]))+
    geom_histogram(bins=50)+
    labs(colnames(ihmo_num)[i])+
    ggtitle(paste0(colnames(ihmo_num)[i]))
  print(histo)
}
dev.off()


### ===== 5.2 INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA ===== ###

## create histograms for all numeric infant outcome phenotypes
pdf(paste0("phenotype_data_distribution/240222_milk_microbiota_infant_outcome_phenotypes_data_distribution_histograms_1_n29.pdf"))
for (i in c(694:699,701:714,716:718,720:725)){
  print(paste0("Creating histogram for ", colnames(ibact_num)[i]))
  histo <- ggplot(ibact_num[!is.na(ibact_num[,i]),], aes(x=ibact_num[!is.na(ibact_num[,i]),i]))+
    geom_histogram(bins=50)+
    labs(colnames(ibact_num)[i])+
    ggtitle(paste0(colnames(ibact_num)[i]))
  print(histo)
}
dev.off()


##### =========================== 6. CREATE SUMMARY STATISTICS FOR CATEGORICAL PHENOTYPES  =========================== #####

### ===== 6.1 FUNCTION TO CREATE SUMMARY STATISTICS ===== ###

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


### ===== 6.2 SUMMARY STATISTICS FOR CATEGORICAL INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs ===== ###

## create summary
factor_ihmo_sum <- create.summary.factors(ihmo_num[,c(70,85,86,90)])

## polish perc
factor_ihmo_sum$perc_answered <- signif(factor_ihmo_sum$perc_answered, 2)
factor_ihmo_sum$perc_missing <- signif(factor_ihmo_sum$perc_missing, 2)

## change colnames
colnames(factor_ihmo_sum)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== W2 ===== ###

## create summary
factor_ihmo_sum_W2 <- create.summary.factors(ihmo_num[ihmo_num$time_point=="0.5_months",c(70,85,86,90)])

## polish perc
factor_ihmo_sum_W2$perc_answered <- signif(factor_ihmo_sum_W2$perc_answered, 2)
factor_ihmo_sum_W2$perc_missing <- signif(factor_ihmo_sum_W2$perc_missing, 2)

## change colnames
colnames(factor_ihmo_sum_W2)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== M1 ===== ###

## create summary
factor_ihmo_sum_M1 <- create.summary.factors(ihmo_num[ihmo_num$time_point=="1_month",c(70,85,86,90)])

## polish perc
factor_ihmo_sum_M1$perc_answered <- signif(factor_ihmo_sum_M1$perc_answered, 2)
factor_ihmo_sum_M1$perc_missing <- signif(factor_ihmo_sum_M1$perc_missing, 2)

## change colnames
colnames(factor_ihmo_sum_M1)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== M2 ===== ###

## create summary
factor_ihmo_sum_M2 <- create.summary.factors(ihmo_num[ihmo_num$time_point=="2_months",c(70,85,86,90)])

## polish perc
factor_ihmo_sum_M2$perc_answered <- signif(factor_ihmo_sum_M2$perc_answered, 2)
factor_ihmo_sum_M2$perc_missing <- signif(factor_ihmo_sum_M2$perc_missing, 2)

## change colnames
colnames(factor_ihmo_sum_M2)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== M3 ===== ###

## create summary
factor_ihmo_sum_M3 <- create.summary.factors(ihmo_num[ihmo_num$time_point=="3_months",c(70,85,86,90)])

## polish perc
factor_ihmo_sum_M3$perc_answered <- signif(factor_ihmo_sum_M3$perc_answered, 2)
factor_ihmo_sum_M3$perc_missing <- signif(factor_ihmo_sum_M3$perc_missing, 2)

## change colnames
colnames(factor_ihmo_sum_M3)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== M6 ===== ###

## create summary
factor_ihmo_sum_M6 <- create.summary.factors(ihmo_num[ihmo_num$time_point=="6_months",c(70,85,86,90)])

## polish perc
factor_ihmo_sum_M6$perc_answered <- signif(factor_ihmo_sum_M6$perc_answered, 2)
factor_ihmo_sum_M6$perc_missing <- signif(factor_ihmo_sum_M6$perc_missing, 2)

## change colnames
colnames(factor_ihmo_sum_M6)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== 6.3 SUMMARY STATISTICS FOR CATEGORICAL INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA ===== ###

## create summary
factor_ibact_sum <- create.summary.factors(ibact_num[,c(700,715,719)])

## polish perc
factor_ibact_sum$perc_answered <- signif(factor_ibact_sum$perc_answered, 2)
factor_ibact_sum$perc_missing <- signif(factor_ibact_sum$perc_missing, 2)

## change colnames
colnames(factor_ibact_sum)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== M1 ===== ###

## create summary
factor_ibact_sum_M1 <- create.summary.factors(ibact_num[ibact_num$time_point=="1_month",c(700,715,719)])

## polish perc
factor_ibact_sum_M1$perc_answered <- signif(factor_ibact_sum_M1$perc_answered, 2)
factor_ibact_sum_M1$perc_missing <- signif(factor_ibact_sum_M1$perc_missing, 2)

## change colnames
colnames(factor_ibact_sum_M1)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== M2 ===== ###

## create summary
factor_ibact_sum_M2 <- create.summary.factors(ibact_num[ibact_num$time_point=="2_months",c(700,715,719)])

## polish perc
factor_ibact_sum_M2$perc_answered <- signif(factor_ibact_sum_M2$perc_answered, 2)
factor_ibact_sum_M2$perc_missing <- signif(factor_ibact_sum_M2$perc_missing, 2)

## change colnames
colnames(factor_ibact_sum_M2)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== M3 ===== ###

## create summary
factor_ibact_sum_M3 <- create.summary.factors(ibact_num[ibact_num$time_point=="3_months",c(700,715,719)])

## polish perc
factor_ibact_sum_M3$perc_answered <- signif(factor_ibact_sum_M3$perc_answered, 2)
factor_ibact_sum_M3$perc_missing <- signif(factor_ibact_sum_M3$perc_missing, 2)

## change colnames
colnames(factor_ibact_sum_M3)[1:3] <- c("Phenotype", "Categories", "n_per_category")


### ===== M6 ===== ###

## create summary
factor_ibact_sum_M6 <- create.summary.factors(ibact_num[ibact_num$time_point=="6_months",c(700,715,719)])

## polish perc
factor_ibact_sum_M6$perc_answered <- signif(factor_ibact_sum_M6$perc_answered, 2)
factor_ibact_sum_M6$perc_missing <- signif(factor_ibact_sum_M6$perc_missing, 2)

## change colnames
colnames(factor_ibact_sum_M6)[1:3] <- c("Phenotype", "Categories", "n_per_category")


##### =========================== 7. CREATE SUMMARY STATISTICS FOR NUMERIC/INTEGER PHENOTYPES  =========================== #####

### ===== 7.1 FUNCTION TO CREATE SUMMARY STATISTICS ===== ###

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


### ===== 7.2 SUMMARY STATISTICS FOR INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs ===== ###

## calculate summary statistics
outcome_sumstats_ihmo <- create.summary.numeric(ihmo_num[,c(64:69,71:84,87:89,91:96)])

## resort columns
outcome_sumstats_ihmo <- outcome_sumstats_ihmo[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ihmo$Note <- NA
# outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo[outcome_sumstats_ihmo$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 64:96){
#   print(colnames(ihmo)[i])
#   print(table(ihmo[,i], useNA="ifany"))
# }


## save phenotype summary statistics table
# write.table(outcome_sumstats_ihmo, file="phenotype_summary_statistics/240222_milk_HMOs_real_numeric_infant_outcome_phenotype_summary_statistics_n29.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== W2 ===== ###

## calculate summary statistics
outcome_sumstats_ihmo_W2 <- create.summary.numeric(ihmo_num[ihmo_num$time_point=="0.5_months",c(64:69,71:84,87:89,91:96)])

## resort columns
outcome_sumstats_ihmo_W2 <- outcome_sumstats_ihmo_W2[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ihmo_W2$Note <- NA
# outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_W2[outcome_sumstats_ihmo_W2$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 64:96){
#   print(colnames(ihmo)[i])
#   print(table(ihmo[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ihmo_W2, file="phenotype_summary_statistics/240311_milk_HMOs_real_numeric_infant_outcome_phenotype_summary_statistics_n29_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M1 ===== ###

## calculate summary statistics
outcome_sumstats_ihmo_M1 <- create.summary.numeric(ihmo_num[ihmo_num$time_point=="1_month",c(64:69,71:84,87:89,91:96)])

## resort columns
outcome_sumstats_ihmo_M1 <- outcome_sumstats_ihmo_M1[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ihmo_M1$Note <- NA
# outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M1[outcome_sumstats_ihmo_M1$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 64:96){
#   print(colnames(ihmo)[i])
#   print(table(ihmo[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ihmo_M1, file="phenotype_summary_statistics/240311_milk_HMOs_real_numeric_infant_outcome_phenotype_summary_statistics_n29_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M2 ===== ###

## calculate summary statistics
outcome_sumstats_ihmo_M2 <- create.summary.numeric(ihmo_num[ihmo_num$time_point=="2_months",c(64:69,71:84,87:89,91:96)])

## resort columns
outcome_sumstats_ihmo_M2 <- outcome_sumstats_ihmo_M2[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ihmo_M2$Note <- NA
# outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M2[outcome_sumstats_ihmo_M2$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 64:96){
#   print(colnames(ihmo)[i])
#   print(table(ihmo[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ihmo_M2, file="phenotype_summary_statistics/240311_milk_HMOs_real_numeric_infant_outcome_phenotype_summary_statistics_n29_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M3 ===== ###

## calculate summary statistics
outcome_sumstats_ihmo_M3 <- create.summary.numeric(ihmo_num[ihmo_num$time_point=="3_months",c(64:69,71:84,87:89,91:96)])

## resort columns
outcome_sumstats_ihmo_M3 <- outcome_sumstats_ihmo_M3[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ihmo_M3$Note <- NA
# outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M3[outcome_sumstats_ihmo_M3$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 64:96){
#   print(colnames(ihmo)[i])
#   print(table(ihmo[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ihmo_M3, file="phenotype_summary_statistics/240311_milk_HMOs_real_numeric_infant_outcome_phenotype_summary_statistics_n29_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M6 ===== ###

## calculate summary statistics
outcome_sumstats_ihmo_M6 <- create.summary.numeric(ihmo_num[ihmo_num$time_point=="6_months",c(64:69,71:84,87:89,91:96)])

## resort columns
outcome_sumstats_ihmo_M6 <- outcome_sumstats_ihmo_M6[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ihmo_M6$Note <- NA
# outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ihmo_M6[outcome_sumstats_ihmo_M6$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 64:96){
#   print(colnames(ihmo)[i])
#   print(table(ihmo[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ihmo_M6, file="phenotype_summary_statistics/240311_milk_HMOs_real_numeric_infant_outcome_phenotype_summary_statistics_n29_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 7.3 SUMMARY STATISTICS FOR INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA ===== ###

## calculate summary statistics
outcome_sumstats_ibact <- create.summary.numeric(ibact_num[,c(694:699,701:714,716:718,720:725)])

## resort columns
outcome_sumstats_ibact <- outcome_sumstats_ibact[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ibact$Note <- NA
# outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact[outcome_sumstats_ibact$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 694:725){
#   print(colnames(ibact)[i])
#   print(table(ibact[,i], useNA="ifany"))
# }


## save phenotype summary statistics table
# write.table(outcome_sumstats_ibact, file="phenotype_summary_statistics/240222_milk_microbiota_numeric_infant_outcome_phenotype_summary_statistics_n29.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M1 ===== ###

## calculate summary statistics
outcome_sumstats_ibact_M1 <- create.summary.numeric(ibact_num[ibact_num$time_point=="1_month",c(694:699,701:714,716:718,720:725)])

## resort columns
outcome_sumstats_ibact_M1 <- outcome_sumstats_ibact_M1[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ibact_M1$Note <- NA
# outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M1[outcome_sumstats_ibact_M1$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 694:725){
#   print(colnames(ibact)[i])
#   print(table(ibact[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ibact_M1, file="phenotype_summary_statistics/240311_milk_microbiota_numeric_infant_outcome_phenotype_summary_statistics_n29_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M2 ===== ###

## calculate summary statistics
outcome_sumstats_ibact_M2 <- create.summary.numeric(ibact_num[ibact_num$time_point=="2_months",c(694:699,701:714,716:718,720:725)])

## resort columns
outcome_sumstats_ibact_M2 <- outcome_sumstats_ibact_M2[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ibact_M2$Note <- NA
# outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M2[outcome_sumstats_ibact_M2$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 694:725){
#   print(colnames(ibact)[i])
#   print(table(ibact[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ibact_M2, file="phenotype_summary_statistics/240311_milk_microbiota_numeric_infant_outcome_phenotype_summary_statistics_n29_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M3 ===== ###

## calculate summary statistics
outcome_sumstats_ibact_M3 <- create.summary.numeric(ibact_num[ibact_num$time_point=="3_months",c(694:699,701:714,716:718,720:725)])

## resort columns
outcome_sumstats_ibact_M3 <- outcome_sumstats_ibact_M3[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ibact_M3$Note <- NA
# outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M3[outcome_sumstats_ibact_M3$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 694:725){
#   print(colnames(ibact)[i])
#   print(table(ibact[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ibact_M3, file="phenotype_summary_statistics/240311_milk_microbiota_numeric_infant_outcome_phenotype_summary_statistics_n29_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M6 ===== ###

## calculate summary statistics
outcome_sumstats_ibact_M6 <- create.summary.numeric(ibact_num[ibact_num$time_point=="6_months",c(694:699,701:714,716:718,720:725)])

## resort columns
outcome_sumstats_ibact_M6 <- outcome_sumstats_ibact_M6[,c(1,9:13,2:8)]

## add a note for categorical phenotypes that were transformed to numeric, indicte which categories the numbers correspond to
outcome_sumstats_ibact_M6$Note <- NA
# outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_cry_crying_more_than_3h_per_day","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_BITSS","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__hard_stools;2__soft_formed_stools;3__mushy_loose_stools;4__watery_stools"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_2_freq_difficult_stool_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_3_freq_spitting_milk_normal_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_5_freq_discomfort_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always;5__always"

outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_6_freq_pain_spitting_milk_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_to_3_times_per_week;4__4_to_6_times_per_week;5__7_or_more_times_per_week"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_7_crying_duration_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_10_minutes_per_day;2__10_to_30_minutes_per_day;3__30_to_60_minutes_per_day;4__1_to_2_hours_per_day;5__more_than_2_hours_per_day"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_8_freq_unstoppable_crying_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_9_freq_crying_after_feeding_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_11_freq_unstoppable_restlessness_baby_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_week;2__1_time_per_week;3__2_times_per_week;4__3_times_per_week;5__4_or_more_times_per_week"

outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_12_freq_gas_per_day_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__0_times_per_day;2__1_time_per_day;3__2_to_3_times_per_day;4__4_to_6_times_per_day;5__7_or_more_times_per_day"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_igsq_13_freq_discomfort_gas_past_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__almost_never;3__sometimes;4__almost_always"
# outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_dis_ROME_colic","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
# outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_dis_ROME_e1_regurgitation","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_dis_ROME_e13_stool_freq_last_month","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__less_than_once_a_week;2__one_to_two_times_a_week;3__between_3_to_7_times_a_week;4__once_a_day;5__two_to_three_times_a_day;6__more_than_three_times_a_day"

outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_dis_ROME_e14_stool_structure","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 2__hard;3__not_hard_or_soft;4__sometimes_hard_sometimes_soft;5__very_soft;6__watery"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_dis_ROME_e25_mucus_stool_last_week","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__never;2__once;3__sometimes;4__almost_always;5__always"
# outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_dis_ROME_e32_sudden_inconsolable_crying_fits","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 0__no;1__yes"
outcome_sumstats_ibact_M6[outcome_sumstats_ibact_M6$Phenotype=="infant_sleep_pattern_problem","Note"] <- "This phenotype was categorical and has been converted to numeric for association with biological data. The numbers correspond to these levels: 1__no_problem;2__small_problem;3__big_problem"

# for (i in 694:725){
#   print(colnames(ibact)[i])
#   print(table(ibact[,i], useNA="ifany"))
# }

## save phenotype summary statistics table
# write.table(outcome_sumstats_ibact_M6, file="phenotype_summary_statistics/240311_milk_microbiota_numeric_infant_outcome_phenotype_summary_statistics_n29_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)


# ### ===== 7.4 SUMMARY STATISTICS FOR MILK HMOs ===== ###
# 
# ## calculate summary statistics
# HMO_sumstats_ihmo <- create.summary.numeric(ihmo_num[,36:63])
# 
# ## resort columns
# HMO_sumstats_ihmo <- HMO_sumstats_ihmo[,c(1,9:13,2:8)]
# 
# ## save phenotype summary statistics table
# write.table(HMO_sumstats_ihmo, file="phenotype_summary_statistics/240222_milk_HMOs_summary_statistics_n28.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 7.5 SUMMARY STATISTICS FOR MILK ALPHA DIVERSITY ===== ###

## calculate summary statistics
alpha_sumstats_ibact <- create.summary.numeric(ibact_num[,67:69])

## resort columns
alpha_sumstats_ibact <- alpha_sumstats_ibact[,c(1,9:13,2:8)]

## save phenotype summary statistics table
# write.table(alpha_sumstats_ibact, file="phenotype_summary_statistics/240222_milk_alpha_diversity_summary_statistics_n3.txt", row.names=F, col.names=T, sep="\t", quote=F)


##### =========================== 8. COMBINE SUMMARY STATISTICS TABLES FOR ALL INFANT OUTCOME PHENOTYPES =========================== #####

### ===== 8.1 COMBINE SUMMARY STATISTICS TABLES WITH CATEGORICAL AND NUMERIC PHENOTYPES FOR INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs (TABLE S35A-D) ===== ###

## final selection of 11 phenotypes
my_sel_phenos <- c("infant_dis_ROME_e1_regurgitation","infant_dis_ROME_e32_sudden_inconsolable_crying_fits",
                        "infant_growth_weight_kg","infant_growth_weight_kg_gain_since_birth",
                        "infant_growth_length_cm","infant_growth_length_cm_gain_since_birth",
                        "infant_cry_prev_week_average_crying_time_per_day_min",
                        "infant_BITSS",
                        "infant_dis_ROME_e13_stool_freq_last_month","infant_dis_ROME_e14_stool_structure","infant_dis_ROME_e25_mucus_stool_last_week")


## add colummns with data_type info
factor_ihmo_sum$data_type <- "factor"
outcome_sumstats_ihmo$data_type <- "numeric"

## combine tables
library(plyr)
ihmo_sum_tmp <- rbind.fill(factor_ihmo_sum, outcome_sumstats_ihmo)

## resort columns
ihmo_sum <- ihmo_sum_tmp[,c(1,9,4:8,2:3,10:17)] #33x17

## save phenotype summary statistics table
# write.table(ihmo_sum, file="phenotype_summary_statistics/240311_milk_HMOs_real_all_infant_outcome_phenotype_summary_statistics_n33.txt", row.names=F, col.names=T, sep="\t", quote=F)


# ### ===== W2 ===== ###
# 
# ## add colummns with data_type info
# factor_ihmo_sum_W2$data_type <- "factor"
# outcome_sumstats_ihmo_W2$data_type <- "numeric"
# 
# ## combine tables
# library(plyr)
# ihmo_sum_tmp_W2 <- rbind.fill(factor_ihmo_sum_W2, outcome_sumstats_ihmo_W2)
# 
# ## add time point column
# ihmo_sum_tmp_W2$time_point <- "0.5_months"
# 
# ## resort columns
# ihmo_sum_W2 <- ihmo_sum_tmp_W2[,c(18,1,9,4:8,2:3,10:17)] #33x18
# 
# ## save phenotype summary statistics table
# # write.table(ihmo_sum_W2, file="phenotype_summary_statistics/240311_milk_HMOs_real_all_infant_outcome_phenotype_summary_statistics_n33_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M1 ===== ###

## add colummns with data_type info
factor_ihmo_sum_M1$data_type <- "factor"
outcome_sumstats_ihmo_M1$data_type <- "numeric"

## combine tables
library(plyr)
ihmo_sum_tmp_M1 <- rbind.fill(factor_ihmo_sum_M1, outcome_sumstats_ihmo_M1)

## add phenotype group column
ihmo_sum_tmp_M1$Phenotype_group <- NA
ihmo_sum_tmp_M1[grep("ROME", ihmo_sum_tmp_M1$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ihmo_sum_tmp_M1[grep("BITSS", ihmo_sum_tmp_M1$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ihmo_sum_tmp_M1[grep("growth", ihmo_sum_tmp_M1$Phenotype), "Phenotype_group"] <- "Infant_growth"
ihmo_sum_tmp_M1[grep("infant_cry_prev_week_average_crying_time_per_day_min", ihmo_sum_tmp_M1$Phenotype), "Phenotype_group"] <- "Infant_crying"

## add time point column
ihmo_sum_tmp_M1$time_point <- "1_month"

## resort columns
ihmo_sum_M1 <- ihmo_sum_tmp_M1[,c(18,1,9,4:8,2:3,10:17,19)] #33x19

## only include selected final phenotypes
ihmo_sum_M1_v2 <- ihmo_sum_M1[ihmo_sum_M1$Phenotype %in% my_sel_phenos,]
ihmo_sum_M1_v2 <- ihmo_sum_M1_v2[order(ihmo_sum_M1_v2$Phenotype_group, ihmo_sum_M1_v2$Phenotype),]

## save phenotype summary statistics table (TABLE S35A)
write.table(ihmo_sum_M1_v2, file="phenotype_summary_statistics/240321_milk_HMOs_real_final_infant_outcome_phenotype_summary_statistics_n11_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M2 ===== ###

## add colummns with data_type info
factor_ihmo_sum_M2$data_type <- "factor"
outcome_sumstats_ihmo_M2$data_type <- "numeric"

## combine tables
library(plyr)
ihmo_sum_tmp_M2 <- rbind.fill(factor_ihmo_sum_M2, outcome_sumstats_ihmo_M2)

## add phenotype group column
ihmo_sum_tmp_M2$Phenotype_group <- NA
ihmo_sum_tmp_M2[grep("ROME", ihmo_sum_tmp_M2$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ihmo_sum_tmp_M2[grep("BITSS", ihmo_sum_tmp_M2$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ihmo_sum_tmp_M2[grep("growth", ihmo_sum_tmp_M2$Phenotype), "Phenotype_group"] <- "Infant_growth"
ihmo_sum_tmp_M2[grep("infant_cry_prev_week_average_crying_time_per_day_min", ihmo_sum_tmp_M2$Phenotype), "Phenotype_group"] <- "Infant_crying"

## add time point column
ihmo_sum_tmp_M2$time_point <- "2_months"

## resort columns
ihmo_sum_M2 <- ihmo_sum_tmp_M2[,c(18,1,9,4:8,2:3,10:17,19)] #33x19

## only include selected final phenotypes
ihmo_sum_M2_v2 <- ihmo_sum_M2[ihmo_sum_M2$Phenotype %in% my_sel_phenos,]
ihmo_sum_M2_v2 <- ihmo_sum_M2_v2[order(ihmo_sum_M2_v2$Phenotype_group, ihmo_sum_M2_v2$Phenotype),]

## save phenotype summary statistics table (TABLE S35B)
write.table(ihmo_sum_M2_v2, file="phenotype_summary_statistics/240321_milk_HMOs_real_final_infant_outcome_phenotype_summary_statistics_n11_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M3 ===== ###

## add colummns with data_type info
factor_ihmo_sum_M3$data_type <- "factor"
outcome_sumstats_ihmo_M3$data_type <- "numeric"

## combine tables
library(plyr)
ihmo_sum_tmp_M3 <- rbind.fill(factor_ihmo_sum_M3, outcome_sumstats_ihmo_M3)

## add phenotype group column
ihmo_sum_tmp_M3$Phenotype_group <- NA
ihmo_sum_tmp_M3[grep("ROME", ihmo_sum_tmp_M3$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ihmo_sum_tmp_M3[grep("BITSS", ihmo_sum_tmp_M3$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ihmo_sum_tmp_M3[grep("growth", ihmo_sum_tmp_M3$Phenotype), "Phenotype_group"] <- "Infant_growth"
ihmo_sum_tmp_M3[grep("infant_cry_prev_week_average_crying_time_per_day_min", ihmo_sum_tmp_M3$Phenotype), "Phenotype_group"] <- "Infant_crying"

## add time point column
ihmo_sum_tmp_M3$time_point <- "3_months"

## resort columns
ihmo_sum_M3 <- ihmo_sum_tmp_M3[,c(18,1,9,4:8,2:3,10:17,19)] #33x19

## only include selected final phenotypes
ihmo_sum_M3_v2 <- ihmo_sum_M3[ihmo_sum_M3$Phenotype %in% my_sel_phenos,]
ihmo_sum_M3_v2 <- ihmo_sum_M3_v2[order(ihmo_sum_M3_v2$Phenotype_group, ihmo_sum_M3_v2$Phenotype),]

## save phenotype summary statistics table (TABLE S35C)
write.table(ihmo_sum_M3_v2, file="phenotype_summary_statistics/240321_milk_HMOs_real_final_infant_outcome_phenotype_summary_statistics_n11_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M6 ===== ###

## add colummns with data_type info
factor_ihmo_sum_M6$data_type <- "factor"
outcome_sumstats_ihmo_M6$data_type <- "numeric"

## combine tables
library(plyr)
ihmo_sum_tmp_M6 <- rbind.fill(factor_ihmo_sum_M6, outcome_sumstats_ihmo_M6)

## add phenotype group column
ihmo_sum_tmp_M6$Phenotype_group <- NA
ihmo_sum_tmp_M6[grep("ROME", ihmo_sum_tmp_M6$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ihmo_sum_tmp_M6[grep("BITSS", ihmo_sum_tmp_M6$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ihmo_sum_tmp_M6[grep("growth", ihmo_sum_tmp_M6$Phenotype), "Phenotype_group"] <- "Infant_growth"
ihmo_sum_tmp_M6[grep("infant_cry_prev_week_average_crying_time_per_day_min", ihmo_sum_tmp_M6$Phenotype), "Phenotype_group"] <- "Infant_crying"

## add time point column
ihmo_sum_tmp_M6$time_point <- "6_months"

## resort columns
ihmo_sum_M6 <- ihmo_sum_tmp_M6[,c(18,1,9,4:8,2:3,10:17,19)] #33x19

## only include selected final phenotypes
ihmo_sum_M6_v2 <- ihmo_sum_M6[ihmo_sum_M6$Phenotype %in% my_sel_phenos,]
ihmo_sum_M6_v2 <- ihmo_sum_M6_v2[order(ihmo_sum_M6_v2$Phenotype_group, ihmo_sum_M6_v2$Phenotype),]

## save phenotype summary statistics table (TABLE S35D)
write.table(ihmo_sum_M6_v2, file="phenotype_summary_statistics/240321_milk_HMOs_real_final_infant_outcome_phenotype_summary_statistics_n11_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 8.2 COMBINE SUMMARY STATISTICS TABLES WITH CATEGORICAL AND NUMERIC PHENOTYPES FOR INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA (TABLE S37A-C) ===== ###

## add colummns with data_type info
factor_ibact_sum$data_type <- "factor"
outcome_sumstats_ibact$data_type <- "numeric"

## combine tables
library(plyr)
ibact_sum_tmp <- rbind.fill(factor_ibact_sum, outcome_sumstats_ibact)

## resort columns
ibact_sum <- ibact_sum_tmp[,c(1,9,4:8,2:3,10:17)] #32x17

## save phenotype summary statistics table
# write.table(ibact_sum, file="phenotype_summary_statistics/240311_milk_microbiota_all_infant_outcome_phenotype_summary_statistics_n32.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M1 ===== ###

## add colummns with data_type info
factor_ibact_sum_M1$data_type <- "factor"
outcome_sumstats_ibact_M1$data_type <- "numeric"

## combine tables
library(plyr)
ibact_sum_tmp_M1 <- rbind.fill(factor_ibact_sum_M1, outcome_sumstats_ibact_M1)

## add phenotype group column
ibact_sum_tmp_M1$Phenotype_group <- NA
ibact_sum_tmp_M1[grep("ROME", ibact_sum_tmp_M1$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ibact_sum_tmp_M1[grep("BITSS", ibact_sum_tmp_M1$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ibact_sum_tmp_M1[grep("growth", ibact_sum_tmp_M1$Phenotype), "Phenotype_group"] <- "Infant_growth"
ibact_sum_tmp_M1[grep("infant_cry_prev_week_average_crying_time_per_day_min", ibact_sum_tmp_M1$Phenotype), "Phenotype_group"] <- "Infant_crying"

## add time point column
ibact_sum_tmp_M1$time_point <- "1_month"

## resort columns
ibact_sum_M1 <- ibact_sum_tmp_M1[,c(18,1,9,4:8,2:3,10:17,19)] #31x19

## only include selected final phenotypes
ibact_sum_M1_v2 <- ibact_sum_M1[ibact_sum_M1$Phenotype %in% my_sel_phenos,]
ibact_sum_M1_v2 <- ibact_sum_M1_v2[order(ibact_sum_M1_v2$Phenotype_group, ibact_sum_M1_v2$Phenotype),]

## save phenotype summary statistics table (TABLE S37A)
write.table(ibact_sum_M1_v2, file="phenotype_summary_statistics/240321_milk_microbiota_final_infant_outcome_phenotype_summary_statistics_n11_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M2 ===== ###

## add colummns with data_type info
factor_ibact_sum_M2$data_type <- "factor"
outcome_sumstats_ibact_M2$data_type <- "numeric"

## combine tables
library(plyr)
ibact_sum_tmp_M2 <- rbind.fill(factor_ibact_sum_M2, outcome_sumstats_ibact_M2)

## add phenotype group column
ibact_sum_tmp_M2$Phenotype_group <- NA
ibact_sum_tmp_M2[grep("ROME", ibact_sum_tmp_M2$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ibact_sum_tmp_M2[grep("BITSS", ibact_sum_tmp_M2$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ibact_sum_tmp_M2[grep("growth", ibact_sum_tmp_M2$Phenotype), "Phenotype_group"] <- "Infant_growth"
ibact_sum_tmp_M2[grep("infant_cry_prev_week_average_crying_time_per_day_min", ibact_sum_tmp_M2$Phenotype), "Phenotype_group"] <- "Infant_crying"

## add time point column
ibact_sum_tmp_M2$time_point <- "2_months"

## resort columns
ibact_sum_M2 <- ibact_sum_tmp_M2[,c(18,1,9,4:8,2:3,10:17,19)] #33x19

## only include selected final phenotypes
ibact_sum_M2_v2 <- ibact_sum_M2[ibact_sum_M2$Phenotype %in% my_sel_phenos,]
ibact_sum_M2_v2 <- ibact_sum_M2_v2[order(ibact_sum_M2_v2$Phenotype_group, ibact_sum_M2_v2$Phenotype),]

## save phenotype summary statistics table
write.table(ibact_sum_M2_v2, file="phenotype_summary_statistics/240321_milk_microbiota_final_infant_outcome_phenotype_summary_statistics_n11_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M3 ===== ###

## add colummns with data_type info
factor_ibact_sum_M3$data_type <- "factor"
outcome_sumstats_ibact_M3$data_type <- "numeric"

## combine tables
library(plyr)
ibact_sum_tmp_M3 <- rbind.fill(factor_ibact_sum_M3, outcome_sumstats_ibact_M3)

## add phenotype group column
ibact_sum_tmp_M3$Phenotype_group <- NA
ibact_sum_tmp_M3[grep("ROME", ibact_sum_tmp_M3$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ibact_sum_tmp_M3[grep("BITSS", ibact_sum_tmp_M3$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ibact_sum_tmp_M3[grep("growth", ibact_sum_tmp_M3$Phenotype), "Phenotype_group"] <- "Infant_growth"
ibact_sum_tmp_M3[grep("infant_cry_prev_week_average_crying_time_per_day_min", ibact_sum_tmp_M3$Phenotype), "Phenotype_group"] <- "Infant_crying"

## add time point column
ibact_sum_tmp_M3$time_point <- "3_months"

## resort columns
ibact_sum_M3 <- ibact_sum_tmp_M3[,c(18,1,9,4:8,2:3,10:17,19)] #33x19

## only include selected final phenotypes
ibact_sum_M3_v2 <- ibact_sum_M3[ibact_sum_M3$Phenotype %in% my_sel_phenos,]
ibact_sum_M3_v2 <- ibact_sum_M3_v2[order(ibact_sum_M3_v2$Phenotype_group, ibact_sum_M3_v2$Phenotype),]

## save phenotype summary statistics table (TABLE S37B)
write.table(ibact_sum_M3_v2, file="phenotype_summary_statistics/240321_milk_microbiota_final_infant_outcome_phenotype_summary_statistics_n11_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M6 ===== ###

## add colummns with data_type info
factor_ibact_sum_M6$data_type <- "factor"
outcome_sumstats_ibact_M6$data_type <- "numeric"

## combine tables
library(plyr)
ibact_sum_tmp_M6 <- rbind.fill(factor_ibact_sum_M6, outcome_sumstats_ibact_M6)

## add phenotype group column
ibact_sum_tmp_M6$Phenotype_group <- NA
ibact_sum_tmp_M6[grep("ROME", ibact_sum_tmp_M6$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ibact_sum_tmp_M6[grep("BITSS", ibact_sum_tmp_M6$Phenotype), "Phenotype_group"] <- "Infant_gastrointestinal_health"
ibact_sum_tmp_M6[grep("growth", ibact_sum_tmp_M6$Phenotype), "Phenotype_group"] <- "Infant_growth"
ibact_sum_tmp_M6[grep("infant_cry_prev_week_average_crying_time_per_day_min", ibact_sum_tmp_M6$Phenotype), "Phenotype_group"] <- "Infant_crying"

## add time point column
ibact_sum_tmp_M6$time_point <- "6_months"

## resort columns
ibact_sum_M6 <- ibact_sum_tmp_M6[,c(18,1,9,4:8,2:3,10:17,19)] #33x19

## only include selected final phenotypes
ibact_sum_M6_v2 <- ibact_sum_M6[ibact_sum_M6$Phenotype %in% my_sel_phenos,]
ibact_sum_M6_v2 <- ibact_sum_M6_v2[order(ibact_sum_M6_v2$Phenotype_group, ibact_sum_M6_v2$Phenotype),]

## save phenotype summary statistics table (TABLE S37C)
write.table(ibact_sum_M6_v2, file="phenotype_summary_statistics/240321_milk_microbiota_final_infant_outcome_phenotype_summary_statistics_n11_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)


#### =========================== 9. CALCULATE MILK RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE =========================== #####

### ===== 9.1 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE ===== ###

## save only columns with absolute bacterial abundances
rownames(ibact_num) <- ibact_num$seq_16S_sample_ID
milkbact <- ibact_num[,70:693] #737x624

## calculate relative bacterial abundances
milkbact_relab <- (milkbact/rowSums(milkbact))
table(rowSums(milkbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
milkbact_relab_tmp <- milkbact_relab[,-624]

## select only bacterial genera with ≥10% prevalence (for this load an earlier file with the milk bacteria with ≥10% prevalence to ensure we select the correct 36 bacteria)
milk_RAprev <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/231215_milk_RA_prev10_summary_statistics.txt", sep="\t", header=T, stringsAsFactors = T)
milk_RAprev_bacteria <- unique(milk_RAprev$bacterium)
milk_RAprev_bacteria <- paste0("g__", milk_RAprev_bacteria)

milkbact_relab_tmp2 <- milkbact_relab_tmp[,colnames(milkbact_relab_tmp) %in% milk_RAprev_bacteria]
length(milkbact_relab_tmp2) #36
# colnames(milkbact_relab_tmp2)

## merge back with metadata
ibact_num_RA <- merge(ibact_num[,c(1:69,694:725)], milkbact_relab_tmp2, by="row.names")
rownames(ibact_num_RA) <- ibact_num_RA$Row.names
ibact_num_RA <- ibact_num_RA[,-1]


### ===== 9.2 SUMMARY STATISTICS FOR MILK RELATIVE BACTERIAL ABUNDANCES ===== ###

## calculate summary statistics
RA_sumstats_ibact <- create.summary.numeric(ibact_num_RA[,102:137])

## resort columns
RA_sumstats_ibact <- RA_sumstats_ibact[,c(1,9:13,2:8)]

## save phenotype summary statistics table
write.table(RA_sumstats_ibact, file="phenotype_summary_statistics/240222_milk_relative_abundances_summary_statistics_n36.txt", row.names=F, col.names=T, sep="\t", quote=F)


##### =========================== 10. INVERSE-RANK TRANSFORM MILK HMOs AND INFANT FAECAL ALPHA DIVERSITY DATA =========================== #####

### ===== 10.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 10.2. INVERSE-RANK TRANSFORM MILK HMOs AND INFANT OUTCOME PHENOTYPES ===== ###

ihmo_num_invr <- ihmo_num
ihmo_num_invr[,c(25,36:63,64:69,71:84,87:89,91:96)] <- as.data.frame(apply(ihmo_num_invr[,c(25,36:63,64:69,71:84,87:89,91:96)], 2, invrank)) #also include mother_birth_gestational_age_weeks for invr transformation
colnames(ihmo_num_invr)[c(25,36:63,64:69,71:84,87:89,91:96)] <- paste0(colnames(ihmo_num_invr)[c(25,36:63,64:69,71:84,87:89,91:96)], "_invr")


### ===== 10.3. INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY MEASURES AND INFANT OUTCOME PHENOTYPES ===== ###

ibact_num_RA_invr <- ibact_num_RA
ibact_num_RA_invr[,c(29,67:69,70:75,77:90,92:94,96:101)] <- as.data.frame(apply(ibact_num_RA_invr[,c(29,67:69,70:75,77:90,92:94,96:101)], 2, invrank)) #also include mother_birth_gestational_age_weeks for invr transformation
colnames(ibact_num_RA_invr)[c(29,67:69,70:75,77:90,92:94,96:101)] <- paste0(colnames(ibact_num_RA_invr)[c(29,67:69,70:75,77:90,92:94,96:101)], "_invr")


##### =========================== 11. HISTOGRAMS FOR INVERSE-RANK-TRANSFORMED DATA =========================== #####

# library(ggplot2)

### ===== 11.1 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK HMOs ===== ###

## create histograms for all numeric infant outcome phenotypes
pdf(paste0("phenotype_data_distribution/240222_milk_HMOs_real_infant_outcome_phenotypes_data_distribution_histograms_2_after_invr_transformation_n29.pdf"))
for (i in c(64:69,71:84,87:89,91:96)){
  print(paste0("Creating histogram for ", colnames(ihmo_num_invr)[i]))
  histo <- ggplot(ihmo_num_invr[!is.na(ihmo_num_invr[,i]),], aes(x=ihmo_num_invr[!is.na(ihmo_num_invr[,i]),i]))+
    geom_histogram(bins=50)+
    labs(colnames(ihmo_num_invr)[i])+
    ggtitle(paste0(colnames(ihmo_num_invr)[i]))
  print(histo)
}
dev.off()


### ===== 11.2 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED INFANT OUTCOME PHENOTYPES IN DATA FRAME WITH MILK MICROBIOTA DATA ===== ###

## create histograms for all numeric infant outcome phenotypes
pdf(paste0("phenotype_data_distribution/240222_milk_microbiota_infant_outcome_phenotypes_data_distribution_histograms_2_after_invr_transformation_n29.pdf"))
for (i in c(70:75,77:90,92:94,96:101)){
  print(paste0("Creating histogram for ", colnames(ibact_num_RA_invr)[i]))
  histo <- ggplot(ibact_num_RA_invr[!is.na(ibact_num_RA_invr[,i]),], aes(x=ibact_num_RA_invr[!is.na(ibact_num_RA_invr[,i]),i]))+
    geom_histogram(bins=50)+
    labs(colnames(ibact_num_RA_invr)[i])+
    ggtitle(paste0(colnames(ibact_num_RA_invr)[i]))
  print(histo)
}
dev.off()


### ===== 11.3 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED MILK HMOs ===== ###

## create histograms for all numeric infant outcome phenotypes
pdf(paste0("HMO_data_distribution/240222_milk_HMOs_real_data_distribution_histograms_after_invr_transformation_n28.pdf"))
for (i in 36:63){
  print(paste0("Creating histogram for ", colnames(ihmo_num_invr)[i]))
  histo <- ggplot(ihmo_num_invr[!is.na(ihmo_num_invr[,i]),], aes(x=ihmo_num_invr[!is.na(ihmo_num_invr[,i]),i]))+
    geom_histogram(bins=50)+
    labs(colnames(ihmo_num_invr)[i])+
    ggtitle(paste0(colnames(ihmo_num_invr)[i]))
  print(histo)
}
dev.off()


### ===== 11.4 HISTOGRAMS FOR INVERSE-RANK TRANSFORMED MILK ALPHA DIVERSITY MEASURES ===== ###

## create histograms for all numeric infant outcome phenotypes
pdf(paste0("microbiota_data_distribution/240222_milk_alpha_diversity_data_distribution_histograms_after_invr_transformation_n3.pdf"))
for (i in 67:69){
  print(paste0("Creating histogram for ", colnames(ibact_num_RA_invr)[i]))
  histo <- ggplot(ibact_num_RA_invr[!is.na(ibact_num_RA_invr[,i]),], aes(x=ibact_num_RA_invr[!is.na(ibact_num_RA_invr[,i]),i]))+
    geom_histogram(bins=50)+
    labs(colnames(ibact_num_RA_invr)[i])+
    ggtitle(paste0(colnames(ibact_num_RA_invr)[i]))
  print(histo)
}
dev.off()


##### =========================== 12. CLR-TRANSFORM MILK RELATIVE ABUNDANCES =========================== #####

### ===== 12.1 FUNCTION FOR CLR-TRANSFORMATION ===== ###

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


### ===== 12.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES ===== ###

ibact_num_RA_invr_clr <- ibact_num_RA_invr
ibact_num_RA_invr_clr[,102:137] <- as.data.frame(do_clr_externalWeighting(ibact_num_RA_invr_clr[,102:137], ibact_num_RA_invr_clr[,102:137]))
colnames(ibact_num_RA_invr_clr)[102:137] <- c(paste0(colnames(ibact_num_RA_invr_clr)[102:137], "_clr"))


##### =========================== 13. HISTOGRAMS FOR CLR-TRANSFORMED DATA =========================== #####

# library(ggplot2)

### ===== 13.1 HISTOGRAMS FOR CLR-TRANSFORMED MILK RELATIVE BACTERIAL ABUNDANCES ===== ###

## create histograms for all numeric infant outcome phenotypes
pdf(paste0("microbiota_data_distribution/240222_milk_relative_abundances_data_distribution_histograms_2_after_invr_transformation_n36.pdf"))
for (i in 102:137){
  print(paste0("Creating histogram for ", colnames(ibact_num_RA_invr_clr)[i]))
  histo <- ggplot(ibact_num_RA_invr_clr[!is.na(ibact_num_RA_invr_clr[,i]),], aes(x=ibact_num_RA_invr_clr[!is.na(ibact_num_RA_invr_clr[,i]),i]))+
    geom_histogram(bins=50)+
    labs(colnames(ibact_num_RA_invr_clr)[i])+
    ggtitle(paste0(colnames(ibact_num_RA_invr_clr)[i]))
  print(histo)
}
dev.off()


##### =========================== 14. MODEL FUNCTION (LMER/LM) =========================== #####

## function to run lmer or lm models with/without correction for covariates
run.models <- function(datasetname, modeltype, covardetails, covariates, inputdata, xcolumn, ycolumn, outputfolder, timepoint){
  
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
      
      print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      if (modeltype == "lmer"){
        
        ## RUN MIXED MODEL
        paste0("Run mixed model")
        
        my.lmer.model <- paste0("my_df[,j] ~ ", covariates, "my_df[,i]")
        print(paste0("Model function: ", my.lmer.model))
        
        m1 <- lmerTest::lmer(as.formula(my.lmer.model), REML=F, data=my_df)
        sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
        sum1$rows <- rownames(sum1)
        print(sum1)
        
        #retrieve estimates and p values from model with phenotype of interest (from sum1)
        levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        stat <- c(stat, rep(paste0(covardetails), length(sum1[grep("my_df",sum1$rows),"rows"])))
        d <- c(d, rep(datasetname, length(sum1[grep("my_df",sum1$rows),"rows"])))
        x <- c(x, rep(colnames(inputdata)[i], length(sum1[grep("my_df",sum1$rows),"rows"])))
        y <- c(y, rep(colnames(inputdata)[j], length(sum1[grep("my_df",sum1$rows),"rows"])))
        n <- c(n, rep(nrow(my_df), length(sum1[grep("my_df",sum1$rows),"rows"])))
      }else{
        
        ##RUN LINEAR MODEL
        paste0("Run linear model")
        
        my.lm.model <- paste0("my_df[,j] ~ ", covariates, "my_df[,i]")
        print(paste0("Model function: ", my.lm.model))
        
        m1 <- lm(as.formula(my.lm.model), data=my_df)
        sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
        sum1$rows <- rownames(sum1)
        print(sum1)
        
        #retrieve estimates and p values from model with phenotype of interest (from sum1)
        levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        stat <- c(stat, rep(paste0(covardetails), length(sum1[grep("my_df",sum1$rows),"rows"])))
        d <- c(d, rep(datasetname, length(sum1[grep("my_df",sum1$rows),"rows"])))
        x <- c(x, rep(colnames(inputdata)[i], length(sum1[grep("my_df",sum1$rows),"rows"])))
        y <- c(y, rep(colnames(inputdata)[j], length(sum1[grep("my_df",sum1$rows),"rows"])))
        n <- c(n, rep(nrow(my_df), length(sum1[grep("my_df",sum1$rows),"rows"])))
      }
      
      ##PLOT DATA
      
      if (class(my_df[,i]) == "numeric" | class(my_df[,i]) == "integer"){
        print("Generating scatterplot")
        scatterplot1 <- ggplot(my_df, aes(x=my_df[,i], y=my_df[,j]))+
          geom_point(alpha=0.5)+
          theme_bw()+
          facet_grid(.~time_point_numeric)+
          geom_smooth(method='lm')+
          labs(x=paste0(colnames(my_df)[i]), y=paste0(colnames(my_df)[j]),
               title=paste0(colnames(my_df)[j], " by ", colnames(my_df)[i], "\n",
                            "n=", nrow(my_df), "\n", 
                            covardetails, "\n",
                            "p=", signif(p[length(p)],3)))
        
        if (p[length(p)]<0.05){
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint, "/nominally_significant/240321_infant_outcomes_", covardetails, "_", colnames(my_df)[j] ,"_by_", colnames(my_df)[i], "_and_time", "_", timepoint, ".pdf"), scatterplot1, dpi=500, width=30, height=15, units="cm")
        }else{
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint,  "/not_significant/240321_infant_outcomes_", covardetails, "_", colnames(my_df)[j] ,"_by_", colnames(my_df)[i], "_and_time", "_", timepoint, ".pdf"), scatterplot1, dpi=500, width=30, height=15, units="cm")
        }
        
      }else{
        print("Generating boxplots")
        boxplot1 <- ggplot(my_df, aes(x=my_df[,i], y=my_df[,j]))+
          geom_jitter(alpha=0.5)+
          geom_boxplot(alpha=0)+
          theme_bw()+
          facet_grid(.~time_point_numeric)+
          labs(x=paste0(colnames(my_df)[i]), y=paste0(colnames(my_df)[j]),
               title=paste0(colnames(my_df)[j], " by ", colnames(my_df)[i], "\n",
                            "n=", nrow(my_df), "\n", 
                            covardetails, "\n",
                            "p=", signif(p[length(p)],3)))
        
        if (p[length(p)]<0.05){
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint,  "/nominally_significant/240321_infant_outcomes_", covardetails, "_", colnames(my_df)[j] ,"_by_", colnames(my_df)[i], "_and_time", "_", timepoint, ".pdf"), boxplot1, dpi=500, width=30, height=15, units="cm")
        }else{
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint,  "/not_significant/240321_infant_outcomes_", covardetails, "_", colnames(my_df)[j] ,"_by_", colnames(my_df)[i], "_and_time", "_", timepoint, ".pdf"), boxplot1, dpi=500, width=30, height=15, units="cm")
        }
      }
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, estimate=e, p=p)
  
  #remove "my_df[, i]" from levels
  res$levels <- substring(res$levels, 11)
  
  #correct for multiple testing
  res$FDR <- p.adjust(res$p, method="BH")
  
  #resort results by FDR and p
  res <- res[order(res$FDR, res$p),]
  
  #save output table
  write.table(res, file=paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint, "/240321_", datasetname, "_", outputfolder, "_", covardetails, "_", timepoint, "_results.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  
  #return output table
  return(res)
}


##### =========================== 15. DETAILS ABOUT ihmo_num_invr DATA FRAME FOR ASSOCIATION OF BASIC PHENOTYPES AND MILK HMOs WITH INFANT OUTCOMES =========================== #####

## Infant data set (ihmo_num_invr)
## 1563x96
## Columns:
## 1:15 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 16:35 = additional maternal/infant phenotypes (check effect for these: c(18:20,22,24:27,29,34))
## 36:63 = invr-transformed milk HMOs
## 64:96 = invr-transformed infant outcome phenotypes


## basic phenotypes (will be used as x in the association models)

for (i in c(18,25:27,29,34)){
  print(paste0(i, "__", colnames(ihmo_num_invr)[i]))
  print(table(ihmo_num_invr[!is.na(ihmo_num_invr[,i]),"time_point"]))
}

## (1) numeric basic phenotype
## -> use numeric basic phenotype as y and numeric/categorical infant outcomes as x
##    lm() model, run per time point
# static, multiple time points, [25] "mother_birth_gestational_age_weeks_invr"  

## (2) categorical basic phenotype
## a) if the infant outcome is numeric:
##    -> use basic phenotype as x and infant outcome as y
##       lm() model, run per time point
## b) if the infant outcome is a categorical, i.e. both x and y are categorical
##    -> create contingency table
##       Fisher's exact test per time point
# [static, multiple time points, [18] "mother_milk_HMO_milk_group"]
# static, multiple time points, [26] "infant_sex"                                                      
# static, multiple time points, [27] "infant_birth_delivery_place"                                     
# static, multiple time points, [29] "infant_birth_delivery_mode"                                      
# dynamic, multiple time points, [34] "infant_ffq_feeding_mode"  


##### =========================== 16. CHECK EFFECTS OF BASIC PHENOTYPES ON INFANT OUTCOMES =========================== #####

### ===== 16.1 RUN MODELS FOR NUMERIC/CATEGORICAL INFANT OUTCOMES ~ NUMERIC BASIC PHENOTYPE ===== ###

## (1) numeric basic phenotype
## -> use numeric basic phenotype as y and numeric/categorical infant outcomes as x
##    lm() model, run per time point
# static, multiple time points, [25] "mother_birth_gestational_age_weeks_invr"  

inff_outcomes_by_GA_results_W2 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_outcomes_by_GA_results_W2",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months",],
                                             xcolumn     = c(64:72,91:92),  #infant outcomes for W2
                                             ycolumn     = c(25),  #basic phenotypes: mother_birth_gestational_age_weeks_invr
                                             outputfolder = "model_results_basic_phenotypes",
                                             timepoint = "W2")

inff_outcomes_by_GA_results_M1 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                          modeltype = "lm",
                                          covardetails = "lm_without_correction_inff_outcomes_by_GA_results_M1",
                                          covariates = "",
                                          inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                          xcolumn     = c(64:72,85:92),  #infant outcomes for M1
                                          ycolumn     = c(25),  #basic phenotypes: mother_birth_gestational_age_weeks_invr
                                          outputfolder = "model_results_basic_phenotypes",
                                          timepoint = "M1")

inff_outcomes_by_GA_results_M2 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_outcomes_by_GA_results_M2",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months",],
                                             xcolumn     = c(64:72,85:92),  #infant outcomes for M2
                                             ycolumn     = c(25),  #basic phenotypes: mother_birth_gestational_age_weeks_invr
                                             outputfolder = "model_results_basic_phenotypes",
                                             timepoint = "M2")

inff_outcomes_by_GA_results_M3 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_outcomes_by_GA_results_M3",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months",],
                                             xcolumn     = c(64:96),  #infant outcomes for M3
                                             ycolumn     = c(25),  #basic phenotypes: mother_birth_gestational_age_weeks_invr
                                             outputfolder = "model_results_basic_phenotypes",
                                             timepoint = "M3")

inff_outcomes_by_GA_results_M6 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_outcomes_by_GA_results_M6",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months",],
                                             xcolumn     = c(64:72,86:90,92:96),  #infant outcomes for M6
                                             ycolumn     = c(25),  #basic phenotypes: mother_birth_gestational_age_weeks_invr
                                             outputfolder = "model_results_basic_phenotypes",
                                             timepoint = "M6")

 
### ===== 16.2 RUN MODELS FOR NUMERIC INFANT OUTCOMES ~ CATEGORICAL BASIC PHENOTYPE ===== ###

## (2) categorical basic phenotype
## a) if the infant outcome is numeric:
##    -> use basic phenotype as x and infant outcome as y
##       lm() model, run per time point
# [static, multiple time points, [18] "mother_milk_HMO_milk_group"]
# static, multiple time points, [26] "infant_sex"                                                      
# static, multiple time points, [27] "infant_birth_delivery_place"                                     
# static, multiple time points, [29] "infant_birth_delivery_mode"                                      
# dynamic, multiple time points, [34] "infant_ffq_feeding_mode"  

inff_numoutcomes_by_catphenos_results_W2 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                                       modeltype = "lm",
                                                       covardetails = "lm_without_correction_inff_numoutcomes_by_catphenos_results_W2",
                                                       covariates = "",
                                                       inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months",],
                                                       xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                                       ycolumn     = c(64:69,71:72,91:92),  #numeric infant outcomes for W2
                                                       outputfolder = "model_results_basic_phenotypes",
                                                       timepoint = "W2")

inff_numoutcomes_by_catphenos_results_M1 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_numoutcomes_by_catphenos_results_M1",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                             xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                             ycolumn     = c(64:69,71:72,87:89,91:92),  #numeric infant outcomes for M1
                                             outputfolder = "model_results_basic_phenotypes",
                                             timepoint = "M1")

inff_numoutcomes_by_catphenos_results_M2 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_numoutcomes_by_catphenos_results_M2",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months",],
                                             xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                             ycolumn     = c(64:69,71:72,87:89,91:92),  #numeric infant outcomes for M2
                                             outputfolder = "model_results_basic_phenotypes",
                                             timepoint = "M2")

inff_numoutcomes_by_catphenos_results_M3 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_numoutcomes_by_catphenos_results_M3",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months",],
                                             xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                             ycolumn     = c(64:69,71:84,87:89,91:96),  #infant outcomes for M3
                                             outputfolder = "model_results_basic_phenotypes",
                                             timepoint = "M3")

inff_numoutcomes_by_catphenos_results_M6 <- run.models(datasetname = "infant_outcomes_by_basic_phenotypes",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_numoutcomes_by_catphenos_results_M6",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months",],
                                             xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                             ycolumn     = c(64:69,71:72,87:89,92:96),  #infant outcomes for M6
                                             outputfolder = "model_results_basic_phenotypes",
                                             timepoint = "M6")


### ===== 16.3 RUN FISHER'S TESTS FOR CATEGORICAL BASIC PHENOTYPE ~ CATEGORICAL INFANT OUTCOMES ===== ###

## (2) categorical basic phenotype
## b) if the infant outcome is a categorical, i.e. both x and y are categorical
##    -> create contingency table
##       Fisher's exact test per time point
# [static, multiple time points, [18] "mother_milk_HMO_milk_group]
# static, multiple time points, [26] "infant_sex"                                                      
# static, multiple time points, [27] "infant_birth_delivery_place"                                     
# static, multiple time points, [29] "infant_birth_delivery_mode"                                      
# dynamic, multiple time points, [34] "infant_ffq_feeding_mode"  

## create matrix

run.fishers.test <- function(datasetname, covardetails, inputdata, xcolumn, ycolumn, outputfolder, timepoint){
  
  x <- c()
  y <- c()
  n <- c()
  e <- c()
  p <- c()
  
  for (i in xcolumn){
    for (j in ycolumn){
      
      ## select only rows without NAs
      print("Creating data frame")
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      x <- c(x, colnames(my_df)[i])
      y <- c(y, colnames(my_df)[j])
      n <- c(n, nrow(my_df))
      
      ## create matrix
      print("Creating matrix")
      my_matrix <- matrix(table(my_df[,i], my_df[,j]), nrow=2, dimnames = list(levels(my_df[,i]), levels(my_df[,j])))
      
      ## run Fisher's test
      print("Running Fisher's test")
      fisher.results <- fisher.test(my_matrix)
      p <- c(p, fisher.results$p.value)
      e <- c(e, fisher.results$estimate)
      
    }
  }
  
  #save results in data frame
  print("Creating results dataframe")
  res <- data.frame(dataset=rep(datasetname, length(p)), statistic=rep("Fishers_test", length(p)), x=x, y=y, n_total=n, estimate=e, p=p)
  
  #correct for multiple testing
  res$FDR <- p.adjust(res$p, method="BH")
  
  #resort results by FDR and p
  res <- res[order(res$FDR, res$p),]
  
  #save output table
  write.table(res, file=paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint, "/240321_", datasetname, "_", outputfolder, "_", covardetails, "_", timepoint, "_results.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  
  #return output table
  return(res)
}


inff_catoutcomes_by_catphenos_results_W2 <- run.fishers.test(datasetname = "infant_outcomes_by_basic_phenotypes",
                                                             covardetails = "Fishers_test_inff_catoutcomes_by_catphenos_results_W2",
                                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months",],
                                                             xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                                             ycolumn     = c(70),  #categorical infant outcomes for W2
                                                             outputfolder = "model_results_basic_phenotypes",
                                                             timepoint = "W2")

inff_catoutcomes_by_catphenos_results_M1 <- run.fishers.test(datasetname = "infant_outcomes_by_basic_phenotypes",
                                                             covardetails = "Fishers_test_inff_catoutcomes_by_catphenos_results_M1",
                                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                                             xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                                             ycolumn     = c(70,85:86,90),  #categorical infant outcomes for M1
                                                             outputfolder = "model_results_basic_phenotypes",
                                                             timepoint = "M1")

inff_catoutcomes_by_catphenos_results_M2 <- run.fishers.test(datasetname = "infant_outcomes_by_basic_phenotypes",
                                                            covardetails = "Fishers_test_inff_catoutcomes_by_catphenos_results_M2",
                                                            inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months",],
                                                            xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                                            ycolumn     = c(70,85:86,90),  #categorical infant outcomes for M2
                                                            outputfolder = "model_results_basic_phenotypes",
                                                            timepoint = "M2")

inff_catoutcomes_by_catphenos_results_M3 <- run.fishers.test(datasetname = "infant_outcomes_by_basic_phenotypes",
                                                             covardetails = "Fishers_test_inff_catoutcomes_by_catphenos_results_M3",
                                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months",],
                                                             xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                                             ycolumn     = c(70,85:86,90),  #categorical infant outcomes for M3
                                                             outputfolder = "model_results_basic_phenotypes",
                                                             timepoint = "M3")

inff_catoutcomes_by_catphenos_results_M6 <- run.fishers.test(datasetname = "infant_outcomes_by_basic_phenotypes",
                                                             covardetails = "Fishers_test_inff_catoutcomes_by_catphenos_results_M6",
                                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months",],
                                                             xcolumn     = c(26,27,29,34),  #categorical basic phenotypes
                                                             ycolumn     = c(70,85:86,90),  #categorical infant outcomes for M6
                                                             outputfolder = "model_results_basic_phenotypes",
                                                             timepoint = "M6")


### ===== 16.4 COMBINE RESULTS PER TIME POINT AND CORRECT FOR MULTIPLE TESTING ===== ###

## add empty levels columns
inff_catoutcomes_by_catphenos_results_W2$levels <- NA
inff_catoutcomes_by_catphenos_results_M1$levels <- NA
inff_catoutcomes_by_catphenos_results_M2$levels <- NA
inff_catoutcomes_by_catphenos_results_M3$levels <- NA
inff_catoutcomes_by_catphenos_results_M6$levels <- NA


### ===== W2 ===== ###

## combine results for W2, 55x9
basic_phenotype_results_W2 <- as.data.frame(rbind(inff_catoutcomes_by_catphenos_results_W2, inff_numoutcomes_by_catphenos_results_W2, inff_outcomes_by_GA_results_W2))

## re-correct for multiple testing
basic_phenotype_results_W2$FDR <- p.adjust(basic_phenotype_results_W2$p, method="BH")

## re-sort results by FDR and p, and resort columns
basic_phenotype_results_W2 <- basic_phenotype_results_W2[order(basic_phenotype_results_W2$FDR, basic_phenotype_results_W2$p), c(1:4,9,5:8)]

## check results
nrow(basic_phenotype_results_W2[basic_phenotype_results_W2$FDR<0.05,]) #4 FDR significant associations (all associations of GA with growth parameters)
nrow(basic_phenotype_results_W2[basic_phenotype_results_W2$p<0.05,])   #9 nominally significant associations

## save combined results table
write.table(basic_phenotype_results_W2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_basic_phenotypes/W2/240222_infant_outcomes_by_basic_phenotypes_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M1 ===== ###

## combine results for M1, 85x9
basic_phenotype_results_M1 <- as.data.frame(rbind(inff_catoutcomes_by_catphenos_results_M1, inff_numoutcomes_by_catphenos_results_M1, inff_outcomes_by_GA_results_M1))

## re-correct for multiple testing
basic_phenotype_results_M1$FDR <- p.adjust(basic_phenotype_results_M1$p, method="BH")

## re-sort results by FDR and p, and resort columns
basic_phenotype_results_M1 <- basic_phenotype_results_M1[order(basic_phenotype_results_M1$FDR, basic_phenotype_results_M1$p), c(1:4,9,5:8)]

## check results
nrow(basic_phenotype_results_M1[basic_phenotype_results_M1$FDR<0.05,]) #10 FDR significant associations (mostly associations of GA and sex with growth parameters)
nrow(basic_phenotype_results_M1[basic_phenotype_results_M1$p<0.05,])   #14 nominally significant associations

## save combined results table
write.table(basic_phenotype_results_M1, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_basic_phenotypes/M1/240222_infant_outcomes_by_basic_phenotypes_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M2 ===== ###

## combine results for M2, 85x9
basic_phenotype_results_M2 <- as.data.frame(rbind(inff_catoutcomes_by_catphenos_results_M2, inff_numoutcomes_by_catphenos_results_M2, inff_outcomes_by_GA_results_M2))

## re-correct for multiple testing
basic_phenotype_results_M2$FDR <- p.adjust(basic_phenotype_results_M2$p, method="BH")

## re-sort results by FDR and p, and resort columns
basic_phenotype_results_M2 <- basic_phenotype_results_M2[order(basic_phenotype_results_M2$FDR, basic_phenotype_results_M2$p), c(1:4,9,5:8)]

## check results
nrow(basic_phenotype_results_M2[basic_phenotype_results_M2$FDR<0.05,]) #15 FDR significant associations (mostly associations of GA and sex with growth parameters)
nrow(basic_phenotype_results_M2[basic_phenotype_results_M2$p<0.05,])   #20 nominally significant associations

## save combined results table
write.table(basic_phenotype_results_M2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_basic_phenotypes/M2/240222_infant_outcomes_by_basic_phenotypes_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M3 ===== ###

## combine results for M3, 85x9
basic_phenotype_results_M3 <- as.data.frame(rbind(inff_catoutcomes_by_catphenos_results_M3, inff_numoutcomes_by_catphenos_results_M3, inff_outcomes_by_GA_results_M3))

## re-correct for multiple testing
basic_phenotype_results_M3$FDR <- p.adjust(basic_phenotype_results_M3$p, method="BH")

## re-sort results by FDR and p, and resort columns
basic_phenotype_results_M3 <- basic_phenotype_results_M3[order(basic_phenotype_results_M3$FDR, basic_phenotype_results_M3$p), c(1:4,9,5:8)]

## check results
nrow(basic_phenotype_results_M3[basic_phenotype_results_M3$FDR<0.05,]) #11 FDR significant associations (mostly associations of GA and sex with growth parameters)
nrow(basic_phenotype_results_M3[basic_phenotype_results_M3$p<0.05,])   #20 nominally significant associations

## save combined results table
write.table(basic_phenotype_results_M3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_basic_phenotypes/M3/240222_infant_outcomes_by_basic_phenotypes_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M6 ===== ###

## combine results for M6, 85x9
basic_phenotype_results_M6 <- as.data.frame(rbind(inff_catoutcomes_by_catphenos_results_M6, inff_numoutcomes_by_catphenos_results_M6, inff_outcomes_by_GA_results_M6))

## re-correct for multiple testing
basic_phenotype_results_M6$FDR <- p.adjust(basic_phenotype_results_M6$p, method="BH")

## re-sort results by FDR and p, and resort columns
basic_phenotype_results_M6 <- basic_phenotype_results_M6[order(basic_phenotype_results_M6$FDR, basic_phenotype_results_M6$p), c(1:4,9,5:8)]

## check results
nrow(basic_phenotype_results_M6[basic_phenotype_results_M6$FDR<0.05,]) # 5 FDR significant associations (mostly associations of GA and sex with growth parameters)
nrow(basic_phenotype_results_M6[basic_phenotype_results_M6$p<0.05,])   #10 nominally significant associations

## save combined results table
write.table(basic_phenotype_results_M6, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_basic_phenotypes/M6/240222_infant_outcomes_by_basic_phenotypes_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 16.5 CHECK WHICH RESULTS ARE CONSISTENT ACROSS TIME POINTS ===== ###

W2_nomsign_associations <- c(paste0(basic_phenotype_results_W2[basic_phenotype_results_W2$p<0.05,"x"], "_", basic_phenotype_results_W2[basic_phenotype_results_W2$p<0.05,"y"]))
M1_nomsign_associations <- c(paste0(basic_phenotype_results_M1[basic_phenotype_results_M1$p<0.05,"x"], "_", basic_phenotype_results_M1[basic_phenotype_results_M1$p<0.05,"y"]))
M2_nomsign_associations <- c(paste0(basic_phenotype_results_M2[basic_phenotype_results_M2$p<0.05,"x"], "_", basic_phenotype_results_M2[basic_phenotype_results_M2$p<0.05,"y"]))
M3_nomsign_associations <- c(paste0(basic_phenotype_results_M3[basic_phenotype_results_M3$p<0.05,"x"], "_", basic_phenotype_results_M3[basic_phenotype_results_M3$p<0.05,"y"]))
M6_nomsign_associations <- c(paste0(basic_phenotype_results_M6[basic_phenotype_results_M6$p<0.05,"x"], "_", basic_phenotype_results_M6[basic_phenotype_results_M6$p<0.05,"y"]))


W2_nomsign_associations
M1_nomsign_associations
M2_nomsign_associations
M3_nomsign_associations
M6_nomsign_associations

M1_M3_nomsign_associations <- intersect(M1_nomsign_associations, M3_nomsign_associations)
M1_M3_nomsign_associations

M1_M2_M3_nomsign_associations <- intersect(M1_M3_nomsign_associations, M2_nomsign_associations)
M1_M2_M3_nomsign_associations

W2_M1_M2_M3_nomsign_associations <- intersect(M1_M2_M3_nomsign_associations, W2_nomsign_associations)
W2_M1_M2_M3_nomsign_associations

W2_M1_M2_M3_M6_nomsign_associations <- intersect(W2_M1_M2_M3_nomsign_associations, M6_nomsign_associations)
W2_M1_M2_M3_M6_nomsign_associations
# [1] "infant_growth_length_cm_invr_mother_birth_gestational_age_weeks_invr"                             
# [2] "infant_growth_head_circumference_cm_gain_since_birth_invr_mother_birth_gestational_age_weeks_invr"
# [3] "infant_sex_infant_growth_head_circumference_cm_invr"                                              
# [4] "infant_sex_infant_growth_length_cm_invr"   

## -> When investigating the effect of milk components on infant growth, correct for infant sex and gestational age in the models!


##### =========================== 17. ASSOCIATIONS OF MILK HMO CONCENTRATIONS WITH INFANT OUTCOME PHENOTYPES =========================== #####

### ===== 17.1 RUN MODELS FOR NUMERIC/CATEGORICAL INFANT OUTCOMES ~ NUMERIC MILK HMO CONCENTRATIONS ===== ###

## -> use numeric milk HMO concentrations as y and numeric/categorical infant outcomes as x
##    lm() model, run per time point
##    no correction factors, except for infant growth outcomes -> for growth, correct for infant sex and gestational age (invr)
##    for each time point: select which infant outcomes to include (have data at that time point)

inff_outcomes_by_HMOs_results_W2 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                               modeltype = "lm",
                                               covardetails = "lm_without_correction_inff_outcomes_by_milk_HMOs_results_W2",
                                               covariates = "",
                                               inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months",],
                                               xcolumn     = c(70:72,91:92),  #infant outcomes for W2 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69)
                                               ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                               outputfolder = "model_results_milk_HMOs",
                                               timepoint = "W2")

inff_outcomes_by_HMOs_results_M1 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_outcomes_by_milk_HMOs_results_M1",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                             xcolumn     = c(70:72,85:92),  #infant outcomes for M1 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69)
                                             ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                             outputfolder = "model_results_milk_HMOs",
                                             timepoint = "M1")

inff_outcomes_by_HMOs_results_M2 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_outcomes_by_milk_HMOs_results_M2",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months",],
                                             xcolumn     = c(70:72,85:92),  #infant outcomes for M2 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69)
                                             ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                             outputfolder = "model_results_milk_HMOs",
                                             timepoint = "M2")

inff_outcomes_by_HMOs_results_M3 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_outcomes_by_milk_HMOs_results_M3",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months",],
                                             xcolumn     = c(70:96),  #infant outcomes for M3 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69)
                                             ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                             outputfolder = "model_results_milk_HMOs",
                                             timepoint = "M3")

inff_outcomes_by_HMOs_results_M6 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                             modeltype = "lm",
                                             covardetails = "lm_without_correction_inff_outcomes_by_milk_HMOs_results_M6",
                                             covariates = "",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months",],
                                             xcolumn     = c(70:72,86:90,92:96),  #infant outcomes for M6 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69)
                                             ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                             outputfolder = "model_results_milk_HMOs",
                                             timepoint = "M6")

## if necessary, re-import results
inff_outcomes_by_HMOs_results_W2 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/W2/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_without_correction_inff_outcomes_by_milk_HMOs_results_W2_W2_results.txt", header=T, sep="\t", stringsAsFactors=T)
inff_outcomes_by_HMOs_results_M1 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M1/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_without_correction_inff_outcomes_by_milk_HMOs_results_M1_M1_results.txt", header=T, sep="\t", stringsAsFactors=T)
inff_outcomes_by_HMOs_results_M2 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M2/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_without_correction_inff_outcomes_by_milk_HMOs_results_M2_M2_results.txt", header=T, sep="\t", stringsAsFactors=T)
inff_outcomes_by_HMOs_results_M3 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M3/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_without_correction_inff_outcomes_by_milk_HMOs_results_M3_M3_results.txt", header=T, sep="\t", stringsAsFactors=T)
inff_outcomes_by_HMOs_results_M6 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M6/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_without_correction_inff_outcomes_by_milk_HMOs_results_M6_M6_results.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 17.2 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ NUMERIC MILK HMO CONCENTRATIONS ===== ###

## -> use numeric milk HMO concentrations as y and numeric/categorical infant outcomes as x
##    lm() model, run per time point
##    no correction factors, except for infant growth outcomes -> for growth, correct for infant sex and gestational age (invr)
##    for each time point: select which infant outcomes to include (have data at that time point)

inff_growth_by_HMOs_results_W2 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                               modeltype = "lm",
                                               covardetails = "lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_W2",
                                               covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                               inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                               xcolumn     = c(64:69),  #numeric infant growth outcomes
                                               ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                               outputfolder = "model_results_milk_HMOs",
                                               timepoint = "W2")

inff_growth_by_HMOs_results_M1 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                               modeltype = "lm",
                                               covardetails = "lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M1",
                                               covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                               inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                               xcolumn     = c(64:69),  #numeric infant growth outcomes
                                               ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                               outputfolder = "model_results_milk_HMOs",
                                               timepoint = "M1")
# inff_growth_by_HMOs_results_M1 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M1/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M1_M1_results.txt", header=T, sep="\t", stringsAsFactors=T)

inff_growth_by_HMOs_results_M2 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                             modeltype = "lm",
                                             covardetails = "lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M2",
                                             covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                             xcolumn     = c(64:69),  #numeric infant growth outcomes
                                             ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                             outputfolder = "model_results_milk_HMOs",
                                             timepoint = "M2")

inff_growth_by_HMOs_results_M3 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                             modeltype = "lm",
                                             covardetails = "lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M3",
                                             covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                             xcolumn     = c(64:69),  #numeric infant growth outcomes
                                             ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                             outputfolder = "model_results_milk_HMOs",
                                             timepoint = "M3")

inff_growth_by_HMOs_results_M6 <- run.models(datasetname = "infant_outcomes_by_milk_HMOs",
                                             modeltype = "lm",
                                             covardetails = "lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M6",
                                             covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                             xcolumn     = c(64:69),  #numeric infant growth outcomes
                                             ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
                                             outputfolder = "model_results_milk_HMOs",
                                             timepoint = "M6")

## if necessary, re-import results
inff_growth_by_HMOs_results_W2 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/W2/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_W2_W2_results.txt", header=T, sep="\t", stringsAsFactors=T)
inff_growth_by_HMOs_results_M1 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M1/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M1_M1_results.txt", header=T, sep="\t", stringsAsFactors=T)
inff_growth_by_HMOs_results_M2 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M2/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M2_M2_results.txt", header=T, sep="\t", stringsAsFactors=T)
inff_growth_by_HMOs_results_M3 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M3/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M3_M3_results.txt", header=T, sep="\t", stringsAsFactors=T)
inff_growth_by_HMOs_results_M6 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/model_results_milk_HMOs/M6/240222_infant_outcomes_by_milk_HMOs_model_results_milk_HMOs_lm_with_correction_for_infant_sex_GA_inff_growth_by_milk_HMOs_results_M6_M6_results.txt", header=T, sep="\t", stringsAsFactors=T)


# ### ===== 17.3 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ NUMERIC MILK HMO CONCENTRATIONS - STRATIFIED BY SEX ===== ###
# 
# ## -> use numeric milk HMO concentrations as y and numeric/categorical infant outcomes as x
# ##    lm() model, run per time point, with corection for gestational age (invr)
# ##    stratify by infant sex
# ##    for each time point: select which infant outcomes to include (have data at that time point)
# 
# ### === BOYS === ###
# 
# boys_growth_by_HMOs_results_W2 <- run.models(datasetname = "boys_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_HMOs_results_W2",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "W2")
# 
# boys_growth_by_HMOs_results_M1 <- run.models(datasetname = "boys_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_HMOs_results_M1",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "M1")
# 
# boys_growth_by_HMOs_results_M2 <- run.models(datasetname = "boys_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_HMOs_results_M2",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "M2")
# 
# boys_growth_by_HMOs_results_M3 <- run.models(datasetname = "boys_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_HMOs_results_M3",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "M3")
# 
# boys_growth_by_HMOs_results_M6 <- run.models(datasetname = "boys_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_HMOs_results_M6",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "M6")
# 
# 
# ### === GIRLS === ###
# 
# girls_growth_by_HMOs_results_W2 <- run.models(datasetname = "girls_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_HMOs_results_W2",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "W2")
# 
# girls_growth_by_HMOs_results_M1 <- run.models(datasetname = "girls_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_HMOs_results_M1",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "M1")
# 
# girls_growth_by_HMOs_results_M2 <- run.models(datasetname = "girls_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_HMOs_results_M2",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "M2")
# 
# girls_growth_by_HMOs_results_M3 <- run.models(datasetname = "girls_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_HMOs_results_M3",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "M3")
# 
# girls_growth_by_HMOs_results_M6 <- run.models(datasetname = "girls_growth_by_milk_HMOs",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_HMOs_results_M6",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(64:69),  #numeric infant growth outcomes
#                                              ycolumn     = c(36:63),  #invr-transformed milk HMO concentrations
#                                              outputfolder = "model_results_milk_HMOs_by_sex",
#                                              timepoint = "M6")


##### =========================== 18. ASSOCIATIONS OF MILK MICROBIOTA (ALPHA DIVERSITY, RELATIVE ABUNDANCES) WITH INFANT OUTCOME PHENOTYPES =========================== #####

for (i in 70:101){
  print(i)
  print(colnames(ibact_num_RA_invr_clr)[i])
  print(table(ibact_num_RA_invr_clr[!is.na(ibact_num_RA_invr_clr[,i]), "time_point"], useNA="ifany"))
}

### ===== 18.1 RUN MODELS FOR NUMERIC/CATEGORICAL INFANT OUTCOMES ~ NUMERIC MILK MICROBIOTA DATA ===== ###

## -> use numeric milk microbiota data as y and numeric/categorical infant outcomes as x
##    lm() model, run per time point
##    no correction factors, except for infant growth outcomes -> for growth, correct for infant sex and gestational age (invr)
##    for each time point: select which infant outcomes to include (have data at that time point)

inff_outcomes_by_milk_microbiota_results_M1 <- run.models(datasetname = "infant_outcomes_by_milk_microbiota",
                                               modeltype = "lm",
                                               covardetails = "lm_with_correction_DNAbatch_mreads_inff_outcomes_by_milk_microbiota_results_M1",
                                               covariates = "DNA_isolation_batch + seq_16S_n_reads_clean + ",
                                               inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="1_month",],
                                               xcolumn     = c(76:78,91:97),  #infant outcomes for M1 (categorical and invr-transformed numeric phenotypes), excl. infant growth (70:75)
                                               ycolumn     = c(67:69,102:137),  #invr-transformed milk allplha diversity and clr-transformed milk bacterial relative abundances
                                               outputfolder = "model_results_milk_microbiota",
                                               timepoint = "M1")

inff_outcomes_by_milk_microbiota_results_M2 <- run.models(datasetname = "infant_outcomes_by_milk_microbiota",
                                               modeltype = "lm",
                                               covardetails = "lm_with_correction_DNAbatch_mreads_inff_outcomes_by_milk_microbiota_results_M2",
                                               covariates = "DNA_isolation_batch + seq_16S_n_reads_clean + ",
                                               inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="2_months",],
                                               xcolumn     = c(76:78,91:97),  #infant outcomes for M2 (categorical and invr-transformed numeric phenotypes), excl. infant growth (70:75)
                                               ycolumn     = c(67:69,102:137),  #invr-transformed milk allplha diversity and clr-transformed milk bacterial relative abundances
                                               outputfolder = "model_results_milk_microbiota",
                                               timepoint = "M2")

inff_outcomes_by_milk_microbiota_results_M3 <- run.models(datasetname = "infant_outcomes_by_milk_microbiota",
                                               modeltype = "lm",
                                               covardetails = "lm_with_correction_DNAbatch_mreads_inff_outcomes_by_milk_microbiota_results_M3",
                                               covariates = "DNA_isolation_batch + seq_16S_n_reads_clean + ",
                                               inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="3_months",],
                                               xcolumn     = c(76:101),  #infant outcomes for M3 (categorical and invr-transformed numeric phenotypes), excl. infant growth (70:75)
                                               ycolumn     = c(67:69,102:137),  #invr-transformed milk allplha diversity and clr-transformed milk bacterial relative abundances
                                               outputfolder = "model_results_milk_microbiota",
                                               timepoint = "M3")

# Note: All M6 samples were in the second milk DNA isolation batch, so we co not correct for the milk DNA isolation batch!
inff_outcomes_by_milk_microbiota_results_M6 <- run.models(datasetname = "infant_outcomes_by_milk_microbiota",
                                                          modeltype = "lm",
                                                          covardetails = "lm_with_correction_mreads_inff_outcomes_by_milk_microbiota_results_M6",
                                                          covariates = "seq_16S_n_reads_clean + ",
                                                          inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="6_months",],
                                                          xcolumn     = c(76:78,91:95,97:101),  #infant outcomes for M6 (categorical and invr-transformed numeric phenotypes), excl. infant growth (70:75)
                                                          ycolumn     = c(67:69,102:137),  #invr-transformed milk allplha diversity and clr-transformed milk bacterial relative abundances
                                                          outputfolder = "model_results_milk_microbiota",
                                                          timepoint = "M6")


### ===== 18.2 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ NUMERIC MILK MICROBIOTA DATA ===== ###

## -> use numeric milk microbiota data as y and numeric/categorical infant outcomes as x
##    lm() model, run per time point
##    no correction factors, except for infant growth outcomes -> for growth, correct for infant sex and gestational age (invr)
##    for each time point: select which infant outcomes to include (have data at that time point)


inff_growth_by_milk_microbiota_results_M1 <- run.models(datasetname = "infant_outcomes_by_milk_microbiota",
                                             modeltype = "lm",
                                             covardetails = "lm_with_correction_for_DNAbatch_mreads_infant_sex_GA_inff_growth_by_milk_microbiota_results_M1",
                                             covariates = "DNA_isolation_batch + seq_16S_n_reads_clean + infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                             inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="1_month" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
                                             xcolumn     = c(70:75),  #numeric infant growth outcomes
                                             ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
                                             outputfolder = "model_results_milk_microbiota",
                                             timepoint = "M1")

inff_growth_by_milk_microbiota_results_M2 <- run.models(datasetname = "infant_outcomes_by_milk_microbiota",
                                             modeltype = "lm",
                                             covardetails = "lm_with_correction_for_DNAbatch_mreads_infant_sex_GA_inff_growth_by_milk_microbiota_results_M2",
                                             covariates = "DNA_isolation_batch + seq_16S_n_reads_clean + infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                             inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="2_months" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
                                             xcolumn     = c(70:75),  #numeric infant growth outcomes
                                             ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
                                             outputfolder = "model_results_milk_microbiota",
                                             timepoint = "M2")

inff_growth_by_milk_microbiota_results_M3 <- run.models(datasetname = "infant_outcomes_by_milk_microbiota",
                                             modeltype = "lm",
                                             covardetails = "lm_with_correction_for_DNAbatch_mreads_infant_sex_GA_inff_growth_by_milk_microbiota_results_M3",
                                             covariates = "DNA_isolation_batch + seq_16S_n_reads_clean + infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                             inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="3_months" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
                                             xcolumn     = c(70:75),  #numeric infant growth outcomes
                                             ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
                                             outputfolder = "model_results_milk_microbiota",
                                             timepoint = "M3")

# Note: All M6 samples were in the second milk DNA isolation batch, so we co not correct for the milk DNA isolation batch!
inff_growth_by_milk_microbiota_results_M6 <- run.models(datasetname = "infant_outcomes_by_milk_microbiota",
                                             modeltype = "lm",
                                             covardetails = "lm_with_correction_for_mreads_infant_sex_GA_inff_growth_by_milk_microbiota_results_M6",
                                             covariates = "seq_16S_n_reads_clean + infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                             inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="6_months" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
                                             xcolumn     = c(70:75),  #numeric infant growth outcomes
                                             ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
                                             outputfolder = "model_results_milk_microbiota",
                                             timepoint = "M6")


# ### ===== 18.3 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ NUMERIC MILK MICROBIOTA DATA - STRATIFIED BY SEX ===== ###
# 
# ## Note: these models need to be updated - they were not corrected for DNA isolation batch and milk sequencing read counts!
# 
# ## -> use numeric milk microbiota data as y and numeric/categorical infant outcomes as x
# ##    lm() model, run per time point, with correction for gestational age (invr)
# ##    stratify by infant sex
# ##    for each time point: select which infant outcomes to include (have data at that time point)
# 
# ### === BOYS === ###
# 
# boys_growth_by_microbiota_results_M1 <- run.models(datasetname = "boys_growth_by_milk_microbiota",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_microbiota_results_M1",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="1_month" & ibact_num_RA_invr_clr$infant_sex=="male" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(70:75),  #numeric infant growth outcomes
#                                              ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
#                                              outputfolder = "model_results_milk_microbiota_by_sex",
#                                              timepoint = "M1")
# 
# boys_growth_by_microbiota_results_M2 <- run.models(datasetname = "boys_growth_by_milk_microbiota",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_microbiota_results_M2",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="2_months" & ibact_num_RA_invr_clr$infant_sex=="male" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(70:75),  #numeric infant growth outcomes
#                                              ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
#                                              outputfolder = "model_results_milk_microbiota_by_sex",
#                                              timepoint = "M2")
# 
# boys_growth_by_microbiota_results_M3 <- run.models(datasetname = "boys_growth_by_milk_microbiota",
#                                              modeltype = "lm",
#                                              covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_microbiota_results_M3",
#                                              covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                              inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="3_months" & ibact_num_RA_invr_clr$infant_sex=="male" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
#                                              xcolumn     = c(70:75),  #numeric infant growth outcomes
#                                              ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
#                                              outputfolder = "model_results_milk_microbiota_by_sex",
#                                              timepoint = "M3")
# 
# boys_growth_by_microbiota_results_M6 <- run.models(datasetname = "boys_growth_by_milk_microbiota",
#                                                    modeltype = "lm",
#                                                    covardetails = "lm_with_correction_for_GA_boys_growth_by_milk_microbiota_results_M6",
#                                                    covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                    inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="6_months" & ibact_num_RA_invr_clr$infant_sex=="male" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
#                                                    xcolumn     = c(70:75),  #numeric infant growth outcomes
#                                                    ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
#                                                    outputfolder = "model_results_milk_microbiota_by_sex",
#                                                    timepoint = "M6")
# 
# 
# ### === GIRLS === ###
# 
# girls_growth_by_microbiota_results_M1 <- run.models(datasetname = "girls_growth_by_milk_microbiota",
#                                                    modeltype = "lm",
#                                                    covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_microbiota_results_M1",
#                                                    covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                    inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="1_month" & ibact_num_RA_invr_clr$infant_sex=="female" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
#                                                    xcolumn     = c(70:75),  #numeric infant growth outcomes
#                                                    ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
#                                                    outputfolder = "model_results_milk_microbiota_by_sex",
#                                                    timepoint = "M1")
# 
# girls_growth_by_microbiota_results_M2 <- run.models(datasetname = "girls_growth_by_milk_microbiota",
#                                                    modeltype = "lm",
#                                                    covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_microbiota_results_M2",
#                                                    covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                    inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="2_months" & ibact_num_RA_invr_clr$infant_sex=="female" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
#                                                    xcolumn     = c(70:75),  #numeric infant growth outcomes
#                                                    ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
#                                                    outputfolder = "model_results_milk_microbiota_by_sex",
#                                                    timepoint = "M2")
# 
# girls_growth_by_microbiota_results_M3 <- run.models(datasetname = "girls_growth_by_milk_microbiota",
#                                                    modeltype = "lm",
#                                                    covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_microbiota_results_M3",
#                                                    covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                    inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="3_months" & ibact_num_RA_invr_clr$infant_sex=="female" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
#                                                    xcolumn     = c(70:75),  #numeric infant growth outcomes
#                                                    ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
#                                                    outputfolder = "model_results_milk_microbiota_by_sex",
#                                                    timepoint = "M3")
# 
# girls_growth_by_microbiota_results_M6 <- run.models(datasetname = "girls_growth_by_milk_microbiota",
#                                                    modeltype = "lm",
#                                                    covardetails = "lm_with_correction_for_GA_girls_growth_by_milk_microbiota_results_M6",
#                                                    covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                    inputdata   = ibact_num_RA_invr_clr[ibact_num_RA_invr_clr$time_point=="6_months" & ibact_num_RA_invr_clr$infant_sex=="female" & !is.na(ibact_num_RA_invr_clr$infant_sex) & !is.na(ibact_num_RA_invr_clr$mother_birth_gestational_age_weeks_invr),],
#                                                    xcolumn     = c(70:75),  #numeric infant growth outcomes
#                                                    ycolumn     = c(67:69,102:137),  #invr-transformed milk alpha diversity and clr-transformed milk bacterial relative abundances
#                                                    outputfolder = "model_results_milk_microbiota_by_sex",
#                                                    timepoint = "M6")


##### =========================== 19. ASSOCIATIONS OF MILK GROUPS WITH INFANT OUTCOME PHENOTYPES =========================== #####

## function to run lmer or lm models with/without correction for covariates
run.models.milkgroups <- function(datasetname, modeltype, covardetails, covariates, inputdata, xcolumn, ycolumn, outputfolder, timepoint){
  
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
      
      print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      # print(table(my_df$mother_milk_HMO_milk_group))
      n_LeplusSeplus <- table(my_df$mother_milk_HMO_milk_group)["Le+Se+"]
      n_LeplusSeminus <- table(my_df$mother_milk_HMO_milk_group)["Le+Se-"]
      n_LeminusSeplus <- table(my_df$mother_milk_HMO_milk_group)["Le-Se+"]
      n_LeminusSeminus <- table(my_df$mother_milk_HMO_milk_group)["Le-Se-"]
      
      if (modeltype == "lmer"){
        
        ## RUN MIXED MODEL
        paste0("Run mixed model")
        
        my.lmer.model <- paste0("my_df[,j] ~ ", covariates, "my_df[,i]")
        print(paste0("Model function: ", my.lmer.model))
        
        m1 <- lmerTest::lmer(as.formula(my.lmer.model), REML=F, data=my_df)
        sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
        sum1$rows <- rownames(sum1)
        print(sum1)
        
        #retrieve estimates and p values from model with phenotype of interest (from sum1)
        levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        stat <- c(stat, rep(paste0(covardetails), length(sum1[grep("my_df",sum1$rows),"rows"])))
        d <- c(d, rep(datasetname, length(sum1[grep("my_df",sum1$rows),"rows"])))
        x <- c(x, rep(colnames(inputdata)[i], length(sum1[grep("my_df",sum1$rows),"rows"])))
        y <- c(y, rep(colnames(inputdata)[j], length(sum1[grep("my_df",sum1$rows),"rows"])))
        n <- c(n, rep(nrow(my_df), length(sum1[grep("my_df",sum1$rows),"rows"])))
      }else{
        
        ##RUN LINEAR MODEL
        paste0("Run linear model")
        
        my.lm.model <- paste0("my_df[,j] ~ ", covariates, "my_df[,i]")
        print(paste0("Model function: ", my.lm.model))
        
        m1 <- lm(as.formula(my.lm.model), data=my_df)
        sum1 <- as.data.frame(as.matrix(summary(m1)$coefficients))
        sum1$rows <- rownames(sum1)
        print(sum1)
        
        #retrieve estimates and p values from model with phenotype of interest (from sum1)
        levels <- c(levels, sum1[grep("my_df",sum1$rows),"rows"]) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        e <- c(e, c(sum1[grep("my_df",sum1$rows), "Estimate"])) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        p <- c(p, c(sum1[grep("my_df",sum1$rows), "Pr(>|t|)"])) #first rows are for covariates, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
        stat <- c(stat, rep(paste0(covardetails), length(sum1[grep("my_df",sum1$rows),"rows"])))
        d <- c(d, rep(datasetname, length(sum1[grep("my_df",sum1$rows),"rows"])))
        x <- c(x, rep(colnames(inputdata)[i], length(sum1[grep("my_df",sum1$rows),"rows"])))
        y <- c(y, rep(colnames(inputdata)[j], length(sum1[grep("my_df",sum1$rows),"rows"])))
        n <- c(n, c(n_LeplusSeplus + n_LeplusSeminus, n_LeplusSeplus + n_LeminusSeplus, n_LeplusSeplus + n_LeminusSeminus))
        # n <- c(n, rep(nrow(my_df), length(sum1[grep("my_df",sum1$rows),"rows"])))
      }
      
      ##PLOT DATA
      
      if (class(my_df[,i]) == "numeric" | class(my_df[,i]) == "integer"){
        print("Generating scatterplot")
        scatterplot1 <- ggplot(my_df, aes(x=my_df[,i], y=my_df[,j]))+
          geom_point(alpha=0.5)+
          theme_bw()+
          facet_grid(.~time_point_numeric)+
          geom_smooth(method='lm')+
          labs(x=paste0(colnames(my_df)[i]), y=paste0(colnames(my_df)[j]),
               title=paste0(colnames(my_df)[j], " by ", colnames(my_df)[i], "\n",
                            "n=", nrow(my_df), "\n",
                            covardetails, "\n",
                            "p=", signif(p[length(p)],3)))

        if (p[length(p)]<0.05){
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint, "/nominally_significant/240321_infant_outcomes_", covardetails, "_", colnames(my_df)[j] ,"_by_", colnames(my_df)[i], "_and_time", "_", timepoint, ".pdf"), scatterplot1, dpi=500, width=30, height=15, units="cm")
        }else{
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint,  "/not_significant/240321_infant_outcomes_", covardetails, "_", colnames(my_df)[j] ,"_by_", colnames(my_df)[i], "_and_time", "_", timepoint, ".pdf"), scatterplot1, dpi=500, width=30, height=15, units="cm")
        }
      
      }else{
        print("Generating boxplots")
        boxplot1 <- ggplot(my_df, aes(x=my_df[,i], y=my_df[,j]))+
          geom_jitter(alpha=0.5)+
          geom_boxplot(alpha=0)+
          theme_bw()+
          facet_grid(.~time_point_numeric)+
          labs(x=paste0(colnames(my_df)[i]), y=paste0(colnames(my_df)[j]),
               title=paste0(colnames(my_df)[j], " by ", colnames(my_df)[i], "\n",
                            "n=", nrow(my_df), "\n",
                            covardetails, "\n",
                            "p=", signif(p[length(p)],3)))

        if (p[length(p)]<0.05){
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint,  "/nominally_significant/240321_infant_outcomes_", covardetails, "_", colnames(my_df)[j] ,"_by_", colnames(my_df)[i], "_and_time", "_", timepoint, ".pdf"), boxplot1, dpi=500, width=30, height=15, units="cm")
        }else{
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint,  "/not_significant/240321_infant_outcomes_", covardetails, "_", colnames(my_df)[j] ,"_by_", colnames(my_df)[i], "_and_time", "_", timepoint, ".pdf"), boxplot1, dpi=500, width=30, height=15, units="cm")
        }
      }
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, estimate=e, p=p)
  
  #remove "my_df[, i]" from levels
  res$levels <- substring(res$levels, 11)
  
  #correct for multiple testing
  res$FDR <- p.adjust(res$p, method="BH")
  
  #resort results by FDR and p
  res <- res[order(res$FDR, res$p),]
  
  #save output table
  write.table(res, file=paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint, "/240321_", datasetname, "_", outputfolder, "_", covardetails, "_", timepoint, "_results.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  
  #return output table
  return(res)
}


for (i in 64:96){
  print(i)
  print(class(ihmo_num_invr[,i]))
  print(colnames(ihmo_num_invr)[i])
  print(table(ihmo_num_invr[!is.na(ihmo_num_invr[,i]), "time_point"], useNA="ifany"))
}

# categorical infant outcome phenotypes:
# 70 infant_cry_crying_more_than_3h_per_day
# 85 infant_dis_ROME_colic
# 86 infant_dis_ROME_e1_regurgitation
# 90 infant_dis_ROME_e32_sudden_inconsolable_crying_fits

### ===== 19.1 RUN MODELS FOR NUMERIC INFANT OUTCOMES ~ MILK GROUPS ===== ###

## -> use numeric infant outcomes as y and categorical milk group as x
##    lm() model, run per time point
##    no correction factors, except for infant growth outcomes -> for growth, correct for infant sex and gestational age (invr)
##    for each time point: select which infant outcomes to include (have data at that time point)

inff_outcomes_by_milkgroup_results_W2 <- run.models.milkgroups(datasetname = "infant_outcomes_by_milkgroup",
                                               modeltype = "lm",
                                               covardetails = "lm_without_correction_inff_outcomes_by_milkgroup_results_W2",
                                               covariates = "",
                                               inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months",],
                                               xcolumn     = c(18),    #categorical maternal milk group
                                               ycolumn     = c(71:72,91:92),  #infant outcomes for W2 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69) and factors (70)
                                               outputfolder = "model_results_milkgroup",
                                               timepoint = "W2")

inff_outcomes_by_milkgroup_results_M1 <- run.models.milkgroups(datasetname = "infant_outcomes_by_milkgroup",
                                               modeltype = "lm",
                                               covardetails = "lm_without_correction_inff_outcomes_by_milkgroup_results_M1",
                                               covariates = "",
                                               inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                               xcolumn     = c(18),    #categorical maternal milk group
                                               ycolumn     = c(71:72,87:89,91:92),  #infant outcomes for M1 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69) and factors (70,85:86,90)
                                               outputfolder = "model_results_milkgroup",
                                               timepoint = "M1")

inff_outcomes_by_milkgroup_results_M2 <- run.models.milkgroups(datasetname = "infant_outcomes_by_milkgroup",
                                               modeltype = "lm",
                                               covardetails = "lm_without_correction_inff_outcomes_by_milkgroup_results_M2",
                                               covariates = "",
                                               inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months",],
                                               xcolumn     = c(18),    #categorical maternal milk group
                                               ycolumn     = c(71:72,87:89,91:92),  #infant outcomes for M2 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69) and factors (70,85:86,90)
                                               outputfolder = "model_results_milkgroup",
                                               timepoint = "M2")

inff_outcomes_by_milkgroup_results_M3 <- run.models.milkgroups(datasetname = "infant_outcomes_by_milkgroup",
                                               modeltype = "lm",
                                               covardetails = "lm_without_correction_inff_outcomes_by_milkgroup_results_M3",
                                               covariates = "",
                                               inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months",],
                                               xcolumn     = c(18),    #categorical maternal milk group
                                               ycolumn     = c(71:84,87:89,91:96),  #infant outcomes for M3 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69) and factors (70,85:86,90)
                                               outputfolder = "model_results_milkgroup",
                                               timepoint = "M3")

inff_outcomes_by_milkgroup_results_M6 <- run.models.milkgroups(datasetname = "infant_outcomes_by_milkgroup",
                                               modeltype = "lm",
                                               covardetails = "lm_without_correction_inff_outcomes_by_milkgroup_results_M6",
                                               covariates = "",
                                               inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months",],
                                               xcolumn     = c(18),    #categorical maternal milk group
                                               ycolumn     = c(71:72,87:89,92:96),  #infant outcomes for M6 (categorical and invr-transformed numeric phenotypes), excl. infant growth (64:69) and factors (70,86,90)
                                               outputfolder = "model_results_milkgroup",
                                               timepoint = "M6")


### ===== 19.2 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ MILK GROUPS ===== ###

## -> use numeric infant growth outcomes as y and categorical milk group as x
##    lm() model, run per time point
##    no correction factors, except for infant growth outcomes -> for growth, correct for infant sex and gestational age (invr)
##    for each time point: select which infant outcomes to include (have data at that time point)

inff_growth_by_milkgroup_results_W2 <- run.models.milkgroups(datasetname = "infant_growth_by_milkgroup",
                                                    modeltype = "lm",
                                                    covardetails = "lm_with_correction_for_infant_sex_and_GA_inff_growth_by_milkgroup_results_W2",
                                                    covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                                    inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                                    xcolumn     = c(18),    #categorical maternal milk group
                                                    ycolumn     = c(64:69),  #infant growth outcomes
                                                    outputfolder = "model_results_milkgroup",
                                                    timepoint = "W2")

inff_growth_by_milkgroup_results_M1 <- run.models.milkgroups(datasetname = "infant_growth_by_milkgroup",
                                                    modeltype = "lm",
                                                    covardetails = "lm_with_correction_for_infant_sex_and_GA_inff_growth_by_milkgroup_results_M1",
                                                    covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                                    inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                                    xcolumn     = c(18),    #categorical maternal milk group
                                                    ycolumn     = c(64:69),  #infant growth outcomes
                                                    outputfolder = "model_results_milkgroup",
                                                    timepoint = "M1")

inff_growth_by_milkgroup_results_M2 <- run.models.milkgroups(datasetname = "infant_growth_by_milkgroup",
                                                    modeltype = "lm",
                                                    covardetails = "lm_with_correction_for_infant_sex_and_GA_inff_growth_by_milkgroup_results_M2",
                                                    covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                                    inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                                    xcolumn     = c(18),    #categorical maternal milk group
                                                    ycolumn     = c(64:69),  #infant growth outcomes
                                                    outputfolder = "model_results_milkgroup",
                                                    timepoint = "M2")

inff_growth_by_milkgroup_results_M3 <- run.models.milkgroups(datasetname = "infant_growth_by_milkgroup",
                                                    modeltype = "lm",
                                                    covardetails = "lm_with_correction_for_infant_sex_and_GA_inff_growth_by_milkgroup_results_M3",
                                                    covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                                    inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                                    xcolumn     = c(18),    #categorical maternal milk group
                                                    ycolumn     = c(64:69),  #infant growth outcomes
                                                    outputfolder = "model_results_milkgroup",
                                                    timepoint = "M3")

inff_growth_by_milkgroup_results_M6 <- run.models.milkgroups(datasetname = "infant_growth_by_milkgroup",
                                                    modeltype = "lm",
                                                    covardetails = "lm_with_correction_for_infant_sex_and_GA_inff_growth_by_milkgroup_results_M6",
                                                    covariates = "infant_sex + mother_birth_gestational_age_weeks_invr + ",
                                                    inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
                                                    xcolumn     = c(18),    #categorical maternal milk group
                                                    ycolumn     = c(64:69),  #infant growth outcomes
                                                    outputfolder = "model_results_milkgroup",
                                                    timepoint = "M6")


# ### ===== 19.3 RUN MODELS FOR NUMERIC INFANT GROWTH OUTCOMES ~ MILK GROUPS - STRATIFIED BY INFANT SEX ===== ###
# 
# ## -> use numeric infant growth outcomes as y and categorical milk group as x
# ##    lm() model, run per time point
# ##    no correction factors, except for infant growth outcomes -> for growth, correct for infant sex and gestational age (invr)
# ##    for each time point: select which infant outcomes to include (have data at that time point)
# 
# ### === BOYS === ###
# 
# boys_growth_by_milkgroup_results_W2 <- run.models(datasetname = "boys_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_boys_growth_by_milkgroup_results_W2",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "W2")
# 
# boys_growth_by_milkgroup_results_M1 <- run.models(datasetname = "boys_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_boys_growth_by_milkgroup_results_M1",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "M1")
# 
# boys_growth_by_milkgroup_results_M2 <- run.models(datasetname = "boys_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_boys_growth_by_milkgroup_results_M2",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "M2")
# 
# boys_growth_by_milkgroup_results_M3 <- run.models(datasetname = "boys_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_boys_growth_by_milkgroup_results_M3",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "M3")
# 
# boys_growth_by_milkgroup_results_M6 <- run.models(datasetname = "boys_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_boys_growth_by_milkgroup_results_M6",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months" & ihmo_num_invr$infant_sex=="male" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "M6")
# 
# 
# ### === GIRLS === ###
# 
# girls_growth_by_milkgroup_results_W2 <- run.models(datasetname = "girls_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_girls_growth_by_milkgroup_results_W2",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "W2")
# 
# girls_growth_by_milkgroup_results_M1 <- run.models(datasetname = "girls_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_girls_growth_by_milkgroup_results_M1",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "M1")
# 
# girls_growth_by_milkgroup_results_M2 <- run.models(datasetname = "girls_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_girls_growth_by_milkgroup_results_M2",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="2_months" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "M2")
# 
# girls_growth_by_milkgroup_results_M3 <- run.models(datasetname = "girls_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_girls_growth_by_milkgroup_results_M3",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="3_months" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "M3")
# 
# girls_growth_by_milkgroup_results_M6 <- run.models(datasetname = "girls_growth_by_milkgroup",
#                                                   modeltype = "lm",
#                                                   covardetails = "lm_with_correction_for_GA_girls_growth_by_milkgroup_results_M6",
#                                                   covariates = "mother_birth_gestational_age_weeks_invr + ",
#                                                   inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="6_months" & ihmo_num_invr$infant_sex=="female" & !is.na(ihmo_num_invr$infant_sex) & !is.na(ihmo_num_invr$mother_birth_gestational_age_weeks_invr),],
#                                                   xcolumn     = c(18),    #categorical maternal milk group
#                                                   ycolumn     = c(64:69),  #infant growth outcomes
#                                                   outputfolder = "model_results_milkgroup_by_sex",
#                                                   timepoint = "M6")


### ===== 19.4 RUN FISHER'S TEST FOR CATEGORICAL INFANT OUTCOMES ~ MILK GROUPS ===== ###

## create matrix
run.fishers.test.by.milkgroup <- function(datasetname, covardetails, inputdata, xcolumn, ycolumn, outputfolder, timepoint){

  x <- c()
  y <- c()
  n <- c()
  p <- c()

  for (i in xcolumn){
    for (j in ycolumn){

      ## select only rows without NAs
      print("Creating data frame")
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]

      x <- c(x, colnames(my_df)[i])
      y <- c(y, colnames(my_df)[j])
      n <- c(n, nrow(my_df))

      ## create matrix
      print("Creating matrix")
      my_matrix <- matrix(table(my_df[,i], my_df[,j]), nrow=4, dimnames = list(levels(my_df[,i]), levels(my_df[,j])))
      print(my_matrix)

      ## run Fisher's test
      print("Running Fisher's test")
      fisher.results <- fisher.test(my_matrix)
      print(fisher.results)
      p <- c(p, fisher.results$p.value)

    }
  }

  #save results in data frame
  print("Creating results dataframe")
  res <- data.frame(dataset=rep(datasetname, length(p)), statistic=rep("Fishers_test", length(p)), x=x, y=y, n_total=n, p=p)

  #correct for multiple testing
  res$FDR <- p.adjust(res$p, method="BH")

  #resort results by FDR and p
  res <- res[order(res$FDR, res$p),]

  #save output table
  write.table(res, file=paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/", outputfolder, "/", timepoint, "/240321_", datasetname, "_", outputfolder, "_", covardetails, "_", timepoint, "_results.txt"), row.names=F, col.names=T, sep="\t", quote=F)

  #return output table
  return(res)
}

for (i in 64:96){
  print(i)
  print(class(ihmo_num_invr[,i]))
  print(colnames(ihmo_num_invr)[i])
  print(table(ihmo_num_invr[!is.na(ihmo_num_invr[,i]), "time_point"], useNA="ifany"))
}

# categorical infant outcome phenotypes:
# 70 infant_cry_crying_more_than_3h_per_day: W2, M1, M2, M3, M6
# 85 infant_dis_ROME_colic: M1, M2, M3
# 86 infant_dis_ROME_e1_regurgitation: M1, M2, M3, M6
# 90 infant_dis_ROME_e32_sudden_inconsolable_crying_fits: M1, M2, M3, M6

# 18 mother_milk_HMO_milk_group

inff_catoutcomes_by_milkgroup_results_W2 <- run.fishers.test.by.milkgroup(datasetname = "infant_outcomes_by_milkgroup",
                                                             covardetails = "Fishers_test_inff_catoutcomes_by_milkgroup_results_W2",
                                                             inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="0.5_months",],
                                                             xcolumn     = c(18),  #categorical milk group
                                                             ycolumn     = c(70),  #categorical infant outcomes for W2
                                                             outputfolder = "model_results_milkgroup",
                                                             timepoint = "W2")

inff_catoutcomes_by_milkgroup_results_M1 <- run.fishers.test.by.milkgroup(datasetname = "infant_outcomes_by_milkgroup",
                                                                          covardetails = "Fishers_test_inff_catoutcomes_by_milkgroup_results_M1",
                                                                          inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                                                          xcolumn     = c(18),  #categorical milk group
                                                                          ycolumn     = c(70,85,86,90),  #categorical infant outcomes for M1
                                                                          outputfolder = "model_results_milkgroup",
                                                                          timepoint = "M1")

inff_catoutcomes_by_milkgroup_results_M2 <- run.fishers.test.by.milkgroup(datasetname = "infant_outcomes_by_milkgroup",
                                                                          covardetails = "Fishers_test_inff_catoutcomes_by_milkgroup_results_M2",
                                                                          inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                                                          xcolumn     = c(18),  #categorical milk group
                                                                          ycolumn     = c(70,85,86,90),  #categorical infant outcomes for M2
                                                                          outputfolder = "model_results_milkgroup",
                                                                          timepoint = "M2")

inff_catoutcomes_by_milkgroup_results_M3 <- run.fishers.test.by.milkgroup(datasetname = "infant_outcomes_by_milkgroup",
                                                                          covardetails = "Fishers_test_inff_catoutcomes_by_milkgroup_results_M3",
                                                                          inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                                                          xcolumn     = c(18),  #categorical milk group
                                                                          ycolumn     = c(70,85,86,90),  #categorical infant outcomes for M3
                                                                          outputfolder = "model_results_milkgroup",
                                                                          timepoint = "M3")

inff_catoutcomes_by_milkgroup_results_M6 <- run.fishers.test.by.milkgroup(datasetname = "infant_outcomes_by_milkgroup",
                                                                          covardetails = "Fishers_test_inff_catoutcomes_by_milkgroup_results_M6",
                                                                          inputdata   = ihmo_num_invr[ihmo_num_invr$time_point=="1_month",],
                                                                          xcolumn     = c(18),  #categorical milk group
                                                                          ycolumn     = c(70,86,90),  #categorical infant outcomes for M6
                                                                          outputfolder = "model_results_milkgroup",
                                                                          timepoint = "M6")


##### =========================== 20. COMBINE RESULTS FROM ASSOCIATIONS OF MILK GROUPS AND MILK HMOs WITH INFANT OUTCOME PHENOTYPES =========================== #####

## done - run models for numeric infant phenotypes by milk HMOs without corrections
## done - run models for numeric infant growth outcomes by milk HMOs with correction for infant sex and GA

## done - run models for numeric infant phenotypes (y) by milk group (x) without corrections
## done - run models for numeric infant growth outcomes (y) by milk group (x) with correction for infant sex and GA
## done - run Fisher's tests for categorical infant phenotypes (y) by milk group (x) without corrections

# [92] "inff_outcomes_by_HMOs_results_M1"           
# [93] "inff_outcomes_by_HMOs_results_M2"           
# [94] "inff_outcomes_by_HMOs_results_M3"           
# [95] "inff_outcomes_by_HMOs_results_M6"           
# [96] "inff_outcomes_by_HMOs_results_W2" 
# [68] "inff_growth_by_HMOs_results_M1"             
# [69] "inff_growth_by_HMOs_results_M2"             
# [70] "inff_growth_by_HMOs_results_M3"             
# [71] "inff_growth_by_HMOs_results_M6"             
# [72] "inff_growth_by_HMOs_results_W2"  
# [101] "inff_outcomes_by_milkgroup_results_M1"      
# [102] "inff_outcomes_by_milkgroup_results_M2"      
# [103] "inff_outcomes_by_milkgroup_results_M3"      
# [104] "inff_outcomes_by_milkgroup_results_M6"      
# [105] "inff_outcomes_by_milkgroup_results_W2"  
# [63] "inff_catoutcomes_by_milkgroup_results_M1"   
# [64] "inff_catoutcomes_by_milkgroup_results_M2"   
# [65] "inff_catoutcomes_by_milkgroup_results_M3"   
# [66] "inff_catoutcomes_by_milkgroup_results_M6"   
# [67] "inff_catoutcomes_by_milkgroup_results_W2"  
# [77] "inff_growth_by_milkgroup_results_M1"        
# [78] "inff_growth_by_milkgroup_results_M2"        
# [79] "inff_growth_by_milkgroup_results_M3"        
# [80] "inff_growth_by_milkgroup_results_M6"        
# [81] "inff_growth_by_milkgroup_results_W2" 


### ===== 20.1 COMBINE RESULTS PER TIME POINT ===== ###

### === W2 === ###

# inff_outcomes_by_milkgroup_results_W2 #12x9; x and y are correct
# inff_catoutcomes_by_milkgroup_results_W2 #1x7; x and y are correct
# inff_growth_by_milkgroup_results_W2 #18x9; x and y are correct
# inff_outcomes_by_HMOs_results_W2 #140x9; swap x and y
# inff_growth_by_HMOs_results_W2 #168x9; swap x and y

## swap x and y
inff_outcomes_by_HMOs_results_W2$new_x <- inff_outcomes_by_HMOs_results_W2$y
inff_outcomes_by_HMOs_results_W2$new_y <- inff_outcomes_by_HMOs_results_W2$x
inff_growth_by_HMOs_results_W2$new_x <- inff_growth_by_HMOs_results_W2$y
inff_growth_by_HMOs_results_W2$new_y <- inff_growth_by_HMOs_results_W2$x

## remove old x and y
inff_outcomes_by_HMOs_results_W2 <- inff_outcomes_by_HMOs_results_W2[,-c(3:4)]
inff_growth_by_HMOs_results_W2 <- inff_growth_by_HMOs_results_W2[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_HMOs_results_W2)[8:9] <- c("x", "y")
colnames(inff_growth_by_HMOs_results_W2)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_HMOs_results_W2 <- as.data.frame(rbind.fill(inff_outcomes_by_milkgroup_results_W2, inff_catoutcomes_by_milkgroup_results_W2, inff_growth_by_milkgroup_results_W2, inff_outcomes_by_HMOs_results_W2, inff_growth_by_HMOs_results_W2))

## check merged data frame
dim(all_outcomes_by_HMOs_results_W2) #339x9
unique(all_outcomes_by_HMOs_results_W2$x) #29
unique(all_outcomes_by_HMOs_results_W2$y) #11

## recalculate FDR
all_outcomes_by_HMOs_results_W2$FDR <- p.adjust(all_outcomes_by_HMOs_results_W2$p, method="BH")

## resort by FDR and p
all_outcomes_by_HMOs_results_W2 <- all_outcomes_by_HMOs_results_W2[order(all_outcomes_by_HMOs_results_W2$FDR, all_outcomes_by_HMOs_results_W2$p),]

## check results
nrow(all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$p<0.05,])   #23 nominally significant results

## save results
write.table(all_outcomes_by_HMOs_results_W2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)
all_outcomes_by_HMOs_results_W2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_W2.txt", header=T, sep="\t")


### === M1 === ###

# inff_outcomes_by_milkgroup_results_M1 #21x9; x and y are correct
# inff_catoutcomes_by_milkgroup_results_M1 #4x7; x and y are correct
# inff_growth_by_milkgroup_results_M1 #18x9; x and y are correct
# inff_outcomes_by_HMOs_results_M1 #308x9; swap x and y
# inff_growth_by_HMOs_results_M1 #168x9; swap x and y

## swap x and y
inff_outcomes_by_HMOs_results_M1$new_x <- inff_outcomes_by_HMOs_results_M1$y
inff_outcomes_by_HMOs_results_M1$new_y <- inff_outcomes_by_HMOs_results_M1$x
inff_growth_by_HMOs_results_M1$new_x <- inff_growth_by_HMOs_results_M1$y
inff_growth_by_HMOs_results_M1$new_y <- inff_growth_by_HMOs_results_M1$x

## remove old x and y
inff_outcomes_by_HMOs_results_M1 <- inff_outcomes_by_HMOs_results_M1[,-c(3:4)]
inff_growth_by_HMOs_results_M1 <- inff_growth_by_HMOs_results_M1[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_HMOs_results_M1)[8:9] <- c("x", "y")
colnames(inff_growth_by_HMOs_results_M1)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_HMOs_results_M1 <- as.data.frame(rbind.fill(inff_outcomes_by_milkgroup_results_M1, inff_catoutcomes_by_milkgroup_results_M1, inff_growth_by_milkgroup_results_M1, inff_outcomes_by_HMOs_results_M1, inff_growth_by_HMOs_results_M1))

## check merged data frame
dim(all_outcomes_by_HMOs_results_M1) #519x9
unique(all_outcomes_by_HMOs_results_M1$x) #29
unique(all_outcomes_by_HMOs_results_M1$y) #17

## recalculate FDR
all_outcomes_by_HMOs_results_M1$FDR <- p.adjust(all_outcomes_by_HMOs_results_M1$p, method="BH")

## resort by FDR and p
all_outcomes_by_HMOs_results_M1 <- all_outcomes_by_HMOs_results_M1[order(all_outcomes_by_HMOs_results_M1$FDR, all_outcomes_by_HMOs_results_M1$p),]

## check results
nrow(all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$p<0.05,])   #55 nominally significant results

## save results
write.table(all_outcomes_by_HMOs_results_M1, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
all_outcomes_by_HMOs_results_M1 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M1.txt", header=T, sep="\t")


### === M2 === ###

# inff_outcomes_by_milkgroup_results_M2 #21x9; x and y are correct
# inff_catoutcomes_by_milkgroup_results_M2 #4x7; x and y are correct
# inff_growth_by_milkgroup_results_M2 #18x9; x and y are correct
# inff_outcomes_by_HMOs_results_M2 #308x9; swap x and y
# inff_growth_by_HMOs_results_M2 #168x9; swap x and y

## swap x and y
inff_outcomes_by_HMOs_results_M2$new_x <- inff_outcomes_by_HMOs_results_M2$y
inff_outcomes_by_HMOs_results_M2$new_y <- inff_outcomes_by_HMOs_results_M2$x
inff_growth_by_HMOs_results_M2$new_x <- inff_growth_by_HMOs_results_M2$y
inff_growth_by_HMOs_results_M2$new_y <- inff_growth_by_HMOs_results_M2$x

## remove old x and y
inff_outcomes_by_HMOs_results_M2 <- inff_outcomes_by_HMOs_results_M2[,-c(3:4)]
inff_growth_by_HMOs_results_M2 <- inff_growth_by_HMOs_results_M2[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_HMOs_results_M2)[8:9] <- c("x", "y")
colnames(inff_growth_by_HMOs_results_M2)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_HMOs_results_M2 <- as.data.frame(rbind.fill(inff_outcomes_by_milkgroup_results_M2, inff_catoutcomes_by_milkgroup_results_M2, inff_growth_by_milkgroup_results_M2, inff_outcomes_by_HMOs_results_M2, inff_growth_by_HMOs_results_M2))

## check merged data frame
dim(all_outcomes_by_HMOs_results_M2) #519x9
unique(all_outcomes_by_HMOs_results_M2$x) #29
unique(all_outcomes_by_HMOs_results_M2$y) #17

## recalculate FDR
all_outcomes_by_HMOs_results_M2$FDR <- p.adjust(all_outcomes_by_HMOs_results_M2$p, method="BH")

## resort by FDR and p
all_outcomes_by_HMOs_results_M2 <- all_outcomes_by_HMOs_results_M2[order(all_outcomes_by_HMOs_results_M2$FDR, all_outcomes_by_HMOs_results_M2$p),]

## check results
nrow(all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$p<0.05,])   #40 nominally significant results

## save results
write.table(all_outcomes_by_HMOs_results_M2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
all_outcomes_by_HMOs_results_M2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M2.txt", header=T, sep="\t")


### === M3 === ###

# inff_outcomes_by_milkgroup_results_M3 #69x9; x and y are correct
# inff_catoutcomes_by_milkgroup_results_M3 #4x7; x and y are correct
# inff_growth_by_milkgroup_results_M3 #18x9; x and y are correct
# inff_outcomes_by_HMOs_results_M3 #756x9; swap x and y
# inff_growth_by_HMOs_results_M3 #168x9; swap x and y

## swap x and y
inff_outcomes_by_HMOs_results_M3$new_x <- inff_outcomes_by_HMOs_results_M3$y
inff_outcomes_by_HMOs_results_M3$new_y <- inff_outcomes_by_HMOs_results_M3$x
inff_growth_by_HMOs_results_M3$new_x <- inff_growth_by_HMOs_results_M3$y
inff_growth_by_HMOs_results_M3$new_y <- inff_growth_by_HMOs_results_M3$x

## remove old x and y
inff_outcomes_by_HMOs_results_M3 <- inff_outcomes_by_HMOs_results_M3[,-c(3:4)]
inff_growth_by_HMOs_results_M3 <- inff_growth_by_HMOs_results_M3[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_HMOs_results_M3)[8:9] <- c("x", "y")
colnames(inff_growth_by_HMOs_results_M3)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_HMOs_results_M3 <- as.data.frame(rbind.fill(inff_outcomes_by_milkgroup_results_M3, inff_catoutcomes_by_milkgroup_results_M3, inff_growth_by_milkgroup_results_M3, inff_outcomes_by_HMOs_results_M3, inff_growth_by_HMOs_results_M3))

## check merged data frame
dim(all_outcomes_by_HMOs_results_M3) #1015x9
unique(all_outcomes_by_HMOs_results_M3$x) #29
unique(all_outcomes_by_HMOs_results_M3$y) #33

## recalculate FDR
all_outcomes_by_HMOs_results_M3$FDR <- p.adjust(all_outcomes_by_HMOs_results_M3$p, method="BH")

## resort by FDR and p
all_outcomes_by_HMOs_results_M3 <- all_outcomes_by_HMOs_results_M3[order(all_outcomes_by_HMOs_results_M3$FDR, all_outcomes_by_HMOs_results_M3$p),]

## check results
nrow(all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$p<0.05,])   #71 nominally significant results

## save results
write.table(all_outcomes_by_HMOs_results_M3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
all_outcomes_by_HMOs_results_M3 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M3.txt", header=T, sep="\t")


### === M6 === ###

# inff_outcomes_by_milkgroup_results_M6 #30x9; x and y are correct
# inff_catoutcomes_by_milkgroup_results_M6 #3x7; x and y are correct
# inff_growth_by_milkgroup_results_M6 #18x9; x and y are correct
# inff_outcomes_by_HMOs_results_M6 #364x9; swap x and y
# inff_growth_by_HMOs_results_M6 #168x9; swap x and y

## swap x and y
inff_outcomes_by_HMOs_results_M6$new_x <- inff_outcomes_by_HMOs_results_M6$y
inff_outcomes_by_HMOs_results_M6$new_y <- inff_outcomes_by_HMOs_results_M6$x
inff_growth_by_HMOs_results_M6$new_x <- inff_growth_by_HMOs_results_M6$y
inff_growth_by_HMOs_results_M6$new_y <- inff_growth_by_HMOs_results_M6$x

## remove old x and y
inff_outcomes_by_HMOs_results_M6 <- inff_outcomes_by_HMOs_results_M6[,-c(3:4)]
inff_growth_by_HMOs_results_M6 <- inff_growth_by_HMOs_results_M6[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_HMOs_results_M6)[8:9] <- c("x", "y")
colnames(inff_growth_by_HMOs_results_M6)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_HMOs_results_M6 <- as.data.frame(rbind.fill(inff_outcomes_by_milkgroup_results_M6, inff_catoutcomes_by_milkgroup_results_M6, inff_growth_by_milkgroup_results_M6, inff_outcomes_by_HMOs_results_M6, inff_growth_by_HMOs_results_M6))

## check merged data frame
dim(all_outcomes_by_HMOs_results_M6) #583x9
unique(all_outcomes_by_HMOs_results_M6$x) #29
unique(all_outcomes_by_HMOs_results_M6$y) #19

## recalculate FDR
all_outcomes_by_HMOs_results_M6$FDR <- p.adjust(all_outcomes_by_HMOs_results_M6$p, method="BH")

## resort by FDR and p
all_outcomes_by_HMOs_results_M6 <- all_outcomes_by_HMOs_results_M6[order(all_outcomes_by_HMOs_results_M6$FDR, all_outcomes_by_HMOs_results_M6$p),]

## check results
nrow(all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$p<0.05,])   #39 nominally significant results

## save results
write.table(all_outcomes_by_HMOs_results_M6, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
all_outcomes_by_HMOs_results_M6 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M6.txt", header=T, sep="\t")


# ### ===== 20.2 CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS ===== ###
# 
# ## check for consistently nominally significant associations at several time points
# nomsign_HMO_results_W2 <- c(paste0(all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$p<0.05, "y"]))
# nomsign_HMO_results_M1 <- c(paste0(all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$p<0.05, "y"]))
# nomsign_HMO_results_M2 <- c(paste0(all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$p<0.05, "y"]))
# nomsign_HMO_results_M3 <- c(paste0(all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$p<0.05, "y"]))
# nomsign_HMO_results_M6 <- c(paste0(all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$p<0.05, "y"]))
# 
# ## consistently nominally significant results at the 2 largest time points (M1, M3)
# nomsign_HMO_results_M1_M3 <- intersect(nomsign_HMO_results_M1, nomsign_HMO_results_M3)
# # [1] "mother_milk_HMO_LNFP_I_ugml_invr_infant_BITSS_invr"                             
# # [2] "mother_milk_HMO_DFLNHa_ugml_invr_infant_BITSS_invr"                             
# # [3] "mother_milk_HMO_LNnFP_V_ugml_invr_infant_dis_atopic_dermatitis_SCORAD_M12_invr" 
# # [4] "mother_milk_HMO_Neut_ugml_invr_infant_dis_ROME_e13_stool_freq_last_month_invr"  
# # [5] "mother_milk_HMO_DFLNHa_ugml_invr_infant_dis_atopic_dermatitis_SCORAD_M12_invr"  
# # [6] "mother_milk_HMO_3FL_ugml_invr_infant_BITSS_invr"                                
# # [7] "mother_milk_HMO_LNT_ugml_invr_infant_dis_ROME_e13_stool_freq_last_month_invr"   
# # [8] "mother_milk_HMO_milk_group_infant_BITSS_invr"                                   
# # [9] "mother_milk_HMO_2FL_ugml_invr_infant_BITSS_invr"                                
# # [10] "mother_milk_HMO_2FL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"    
# # [11] "mother_milk_HMO_DFLNHa_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr" 
# # [12] "mother_milk_HMO_3FL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"    
# # [13] "mother_milk_HMO_3F3SL_ugml_invr_infant_BITSS_invr"                              
# # [14] "mother_milk_HMO_LNFP_I_ugml_invr_infant_dis_atopic_dermatitis_SCORAD_M12_invr"  
# # [15] "mother_milk_HMO_3F3SL_ugml_invr_infant_dis_atopic_dermatitis_SCORAD_M12_invr"   
# # [16] "mother_milk_HMO_3F3SL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"  
# # [17] "mother_milk_HMO_LNFP_II_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"
# # [18] "mother_milk_HMO_A_tetra_ugml_invr_infant_BITSS_invr"                            
# # [19] "mother_milk_HMO_3F3SL_ugml_invr_infant_growth_length_cm_invr"
# 
# ## consistently nominally significant results at the 3 largest time points (M1, M2, M3)
# nomsign_HMO_results_M1_M2_M3 <- intersect(nomsign_HMO_results_M1_M3, nomsign_HMO_results_M2)
# # [6] "mother_milk_HMO_3F3SL_ugml_invr_infant_growth_length_cm_invr"
# all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_W2$y=="infant_growth_length_cm_invr",]
# all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_M1$y=="infant_growth_length_cm_invr",]
# all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_M2$y=="infant_growth_length_cm_invr",]
# all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_M3$y=="infant_growth_length_cm_invr",]
# all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_M6$y=="infant_growth_length_cm_invr",]
# # -> more 3F3SL is !!inconsistently!! associated with increased infant weight gain
# 
# # [4] "mother_milk_HMO_3F3SL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"
# all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_W2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_M1$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_M2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_M3$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$x=="mother_milk_HMO_3F3SL_ugml_invr" & all_outcomes_by_HMOs_results_M6$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# # -> more 3F3SL is consistently associated with increased infant weight gain
# 
# # [3] "mother_milk_HMO_3FL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr" 
# all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$x=="mother_milk_HMO_3FL_ugml_invr" & all_outcomes_by_HMOs_results_W2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$x=="mother_milk_HMO_3FL_ugml_invr" & all_outcomes_by_HMOs_results_M1$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$x=="mother_milk_HMO_3FL_ugml_invr" & all_outcomes_by_HMOs_results_M2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$x=="mother_milk_HMO_3FL_ugml_invr" & all_outcomes_by_HMOs_results_M3$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$x=="mother_milk_HMO_3FL_ugml_invr" & all_outcomes_by_HMOs_results_M6$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# # -> more 3FL is consistently associated with increased infant weight gain
# 
# # [5] "mother_milk_HMO_LNFP_II_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"
# all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_W2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_M1$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_M2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_M3$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_M6$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# # -> more LNFP II is consistently associated with increased infant weight gain
# 
# # [2] "mother_milk_HMO_DFLNHa_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr" 
# all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_W2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_M1$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_M2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_M3$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_M6$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# # -> more DFLNHa is consistently associated with decreased infant weight gain
# 
# # [1] "mother_milk_HMO_DFLNHa_ugml_invr_infant_dis_atopic_dermatitis_SCORAD_M12_invr"  
# all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_W2$y=="infant_dis_atopic_dermatitis_SCORAD_M12_invr",]
# all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_M1$y=="infant_dis_atopic_dermatitis_SCORAD_M12_invr",]
# all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_M2$y=="infant_dis_atopic_dermatitis_SCORAD_M12_invr",]
# all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_M3$y=="infant_dis_atopic_dermatitis_SCORAD_M12_invr",]
# all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$x=="mother_milk_HMO_DFLNHa_ugml_invr" & all_outcomes_by_HMOs_results_M6$y=="infant_dis_atopic_dermatitis_SCORAD_M12_invr",]
# # -> more DFLNHa is consistently associated with increased SCORAD scores at M12 (more atopic dermatitis)   
# 
# ## consistently nominally significant results at the 4 largest time points (W2, M1, M2, M3)
# nomsign_HMO_results_W2_M1_M2_M3 <- intersect(nomsign_HMO_results_M1_M2_M3, nomsign_HMO_results_W2)
# #mother_milk_HMO_LNFP_II_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr is nominally significant at W2-M3
# all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_W2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_M1$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_M2$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_M3$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$x=="mother_milk_HMO_LNFP_II_ugml_invr" & all_outcomes_by_HMOs_results_M6$y=="infant_growth_weight_kg_gain_since_birth_invr",]
# # -> more LNFP II is consistently associated with increased infant weight gain
# 
# ## consistently nominally significant results at all time points
# nomsign_HMO_results_W2_M1_M2_M3_M6 <- intersect(nomsign_HMO_results_W2_M1_M2_M3, nomsign_HMO_results_M6)
# #no nominally significant associations are found at all time points
# 
# ## !! Note: Keep in mind that not all phenotypes were tested at all / ≥1 time point(s)!!


### ===== 20.3 MORE FOCUSED ANALYSIS: ONLY USE SELECTED TIME POINTS AND PHENOTYPES AND RECALCULATE FDR (TABLE S36A-D) ===== ###

## if necessary, re-import results:
all_outcomes_by_HMOs_results_W2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_W2.txt", header=T, sep="\t")
all_outcomes_by_HMOs_results_M1 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M1.txt", header=T, sep="\t")
all_outcomes_by_HMOs_results_M2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M2.txt", header=T, sep="\t")
all_outcomes_by_HMOs_results_M3 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M3.txt", header=T, sep="\t")
all_outcomes_by_HMOs_results_M6 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_HMOs_M6.txt", header=T, sep="\t")


## final selection of 11 phenotypes
sel_phenos_HMOs <- c("infant_dis_ROME_e1_regurgitation","infant_dis_ROME_e32_sudden_inconsolable_crying_fits",
                     "infant_growth_weight_kg_invr","infant_growth_weight_kg_gain_since_birth_invr",
                     "infant_growth_length_cm_invr","infant_growth_length_cm_gain_since_birth_invr",
                     "infant_cry_prev_week_average_crying_time_per_day_min_invr",
                     "infant_BITSS_invr",
                     "infant_dis_ROME_e13_stool_freq_last_month_invr","infant_dis_ROME_e14_stool_structure_invr","infant_dis_ROME_e25_mucus_stool_last_week_invr")


### ===== W2 ===== ###

## select data of interest
all_outcomes_by_HMOs_results_W2_v2 <- all_outcomes_by_HMOs_results_W2[all_outcomes_by_HMOs_results_W2$y %in% sel_phenos_HMOs,]

## recalculate FDR
all_outcomes_by_HMOs_results_W2_v2$FDR <- p.adjust(all_outcomes_by_HMOs_results_W2_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_HMOs_results_W2_v2 <- all_outcomes_by_HMOs_results_W2_v2[order(all_outcomes_by_HMOs_results_W2_v2$FDR, all_outcomes_by_HMOs_results_W2_v2$p),]

## add column for phenotype group
all_outcomes_by_HMOs_results_W2_v2$Phenotype_group <- NA
all_outcomes_by_HMOs_results_W2_v2[grep("ROME", all_outcomes_by_HMOs_results_W2_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_W2_v2[grep("BITSS", all_outcomes_by_HMOs_results_W2_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_W2_v2[grep("growth", all_outcomes_by_HMOs_results_W2_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_HMOs_results_W2_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_HMOs_results_W2_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column for reference category and ensure correct n
all_outcomes_by_HMOs_results_W2_v2[all_outcomes_by_HMOs_results_W2_v2$levels=="" & !is.na(all_outcomes_by_HMOs_results_W2_v2$levels), "levels"] <- NA
all_outcomes_by_HMOs_results_W2_v2$Reference_category <- NA
all_outcomes_by_HMOs_results_W2_v2[!is.na(all_outcomes_by_HMOs_results_W2_v2$levels), "Reference_category"] <- "Le+Se+"

## resort columns
all_outcomes_by_HMOs_results_W2_v3 <- all_outcomes_by_HMOs_results_W2_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_HMOs_results_W2_v3)

# write.table(all_outcomes_by_HMOs_results_W2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_HMOs_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_W2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_HMOs_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_W2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_HMOs_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_W2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_HMOs_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_HMOs_results_W2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_HMOs_W2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M1 (TABLE S36A) ===== ###

## select data of interest
all_outcomes_by_HMOs_results_M1_v2 <- all_outcomes_by_HMOs_results_M1[all_outcomes_by_HMOs_results_M1$y %in% sel_phenos_HMOs,]

## recalculate FDR
all_outcomes_by_HMOs_results_M1_v2$FDR <- p.adjust(all_outcomes_by_HMOs_results_M1_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_HMOs_results_M1_v2 <- all_outcomes_by_HMOs_results_M1_v2[order(all_outcomes_by_HMOs_results_M1_v2$FDR, all_outcomes_by_HMOs_results_M1_v2$p),]

## add column for phenotype group
all_outcomes_by_HMOs_results_M1_v2$Phenotype_group <- NA
all_outcomes_by_HMOs_results_M1_v2[grep("ROME", all_outcomes_by_HMOs_results_M1_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_M1_v2[grep("BITSS", all_outcomes_by_HMOs_results_M1_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_M1_v2[grep("growth", all_outcomes_by_HMOs_results_M1_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_HMOs_results_M1_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_HMOs_results_M1_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column for reference category and ensure correct n
all_outcomes_by_HMOs_results_M1_v2[all_outcomes_by_HMOs_results_M1_v2$levels=="" & !is.na(all_outcomes_by_HMOs_results_M1_v2$levels), "levels"] <- NA
all_outcomes_by_HMOs_results_M1_v2$Reference_category <- NA
all_outcomes_by_HMOs_results_M1_v2[!is.na(all_outcomes_by_HMOs_results_M1_v2$levels) & all_outcomes_by_HMOs_results_M1_v2$levels!="yes", "Reference_category"] <- "Le+Se+"
all_outcomes_by_HMOs_results_M1_v2[!is.na(all_outcomes_by_HMOs_results_M1_v2$levels) & all_outcomes_by_HMOs_results_M1_v2$levels=="yes", "Reference_category"] <- "no"

## resort columns
all_outcomes_by_HMOs_results_M1_v3 <- all_outcomes_by_HMOs_results_M1_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_HMOs_results_M1_v3)

# write.table(all_outcomes_by_HMOs_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_HMOs_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_HMOs_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_HMOs_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_HMOs_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_HMOs_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_HMOs_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M2 (TABLE S36B) ===== ###

## select data of interest
all_outcomes_by_HMOs_results_M2_v2 <- all_outcomes_by_HMOs_results_M2[all_outcomes_by_HMOs_results_M2$y %in% sel_phenos_HMOs,]

## recalculate FDR
all_outcomes_by_HMOs_results_M2_v2$FDR <- p.adjust(all_outcomes_by_HMOs_results_M2_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_HMOs_results_M2_v2 <- all_outcomes_by_HMOs_results_M2_v2[order(all_outcomes_by_HMOs_results_M2_v2$FDR, all_outcomes_by_HMOs_results_M2_v2$p),]

## add column for phenotype group
all_outcomes_by_HMOs_results_M2_v2$Phenotype_group <- NA
all_outcomes_by_HMOs_results_M2_v2[grep("ROME", all_outcomes_by_HMOs_results_M2_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_M2_v2[grep("BITSS", all_outcomes_by_HMOs_results_M2_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_M2_v2[grep("growth", all_outcomes_by_HMOs_results_M2_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_HMOs_results_M2_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_HMOs_results_M2_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column for reference category and ensure correct n
all_outcomes_by_HMOs_results_M2_v2[all_outcomes_by_HMOs_results_M2_v2$levels=="" & !is.na(all_outcomes_by_HMOs_results_M2_v2$levels), "levels"] <- NA
all_outcomes_by_HMOs_results_M2_v2$Reference_category <- NA
all_outcomes_by_HMOs_results_M2_v2[!is.na(all_outcomes_by_HMOs_results_M2_v2$levels) & all_outcomes_by_HMOs_results_M2_v2$levels!="yes", "Reference_category"] <- "Le+Se+"
all_outcomes_by_HMOs_results_M2_v2[!is.na(all_outcomes_by_HMOs_results_M2_v2$levels) & all_outcomes_by_HMOs_results_M2_v2$levels=="yes", "Reference_category"] <- "no"

## resort columns
all_outcomes_by_HMOs_results_M2_v3 <- all_outcomes_by_HMOs_results_M2_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_HMOs_results_M2_v3)

# write.table(all_outcomes_by_HMOs_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_HMOs_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_HMOs_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_HMOs_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_HMOs_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_HMOs_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_HMOs_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M3 (TABLE S36C) ===== ###

## select data of interest
all_outcomes_by_HMOs_results_M3_v2 <- all_outcomes_by_HMOs_results_M3[all_outcomes_by_HMOs_results_M3$y %in% sel_phenos_HMOs,]

## recalculate FDR
all_outcomes_by_HMOs_results_M3_v2$FDR <- p.adjust(all_outcomes_by_HMOs_results_M3_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_HMOs_results_M3_v2 <- all_outcomes_by_HMOs_results_M3_v2[order(all_outcomes_by_HMOs_results_M3_v2$FDR, all_outcomes_by_HMOs_results_M3_v2$p),]

## add column for phenotype group
all_outcomes_by_HMOs_results_M3_v2$Phenotype_group <- NA
all_outcomes_by_HMOs_results_M3_v2[grep("ROME", all_outcomes_by_HMOs_results_M3_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_M3_v2[grep("BITSS", all_outcomes_by_HMOs_results_M3_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_M3_v2[grep("growth", all_outcomes_by_HMOs_results_M3_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_HMOs_results_M3_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_HMOs_results_M3_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column for reference category and ensure correct n
all_outcomes_by_HMOs_results_M3_v2[all_outcomes_by_HMOs_results_M3_v2$levels=="" & !is.na(all_outcomes_by_HMOs_results_M3_v2$levels), "levels"] <- NA
all_outcomes_by_HMOs_results_M3_v2$Reference_category <- NA
all_outcomes_by_HMOs_results_M3_v2[!is.na(all_outcomes_by_HMOs_results_M3_v2$levels) & all_outcomes_by_HMOs_results_M3_v2$levels!="yes", "Reference_category"] <- "Le+Se+"
all_outcomes_by_HMOs_results_M3_v2[!is.na(all_outcomes_by_HMOs_results_M3_v2$levels) & all_outcomes_by_HMOs_results_M3_v2$levels=="yes", "Reference_category"] <- "no"

## resort columns
all_outcomes_by_HMOs_results_M3_v3 <- all_outcomes_by_HMOs_results_M3_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_HMOs_results_M3_v3)

# write.table(all_outcomes_by_HMOs_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_HMOs_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_HMOs_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_HMOs_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_HMOs_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_HMOs_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_HMOs_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M6 (TABLE S36D) ===== ###

## select data of interest
all_outcomes_by_HMOs_results_M6_v2 <- all_outcomes_by_HMOs_results_M6[all_outcomes_by_HMOs_results_M6$y %in% sel_phenos_HMOs,]

## recalculate FDR
all_outcomes_by_HMOs_results_M6_v2$FDR <- p.adjust(all_outcomes_by_HMOs_results_M6_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_HMOs_results_M6_v2 <- all_outcomes_by_HMOs_results_M6_v2[order(all_outcomes_by_HMOs_results_M6_v2$FDR, all_outcomes_by_HMOs_results_M6_v2$p),]

## add column for phenotype group
all_outcomes_by_HMOs_results_M6_v2$Phenotype_group <- NA
all_outcomes_by_HMOs_results_M6_v2[grep("ROME", all_outcomes_by_HMOs_results_M6_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_M6_v2[grep("BITSS", all_outcomes_by_HMOs_results_M6_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_HMOs_results_M6_v2[grep("growth", all_outcomes_by_HMOs_results_M6_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_HMOs_results_M6_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_HMOs_results_M6_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column for reference category and ensure correct n
all_outcomes_by_HMOs_results_M6_v2[all_outcomes_by_HMOs_results_M6_v2$levels=="" & !is.na(all_outcomes_by_HMOs_results_M6_v2$levels), "levels"] <- NA
all_outcomes_by_HMOs_results_M6_v2$Reference_category <- NA
all_outcomes_by_HMOs_results_M6_v2[!is.na(all_outcomes_by_HMOs_results_M6_v2$levels) & all_outcomes_by_HMOs_results_M6_v2$levels!="yes", "Reference_category"] <- "Le+Se+"
all_outcomes_by_HMOs_results_M6_v2[!is.na(all_outcomes_by_HMOs_results_M6_v2$levels) & all_outcomes_by_HMOs_results_M6_v2$levels=="yes", "Reference_category"] <- "no"

## resort columns
all_outcomes_by_HMOs_results_M6_v3 <- all_outcomes_by_HMOs_results_M6_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_HMOs_results_M6_v3)

# write.table(all_outcomes_by_HMOs_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_HMOs_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_HMOs_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_HMOs_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_HMOs_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_HMOs_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_HMOs_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_HMOs_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 20.4 MORE FOCUSED ANALYSIS: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS ===== ###

## check for consistently nominally significant associations at several time points
nomsign_HMO_results_W2 <- c(paste0(all_outcomes_by_HMOs_results_W2_v3[all_outcomes_by_HMOs_results_W2_v3$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_W2_v3[all_outcomes_by_HMOs_results_W2_v3$p<0.05, "y"]))
nomsign_HMO_results_M1 <- c(paste0(all_outcomes_by_HMOs_results_M1_v3[all_outcomes_by_HMOs_results_M1_v3$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_M1_v3[all_outcomes_by_HMOs_results_M1_v3$p<0.05, "y"]))
nomsign_HMO_results_M2 <- c(paste0(all_outcomes_by_HMOs_results_M2_v3[all_outcomes_by_HMOs_results_M2_v3$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_M2_v3[all_outcomes_by_HMOs_results_M2_v3$p<0.05, "y"]))
nomsign_HMO_results_M3 <- c(paste0(all_outcomes_by_HMOs_results_M3_v3[all_outcomes_by_HMOs_results_M3_v3$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_M3_v3[all_outcomes_by_HMOs_results_M3_v3$p<0.05, "y"]))
nomsign_HMO_results_M6 <- c(paste0(all_outcomes_by_HMOs_results_M6_v3[all_outcomes_by_HMOs_results_M6_v3$p<0.05, "x"], "_", all_outcomes_by_HMOs_results_M6_v3[all_outcomes_by_HMOs_results_M6_v3$p<0.05, "y"]))

## consistently nominally significant results at M1, M3
nomsign_HMO_results_M1_M3 <- intersect(nomsign_HMO_results_M1, nomsign_HMO_results_M3)
# [1] "mother_milk_HMO_LNFP_I_ugml_invr_infant_BITSS_invr"                             
# [2] "mother_milk_HMO_DFLNHa_ugml_invr_infant_BITSS_invr"                             
# [3] "mother_milk_HMO_Neut_ugml_invr_infant_dis_ROME_e13_stool_freq_last_month_invr"  
# [4] "mother_milk_HMO_3FL_ugml_invr_infant_BITSS_invr"                                
# [5] "mother_milk_HMO_LNT_ugml_invr_infant_dis_ROME_e13_stool_freq_last_month_invr"   
# [6] "mother_milk_HMO_milk_group_infant_BITSS_invr"                                   
# [7] "mother_milk_HMO_2FL_ugml_invr_infant_BITSS_invr"                                
# [8] "mother_milk_HMO_2FL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"    
# [9] "mother_milk_HMO_DFLNHa_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr" 
# [10] "mother_milk_HMO_3FL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"    
# [11] "mother_milk_HMO_3F3SL_ugml_invr_infant_BITSS_invr"                              
# [12] "mother_milk_HMO_3F3SL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"  
# [13] "mother_milk_HMO_LNFP_II_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"
# [14] "mother_milk_HMO_A_tetra_ugml_invr_infant_BITSS_invr"                            
# [15] "mother_milk_HMO_3F3SL_ugml_invr_infant_growth_length_cm_invr" 

## consistently nominally significant results at M1, M2, M3
nomsign_HMO_results_M1_M2_M3 <- intersect(nomsign_HMO_results_M1_M3, nomsign_HMO_results_M2)
# [1] "mother_milk_HMO_DFLNHa_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr" 
# [2] "mother_milk_HMO_3FL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"    
# [3] "mother_milk_HMO_3F3SL_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"  
# [4] "mother_milk_HMO_LNFP_II_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"
# [5] "mother_milk_HMO_3F3SL_ugml_invr_infant_growth_length_cm_invr"   

## consistently nominally significant results at W2, M1, M2, M3
nomsign_HMO_results_W2_M1_M2_M3 <- intersect(nomsign_HMO_results_M1_M2_M3, nomsign_HMO_results_W2)
# [1] "mother_milk_HMO_LNFP_II_ugml_invr_infant_growth_weight_kg_gain_since_birth_invr"

## consistently nominally significant results at M1, M3, M6
nomsign_HMO_results_M1_M3_M6 <- intersect(nomsign_HMO_results_M1_M3, nomsign_HMO_results_M6)
# [1] "mother_milk_HMO_A_tetra_ugml_invr_infant_BITSS_invr"

## consistently nominally significant results at M1, M2, M3, M6
nomsign_HMO_results_M1_M2_M3_M6 <- intersect(nomsign_HMO_results_M1_M3_M6, nomsign_HMO_results_M2)
# none


##### =========================== 21. COMBINE RESULTS FROM ASSOCIATIONS OF MILK MILK MICROBIOTA DATA WITH INFANT OUTCOME PHENOTYPES =========================== #####

## combine and then split by alpha diversity and relative abundances
# "inff_outcomes_by_milk_microbiota_results_M1"
# "inff_growth_by_milk_microbiota_results_M1"  
# 
# "inff_outcomes_by_milk_microbiota_results_M2"
# "inff_growth_by_milk_microbiota_results_M2"  
# 
# "inff_outcomes_by_milk_microbiota_results_M3"
# "inff_growth_by_milk_microbiota_results_M3" 
# 
# "inff_outcomes_by_milk_microbiota_results_M6"
# "inff_growth_by_milk_microbiota_results_M6"  


### ===== 21.1 COMBINE RESULTS PER TIME POINT ===== ###

### === M1 === ###

# "inff_outcomes_by_milk_microbiota_results_M1" #390x9; swap x and y
# "inff_growth_by_milk_microbiota_results_M1" #234x9; swap x and y

## swap x and y
inff_outcomes_by_milk_microbiota_results_M1$new_x <- inff_outcomes_by_milk_microbiota_results_M1$y
inff_outcomes_by_milk_microbiota_results_M1$new_y <- inff_outcomes_by_milk_microbiota_results_M1$x
inff_growth_by_milk_microbiota_results_M1$new_x <- inff_growth_by_milk_microbiota_results_M1$y
inff_growth_by_milk_microbiota_results_M1$new_y <- inff_growth_by_milk_microbiota_results_M1$x

## remove old x and y
inff_outcomes_by_milk_microbiota_results_M1 <- inff_outcomes_by_milk_microbiota_results_M1[,-c(3:4)]
inff_growth_by_milk_microbiota_results_M1 <- inff_growth_by_milk_microbiota_results_M1[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_milk_microbiota_results_M1)[8:9] <- c("x", "y")
colnames(inff_growth_by_milk_microbiota_results_M1)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_microbiota_results_M1 <- as.data.frame(rbind.fill(inff_outcomes_by_milk_microbiota_results_M1, inff_growth_by_milk_microbiota_results_M1))

## resort columns
all_outcomes_by_microbiota_results_M1 <- all_outcomes_by_microbiota_results_M1[,c(1:2,8:9,3:7)]

## check merged data frame
dim(all_outcomes_by_microbiota_results_M1) #624
unique(all_outcomes_by_microbiota_results_M1$x) #39 (3+36)
unique(all_outcomes_by_microbiota_results_M1$y) #16

## split data frame into 1 for alpha diversity and 1 for relative abundances
all_outcomes_by_alpha_diversity_results_M1 <- all_outcomes_by_microbiota_results_M1[grep("alpha_div", all_outcomes_by_microbiota_results_M1$x),] #48x9
all_outcomes_by_relative_abundances_results_M1 <- all_outcomes_by_microbiota_results_M1[grep("g__", all_outcomes_by_microbiota_results_M1$x),] #576x9

## recalculate FDR
all_outcomes_by_alpha_diversity_results_M1$FDR <- p.adjust(all_outcomes_by_alpha_diversity_results_M1$p, method="BH")
all_outcomes_by_relative_abundances_results_M1$FDR <- p.adjust(all_outcomes_by_relative_abundances_results_M1$p, method="BH")

## resort by FDR and p
all_outcomes_by_alpha_diversity_results_M1 <- all_outcomes_by_alpha_diversity_results_M1[order(all_outcomes_by_alpha_diversity_results_M1$FDR, all_outcomes_by_alpha_diversity_results_M1$p),]
all_outcomes_by_relative_abundances_results_M1 <- all_outcomes_by_relative_abundances_results_M1[order(all_outcomes_by_relative_abundances_results_M1$FDR, all_outcomes_by_relative_abundances_results_M1$p),]

## check results
nrow(all_outcomes_by_alpha_diversity_results_M1[all_outcomes_by_alpha_diversity_results_M1$FDR<0.05,]) #0 FDR significant results
nrow(all_outcomes_by_alpha_diversity_results_M1[all_outcomes_by_alpha_diversity_results_M1$p<0.05,])   #1 nominally significant results

nrow(all_outcomes_by_relative_abundances_results_M1[all_outcomes_by_relative_abundances_results_M1$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_relative_abundances_results_M1[all_outcomes_by_relative_abundances_results_M1$p<0.05,])   #29 nominally significant results

## save results
write.table(all_outcomes_by_alpha_diversity_results_M1, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_relative_abundances_results_M1, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# all_outcomes_by_alpha_diversity_results_M1 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M1.txt", header=T, sep="\t")
# all_outcomes_by_relative_abundances_results_M1 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M1.txt", header=T, sep="\t")


### === M2 === ###

# "inff_outcomes_by_milk_microbiota_results_M2" #390x9; swap x and y
# "inff_growth_by_milk_microbiota_results_M2" #234x9; swap x and y

## swap x and y
inff_outcomes_by_milk_microbiota_results_M2$new_x <- inff_outcomes_by_milk_microbiota_results_M2$y
inff_outcomes_by_milk_microbiota_results_M2$new_y <- inff_outcomes_by_milk_microbiota_results_M2$x
inff_growth_by_milk_microbiota_results_M2$new_x <- inff_growth_by_milk_microbiota_results_M2$y
inff_growth_by_milk_microbiota_results_M2$new_y <- inff_growth_by_milk_microbiota_results_M2$x

## remove old x and y
inff_outcomes_by_milk_microbiota_results_M2 <- inff_outcomes_by_milk_microbiota_results_M2[,-c(3:4)]
inff_growth_by_milk_microbiota_results_M2 <- inff_growth_by_milk_microbiota_results_M2[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_milk_microbiota_results_M2)[8:9] <- c("x", "y")
colnames(inff_growth_by_milk_microbiota_results_M2)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_microbiota_results_M2 <- as.data.frame(rbind.fill(inff_outcomes_by_milk_microbiota_results_M2, inff_growth_by_milk_microbiota_results_M2))

## resort columns
all_outcomes_by_microbiota_results_M2 <- all_outcomes_by_microbiota_results_M2[,c(1:2,8:9,3:7)]

## check merged data frame
dim(all_outcomes_by_microbiota_results_M2) #624
unique(all_outcomes_by_microbiota_results_M2$x) #39 (3+36)
unique(all_outcomes_by_microbiota_results_M2$y) #16

## split data frame into 1 for alpha diversity and 1 for relative abundances
all_outcomes_by_alpha_diversity_results_M2 <- all_outcomes_by_microbiota_results_M2[grep("alpha_div", all_outcomes_by_microbiota_results_M2$x),] #48x9
all_outcomes_by_relative_abundances_results_M2 <- all_outcomes_by_microbiota_results_M2[grep("g__", all_outcomes_by_microbiota_results_M2$x),] #576x9

## recalculate FDR
all_outcomes_by_alpha_diversity_results_M2$FDR <- p.adjust(all_outcomes_by_alpha_diversity_results_M2$p, method="BH")
all_outcomes_by_relative_abundances_results_M2$FDR <- p.adjust(all_outcomes_by_relative_abundances_results_M2$p, method="BH")

## resort by FDR and p
all_outcomes_by_alpha_diversity_results_M2 <- all_outcomes_by_alpha_diversity_results_M2[order(all_outcomes_by_alpha_diversity_results_M2$FDR, all_outcomes_by_alpha_diversity_results_M2$p),]
all_outcomes_by_relative_abundances_results_M2 <- all_outcomes_by_relative_abundances_results_M2[order(all_outcomes_by_relative_abundances_results_M2$FDR, all_outcomes_by_relative_abundances_results_M2$p),]

## check results
nrow(all_outcomes_by_alpha_diversity_results_M2[all_outcomes_by_alpha_diversity_results_M2$FDR<0.05,]) #0 FDR significant results
nrow(all_outcomes_by_alpha_diversity_results_M2[all_outcomes_by_alpha_diversity_results_M2$p<0.05,])   #2 nominally significant results

nrow(all_outcomes_by_relative_abundances_results_M2[all_outcomes_by_relative_abundances_results_M2$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_relative_abundances_results_M2[all_outcomes_by_relative_abundances_results_M2$p<0.05,])   #32 nominally significant results

## save results
write.table(all_outcomes_by_alpha_diversity_results_M2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_relative_abundances_results_M2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# all_outcomes_by_alpha_diversity_results_M2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M2.txt", header=T, sep="\t")
# all_outcomes_by_relative_abundances_results_M2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M2.txt", header=T, sep="\t")


### === M3 === ###

# "inff_outcomes_by_milk_microbiota_results_M3" #1014x9; swap x and y
# "inff_growth_by_milk_microbiota_results_M3" #234x9; swap x and y

## swap x and y
inff_outcomes_by_milk_microbiota_results_M3$new_x <- inff_outcomes_by_milk_microbiota_results_M3$y
inff_outcomes_by_milk_microbiota_results_M3$new_y <- inff_outcomes_by_milk_microbiota_results_M3$x
inff_growth_by_milk_microbiota_results_M3$new_x <- inff_growth_by_milk_microbiota_results_M3$y
inff_growth_by_milk_microbiota_results_M3$new_y <- inff_growth_by_milk_microbiota_results_M3$x

## remove old x and y
inff_outcomes_by_milk_microbiota_results_M3 <- inff_outcomes_by_milk_microbiota_results_M3[,-c(3:4)]
inff_growth_by_milk_microbiota_results_M3 <- inff_growth_by_milk_microbiota_results_M3[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_milk_microbiota_results_M3)[8:9] <- c("x", "y")
colnames(inff_growth_by_milk_microbiota_results_M3)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_microbiota_results_M3 <- as.data.frame(rbind.fill(inff_outcomes_by_milk_microbiota_results_M3, inff_growth_by_milk_microbiota_results_M3))

## resort columns
all_outcomes_by_microbiota_results_M3 <- all_outcomes_by_microbiota_results_M3[,c(1:2,8:9,3:7)]

## check merged data frame
dim(all_outcomes_by_microbiota_results_M3) #1248
unique(all_outcomes_by_microbiota_results_M3$x) #39 (3+36)
unique(all_outcomes_by_microbiota_results_M3$y) #32

## split data frame into 1 for alpha diversity and 1 for relative abundances
all_outcomes_by_alpha_diversity_results_M3 <- all_outcomes_by_microbiota_results_M3[grep("alpha_div", all_outcomes_by_microbiota_results_M3$x),] #96x9
all_outcomes_by_relative_abundances_results_M3 <- all_outcomes_by_microbiota_results_M3[grep("g__", all_outcomes_by_microbiota_results_M3$x),] #1152x9

## recalculate FDR
all_outcomes_by_alpha_diversity_results_M3$FDR <- p.adjust(all_outcomes_by_alpha_diversity_results_M3$p, method="BH")
all_outcomes_by_relative_abundances_results_M3$FDR <- p.adjust(all_outcomes_by_relative_abundances_results_M3$p, method="BH")

## resort by FDR and p
all_outcomes_by_alpha_diversity_results_M3 <- all_outcomes_by_alpha_diversity_results_M3[order(all_outcomes_by_alpha_diversity_results_M3$FDR, all_outcomes_by_alpha_diversity_results_M3$p),]
all_outcomes_by_relative_abundances_results_M3 <- all_outcomes_by_relative_abundances_results_M3[order(all_outcomes_by_relative_abundances_results_M3$FDR, all_outcomes_by_relative_abundances_results_M3$p),]

## check results
nrow(all_outcomes_by_alpha_diversity_results_M3[all_outcomes_by_alpha_diversity_results_M3$FDR<0.05,]) #0 FDR significant results
nrow(all_outcomes_by_alpha_diversity_results_M3[all_outcomes_by_alpha_diversity_results_M3$p<0.05,])   #5 nominally significant results

nrow(all_outcomes_by_relative_abundances_results_M3[all_outcomes_by_relative_abundances_results_M3$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_relative_abundances_results_M3[all_outcomes_by_relative_abundances_results_M3$p<0.05,])   #73 nominally significant results

## save results
write.table(all_outcomes_by_alpha_diversity_results_M3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_relative_abundances_results_M3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# all_outcomes_by_alpha_diversity_results_M3 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M3.txt", header=T, sep="\t")
# all_outcomes_by_relative_abundances_results_M3 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M3.txt", header=T, sep="\t")


### === M6 === ###

# "inff_outcomes_by_milk_microbiota_results_M6" #507x9; swap x and y
# "inff_growth_by_milk_microbiota_results_M6" #234x9; swap x and y

## swap x and y
inff_outcomes_by_milk_microbiota_results_M6$new_x <- inff_outcomes_by_milk_microbiota_results_M6$y
inff_outcomes_by_milk_microbiota_results_M6$new_y <- inff_outcomes_by_milk_microbiota_results_M6$x
inff_growth_by_milk_microbiota_results_M6$new_x <- inff_growth_by_milk_microbiota_results_M6$y
inff_growth_by_milk_microbiota_results_M6$new_y <- inff_growth_by_milk_microbiota_results_M6$x

## remove old x and y
inff_outcomes_by_milk_microbiota_results_M6 <- inff_outcomes_by_milk_microbiota_results_M6[,-c(3:4)]
inff_growth_by_milk_microbiota_results_M6 <- inff_growth_by_milk_microbiota_results_M6[,-c(3:4)]

## rename new_x and new_y to x and y
colnames(inff_outcomes_by_milk_microbiota_results_M6)[8:9] <- c("x", "y")
colnames(inff_growth_by_milk_microbiota_results_M6)[8:9] <- c("x", "y")

## merge result data frames
all_outcomes_by_microbiota_results_M6 <- as.data.frame(rbind.fill(inff_outcomes_by_milk_microbiota_results_M6, inff_growth_by_milk_microbiota_results_M6))

## resort columns
all_outcomes_by_microbiota_results_M6 <- all_outcomes_by_microbiota_results_M6[,c(1:2,8:9,3:7)]

## check merged data frame
dim(all_outcomes_by_microbiota_results_M6) #741x9
unique(all_outcomes_by_microbiota_results_M6$x) #39 (3+36)
unique(all_outcomes_by_microbiota_results_M6$y) #19

## split data frame into 1 for alpha diversity and 1 for relative abundances
all_outcomes_by_alpha_diversity_results_M6 <- all_outcomes_by_microbiota_results_M6[grep("alpha_div", all_outcomes_by_microbiota_results_M6$x),] #57x9
all_outcomes_by_relative_abundances_results_M6 <- all_outcomes_by_microbiota_results_M6[grep("g__", all_outcomes_by_microbiota_results_M6$x),] #684x9

## recalculate FDR
all_outcomes_by_alpha_diversity_results_M6$FDR <- p.adjust(all_outcomes_by_alpha_diversity_results_M6$p, method="BH")
all_outcomes_by_relative_abundances_results_M6$FDR <- p.adjust(all_outcomes_by_relative_abundances_results_M6$p, method="BH")

## resort by FDR and p
all_outcomes_by_alpha_diversity_results_M6 <- all_outcomes_by_alpha_diversity_results_M6[order(all_outcomes_by_alpha_diversity_results_M6$FDR, all_outcomes_by_alpha_diversity_results_M6$p),]
all_outcomes_by_relative_abundances_results_M6 <- all_outcomes_by_relative_abundances_results_M6[order(all_outcomes_by_relative_abundances_results_M6$FDR, all_outcomes_by_relative_abundances_results_M6$p),]

## check results
nrow(all_outcomes_by_alpha_diversity_results_M6[all_outcomes_by_alpha_diversity_results_M6$FDR<0.05,]) #0 FDR significant results
nrow(all_outcomes_by_alpha_diversity_results_M6[all_outcomes_by_alpha_diversity_results_M6$p<0.05,])   #3 nominally significant results

nrow(all_outcomes_by_relative_abundances_results_M6[all_outcomes_by_relative_abundances_results_M6$FDR<0.05,]) # 0 FDR significant results
nrow(all_outcomes_by_relative_abundances_results_M6[all_outcomes_by_relative_abundances_results_M6$p<0.05,])   #32 nominally significant results

## save results
write.table(all_outcomes_by_alpha_diversity_results_M6, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_relative_abundances_results_M6, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# all_outcomes_by_alpha_diversity_results_M6 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M6.txt", header=T, sep="\t")
# all_outcomes_by_relative_abundances_results_M6 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M6.txt", header=T, sep="\t")


# ### ===== 21.2 MILK ALPHA DIVERSITY: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS ===== ###
# 
# ## check for consistently nominally significant associations at several time points
# nomsign_alpha_results_M1 <- c(paste0(all_outcomes_by_alpha_diversity_results_M1[all_outcomes_by_alpha_diversity_results_M1$p<0.05, "x"], "_", all_outcomes_by_alpha_diversity_results_M1[all_outcomes_by_alpha_diversity_results_M1$p<0.05, "y"]))
# nomsign_alpha_results_M2 <- c(paste0(all_outcomes_by_alpha_diversity_results_M2[all_outcomes_by_alpha_diversity_results_M2$p<0.05, "x"], "_", all_outcomes_by_alpha_diversity_results_M2[all_outcomes_by_alpha_diversity_results_M2$p<0.05, "y"]))
# nomsign_alpha_results_M3 <- c(paste0(all_outcomes_by_alpha_diversity_results_M3[all_outcomes_by_alpha_diversity_results_M3$p<0.05, "x"], "_", all_outcomes_by_alpha_diversity_results_M3[all_outcomes_by_alpha_diversity_results_M3$p<0.05, "y"]))
# nomsign_alpha_results_M6 <- c(paste0(all_outcomes_by_alpha_diversity_results_M6[all_outcomes_by_alpha_diversity_results_M6$p<0.05, "x"], "_", all_outcomes_by_alpha_diversity_results_M6[all_outcomes_by_alpha_diversity_results_M6$p<0.05, "y"]))
# 
# # nomsign_alpha_results_M1
# # [1] "_"
# # 
# # > nomsign_alpha_results_M2
# # [1] "alpha_div_simpson_invr_infant_dis_atopic_dermatitis_SCORAD_M12_invr"
# # [2] "alpha_div_shannon_invr_infant_dis_atopic_dermatitis_SCORAD_M12_invr"
# # 
# # > nomsign_alpha_results_M3
# # [1] "alpha_div_simpson_invr_infant_cry_prev_week_average_crying_time_per_day_min_invr"         
# # [2] "alpha_div_genera_richness_invr_infant_igsq_3_freq_spitting_milk_normal_day_past_week_invr"
# # [3] "alpha_div_shannon_invr_infant_cry_prev_week_average_crying_time_per_day_min_invr"         
# # [4] "alpha_div_shannon_invr_infant_dis_ROME_e14_stool_structure_invr"                          
# # [5] "alpha_div_simpson_invr_infant_dis_ROME_e14_stool_structure_invr" 
# # 
# # > nomsign_alpha_results_M6
# # [1] "alpha_div_genera_richness_invr_infant_sleep_day_hours_invr"     
# # [2] "alpha_div_simpson_invr_infant_growth_head_circumference_cm_invr"
# # [3] "alpha_div_shannon_invr_infant_BITSS_invr"
# 
# # -> No nominally signifciant associations were observed at ≥1 time pint - or they were not tested at ≥1 time point.


# ### ===== 21.3 MILK RELATIVE ABUNDANCES: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS ===== ###
# 
# ## check for consistently nominally significant associations at several time points
# nomsign_RA_results_M1 <- c(paste0(all_outcomes_by_relative_abundances_results_M1[all_outcomes_by_relative_abundances_results_M1$p<0.05, "x"], "_", all_outcomes_by_relative_abundances_results_M1[all_outcomes_by_relative_abundances_results_M1$p<0.05, "y"]))
# nomsign_RA_results_M2 <- c(paste0(all_outcomes_by_relative_abundances_results_M2[all_outcomes_by_relative_abundances_results_M2$p<0.05, "x"], "_", all_outcomes_by_relative_abundances_results_M2[all_outcomes_by_relative_abundances_results_M2$p<0.05, "y"]))
# nomsign_RA_results_M3 <- c(paste0(all_outcomes_by_relative_abundances_results_M3[all_outcomes_by_relative_abundances_results_M3$p<0.05, "x"], "_", all_outcomes_by_relative_abundances_results_M3[all_outcomes_by_relative_abundances_results_M3$p<0.05, "y"]))
# nomsign_RA_results_M6 <- c(paste0(all_outcomes_by_relative_abundances_results_M6[all_outcomes_by_relative_abundances_results_M6$p<0.05, "x"], "_", all_outcomes_by_relative_abundances_results_M6[all_outcomes_by_relative_abundances_results_M6$p<0.05, "y"]))
# 
# 
# ## consistently nominally significant results at the 2 largest time points (M1, M3)
# nomsign_RA_results_M1_M3 <- intersect(nomsign_RA_results_M1, nomsign_RA_results_M3)
# nomsign_RA_results_M1_M3
# # [1] "g__Escherichia_Shigella_clr_infant_dis_ROME_e25_mucus_stool_last_week_invr" 
# all_outcomes_by_relative_abundances_results_M1[all_outcomes_by_relative_abundances_results_M1$x=="g__Escherichia_Shigella_clr" & all_outcomes_by_relative_abundances_results_M1$y=="infant_dis_ROME_e25_mucus_stool_last_week_invr",]
# all_outcomes_by_relative_abundances_results_M2[all_outcomes_by_relative_abundances_results_M2$x=="g__Escherichia_Shigella_clr" & all_outcomes_by_relative_abundances_results_M2$y=="infant_dis_ROME_e25_mucus_stool_last_week_invr",]
# all_outcomes_by_relative_abundances_results_M3[all_outcomes_by_relative_abundances_results_M3$x=="g__Escherichia_Shigella_clr" & all_outcomes_by_relative_abundances_results_M3$y=="infant_dis_ROME_e25_mucus_stool_last_week_invr",]
# all_outcomes_by_relative_abundances_results_M6[all_outcomes_by_relative_abundances_results_M6$x=="g__Escherichia_Shigella_clr" & all_outcomes_by_relative_abundances_results_M6$y=="infant_dis_ROME_e25_mucus_stool_last_week_invr",]
# # -> higher RA of Escherichia-Shigella were associated with mucus stool
# 
# # [2] "g__Actinomyces_clr_infant_cry_prev_week_average_crying_time_per_day_min_invr"
# all_outcomes_by_relative_abundances_results_M1[all_outcomes_by_relative_abundances_results_M1$x=="g__Actinomyces_clr" & all_outcomes_by_relative_abundances_results_M1$y=="infant_cry_prev_week_average_crying_time_per_day_min_invr",]
# all_outcomes_by_relative_abundances_results_M2[all_outcomes_by_relative_abundances_results_M2$x=="g__Actinomyces_clr" & all_outcomes_by_relative_abundances_results_M2$y=="infant_cry_prev_week_average_crying_time_per_day_min_invr",]
# all_outcomes_by_relative_abundances_results_M3[all_outcomes_by_relative_abundances_results_M3$x=="g__Actinomyces_clr" & all_outcomes_by_relative_abundances_results_M3$y=="infant_cry_prev_week_average_crying_time_per_day_min_invr",]
# all_outcomes_by_relative_abundances_results_M6[all_outcomes_by_relative_abundances_results_M6$x=="g__Actinomyces_clr" & all_outcomes_by_relative_abundances_results_M6$y=="infant_cry_prev_week_average_crying_time_per_day_min_invr",]
# # -> Actinomyces RA was inconsistently (higher/lower) associated with more infant crying time
# 
# # [3] "g__Enhydrobacter_clr_infant_cry_crying_more_than_3h_per_day" 
# all_outcomes_by_relative_abundances_results_M1[all_outcomes_by_relative_abundances_results_M1$x=="g__Enhydrobacter_clr" & all_outcomes_by_relative_abundances_results_M1$y=="infant_cry_crying_more_than_3h_per_day",]
# all_outcomes_by_relative_abundances_results_M2[all_outcomes_by_relative_abundances_results_M2$x=="g__Enhydrobacter_clr" & all_outcomes_by_relative_abundances_results_M2$y=="infant_cry_crying_more_than_3h_per_day",]
# all_outcomes_by_relative_abundances_results_M3[all_outcomes_by_relative_abundances_results_M3$x=="g__Enhydrobacter_clr" & all_outcomes_by_relative_abundances_results_M3$y=="infant_cry_crying_more_than_3h_per_day",]
# all_outcomes_by_relative_abundances_results_M6[all_outcomes_by_relative_abundances_results_M6$x=="g__Enhydrobacter_clr" & all_outcomes_by_relative_abundances_results_M6$y=="infant_cry_crying_more_than_3h_per_day",]
# # -> Enhydrobacter RA was inconsistently (lower/higher) associated with more infant crying
# 
# 
# ## consistently nominally significant results at the 3 largest time points (M1, M3, M6)
# nomsign_RA_results_M1_M3_M6 <- intersect(nomsign_RA_results_M1_M3, nomsign_RA_results_M6)
# nomsign_RA_results_M1_M3_M6
# #no nominally significant associations are found at M1, M3 and M6
# 
# ## consistently nominally significant results at the 3 first time points (M1, M2, M3)
# nomsign_RA_results_M1_M2_M3 <- intersect(nomsign_RA_results_M1_M3, nomsign_RA_results_M2)
# nomsign_RA_results_M1_M2_M3
# #no nominally significant associations are found at M1, M2 and M3


### ===== 21.4 MILK ALPHA DIVERSITY: MORE FOCUSED ANALYSIS: ONLY USE SELECTED TIME POINTS AND PHENOTYPES AND RECALCULATE FDR (TABLE S38A-C) ===== ###

## if necessary, re-import results files
all_outcomes_by_alpha_diversity_results_M1 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M1.txt", header=T, sep="\t")
all_outcomes_by_relative_abundances_results_M1 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M1.txt", header=T, sep="\t")
all_outcomes_by_alpha_diversity_results_M2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M2.txt", header=T, sep="\t")
all_outcomes_by_relative_abundances_results_M2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M2.txt", header=T, sep="\t")
all_outcomes_by_alpha_diversity_results_M3 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M3.txt", header=T, sep="\t")
all_outcomes_by_relative_abundances_results_M3 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M3.txt", header=T, sep="\t")
all_outcomes_by_alpha_diversity_results_M6 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M6.txt", header=T, sep="\t")
all_outcomes_by_relative_abundances_results_M6 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M6.txt", header=T, sep="\t")

## final selection of 11 phenotypes
sel_phenos_microbes <- c("infant_dis_ROME_e1_regurgitation","infant_dis_ROME_e32_sudden_inconsolable_crying_fits",
                     "infant_growth_weight_kg_invr","infant_growth_weight_kg_gain_since_birth_invr",
                     "infant_growth_length_cm_invr","infant_growth_length_cm_gain_since_birth_invr",
                     "infant_cry_prev_week_average_crying_time_per_day_min_invr",
                     "infant_BITSS_invr",
                     "infant_dis_ROME_e13_stool_freq_last_month_invr","infant_dis_ROME_e14_stool_structure_invr","infant_dis_ROME_e25_mucus_stool_last_week_invr")


### ===== M1 (TABLE S38A) ===== ###

## select data of interest
all_outcomes_by_alpha_diversity_results_M1_v2 <- all_outcomes_by_alpha_diversity_results_M1[all_outcomes_by_alpha_diversity_results_M1$y %in% sel_phenos_microbes,]

## exclude results for simpson
all_outcomes_by_alpha_diversity_results_M1_v2 <- all_outcomes_by_alpha_diversity_results_M1_v2[all_outcomes_by_alpha_diversity_results_M1_v2$x!="alpha_div_simpson_invr",]

## recalculate FDR
all_outcomes_by_alpha_diversity_results_M1_v2$FDR <- p.adjust(all_outcomes_by_alpha_diversity_results_M1_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_alpha_diversity_results_M1_v2 <- all_outcomes_by_alpha_diversity_results_M1_v2[order(all_outcomes_by_alpha_diversity_results_M1_v2$FDR, all_outcomes_by_alpha_diversity_results_M1_v2$p),]

## add column showing phenotype group
all_outcomes_by_alpha_diversity_results_M1_v2$Phenotype_group <- NA
all_outcomes_by_alpha_diversity_results_M1_v2[grep("ROME", all_outcomes_by_alpha_diversity_results_M1_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_alpha_diversity_results_M1_v2[grep("BITSS", all_outcomes_by_alpha_diversity_results_M1_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_alpha_diversity_results_M1_v2[grep("growth", all_outcomes_by_alpha_diversity_results_M1_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_alpha_diversity_results_M1_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_alpha_diversity_results_M1_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column showing reference level
all_outcomes_by_alpha_diversity_results_M1_v2[all_outcomes_by_alpha_diversity_results_M1_v2$levels=="" & !is.na(all_outcomes_by_alpha_diversity_results_M1_v2$levels),"levels"] <- NA
all_outcomes_by_alpha_diversity_results_M1_v2$Reference_level <- NA
all_outcomes_by_alpha_diversity_results_M1_v2[all_outcomes_by_alpha_diversity_results_M1_v2$levels=="yes" &!is.na(all_outcomes_by_alpha_diversity_results_M1_v2$levels), "Reference_level"] <- "no"

## resort columns
all_outcomes_by_alpha_diversity_results_M1_v3 <- all_outcomes_by_alpha_diversity_results_M1_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_alpha_diversity_results_M1_v3)

# write.table(all_outcomes_by_alpha_diversity_results_M1_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_alpha_diversity_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M1_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_alpha_diversity_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M1_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_alpha_diversity_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M1_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_alpha_diversity_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_alpha_diversity_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_alpha_diversity_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)

# -> 1 nominally significant association between milk alpha diversity and infant outcomes at M1.
# alpha_div_genera_richness_invr ~ infant_dis_ROME_e14_stool_structure_invr (increase)


### ===== M2 ===== ###

## select data of interest
all_outcomes_by_alpha_diversity_results_M2_v2 <- all_outcomes_by_alpha_diversity_results_M2[all_outcomes_by_alpha_diversity_results_M2$y %in% sel_phenos_microbes,]

## exclude results for simpson
all_outcomes_by_alpha_diversity_results_M2_v2 <- all_outcomes_by_alpha_diversity_results_M2_v2[all_outcomes_by_alpha_diversity_results_M2_v2$x!="alpha_div_simpson_invr",]

## recalculate FDR
all_outcomes_by_alpha_diversity_results_M2_v2$FDR <- p.adjust(all_outcomes_by_alpha_diversity_results_M2_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_alpha_diversity_results_M2_v2 <- all_outcomes_by_alpha_diversity_results_M2_v2[order(all_outcomes_by_alpha_diversity_results_M2_v2$FDR, all_outcomes_by_alpha_diversity_results_M2_v2$p),]

## add column showing phenotype group
all_outcomes_by_alpha_diversity_results_M2_v2$Phenotype_group <- NA
all_outcomes_by_alpha_diversity_results_M2_v2[grep("ROME", all_outcomes_by_alpha_diversity_results_M2_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_alpha_diversity_results_M2_v2[grep("BITSS", all_outcomes_by_alpha_diversity_results_M2_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_alpha_diversity_results_M2_v2[grep("growth", all_outcomes_by_alpha_diversity_results_M2_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_alpha_diversity_results_M2_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_alpha_diversity_results_M2_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column showing reference level
all_outcomes_by_alpha_diversity_results_M2_v2[all_outcomes_by_alpha_diversity_results_M2_v2$levels=="" & !is.na(all_outcomes_by_alpha_diversity_results_M2_v2$levels),"levels"] <- NA
all_outcomes_by_alpha_diversity_results_M2_v2$Reference_level <- NA
all_outcomes_by_alpha_diversity_results_M2_v2[all_outcomes_by_alpha_diversity_results_M2_v2$levels=="yes" &!is.na(all_outcomes_by_alpha_diversity_results_M2_v2$levels), "Reference_level"] <- "no"

## resort columns
all_outcomes_by_alpha_diversity_results_M2_v3 <- all_outcomes_by_alpha_diversity_results_M2_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_alpha_diversity_results_M2_v3)

# write.table(all_outcomes_by_alpha_diversity_results_M2_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_alpha_diversity_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M2_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_alpha_diversity_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M2_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_alpha_diversity_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M2_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_alpha_diversity_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_alpha_diversity_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_alpha_diversity_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)

# -> No nominally significant associations between milk alpha diversity and infant outcomes at M2.


### ===== M3 (TABLE S38B) ===== ###

## select data of interest
all_outcomes_by_alpha_diversity_results_M3_v2 <- all_outcomes_by_alpha_diversity_results_M3[all_outcomes_by_alpha_diversity_results_M3$y %in% sel_phenos_microbes,]

## exclude results for simpson
all_outcomes_by_alpha_diversity_results_M3_v2 <- all_outcomes_by_alpha_diversity_results_M3_v2[all_outcomes_by_alpha_diversity_results_M3_v2$x!="alpha_div_simpson_invr",]

## recalculate FDR
all_outcomes_by_alpha_diversity_results_M3_v2$FDR <- p.adjust(all_outcomes_by_alpha_diversity_results_M3_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_alpha_diversity_results_M3_v2 <- all_outcomes_by_alpha_diversity_results_M3_v2[order(all_outcomes_by_alpha_diversity_results_M3_v2$FDR, all_outcomes_by_alpha_diversity_results_M3_v2$p),]

## add column showing phenotype group
all_outcomes_by_alpha_diversity_results_M3_v2$Phenotype_group <- NA
all_outcomes_by_alpha_diversity_results_M3_v2[grep("ROME", all_outcomes_by_alpha_diversity_results_M3_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_alpha_diversity_results_M3_v2[grep("BITSS", all_outcomes_by_alpha_diversity_results_M3_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_alpha_diversity_results_M3_v2[grep("growth", all_outcomes_by_alpha_diversity_results_M3_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_alpha_diversity_results_M3_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_alpha_diversity_results_M3_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column showing reference level
all_outcomes_by_alpha_diversity_results_M3_v2[all_outcomes_by_alpha_diversity_results_M3_v2$levels=="" & !is.na(all_outcomes_by_alpha_diversity_results_M3_v2$levels),"levels"] <- NA
all_outcomes_by_alpha_diversity_results_M3_v2$Reference_level <- NA
all_outcomes_by_alpha_diversity_results_M3_v2[all_outcomes_by_alpha_diversity_results_M3_v2$levels=="yes" &!is.na(all_outcomes_by_alpha_diversity_results_M3_v2$levels), "Reference_level"] <- "no"

## resort columns
all_outcomes_by_alpha_diversity_results_M3_v3 <- all_outcomes_by_alpha_diversity_results_M3_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_alpha_diversity_results_M3_v3)

# write.table(all_outcomes_by_alpha_diversity_results_M3_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_alpha_diversity_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M3_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_alpha_diversity_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M3_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_alpha_diversity_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M3_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_alpha_diversity_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_alpha_diversity_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_alpha_diversity_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)

# -> 2 nominally significant associations between milk alpha diversity and infant outcomes at M3.
# Shannon ~ crying time (increase)
# Shannon ~ ROME e14 stool structure (decrease)


### ===== M6 (TABLE S38C) ===== ###

## select data of interest
all_outcomes_by_alpha_diversity_results_M6_v2 <- all_outcomes_by_alpha_diversity_results_M6[all_outcomes_by_alpha_diversity_results_M6$y %in% sel_phenos_microbes,]

## exclude results for simpson
all_outcomes_by_alpha_diversity_results_M6_v2 <- all_outcomes_by_alpha_diversity_results_M6_v2[all_outcomes_by_alpha_diversity_results_M6_v2$x!="alpha_div_simpson_invr",]

## recalculate FDR
all_outcomes_by_alpha_diversity_results_M6_v2$FDR <- p.adjust(all_outcomes_by_alpha_diversity_results_M6_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_alpha_diversity_results_M6_v2 <- all_outcomes_by_alpha_diversity_results_M6_v2[order(all_outcomes_by_alpha_diversity_results_M6_v2$FDR, all_outcomes_by_alpha_diversity_results_M6_v2$p),]

## add column showing phenotype group
all_outcomes_by_alpha_diversity_results_M6_v2$Phenotype_group <- NA
all_outcomes_by_alpha_diversity_results_M6_v2[grep("ROME", all_outcomes_by_alpha_diversity_results_M6_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_alpha_diversity_results_M6_v2[grep("BITSS", all_outcomes_by_alpha_diversity_results_M6_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_alpha_diversity_results_M6_v2[grep("growth", all_outcomes_by_alpha_diversity_results_M6_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_alpha_diversity_results_M6_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_alpha_diversity_results_M6_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column showing reference level
all_outcomes_by_alpha_diversity_results_M6_v2[all_outcomes_by_alpha_diversity_results_M6_v2$levels=="" & !is.na(all_outcomes_by_alpha_diversity_results_M6_v2$levels),"levels"] <- NA
all_outcomes_by_alpha_diversity_results_M6_v2$Reference_level <- NA
all_outcomes_by_alpha_diversity_results_M6_v2[all_outcomes_by_alpha_diversity_results_M6_v2$levels=="yes" &!is.na(all_outcomes_by_alpha_diversity_results_M6_v2$levels), "Reference_level"] <- "no"

## resort columns
all_outcomes_by_alpha_diversity_results_M6_v3 <- all_outcomes_by_alpha_diversity_results_M6_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_alpha_diversity_results_M6_v3)

# write.table(all_outcomes_by_alpha_diversity_results_M6_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_alpha_diversity_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M6_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_alpha_diversity_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M6_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_alpha_diversity_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_alpha_diversity_results_M6_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_alpha_diversity_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_alpha_diversity_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_alpha_diversity_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)

# -> 1 nominally significant associations between milk alpha diversity and infant outcomes at M6.
# Shannon ~ BITTS (decrease)


### ===== 21.5 MILK ALPHA DIVERSITY: MORE FOCUSED ANALYSIS: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS ===== ###

## check for consistently nominally significant associations at several time points
nomsign_alpha_results_M1 <- c(paste0(all_outcomes_by_alpha_diversity_results_M1_v3[all_outcomes_by_alpha_diversity_results_M1_v3$p<0.05, "x"], "_", all_outcomes_by_alpha_diversity_results_M1_v3[all_outcomes_by_alpha_diversity_results_M1_v3$p<0.05, "y"]))
# nomsign_alpha_results_M2 <- c(paste0(all_outcomes_by_alpha_diversity_results_M2_v3[all_outcomes_by_alpha_diversity_results_M2_v3$p<0.05, "x"], "_", all_outcomes_by_alpha_diversity_results_M2_v3[all_outcomes_by_alpha_diversity_results_M2_v3$p<0.05, "y"]))
nomsign_alpha_results_M3 <- c(paste0(all_outcomes_by_alpha_diversity_results_M3_v3[all_outcomes_by_alpha_diversity_results_M3_v3$p<0.05, "x"], "_", all_outcomes_by_alpha_diversity_results_M3_v3[all_outcomes_by_alpha_diversity_results_M3_v3$p<0.05, "y"]))
nomsign_alpha_results_M6 <- c(paste0(all_outcomes_by_alpha_diversity_results_M6_v3[all_outcomes_by_alpha_diversity_results_M6_v3$p<0.05, "x"], "_", all_outcomes_by_alpha_diversity_results_M6_v3[all_outcomes_by_alpha_diversity_results_M6_v3$p<0.05, "y"]))


## consistently nominally significant results at M1, M3
nomsign_alpha_results_M1_M3 <- intersect(nomsign_alpha_results_M1, nomsign_alpha_results_M3)
# none

## consistently nominally significant results at M1, M6
nomsign_alpha_results_M1_M6 <- intersect(nomsign_alpha_results_M1, nomsign_alpha_results_M6)
# none

## consistently nominally significant results at M3, M6
nomsign_alpha_results_M3_M6 <- intersect(nomsign_alpha_results_M3, nomsign_alpha_results_M6)
# none


### ===== 21.6 MILK RELATIVE ABUNDANCES: MORE FOCUSED ANALYSIS: ONLY USE SELECTED TIME POINTS AND PHENOTYPES AND RECALCULATE FDR (TABLE S39A-C) ===== ###

## if necessary, re-import results files
all_outcomes_by_alpha_diversity_results_M1 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M1.txt", header=T, sep="\t")
all_outcomes_by_relative_abundances_results_M1 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M1.txt", header=T, sep="\t")
all_outcomes_by_alpha_diversity_results_M2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M2.txt", header=T, sep="\t")
all_outcomes_by_relative_abundances_results_M2 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M2.txt", header=T, sep="\t")
all_outcomes_by_alpha_diversity_results_M3 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M3.txt", header=T, sep="\t")
all_outcomes_by_relative_abundances_results_M3 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M3.txt", header=T, sep="\t")
all_outcomes_by_alpha_diversity_results_M6 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_alpha_diversity_M6.txt", header=T, sep="\t")
all_outcomes_by_relative_abundances_results_M6 <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_infant_outcomes_by_milk_relative_abundances_M6.txt", header=T, sep="\t")

## final selection of 11 phenotypes
sel_phenos_microbes <- c("infant_dis_ROME_e1_regurgitation","infant_dis_ROME_e32_sudden_inconsolable_crying_fits",
                     "infant_growth_weight_kg_invr","infant_growth_weight_kg_gain_since_birth_invr",
                     "infant_growth_length_cm_invr","infant_growth_length_cm_gain_since_birth_invr",
                     "infant_cry_prev_week_average_crying_time_per_day_min_invr",
                     "infant_BITSS_invr",
                     "infant_dis_ROME_e13_stool_freq_last_month_invr","infant_dis_ROME_e14_stool_structure_invr","infant_dis_ROME_e25_mucus_stool_last_week_invr")


### ===== M1 (TABLE S39A) ===== ###

## select data of interest
all_outcomes_by_relative_abundances_results_M1_v2 <- all_outcomes_by_relative_abundances_results_M1[all_outcomes_by_relative_abundances_results_M1$y %in% sel_phenos_microbes,]

## recalculate FDR
all_outcomes_by_relative_abundances_results_M1_v2$FDR <- p.adjust(all_outcomes_by_relative_abundances_results_M1_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_relative_abundances_results_M1_v2 <- all_outcomes_by_relative_abundances_results_M1_v2[order(all_outcomes_by_relative_abundances_results_M1_v2$FDR, all_outcomes_by_relative_abundances_results_M1_v2$p),]

## add column showing phenotype group
all_outcomes_by_relative_abundances_results_M1_v2$Phenotype_group <- NA
all_outcomes_by_relative_abundances_results_M1_v2[grep("ROME", all_outcomes_by_relative_abundances_results_M1_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_relative_abundances_results_M1_v2[grep("BITSS", all_outcomes_by_relative_abundances_results_M1_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_relative_abundances_results_M1_v2[grep("growth", all_outcomes_by_relative_abundances_results_M1_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_relative_abundances_results_M1_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_relative_abundances_results_M1_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column showing reference level
all_outcomes_by_relative_abundances_results_M1_v2[all_outcomes_by_relative_abundances_results_M1_v2$levels=="" & !is.na(all_outcomes_by_relative_abundances_results_M1_v2$levels),"levels"] <- NA
all_outcomes_by_relative_abundances_results_M1_v2$Reference_level <- NA
all_outcomes_by_relative_abundances_results_M1_v2[all_outcomes_by_relative_abundances_results_M1_v2$levels=="yes" &!is.na(all_outcomes_by_relative_abundances_results_M1_v2$levels), "Reference_level"] <- "no"

## resort columns
all_outcomes_by_relative_abundances_results_M1_v3 <- all_outcomes_by_relative_abundances_results_M1_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_relative_abundances_results_M1_v3)

# write.table(all_outcomes_by_relative_abundances_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_relative_abundances_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_relative_abundances_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_relative_abundances_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_relative_abundances_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_relative_abundances_results_M1_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_relative_abundances_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M2 ===== ###

## select data of interest
all_outcomes_by_relative_abundances_results_M2_v2 <- all_outcomes_by_relative_abundances_results_M2[all_outcomes_by_relative_abundances_results_M2$y %in% sel_phenos_microbes,]

## recalculate FDR
all_outcomes_by_relative_abundances_results_M2_v2$FDR <- p.adjust(all_outcomes_by_relative_abundances_results_M2_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_relative_abundances_results_M2_v2 <- all_outcomes_by_relative_abundances_results_M2_v2[order(all_outcomes_by_relative_abundances_results_M2_v2$FDR, all_outcomes_by_relative_abundances_results_M2_v2$p),]

## add column showing phenotype group
all_outcomes_by_relative_abundances_results_M2_v2$Phenotype_group <- NA
all_outcomes_by_relative_abundances_results_M2_v2[grep("ROME", all_outcomes_by_relative_abundances_results_M2_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_relative_abundances_results_M2_v2[grep("BITSS", all_outcomes_by_relative_abundances_results_M2_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_relative_abundances_results_M2_v2[grep("growth", all_outcomes_by_relative_abundances_results_M2_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_relative_abundances_results_M2_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_relative_abundances_results_M2_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column showing reference level
all_outcomes_by_relative_abundances_results_M2_v2[all_outcomes_by_relative_abundances_results_M2_v2$levels=="" & !is.na(all_outcomes_by_relative_abundances_results_M2_v2$levels),"levels"] <- NA
all_outcomes_by_relative_abundances_results_M2_v2$Reference_level <- NA
all_outcomes_by_relative_abundances_results_M2_v2[all_outcomes_by_relative_abundances_results_M2_v2$levels=="yes" &!is.na(all_outcomes_by_relative_abundances_results_M2_v2$levels), "Reference_level"] <- "no"

## resort columns
all_outcomes_by_relative_abundances_results_M2_v3 <- all_outcomes_by_relative_abundances_results_M2_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_relative_abundances_results_M2_v3)

# write.table(all_outcomes_by_relative_abundances_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_relative_abundances_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_relative_abundances_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_relative_abundances_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_relative_abundances_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_relative_abundances_results_M2_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_relative_abundances_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M3 (TABLE S39B) ===== ###

## select data of interest
all_outcomes_by_relative_abundances_results_M3_v2 <- all_outcomes_by_relative_abundances_results_M3[all_outcomes_by_relative_abundances_results_M3$y %in% sel_phenos_microbes,]

## recalculate FDR
all_outcomes_by_relative_abundances_results_M3_v2$FDR <- p.adjust(all_outcomes_by_relative_abundances_results_M3_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_relative_abundances_results_M3_v2 <- all_outcomes_by_relative_abundances_results_M3_v2[order(all_outcomes_by_relative_abundances_results_M3_v2$FDR, all_outcomes_by_relative_abundances_results_M3_v2$p),]

## add column showing phenotype group
all_outcomes_by_relative_abundances_results_M3_v2$Phenotype_group <- NA
all_outcomes_by_relative_abundances_results_M3_v2[grep("ROME", all_outcomes_by_relative_abundances_results_M3_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_relative_abundances_results_M3_v2[grep("BITSS", all_outcomes_by_relative_abundances_results_M3_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_relative_abundances_results_M3_v2[grep("growth", all_outcomes_by_relative_abundances_results_M3_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_relative_abundances_results_M3_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_relative_abundances_results_M3_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column showing reference level
all_outcomes_by_relative_abundances_results_M3_v2[all_outcomes_by_relative_abundances_results_M3_v2$levels=="" & !is.na(all_outcomes_by_relative_abundances_results_M3_v2$levels),"levels"] <- NA
all_outcomes_by_relative_abundances_results_M3_v2$Reference_level <- NA
all_outcomes_by_relative_abundances_results_M3_v2[all_outcomes_by_relative_abundances_results_M3_v2$levels=="yes" &!is.na(all_outcomes_by_relative_abundances_results_M3_v2$levels), "Reference_level"] <- "no"

## resort columns
all_outcomes_by_relative_abundances_results_M3_v3 <- all_outcomes_by_relative_abundances_results_M3_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_relative_abundances_results_M3_v3)

# write.table(all_outcomes_by_relative_abundances_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_relative_abundances_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_relative_abundances_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_relative_abundances_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_relative_abundances_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_relative_abundances_results_M3_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_relative_abundances_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== M6 (TABLE S39C) ===== ###

## select data of interest
all_outcomes_by_relative_abundances_results_M6_v2 <- all_outcomes_by_relative_abundances_results_M6[all_outcomes_by_relative_abundances_results_M6$y %in% sel_phenos_microbes,]

## recalculate FDR
all_outcomes_by_relative_abundances_results_M6_v2$FDR <- p.adjust(all_outcomes_by_relative_abundances_results_M6_v2$p, method="BH")

## sort results by FDR and p
all_outcomes_by_relative_abundances_results_M6_v2 <- all_outcomes_by_relative_abundances_results_M6_v2[order(all_outcomes_by_relative_abundances_results_M6_v2$FDR, all_outcomes_by_relative_abundances_results_M6_v2$p),]

## add column showing phenotype group
all_outcomes_by_relative_abundances_results_M6_v2$Phenotype_group <- NA
all_outcomes_by_relative_abundances_results_M6_v2[grep("ROME", all_outcomes_by_relative_abundances_results_M6_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_relative_abundances_results_M6_v2[grep("BITSS", all_outcomes_by_relative_abundances_results_M6_v2$y), "Phenotype_group"] <- "Infant_gastrointestinal_health"
all_outcomes_by_relative_abundances_results_M6_v2[grep("growth", all_outcomes_by_relative_abundances_results_M6_v2$y), "Phenotype_group"] <- "Infant_growth"
all_outcomes_by_relative_abundances_results_M6_v2[grep("infant_cry_prev_week_average_crying_time_per_day_min_invr", all_outcomes_by_relative_abundances_results_M6_v2$y), "Phenotype_group"] <- "Infant_crying"

## add column showing reference level
all_outcomes_by_relative_abundances_results_M6_v2[all_outcomes_by_relative_abundances_results_M6_v2$levels=="" & !is.na(all_outcomes_by_relative_abundances_results_M6_v2$levels),"levels"] <- NA
all_outcomes_by_relative_abundances_results_M6_v2$Reference_level <- NA
all_outcomes_by_relative_abundances_results_M6_v2[all_outcomes_by_relative_abundances_results_M6_v2$levels=="yes" &!is.na(all_outcomes_by_relative_abundances_results_M6_v2$levels), "Reference_level"] <- "no"

## resort columns
all_outcomes_by_relative_abundances_results_M6_v3 <- all_outcomes_by_relative_abundances_results_M6_v2[,c(1,10,3:5,11,6:9,2)]
head(all_outcomes_by_relative_abundances_results_M6_v3)

# write.table(all_outcomes_by_relative_abundances_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_14_infant_outcomes_by_milk_relative_abundances_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_11_infant_outcomes_by_milk_relative_abundances_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_8_infant_outcomes_by_milk_relative_abundances_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
# write.table(all_outcomes_by_relative_abundances_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240311_3_infant_outcomes_by_milk_relative_abundances_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(all_outcomes_by_relative_abundances_results_M6_v3, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/9_infant_outcomes/240321_final_infant_outcomes_by_milk_relative_abundances_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 21.7 MILK RELATIVE ABUNDANCES: MORE FOCUSED ANALYSIS: CHECK FOR CONSISTENTLY NOMINALLY SIGNIFICANT ASSOCIATIONS AT SEVERAL TIME POINTS ===== ###

## check for consistently nominally significant associations at several time points
nomsign_RA_results_M1 <- c(paste0(all_outcomes_by_relative_abundances_results_M1_v3[all_outcomes_by_relative_abundances_results_M1_v3$p<0.05, "x"], "_", all_outcomes_by_relative_abundances_results_M1_v3[all_outcomes_by_relative_abundances_results_M1_v3$p<0.05, "y"]))
nomsign_RA_results_M2 <- c(paste0(all_outcomes_by_relative_abundances_results_M2_v3[all_outcomes_by_relative_abundances_results_M2_v3$p<0.05, "x"], "_", all_outcomes_by_relative_abundances_results_M2_v3[all_outcomes_by_relative_abundances_results_M2_v3$p<0.05, "y"]))
nomsign_RA_results_M3 <- c(paste0(all_outcomes_by_relative_abundances_results_M3_v3[all_outcomes_by_relative_abundances_results_M3_v3$p<0.05, "x"], "_", all_outcomes_by_relative_abundances_results_M3_v3[all_outcomes_by_relative_abundances_results_M3_v3$p<0.05, "y"]))
nomsign_RA_results_M6 <- c(paste0(all_outcomes_by_relative_abundances_results_M6_v3[all_outcomes_by_relative_abundances_results_M6_v3$p<0.05, "x"], "_", all_outcomes_by_relative_abundances_results_M6_v3[all_outcomes_by_relative_abundances_results_M6_v3$p<0.05, "y"]))


## consistently nominally significant results at M1, M3
nomsign_RA_results_M1_M3 <- intersect(nomsign_RA_results_M1, nomsign_RA_results_M3)
# [1] "g__Escherichia_Shigella_clr_infant_dis_ROME_e25_mucus_stool_last_week_invr"

## consistently nominally significant results at M1, M2, M3
nomsign_RA_results_M1_M2_M3 <- intersect(nomsign_RA_results_M1_M3, nomsign_RA_results_M2)
# none (but M2 is very small)

## consistently nominally significant results at M1, M3, M6
nomsign_RA_results_M1_M3_M6 <- intersect(nomsign_RA_results_M1_M3, nomsign_RA_results_M6)
# none




