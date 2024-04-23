################################################################################################################
### CREATE BASIC COHORT SUMMARY STATISTICS FOR MICROBIOTA DATA #################################################
################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.0.3 (2020-10-10):     R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE
# 0. DATA IMPORT
# 1. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
# 2. BASIC INFORMATION / STATISTICS
# 3. READ COUNT SUMMARY STATISTICS
# 4. CREATE BASIC COHORT SUMMARY STATISTICS


##### =========================== 0. DATA IMPORT =========================== #####

## import files, phenotypes+microbiota data on genus level
df <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/231129_milk_composition_paper_sel_phenotypes_sel_microbiota_data_genera_n1515.rds")
#1515x1154

# Notes:
# This file contains the following:
# - Phenotypes linked to milk, maternal faecal and infant faecal microbiota data.
#   The file also includes real measured levels of 24 single HMOs and 4 grouped HMOs in µg/ml and HMO-based Le and Se status and milk groups as phenotypes.
#   !! Note that maternal phenotype and HMO data was duplicated for twins. This is indicated by ";1" (single baby/first baby) and ";2" (twin baby/duplicated data) in the mother_ID and mother_sample_ID.
#      -> For association with maternal phenotypes (i.e. non-duplicated data), grep only ";1"
#      -> For association with infant phenotypes, use both ";1" and ";2"


# ##### =========================== 1. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####
# 
# ## for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
# df$mother_milk_collection_season <- factor(df$mother_milk_collection_season, levels = c("Spring", "Summer", "Autumn", "Winter"))
# df$mother_milk_collection_month <- factor(df$mother_milk_collection_month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",  "Oct",  "Nov", "Dec"))
# df$mother_milk_collection_notes <- factor(df$mother_milk_collection_notes, levels = c("No_records_of_incorrect_sample_handling", levels(df$mother_milk_collection_notes)[c(3,5,7,6,2,1)]))
# df$mother_milk_collection_breasts_sampled <- factor(df$mother_milk_collection_breasts_sampled, levels=c("one_breast", "both_breasts"))
# df$mother_genetics_blood_group_genotype <- factor(df$mother_genetics_blood_group_genotype, levels = c("O/O", "A/O", "A/A", "B/O", "B/B", "A/B"))
# df$mother_blood_group <- factor(df$mother_blood_group, levels = c("O", "A", "B", "AB"))
# df$mother_exp_living_situation <- factor(df$mother_exp_living_situation, levels = c("big_house", "small_house", "flat", "farm", "other"))
# df$mother_birth_gestational_age_categories <- factor(df$mother_birth_gestational_age_categories, levels = c("term", "preterm"))
# df$mother_breastpump_brand <- factor(df$mother_breastpump_brand, levels = c("Philips", "Medela", "Ardo", "Other"))
# df$infant_birth_delivery_mode <- factor(df$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
# df$infant_birth_delivery_mode_detailed <- factor(df$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
#                                                                                                       "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
# df$infant_ffq_feeding_type_delivery <- factor(df$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))
# df$mother_milk_HMO_milk_group <- factor(df$mother_milk_HMO_milk_group, levels = c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 2. BASIC INFORMATION / STATISTICS =========================== #####

## check sample numbers
table(df$sample_origin_type, useNA="ifany")
# 737 human milk (incl. 12 duplicated rows with maternal data for infant twins)
# 552 infant faeces
# 172 maternal faeces (incl. 2 duplicated rows with maternal data for infant twins)
# 7 isolation negative controls
# 25 lib prep negative controls
# 1 isolation positive control
# 21 lib prep positive controls

## split data frame by sample type
milk <- df[df$sample_origin_type=="mother_human_milk",]
inff <- df[df$sample_origin_type=="infant_faeces",]
matf <- df[df$sample_origin_type=="mother_faeces",]

## for maternal data frames: select every mother only 1x (i.e. remove duplication of data generated for twin babies)
smilk <- milk[grep(";1", milk$mother_sample_ID),] #725x1154
smatf <- matf[grep(";1", matf$mother_sample_ID),] #170x1154


## number of total measured samples per sample type
length(unique(smilk$mother_sample_ID)) #725
length(unique(inff$mother_sample_ID)) #552
length(unique(smatf$mother_sample_ID)) #170

## number of samples by time point per sample type
table(smilk$time_point)
# 1_month 2_months 3_months 6_months 
#     300       50      268      107

table(inff$time_point)
# 1_month 2_months 3_months 6_months 
#     229       38      195       90

table(smatf$time_point)
# 1_month 2_months 3_months 6_months 
#       0        0      170        0


## number of unique participations per sample type (this includes second pregnancies as a unique participation)
length(unique(smilk$mother_ID)) #336 unique mother-participations
length(unique(inff$infant_ID)) #267 unique infants
length(unique(smatf$mother_ID)) #170 unique mother-participations

## check number of mothers who have data for first and second NEXT pregnancy, per sample type
smilk_mother_ID_simple <- gsub("_1_", "_", smilk$mother_ID)
smilk_mother_ID_simple <- as.factor(as.character(gsub("_2_", "_", smilk_mother_ID_simple)))
length(unique(smilk_mother_ID_simple)) #324 -> 336-324=12 -> 12 mothers have milk microbiota data for first and second NEXT pregnancy

smatf_mother_ID_simple <- gsub("_1_", "_", smatf$mother_ID)
smatf_mother_ID_simple <- as.factor(as.character(gsub("_2_", "_", smatf_mother_ID_simple)))
length(unique(smatf_mother_ID_simple)) #166 -> 170-166=4 -> 4 mothers have milk microbiota data for first and second NEXT pregnancy (all from M3)


## check number of twins
milktwins <- milk[milk$mother_type_pregnancy=="twin_pregnancy",]
milktwins1 <- milktwins[grep("Infant1", milktwins$infant_ID),]
milktwins2 <- milktwins[grep("Infant2", milktwins$infant_ID),]
length(unique(milktwins1$infant_ID)) #5
length(unique(milktwins2$infant_ID)) #5
# -> 5 mothers with twins have milk microbiota data

infftwins <- inff[inff$mother_type_pregnancy=="twin_pregnancy",]
infftwins1 <- infftwins[grep("Infant1", infftwins$infant_ID),]
infftwins2 <- infftwins[grep("Infant2", infftwins$infant_ID),]
length(unique(infftwins1$infant_ID)) #3
length(unique(infftwins2$infant_ID)) #2
# -> for 5 twin babies, we have microbiota data (3x Infant1 and 2x Infant2)


fam <- milk[,c(7,8,11)]
fam <- fam[!duplicated(fam),]
table(fam$mother_type_pregnancy, useNA="ifany") #331 single pregnancies + 10 twin pregnancies (this counts every twin pregnancy 2x)
length(unique(fam$infant_ID)) #341 unique infants
length(fam[grep("Infant1", fam$infant_ID),"infant_ID"]) #336 unique Infant1
length(fam[grep("Infant2", fam$infant_ID),"infant_ID"]) #5 unique Infant2


##### =========================== 3. READ COUNT SUMMARY STATISTICS  =========================== #####

### CLEAN READ COUNT FOR 725 MILK SAMPLES

summary(smilk$seq_16S_n_reads_clean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6761   17334   25335   26573   34828   70876

sd(smilk$seq_16S_n_reads_clean)
# 11839.5


### CLEAN READ COUNT FOR 552 INFANT FAECAL SAMPLES

summary(inff$seq_16S_n_reads_clean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 6836   13870   18113   20394   26313   56027

sd(inff$seq_16S_n_reads_clean)
# 8472.07


### CLEAN READ COUNT FOR 170 MATERNAL FAECAL SAMPLES

summary(smatf$seq_16S_n_reads_clean)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 7595   15694   19408   20779   26612   43020

sd(smatf$seq_16S_n_reads_clean)
# 7664.265


##### =========================== 4. CHECK BASIC COHORT PHENOTYPE SUMMARY STATISTICS  =========================== #####

## save maternal cross-sectional phenotypes for basic summary statistics separately, exclude duplicated data for twins
mc <- smilk[,c("mother_ID", "mother_age_birth","mother_birth_gestational_age_weeks")]
colnames(mc)

mc2 <- mc[!duplicated(mc),] #336 unique mother-participations


for (i in 1:ncol(mc2)){
  print(colnames(mc2)[i])
  print(table(mc2[,i], useNA="ifany"))
}


summary(mc2$mother_age_birth)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   24.00   30.00   32.00   32.56   35.00   43.00       1 

sd(mc2$mother_age_birth, na.rm=T)
#3.722552


summary(mc2$mother_birth_gestational_age_weeks)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#   32.29   38.86   39.71   39.58   40.71   42.14      69

sd(mc2$mother_birth_gestational_age_weeks, na.rm=T)
#1.570084



## save infant cross-sectional phenotypes for basic summary statistics separately
ic <- smilk[,c("infant_ID", "infant_birthweight_g")]
colnames(ic)

ic2 <- ic[!duplicated(ic),] #336 unique infant-participations


for (i in 1:ncol(ic2)){
  print(colnames(ic2)[i])
  print(table(ic2[,i], useNA="ifany"))
}


summary(ic2$infant_birthweight_g)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#    830    3281    3630    3609    3964    5055      74
    
sd(ic2$infant_birthweight_g, na.rm=T)
#545.8957


## save dynamic phenotypes for basic summary statistics separately
id <- smilk[,c("infant_sample_ID", "time_point", "infant_ffq_feeding_mode")]

id_M1 <- id[id$time_point=="1_month",]    #300
id_M2 <- id[id$time_point=="2_months",]   #50
id_M3 <- id[id$time_point=="3_months",]   #268
id_M6 <- id[id$time_point=="6_months",]   #107

table(id_M1$infant_ffq_feeding_mode, useNA="ifany")
# breastfeeding mixed_feeding          <NA> 
#     193 (80%)      49 (20%)           58

table(id_M2$infant_ffq_feeding_mode, useNA="ifany") 
# breastfeeding mixed_feeding          <NA> 
#      20 (53%)      18 (47%)           12
# -> the higher % of mixed feeding at M2 is probably because these mothers stopped breastfeeding at M3 and thus did not provide M3 samples for analysis and
#    as we included M2 samples only if M3 samples were not available, we possibly selectively included M2 samples for mothers who fed more formula at M3

table(id_M3$infant_ffq_feeding_mode, useNA="ifany")
# breastfeeding mixed_feeding          <NA> 
#     161 (79%)      42 (21%)           65

table(id_M6$infant_ffq_feeding_mode, useNA="ifany")
# breastfeeding mixed_feeding          <NA> 
#      54 (67%)      27 (33%)           26


## check number of mothers with info on subclinical mastitis (mineral measurements)
table(smilk$mother_dis_subclinical_mastitis, useNA="ifany")
#  no  yes <NA> 
# 347   16  362
## -> 347+16=363 samples had milk microbiota and mineral measurements

scm <- smilk[!is.na(smilk$mother_dis_subclinical_mastitis),]
table(scm$time_point, useNA="ifany")
# 1_month 2_months 3_months 6_months 
#     204       36      123        0







