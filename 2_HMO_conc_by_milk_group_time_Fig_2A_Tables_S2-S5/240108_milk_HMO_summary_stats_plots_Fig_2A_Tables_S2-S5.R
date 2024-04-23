################################################################################################################
### MILK HMO CONCENTRATIONS BY MILK GROUP AND TIME (FIGURE 2) ##################################################
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
# 2. ENSURE CORRECT DATA STRUCTURE
# 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
# 4. CHECK BASIC INFORMATION / STATISTICS
# 5. CHECK MOTHER DISTRIBUTION AMONG MILK GROUPS
# 6. CHECK ASSOCIATION BETWEEN FUT GENOTYPES AND HMO-LEVEL-BASED Le AND Se STATUS
#     6.1 CHECK ASSOCIATION BETWEEN FUT2/FUT3 GENOTYPES AND HMO-BASED LEWIS AND SECRETOR STATUS (FISHER'S TESTS) (TABLE S2)
#     6.2 CONTINGENCY TABLES FOR FUT3 GENOTYPE AND LEWIS STATUS
#     6.3 CONTINGENCY TABLES FOR FUT2 GENOTYPE AND SECRETOR STATUS
# 7. CHECK HMO DATA DISTRIBUTION
# 8. STACKED BARPLOTS FOR MEAN SINGLE HMO ABSOLUTE ABUNDANCES BY MILK GROUP AND TIME POINT (FIGURE 2A) 
# 9. COMBINE MEAN SINGLE HMO ABSOLUTE ABUNDANCE BARPLOTS IN 1 FIGURE (FIGURE 2A) 
# 10. SUMMARY STATISTICS FOR MEASURED MILK HMO LEVELS (TABLE S3)
# 11. MODELS FOR ASSOCIATION OF TIME POINT AND MATERNAL MILK GROUP WITH HMOs (TABLES S4, S5) 
#      11.1 INVERSE-RANK TRANSFORM HMO LEVELS
#      11.2 CHECK THAT REFERENCE CATEGORIES FOR MILK GROUP (Le+Se+) AND TIME POINT (0.5_months) ARE CORRECTLY ASSIGNED
#      11.3 HISTOGRAMS FOR INVERSE_RANK TRANSFORMED HMO LEVELS
#      11.4 MODEL TO TEST FOR EFFECT OF time_point_numeric ON HMOs WITH CORRECTION FOR milk_group (factor)
#      11.5 MODEL TO TEST FOR EFFECT OF milk_group (factor) ON HMOs WITH CORRECTION FOR time_point_numeric
#      11.6 COMBINE RESULTS FOR EFFECTS OF TIME POINT AND MILK GROUP ON HMOs, CORRECT FOR MULTIPLE TESTING AND SAVE RESULTS                                           


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/2_HMO_conc_by_milk_group_and_time/")

## import file
h <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/phenotypes/240108_milk_composition_paper_sel_phenotypes_hmo_data_n1563.txt", header=T, sep="\t", stringsAsFactors=T) #1563x512

# Notes:
# This file contains the following:
# - Phenotypes linked to measured HMO data.
#   !! Note that maternal phenotype and HMO data was duplicated for twins. This is indicated by ";1" (single baby/first baby) and ";2" (twin baby/duplicated data) in the mother_ID and mother_sample_ID.
#      -> For association with maternal phenotypes (i.e. non-duplicated data), grep only ";1"
#      -> For association with infant phenotypes, use both ";1" and ";2"
# - Real measured levels of 24 single HMOs and 4 grouped HMOs in µg/ml and HMO-based Le and Se status and milk groups
# - A quality control sample was included in all HMO batches and if there was less than <15% variation, the quality of the run was considered good. No need to correct for a batch effect from HMO measurements.


##### =========================== 1. SELECT DATA OF INTEREST =========================== #####

## select every mother only 1x (i.e. remove duplication of data generated for twin babies)
hmo <- h[grep(";1", h$mother_sample_ID),] #1542x512


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

## check that data structure for relevant columns is correct
str(hmo[,1:28])
str(hmo[,481:512])
# -> data structure in real HMO data file is correct.


##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

## for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
hmo$mother_milk_collection_season <- factor(hmo$mother_milk_collection_season, levels = c("Spring", "Summer", "Autumn", "Winter"))
hmo$mother_milk_collection_month <- factor(hmo$mother_milk_collection_month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",  "Oct",  "Nov", "Dec"))
hmo$mother_milk_collection_notes <- factor(hmo$mother_milk_collection_notes, levels = c("No_records_of_incorrect_sample_handling", levels(hmo$mother_milk_collection_notes)[c(3,5,7,6,2,1)]))
hmo$mother_milk_collection_breasts_sampled <- factor(hmo$mother_milk_collection_breasts_sampled, levels=c("one_breast", "both_breasts"))
hmo$mother_genetics_blood_group_genotype <- factor(hmo$mother_genetics_blood_group_genotype, levels = c("O/O", "A/O", "A/A", "B/O", "B/B", "A/B"))
hmo$mother_blood_group <- factor(hmo$mother_blood_group, levels = c("O", "A", "B", "AB"))
hmo$mother_exp_living_situation <- factor(hmo$mother_exp_living_situation, levels = c("big_house", "small_house", "flat", "farm", "other"))
hmo$mother_birth_gestational_age_categories <- factor(hmo$mother_birth_gestational_age_categories, levels = c("term", "preterm"))
hmo$mother_breastpump_brand <- factor(hmo$mother_breastpump_brand, levels = c("Philips", "Medela", "Ardo", "Other"))
hmo$infant_birth_delivery_mode <- factor(hmo$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
hmo$infant_birth_delivery_mode_detailed <- factor(hmo$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                      "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
hmo$infant_ffq_feeding_type_delivery <- factor(hmo$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))
hmo$mother_milk_HMO_milk_group <- factor(hmo$mother_milk_HMO_milk_group, levels = c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 4. CHECK BASIC INFORMATION / STATISTICS =========================== #####

## number of total measured samples 
length(unique(hmo$mother_sample_ID)) #1542

## number of samples by time point
table(hmo$time_point)
# 0.5_months    1_month   2_months   3_months   6_months 
#        214        419        340        433        136

## number of mothers
length(unique(hmo$mother_ID)) #524 mother-participations

hmo$mother_ID_simple <- gsub("_1_", "_", hmo$mother_ID)
hmo$mother_ID_simple <- as.factor(as.character(gsub("_2_", "_", hmo$mother_ID_simple)))
hmo <- hmo[,c(1:6,513,7:512)]
length(unique(hmo$mother_ID_simple)) #500 -> 524-500=24 -> 24 mothers have HMO measurements for first and second NEXT pregnancy


##### =========================== 5. CHECK MOTHER DISTRIBUTION AMONG MILK GROUPS =========================== #####

## check the number of mother-participations and unique mothers per milk group
length(unique(hmo[hmo$mother_milk_HMO_milk_group=="Le-Se-","mother_ID"])) #11 mother-participations
length(unique(hmo[hmo$mother_milk_HMO_milk_group=="Le-Se-","mother_ID_simple"])) #11 unique mothers
length(unique(hmo[hmo$mother_milk_HMO_milk_group=="Le-Se+","mother_ID"])) #32 mother-participations
length(unique(hmo[hmo$mother_milk_HMO_milk_group=="Le-Se+","mother_ID_simple"])) #32 unique mothers
length(unique(hmo[hmo$mother_milk_HMO_milk_group=="Le+Se-","mother_ID"])) #113 mother-participations
length(unique(hmo[hmo$mother_milk_HMO_milk_group=="Le+Se-","mother_ID_simple"])) #105 unique mothers
length(unique(hmo[hmo$mother_milk_HMO_milk_group=="Le+Se+","mother_ID"])) #368 mother-participations
length(unique(hmo[hmo$mother_milk_HMO_milk_group=="Le+Se+","mother_ID_simple"])) #352 unique mothers

## show distribution (n) of mother-participations across milk groups (total n=524)
table(hmo[!duplicated(hmo$mother_ID),"mother_milk_HMO_Le"], hmo[!duplicated(hmo$mother_ID),"mother_milk_HMO_Se"])
#     Se- Se+
# Le-  11  32
# Le+ 113 368

## show distribution (%) of mothers across milk groups (total n=524)
# prop.table(table(hmo[!duplicated(hmo$mother_ID),"mother_milk_HMO_Le"], hmo[!duplicated(hmo$mother_ID),"mother_milk_HMO_Se"]))*100
round(prop.table(table(hmo[!duplicated(hmo$mother_ID),"mother_milk_HMO_Le"], hmo[!duplicated(hmo$mother_ID),"mother_milk_HMO_Se"]))*100,1)
#      Se-  Se+
# Le-  2.1  6.1
# Le+ 21.6 70.2


## show distribution (n) of unique mothers across milk groups (total n=500)
table(hmo[!duplicated(hmo$mother_ID_simple),"mother_milk_HMO_Le"], hmo[!duplicated(hmo$mother_ID_simple),"mother_milk_HMO_Se"])
#     Se- Se+
# Le-  11  32
# Le+ 105 352

## show distribution (%) of mothers across milk groups (total n=500)
prop.table(table(hmo[!duplicated(hmo$mother_ID_simple),"mother_milk_HMO_Le"], hmo[!duplicated(hmo$mother_ID_simple),"mother_milk_HMO_Se"]))*100
#      Se-  Se+
# Le-  2.2  6.4
# Le+ 21.0 70.4


##### =========================== 6. CHECK ASSOCIATION BETWEEN FUT GENOTYPES AND HMO-LEVEL-BASED Le AND Se STATUS =========================== #####

## check for how many unique mothers genetics data is available
nrow(hmo[!duplicated(hmo$mother_ID_simple) & !is.na(hmo$mother_genetics_FUT3_C314T_rs778986),])
# -> 229/500 (46%) unique mothers have genetic data

## save list of mother IDs for which genotypes are available
# write.table(hmo[!duplicated(hmo$mother_ID_simple) & !is.na(hmo$mother_genetics_FUT3_C314T_rs778986),1:15], file="230821_sample_IDs_with_available_genotypes_n229.txt", row.names=F, col.names=T, sep="\t", quote=F)


## check for how many mother-participations genetics data is available
nrow(hmo[!duplicated(hmo$mother_ID) & !is.na(hmo$mother_genetics_FUT3_C314T_rs778986),])
# -> 249/524 (48%) unique mother-participations have genetic data, 249-229=20 mothers who have genetic info available also participated with a second pregnancy


### ===== 6.1 CHECK ASSOCIATION BETWEEN FUT2/FUT3 GENOTYPES AND HMO-BASED LEWIS AND SECRETOR STATUS (FISHER'S TESTS) (TABLE S2) ===== ###

## function to associate FUT2/FUT3 genotypes with Lewis/Secretor status
associate.FUT.genotypes.with.milk.groups <- function(inputdata, x_columns, y_columns){
  x <- c()
  y <- c()
  n <- c()
  p_value <- c()
  
  for (i in x_columns){
    for (j in y_columns){
      print(paste0("Associating ", colnames(inputdata)[i], " with ", colnames(inputdata)[j]))
      
      x <- c(x, paste0(colnames(inputdata)[i]))
      y <- c(y, paste0(colnames(inputdata)[j]))
      n <- c(n, nrow(inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]))
      
      #Fisher's test
      fres <- fisher.test(inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),i], inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),j])
      p_value <- c(p_value, fres$p.value)
    }
  }
  
  results <- data.frame(x=x, y=y, n=n, p_value=p_value)
  
  return(results)
}

## check data distribution in genotype columns
for (i in 23:29){
  print(colnames(hmo)[i])
  print(table(hmo[!duplicated(hmo$mother_ID_simple),i], useNA="ifany"))
}

## FUT3 and Lewis status
# There are 3 FUT3 SNPs and all have sufficient variation to run Fisher's tests.
FUT3_results <- associate.FUT.genotypes.with.milk.groups(hmo[!duplicated(hmo$mother_ID_simple),c(7,23:25,484)], x_columns=c(2:4), y_columns=c(5))

## FUT2 and secretor status
# There is too little variation in the FUT2 SNPs rs200157007 (all unique mothers have C/C genotype) and rs1047781 (228 unique mothers have A/A genotype and only 1 mother has T/A genotype (weak/non-secretor)) to run Fisher's tests for these SNPs.
# Only associate the main FUT2 SNP (rs601338) with the HMO-based secretor status.
FUT2_results <- associate.FUT.genotypes.with.milk.groups(hmo[!duplicated(hmo$mother_ID_simple),c(7,26,485)], x_columns=c(2), y_columns=c(3))

## combine results in one table and adjust for multiple testing
FUT_results <- rbind(FUT3_results, FUT2_results)
FUT_results$BH_adj_p_value <- p.adjust(FUT_results$p_value, method="BH")
FUT_results <- FUT_results[order(FUT_results$BH_adj_p_value),]
FUT_results

## save results
# write.table(FUT_results, file="240108_association_FUT2_FUT2_Le_Se_status_n229.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 6.2 CONTINGENCY TABLES FOR FUT3 GENOTYPE AND LEWIS STATUS ===== ###

gen_n229 <- hmo[!duplicated(hmo$mother_ID_simple) & !is.na(hmo$mother_genetics_FUT3_C314T_rs778986),c(7,23:28,484:485)]

table(gen_n229$mother_genetics_FUT3_C314T_rs778986, gen_n229$mother_milk_HMO_Le)
#     Le- Le+
# A/A  11   0
# A/G  10  55
# G/G   3 150

table(gen_n229$mother_genetics_FUT3_T202C_rs812936, gen_n229$mother_milk_HMO_Le)
#     Le- Le+
# A/A   3 146
# G/A  10  59
# G/G  11   0

table(gen_n229$mother_genetics_FUT3_T59G_rs28362459, gen_n229$mother_milk_HMO_Le)
#     Le- Le+
# A/A  12 173
# C/A   9  32
# C/C   3   0


### ===== 6.3 CONTINGENCY TABLES FOR FUT2 GENOTYPE AND SECRETOR STATUS ===== ###

table(gen_n229$mother_genetics_FUT2_G_to_A_rs601338, gen_n229$mother_milk_HMO_Se)
#     Se- Se+
# A/A  51   0
# A/G   1 111
# G/G   0  66

table(gen_n229$mother_genetics_FUT2_C_to_T_rs200157007, gen_n229$mother_milk_HMO_Se)
#     Se- Se+
# C/C  52 177

table(gen_n229$mother_genetics_FUT2_A_to_T_rs1047781, gen_n229$mother_milk_HMO_Se)
#     Se- Se+
# A/A  51 177
# T/A   1   0


### check the individual with mismatches in main FUT2 SNP - Se status (FAM1182_1_Mother)
hmo[!is.na(hmo$mother_genetics_FUT2_G_to_A_rs601338) & hmo$mother_genetics_FUT2_G_to_A_rs601338=="A/G" & hmo$mother_milk_HMO_Se=="Se-",c(1,15,23:28,482:485,486,499,504,501)]
# this is the mother with the T/A genotype for SNP rs1047781 so even though the main FUT2 SNP makes her a secretor, she still is a weak/non-secretor


##### =========================== 7. CHECK HMO DATA DISTRIBUTION  =========================== #####

## work with 524 unique mother-participations, i.e. see the 24 mothers who participated with 2 NEXT pregnancies as individual samples
## check distribution of single and grouped HMO levels per milk group, pooled and separately for all time points

library(ggplot2)

## make histograms for each milk group
check.distribution <- function(inputdata, milk_group, time_point, hmo_columns){
  print(paste0("Creating histograms for ", milk_group, " at ", time_point, " time point(s)"))
  
  pdf(paste0("HMO_data_distribution/240108_HMO_levels_histograms_", milk_group, "_", time_point, ".pdf"))
  
  for (i in c(hmo_columns)){
    hmohist <- ggplot(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group,],
                      aes(x=as.numeric(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group, i]))) + geom_histogram() + labs(x=colnames(inputdata)[i])
    print(hmohist)
  }
  
  dev.off()
}

# check HMO data pooled from all time points
check.distribution(hmo, milk_group="Le-Se-", time_point="all", hmo_columns=c(486:513))
check.distribution(hmo, milk_group="Le-Se+", time_point="all", hmo_columns=c(486:513))
check.distribution(hmo, milk_group="Le+Se-", time_point="all", hmo_columns=c(486:513))
check.distribution(hmo, milk_group="Le+Se+", time_point="all", hmo_columns=c(486:513))

# W2
check.distribution(hmo[hmo$time_point=="0.5_months",], milk_group="Le-Se-", time_point="W2", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="0.5_months",], milk_group="Le-Se+", time_point="W2", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="0.5_months",], milk_group="Le+Se-", time_point="W2", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="0.5_months",], milk_group="Le+Se+", time_point="W2", hmo_columns=c(486:513))

# M1
check.distribution(hmo[hmo$time_point=="1_month",], milk_group="Le-Se-", time_point="M1", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="1_month",], milk_group="Le-Se+", time_point="M1", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="1_month",], milk_group="Le+Se-", time_point="M1", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="1_month",], milk_group="Le+Se+", time_point="M1", hmo_columns=c(486:513))

# M2
check.distribution(hmo[hmo$time_point=="2_months",], milk_group="Le-Se-", time_point="M2", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="2_months",], milk_group="Le-Se+", time_point="M2", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="2_months",], milk_group="Le+Se-", time_point="M2", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="2_months",], milk_group="Le+Se+", time_point="M2", hmo_columns=c(486:513))

# M3
check.distribution(hmo[hmo$time_point=="3_months",], milk_group="Le-Se-", time_point="M3", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="3_months",], milk_group="Le-Se+", time_point="M3", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="3_months",], milk_group="Le+Se-", time_point="M3", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="3_months",], milk_group="Le+Se+", time_point="M3", hmo_columns=c(486:513))

# M6
check.distribution(hmo[hmo$time_point=="6_months",], milk_group="Le-Se-", time_point="M6", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="6_months",], milk_group="Le-Se+", time_point="M6", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="6_months",], milk_group="Le+Se-", time_point="M6", hmo_columns=c(486:513))
check.distribution(hmo[hmo$time_point=="6_months",], milk_group="Le+Se+", time_point="M6", hmo_columns=c(486:513))

## -> most look fairly normal, some a bit left-skewed and for certain HMOs there is a zero-inflation
## -> work with means for now


##### =========================== 8. STACKED BARPLOTS FOR MEAN SINGLE HMO ABSOLUTE ABUNDANCES BY MILK GROUP AND TIME POINT (FIGURE 2A) =========================== #####

## work with 524 unique mother-participations, i.e. see the 24 mothers who participated with 2 NEXT pregnancies as individual samples

## select only relevant columns (first info columns and single HMO levels, exclude grouped HMO columns) and separate data by milk group
single_hmo_Le0Se0 <- hmo[hmo$mother_milk_HMO_milk_group=="Le-Se-",c(1,9,483,486:509)] #  30 samples
single_hmo_Le0Se1 <- hmo[hmo$mother_milk_HMO_milk_group=="Le-Se+",c(1,9,483,486:509)] #  91 samples
single_hmo_Le1Se0 <- hmo[hmo$mother_milk_HMO_milk_group=="Le+Se-",c(1,9,483,486:509)] # 345 samples
single_hmo_Le1Se1 <- hmo[hmo$mother_milk_HMO_milk_group=="Le+Se+",c(1,9,483,486:509)] #1076 samples

## function to calculate means for each time point and milk group and change data frame to long format for plotting
library(reshape)

calculate.single.hmo.means <- function(inputdata, hmo_columns, milk_group){
  print(paste0("Calculating means and SD for ", milk_group))
  
  mean_W2 <- as.data.frame(apply(inputdata[inputdata$time_point=="0.5_months",hmo_columns], 2, mean))
  mean_M1 <- as.data.frame(apply(inputdata[inputdata$time_point=="1_month",hmo_columns], 2, mean))
  mean_M2 <- as.data.frame(apply(inputdata[inputdata$time_point=="2_months",hmo_columns], 2, mean))
  mean_M3 <- as.data.frame(apply(inputdata[inputdata$time_point=="3_months",hmo_columns], 2, mean))
  mean_M6 <- as.data.frame(apply(inputdata[inputdata$time_point=="6_months",hmo_columns], 2, mean))
  
  print(paste0("Combining means"))
  my_means <- rbind(as.data.frame(t(mean_W2)),
                    as.data.frame(t(mean_M1)),
                    as.data.frame(t(mean_M2)),
                    as.data.frame(t(mean_M3)),
                    as.data.frame(t(mean_M6)))
  
  print(paste0("Adding columns for milk group and time point"))
  my_means$milk_group <- c(rep(milk_group,5))
  my_means$time_point <- c("0.5_months","1_month","2_months","3_months","6_months")
  rownames(my_means) <- 1:5
  
  print(paste0("Changing data frame structure"))
  my_means_long <- melt(my_means, id.vars=colnames(my_means)[(length(hmo_columns)+1):ncol(my_means)], variable_name="HMO")
  my_means_long$time_point <- factor(my_means_long$time_point, levels=c("0.5_months","1_month","2_months","3_months","6_months"))
  my_means_long$HMO <- gsub("mother_milk_HMO_", "", my_means_long$HMO)
  my_means_long$HMO <- gsub("_ugml", "", my_means_long$HMO)
  
  return(my_means_long)
}

Le0Se0_means_long <- calculate.single.hmo.means(single_hmo_Le0Se0, hmo_columns = c(4:27), milk_group = "Le-Se-")
Le0Se1_means_long <- calculate.single.hmo.means(single_hmo_Le0Se1, hmo_columns = c(4:27), milk_group = "Le-Se+")
Le1Se0_means_long <- calculate.single.hmo.means(single_hmo_Le1Se0, hmo_columns = c(4:27), milk_group = "Le+Se-")
Le1Se1_means_long <- calculate.single.hmo.means(single_hmo_Le1Se1, hmo_columns = c(4:27), milk_group = "Le+Se+")

# define x_axis_labels
x_axis_labels_single_Le0Se0 <- c(paste0("0.5_months_n=", nrow(single_hmo_Le0Se0[single_hmo_Le0Se0$time_point=="0.5_months",])),
                                 paste0("1_month_n=", nrow(single_hmo_Le0Se0[single_hmo_Le0Se0$time_point=="1_month",])),
                                 paste0("2_months_n=", nrow(single_hmo_Le0Se0[single_hmo_Le0Se0$time_point=="2_months",])),
                                 paste0("3_months_n=", nrow(single_hmo_Le0Se0[single_hmo_Le0Se0$time_point=="3_months",])),
                                 paste0("6_months_n=", nrow(single_hmo_Le0Se0[single_hmo_Le0Se0$time_point=="6_months",])))
x_axis_labels_single_Le0Se1 <- c(paste0("0.5_months_n=", nrow(single_hmo_Le0Se1[single_hmo_Le0Se1$time_point=="0.5_months",])),
                                 paste0("1_month_n=", nrow(single_hmo_Le0Se1[single_hmo_Le0Se1$time_point=="1_month",])),
                                 paste0("2_months_n=", nrow(single_hmo_Le0Se1[single_hmo_Le0Se1$time_point=="2_months",])),
                                 paste0("3_months_n=", nrow(single_hmo_Le0Se1[single_hmo_Le0Se1$time_point=="3_months",])),
                                 paste0("6_months_n=", nrow(single_hmo_Le0Se1[single_hmo_Le0Se1$time_point=="6_months",])))
x_axis_labels_single_Le1Se0 <- c(paste0("0.5_months_n=", nrow(single_hmo_Le1Se0[single_hmo_Le1Se0$time_point=="0.5_months",])),
                                 paste0("1_month_n=", nrow(single_hmo_Le1Se0[single_hmo_Le1Se0$time_point=="1_month",])),
                                 paste0("2_months_n=", nrow(single_hmo_Le1Se0[single_hmo_Le1Se0$time_point=="2_months",])),
                                 paste0("3_months_n=", nrow(single_hmo_Le1Se0[single_hmo_Le1Se0$time_point=="3_months",])),
                                 paste0("6_months_n=", nrow(single_hmo_Le1Se0[single_hmo_Le1Se0$time_point=="6_months",])))
x_axis_labels_single_Le1Se1 <- c(paste0("0.5_months_n=", nrow(single_hmo_Le1Se1[single_hmo_Le1Se1$time_point=="0.5_months",])),
                                 paste0("1_month_n=", nrow(single_hmo_Le1Se1[single_hmo_Le1Se1$time_point=="1_month",])),
                                 paste0("2_months_n=", nrow(single_hmo_Le1Se1[single_hmo_Le1Se1$time_point=="2_months",])),
                                 paste0("3_months_n=", nrow(single_hmo_Le1Se1[single_hmo_Le1Se1$time_point=="3_months",])),
                                 paste0("6_months_n=", nrow(single_hmo_Le1Se1[single_hmo_Le1Se1$time_point=="6_months",])))

## change the order of HMOs to show 1. neutral, 2. fucosylated and 3. sialylated HMOs (within groups sort them alphabetically)
# convert columns from character to factor
for (i in c(1,3)){Le0Se0_means_long[,i] <- as.factor(as.character(Le0Se0_means_long[,i]))}
for (i in c(1,3)){Le0Se1_means_long[,i] <- as.factor(as.character(Le0Se1_means_long[,i]))}
for (i in c(1,3)){Le1Se0_means_long[,i] <- as.factor(as.character(Le1Se0_means_long[,i]))}
for (i in c(1,3)){Le1Se1_means_long[,i] <- as.factor(as.character(Le1Se1_means_long[,i]))}

# change the order of the HMO levels
Le0Se0_means_long$HMO<-factor(Le0Se0_means_long$HMO,
                              levels=c("3GL","6GL","LNH","LNnT","LNT", #neutral HMOs
                                       "2FL","3FL","A_tetra","DFLNHa","LDFT", #neutral+fucosylated HMOs
                                       "LNDH_I","LNFP_I","LNFP_II","LNFP_III","LNFP_V",
                                       "LNnDFH","LNnFP_V","MFLNH_III",
                                       "3F3SL",
                                       "3SL","6SL","DSLNT","LSTb","LSTc")) #sialylated HMOs #note that 3F3'SL is both fucosylated and sialylated
Le0Se1_means_long$HMO<-factor(Le0Se1_means_long$HMO,
                              levels=c("3GL","6GL","LNH","LNnT","LNT", #neutral HMOs
                                       "2FL","3FL","A_tetra","DFLNHa","LDFT", #neutral+fucosylated HMOs
                                       "LNDH_I","LNFP_I","LNFP_II","LNFP_III","LNFP_V",
                                       "LNnDFH","LNnFP_V","MFLNH_III",
                                       "3F3SL",
                                       "3SL","6SL","DSLNT","LSTb","LSTc")) #sialylated HMOs #note that 3F3'SL is both fucosylated and sialylated
Le1Se0_means_long$HMO<-factor(Le1Se0_means_long$HMO,
                              levels=c("3GL","6GL","LNH","LNnT","LNT", #neutral HMOs
                                       "2FL","3FL","A_tetra","DFLNHa","LDFT", #neutral+fucosylated HMOs
                                       "LNDH_I","LNFP_I","LNFP_II","LNFP_III","LNFP_V",
                                       "LNnDFH","LNnFP_V","MFLNH_III",
                                       "3F3SL",
                                       "3SL","6SL","DSLNT","LSTb","LSTc")) #sialylated HMOs #note that 3F3'SL is both fucosylated and sialylated
Le1Se1_means_long$HMO<-factor(Le1Se1_means_long$HMO,
                              levels=c("3GL","6GL","LNH","LNnT","LNT", #neutral HMOs
                                       "2FL","3FL","A_tetra","DFLNHa","LDFT", #neutral+fucosylated HMOs
                                       "LNDH_I","LNFP_I","LNFP_II","LNFP_III","LNFP_V",
                                       "LNnDFH","LNnFP_V","MFLNH_III",
                                       "3F3SL",
                                       "3SL","6SL","DSLNT","LSTb","LSTc")) #sialylated HMOs #note that 3F3'SL is both fucosylated and sialylated

## make stacked barplots for HMO mean absolute abundances by mlk group and time point
make.hmo.plots <- function(inputdata, x_axis_labels, milk_group){
  print(paste0("Creating plot for mean single HMOs for ", milk_group))
  
  p1 <- ggplot(inputdata, aes(x=time_point, y=value, fill=HMO))+
    geom_bar(position="stack", stat="identity")+
    theme_bw()+
    theme(panel.grid = element_blank(),
          text = element_text(size=20),
          axis.text.x = element_text(angle = 90))+
    labs(x="", y="Mean HMO abundance [µg/ml]", title=c(paste0("Milk group:\n", milk_group)))+
    scale_x_discrete(labels=x_axis_labels)+
    scale_y_continuous(limits=c(0,12000), breaks=seq(from=0,to=12000,by=2000))+
    scale_fill_manual(values=c("palegreen1","#CC79A7","lightgoldenrod1","orangered","#E69F00",
                               "blue","turquoise","gold","navy","thistle2",
                               "seagreen3","darkgreen","deeppink3","burlywood","chocolate4",
                               "olivedrab2","darkcyan","darkorange1","mistyrose","#0072B2",
                               "paleturquoise2","firebrick3","#56B4E9","plum1"))
  
  ggsave(paste0("plots_HMO_levels_by_milk_group_and_time/240108_stacked_barplot_single_HMOs_", milk_group, "_by_time.pdf"), p1, dpi=500, width=30, height=20, units="cm", useDingbats=F)
  return(p1)
}

# plots
Le0Se0_single_means_plot <- make.hmo.plots(Le0Se0_means_long, x_axis_labels=x_axis_labels_single_Le0Se0, milk_group="Le-Se-")
Le0Se1_single_means_plot <- make.hmo.plots(Le0Se1_means_long, x_axis_labels=x_axis_labels_single_Le0Se1, milk_group="Le-Se+")
Le1Se0_single_means_plot <- make.hmo.plots(Le1Se0_means_long, x_axis_labels=x_axis_labels_single_Le1Se0, milk_group="Le+Se-")
Le1Se1_single_means_plot <- make.hmo.plots(Le1Se1_means_long, x_axis_labels=x_axis_labels_single_Le1Se1, milk_group="Le+Se+")


##### =========================== 9. COMBINE MEAN SINGLE HMO ABSOLUTE ABUNDANCE BARPLOTS IN 1 FIGURE (FIGURE 2A) =========================== #####

library(ggpubr)
library(cowplot)

## save legend separately
single_legend <- as_ggplot(get_legend(Le0Se0_single_means_plot))

## save plots without legend
Le1Se1_single_means_plot_pure <- Le1Se1_single_means_plot + theme(legend.position = "none", axis.text.x=element_text(angle=90))
Le1Se0_single_means_plot_pure <- Le1Se0_single_means_plot + theme(legend.position = "none", axis.text.x=element_text(angle=90)) + labs(y="")
Le0Se1_single_means_plot_pure <- Le0Se1_single_means_plot + theme(legend.position = "none", axis.text.x=element_text(angle=90))
Le0Se0_single_means_plot_pure <- Le0Se0_single_means_plot + theme(legend.position = "none", axis.text.x=element_text(angle=90)) + labs(y="")

## combine plots
HMO_plots <- ggdraw(plot_grid(Le1Se1_single_means_plot_pure, Le1Se0_single_means_plot_pure, Le0Se1_single_means_plot_pure,  Le0Se0_single_means_plot_pure, single_legend,
                              nrow=1, align="h"))

## save plot
ggsave("plots_HMO_levels_by_milk_group_and_time/240108_barplots_HMOs_single_all_milk_groups.pdf", HMO_plots, dpi=500, width=40, height=20, units="cm", useDingbats=F)


##### =========================== 10. SUMMARY STATISTICS FOR MEASURED MILK HMO LEVELS (TABLE S3) =========================== #####

## work with 524 unique mother-participations, i.e. see the 24 mothers who participated with 2 NEXT pregnancies as individual samples

## function to retrieve summary statistics
get.HMO.sum.stats <- function(inputdata, milk_group, time_point){
  print(paste0("Generating summary statistics for ", milk_group, " at ", time_point))
  
  #summary statistics
  SumStats <- as.data.frame(apply(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group, c(486:513)], 2, summary))
  
  #add SD
  SumSD <- c()
  for (i in c(486:513)){SumSD <- c(SumSD, sd(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group, i]))}
  SumStats[7,] <- SumSD
  
  #add n
  SumStats[8,] <- rep(nrow(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group,]), 4)
  
  #clean summary statistics file
  SumStats$statistic <- c(rownames(SumStats)[1:6], "SD", "n")
  rownames(SumStats) <- 1:8
  
  SumStats <- SumStats[c(4,7,3,1,6,2,5,8),c(29,25,27,26,28, #resort to show HMOs in the same order as in other tables and figures
                                            4,5,21,9,8,
                                            1,2,6,24,3,
                                            19,14,16,15,12,
                                            17,13,22,
                                              11,
                                              7,10,23,18,20)]
  
  #transpose and add time point and milk group to create supp. table
  rownames(SumStats) <- SumStats$statistic
  SumStats <- SumStats[,-1]
  SumStats2 <- as.data.frame(t(SumStats))
  SumStats2$time_point <- c(rep(time_point,nrow(SumStats2)))
  SumStats2$milk_group <- c(rep(milk_group,nrow(SumStats2)))
  SumStats2$HMO <- rownames(SumStats2)
  SumStats2$HMO <- gsub("mother_milk_HMO_", "", SumStats2$HMO)
  SumStats2$HMO <- gsub("_ugml", "", SumStats2$HMO)
  rownames(SumStats2) <- 1:nrow(SumStats2)
  SumStats2 <- SumStats2[,c(9:11,1:8)]

  #save summary statistics
  # write.table(SumStats2, file=paste0("HMO_data_summary_statistics/240108_HMO_summary_statistics_", milk_group, "_", time_point, ".txt"), sep="\t", row.names=F, quote=F)

  return(SumStats2)
}

# for all time points together
HMO_sum_Le0Se0_all <- get.HMO.sum.stats(hmo, milk_group="Le-Se-", time_point="all")
HMO_sum_Le0Se1_all <- get.HMO.sum.stats(hmo, milk_group="Le-Se+", time_point="all")
HMO_sum_Le1Se0_all <- get.HMO.sum.stats(hmo, milk_group="Le+Se-", time_point="all")
HMO_sum_Le1Se1_all <- get.HMO.sum.stats(hmo, milk_group="Le+Se+", time_point="all")

# W2
HMO_sum_Le0Se0_W2 <- get.HMO.sum.stats(hmo[hmo$time_point=="0.5_months",], milk_group="Le-Se-", time_point="W2")
HMO_sum_Le0Se1_W2 <- get.HMO.sum.stats(hmo[hmo$time_point=="0.5_months",], milk_group="Le-Se+", time_point="W2")
HMO_sum_Le1Se0_W2 <- get.HMO.sum.stats(hmo[hmo$time_point=="0.5_months",], milk_group="Le+Se-", time_point="W2")
HMO_sum_Le1Se1_W2 <- get.HMO.sum.stats(hmo[hmo$time_point=="0.5_months",], milk_group="Le+Se+", time_point="W2")

# M1
HMO_sum_Le0Se0_M1 <- get.HMO.sum.stats(hmo[hmo$time_point=="1_month",], milk_group="Le-Se-", time_point="M1")
HMO_sum_Le0Se1_M1 <- get.HMO.sum.stats(hmo[hmo$time_point=="1_month",], milk_group="Le-Se+", time_point="M1")
HMO_sum_Le1Se0_M1 <- get.HMO.sum.stats(hmo[hmo$time_point=="1_month",], milk_group="Le+Se-", time_point="M1")
HMO_sum_Le1Se1_M1 <- get.HMO.sum.stats(hmo[hmo$time_point=="1_month",], milk_group="Le+Se+", time_point="M1")

# M2
HMO_sum_Le0Se0_M2 <- get.HMO.sum.stats(hmo[hmo$time_point=="2_months",], milk_group="Le-Se-", time_point="M2")
HMO_sum_Le0Se1_M2 <- get.HMO.sum.stats(hmo[hmo$time_point=="2_months",], milk_group="Le-Se+", time_point="M2")
HMO_sum_Le1Se0_M2 <- get.HMO.sum.stats(hmo[hmo$time_point=="2_months",], milk_group="Le+Se-", time_point="M2")
HMO_sum_Le1Se1_M2 <- get.HMO.sum.stats(hmo[hmo$time_point=="2_months",], milk_group="Le+Se+", time_point="M2")

# M3
HMO_sum_Le0Se0_M3 <- get.HMO.sum.stats(hmo[hmo$time_point=="3_months",], milk_group="Le-Se-", time_point="M3")
HMO_sum_Le0Se1_M3 <- get.HMO.sum.stats(hmo[hmo$time_point=="3_months",], milk_group="Le-Se+", time_point="M3")
HMO_sum_Le1Se0_M3 <- get.HMO.sum.stats(hmo[hmo$time_point=="3_months",], milk_group="Le+Se-", time_point="M3")
HMO_sum_Le1Se1_M3 <- get.HMO.sum.stats(hmo[hmo$time_point=="3_months",], milk_group="Le+Se+", time_point="M3")

# M6
HMO_sum_Le0Se0_M6 <- get.HMO.sum.stats(hmo[hmo$time_point=="6_months",], milk_group="Le-Se-", time_point="M6")
HMO_sum_Le0Se1_M6 <- get.HMO.sum.stats(hmo[hmo$time_point=="6_months",], milk_group="Le-Se+", time_point="M6")
HMO_sum_Le1Se0_M6 <- get.HMO.sum.stats(hmo[hmo$time_point=="6_months",], milk_group="Le+Se-", time_point="M6")
HMO_sum_Le1Se1_M6 <- get.HMO.sum.stats(hmo[hmo$time_point=="6_months",], milk_group="Le+Se+", time_point="M6")


## combine tables with summary statistics from different time points and milk groups into 1 supp. table
HMO_sums <- rbind(HMO_sum_Le1Se1_W2, HMO_sum_Le1Se0_W2, HMO_sum_Le0Se1_W2, HMO_sum_Le0Se0_W2,
                  HMO_sum_Le1Se1_M1, HMO_sum_Le1Se0_M1, HMO_sum_Le0Se1_M1, HMO_sum_Le0Se0_M1,
                  HMO_sum_Le1Se1_M2, HMO_sum_Le1Se0_M2, HMO_sum_Le0Se1_M2, HMO_sum_Le0Se0_M2,
                  HMO_sum_Le1Se1_M3, HMO_sum_Le1Se0_M3, HMO_sum_Le0Se1_M3, HMO_sum_Le0Se0_M3,
                  HMO_sum_Le1Se1_M6, HMO_sum_Le1Se0_M6, HMO_sum_Le0Se1_M6, HMO_sum_Le0Se0_M6)
for (i in 1:3){HMO_sums[,i] <- as.factor(as.character(HMO_sums[,i]))}
HMO_sums <- HMO_sums[,c(1:6,9:10,7:8,11)]

# change time points annotation
levels(HMO_sums$time_point) <- c("1_month", "2_months", "3_months", "6_months", "0.5_months")

# save supp. table
# write.table(HMO_sums, file="HMO_data_summary_statistics/240108_HMO_summary_statistics_all_milk_groups_W2-M6.txt", sep="\t", row.names=F, quote=F)


##### =========================== 11. MODELS FOR ASSOCIATION OF TIME POINT AND MATERNAL MILK GROUP WITH HMOs (TABLES S4, S5) =========================== #####

## Use data frame hmo (this correctly contains every mothers data 1x, the duplicated data for twins has already been removed above)
# location for saving output: /groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/2_HMO_conc_by_milk_group_and_time/model_results_HMO_levels_by_milk_group_and_time/

### ===== 11.1 INVERSE-RANK TRANSFORM HMO LEVELS

### function for inverse rank tranformation
invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}

hmoinvr <- hmo
hmoinvr[,486:513] <- as.data.frame(apply(hmoinvr[,486:513], 2, invrank))
for (i in 486:513){colnames(hmoinvr)[i] <- paste0(colnames(hmoinvr)[i], "_invr")}


### ===== 11.2 CHECK THAT REFERENCE CATEGORIES FOR MILK GROUP (Le+Se+) AND TIME POINT (0.5_months) ARE CORRECTLY ASSIGNED

### change reference category before running models
levels(hmoinvr$time_point) #time points are correctly ordered
levels(hmoinvr$mother_milk_HMO_milk_group) #milk groups are correctly ordered


### ===== 11.3 HISTOGRAMS FOR INVERSE_RANK TRANSFORMED HMO LEVELS

## make histograms for each milk group
check.distribution.invr <- function(inputdata, milk_group, time_point, hmo_columns){
  print(paste0("Creating histograms for ", milk_group, " at ", time_point, " time point(s)"))
  
  pdf(paste0("HMO_data_distribution/240108_HMO_invr_levels_histograms_", milk_group, "_", time_point, ".pdf"))
  
  for (i in c(hmo_columns)){
    hmohist <- ggplot(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group,],
                      aes(x=as.numeric(inputdata[inputdata$mother_milk_HMO_milk_group==milk_group, i]))) + geom_histogram() + labs(x=colnames(inputdata)[i])
    print(hmohist)
  }
  
  dev.off()
}

# check HMO data pooled from all time points
check.distribution.invr(hmoinvr, milk_group="Le+Se+", time_point="all", hmo_columns=c(486:513))
check.distribution.invr(hmoinvr, milk_group="Le+Se-", time_point="all", hmo_columns=c(486:513))
check.distribution.invr(hmoinvr, milk_group="Le-Se+", time_point="all", hmo_columns=c(486:513))
check.distribution.invr(hmoinvr, milk_group="Le-Se-", time_point="all", hmo_columns=c(486:513))


### ===== 11.4 MODEL TO TEST FOR EFFECT OF time_point_numeric ON HMOs WITH CORRECTION FOR milk_group (factor)

## model function
run.lmer.cor.milkgroup <- function(datasetname, inputdata, xcolumn, ycolumn){
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
      print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #first mixed model with correction for mother_milk_HMO_milk_group as a fixed effect and without the phenotype of interest (i)
      m0 <- lmerTest::lmer(my_df[,j] ~ mother_milk_HMO_milk_group + (1|mother_ID), REML=F, data=my_df)
      sum0 <- summary(m0)
      
      #second mixed model with correction for mother_milk_HMO_milk_group as a fixed effect and WITH the phenotype of interest (i)
      m1 <- lmerTest::lmer(my_df[,j] ~ mother_milk_HMO_milk_group + (1|mother_ID) + my_df[,i], REML=F, data=my_df)
      sum1 <- summary(m1)
      
      #compare models and save p value from model comparison
      an1 <- anova(m0, m1)
      
      p_models <- c(p_models, an1["m1", "Pr(>Chisq)"])
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      levels <- c(levels, paste(collapse=";",levels(my_df[,i])))
      n_levels <- c(n_levels, paste(collapse=";",table(my_df[,i])))
      # e <- c(e, paste(collapse=";",c(0,signif(digits=3,sum1$coef[5:nrow(sum1$coef), "Estimate"])))) #first 4 rows are for 4 milk groups, ensure that this picks from row 5 on as these are the results for i
      e <- c(e, paste(signif(digits=3,sum1$coef[5:nrow(sum1$coef), "Estimate"]))) #first 4 rows are for 4 milk groups, ensure that this picks from row 5 on as these are the results for i
      
      stat <- c(stat, paste0("lmer_correction_for_milkgroup_ID"))
      
      d <- c(d, datasetname)
      x <- c(x, colnames(inputdata)[i])
      y <- c(y, colnames(inputdata)[j])
      n <- c(n, nrow(my_df))
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, n_levels=n_levels, estimate=e, p_models=p_models)
  
  #add FDR correction
  res$FDR <- p.adjust(res$p_models, method="BH")
  
  #sort results by FDR and p value
  res2 <- res[order(res$FDR, res$p_models),]
  
  return(res2)
}

## run model
lmer.results.by.time <- run.lmer.cor.milkgroup(datasetname="mother",
                                               inputdata=hmoinvr,
                                               xcolumn=c(10), #test the effect of time point (numeric!) ..
                                               ycolumn=c(486:513)) #..on milk HMOs

# ## check that models return correct p values: HMO with milk group and time point
# m0 <- lmerTest::lmer(mother_milk_HMO_3SL_ugml_invr ~ mother_milk_HMO_milk_group + (1|mother_ID), REML=F, data=hmoinvr)
# m1 <- lmerTest::lmer(mother_milk_HMO_3SL_ugml_invr ~ mother_milk_HMO_milk_group + (1|mother_ID) + time_point_numeric, REML=F, data=hmoinvr)
# anova(m0,m1)
# 
# m0 <- lmerTest::lmer(mother_milk_HMO_LNnDFH_ugml_invr ~ mother_milk_HMO_milk_group + (1|mother_ID), REML=F, data=hmoinvr)
# m1 <- lmerTest::lmer(mother_milk_HMO_LNnDFH_ugml_invr ~ mother_milk_HMO_milk_group + (1|mother_ID) + time_point_numeric, REML=F, data=hmoinvr)
# anova(m0,m1)
# 
# m0 <- lmerTest::lmer(mother_milk_HMO_2FL_ugml_invr ~ mother_milk_HMO_milk_group + (1|mother_ID), REML=F, data=hmoinvr)
# m1 <- lmerTest::lmer(mother_milk_HMO_2FL_ugml_invr ~ mother_milk_HMO_milk_group + (1|mother_ID) + time_point_numeric, REML=F, data=hmoinvr)
# anova(m0,m1)


## save results table
# write.table(lmer.results.by.time, file="model_results_HMO_levels_by_milk_group_and_time/240108_lmer.results.by.time.txt", sep="\t", row.names=F, quote=F)


### ===== 11.5 MODEL TO TEST FOR EFFECT OF milk_group (factor) ON HMOs WITH CORRECTION FOR time_point_numeric

## model function
run.lmer.cor.time <- function(datasetname, inputdata, xcolumn, ycolumn){
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
      print(paste0("Start running model for ", colnames(inputdata)[i], " and ", colnames(inputdata)[j]))
      
      my_df <- inputdata[!is.na(inputdata[,i]) & !is.na(inputdata[,j]),]
      
      #first mixed model with correction for time_point as a fixed effect and without the phenotype of interest (i)
      m0 <- lmerTest::lmer(my_df[,j] ~ time_point_numeric + (1|mother_ID), REML=F, data=my_df)
      sum0 <- summary(m0)
      
      #second mixed model with correction for time_point as a fixed effect and WITH the phenotype of interest (i)
      m1 <- lmerTest::lmer(my_df[,j] ~ time_point_numeric + (1|mother_ID) + my_df[,i], REML=F, data=my_df)
      sum1 <- summary(m1)
      
      #compare models and save p value from model comparison
      an1 <- anova(m0, m1)
      
      p_models <- c(p_models, an1["m1", "Pr(>Chisq)"])
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      levels <- c(levels, paste(collapse=";",levels(my_df[,i])))
      n_levels <- c(n_levels, paste(collapse=";",table(my_df[,i])))
      e <- c(e, paste(collapse=";",c(0,signif(digits=3,sum1$coef[3:nrow(sum1$coef), "Estimate"])))) #first rows are for the intercept and numeric time point, ensure that this picks from row 3 on as these are the results for i
      
      stat <- c(stat, paste0("lmer_correction_for_time_ID"))
      
      d <- c(d, datasetname)
      x <- c(x, colnames(inputdata)[i])
      y <- c(y, colnames(inputdata)[j])
      n <- c(n, nrow(my_df))
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, levels=levels, n_levels=n_levels, estimate=e, p_models=p_models)
  
  #add FDR correction
  res$FDR <- p.adjust(res$p_models, method="BH")
  
  #sort results by FDR and p value
  res2 <- res[order(res$FDR, res$p_models),]
  
  return(res2)
}

## run model
lmer.results.by.milkgroup <- run.lmer.cor.time(datasetname="mother",
                                               inputdata=hmoinvr,
                                               xcolumn=c(483), #test the effect of mother_milk_HMO_milk_group (factor)..
                                               ycolumn=c(486:513)) #..on milk HMOs


# ## check that models return correct p values: HMO with milk group and time point
# m0 <- lmerTest::lmer(mother_milk_HMO_3SL_ugml_invr ~ time_point_numeric + (1|mother_ID), REML=F, data=hmoinvr)
# m1 <- lmerTest::lmer(mother_milk_HMO_3SL_ugml_invr ~ time_point_numeric + (1|mother_ID) + mother_milk_HMO_milk_group, REML=F, data=hmoinvr)
# anova(m0,m1)
# 
# m0 <- lmerTest::lmer(mother_milk_HMO_LNnDFH_ugml_invr ~ time_point_numeric + (1|mother_ID), REML=F, data=hmoinvr)
# m1 <- lmerTest::lmer(mother_milk_HMO_LNnDFH_ugml_invr ~ time_point_numeric + (1|mother_ID) + mother_milk_HMO_milk_group, REML=F, data=hmoinvr)
# anova(m0,m1)
# 
# m0 <- lmerTest::lmer(mother_milk_HMO_2FL_ugml_invr ~ time_point_numeric + (1|mother_ID), REML=F, data=hmoinvr)
# m1 <- lmerTest::lmer(mother_milk_HMO_2FL_ugml_invr ~ time_point_numeric + (1|mother_ID) + mother_milk_HMO_milk_group, REML=F, data=hmoinvr)
# anova(m0,m1)

## save results table
# write.table(lmer.results.by.milkgroup, file="model_results_HMO_levels_by_milk_group_and_time/240108_lmer.results.by.milkgroup.txt", sep="\t", row.names=F, quote=F)


## ===== 11.6 COMBINE RESULTS FOR EFFECTS OF TIME POINT AND MILK GROUP ON HMOs, CORRECT FOR MULTIPLE TESTING AND SAVE RESULTS
# 
# ## combine results from the 2 models
# lmer.results.by.time.milk.group <- as.data.frame(rbind(lmer.results.by.time,lmer.results.by.milkgroup))
# 
# ## correct for multiple testing
# lmer.results.by.time.milk.group$FDR <- p.adjust(lmer.results.by.time.milk.group$p_models, method="BH")
# 
# ## resort results by FDR
# lmer.results.by.time.milk.group <- lmer.results.by.time.milk.group[order(lmer.results.by.time.milk.group$FDR, lmer.results.by.time.milk.group$p_models),]
# 
# ## save output
# write.table(lmer.results.by.time.milk.group, file="model_results_HMO_levels_by_milk_group_and_time/240108_lmer_results_by_time_milk_group.txt", sep="\t", row.names=F, quote=F)
# write.table(lmer.results.by.time.milk.group[lmer.results.by.time.milk.group$FDR<0.05,], file="model_results_HMO_levels_by_milk_group_and_time/240108_lmer_results_by_time_milk_group_FDR_below_0.05.txt", sep="\t", row.names=F, quote=F)

