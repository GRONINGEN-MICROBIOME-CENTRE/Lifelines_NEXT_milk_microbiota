####################################################################################################################################
### ASSOCIATE MILKGROUP AND MILK HMOs WITH INFANT FAECAL MICROBIOTA DATA (BETA DIVERSITY) FOR MILK COMPOSITION PAPER (TABLE S32) ###
####################################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessions type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 0. DATA IMPORT
# 
# 1. SELECT DATA OF INTEREST: INFANT FAECAL MICROBIOTA DATA AND MATERNAL MILK HMO DATA
# 
# 2. ENSURE CORRECT DATA STRUCTURE
# 
# 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
# 
# 4. INVERSE-RANK TRANSFORM NUMERIC MILK HMO DATA AND NUMERIC PHENOTYPES
# 4.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 4.2 INVERSE-RANK TRANSFORM MILK HMO DATA AND NUMERIC PHENOTYPES
# 
# 5. CALCULATE INFANT FAECAL RELATIVE ABUNDANCES
# 
# 6. CLR-TRANSFORM INFANT FAECAL RELATIVE ABUNDANCES
# 6.1 FUNCTION FOR CLR-TRANSFORMATION
# 6.2 CLR-TRANSFORM INFANT FAECAL RELATIVE ABUNDANCES
# 
# 7. FUNCTION FOR ADONIS
# 
# 8. ADONIS: CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AND DECIDE ON CORRECTION FACTORS
# 8.1 CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AT M1
# 8.2 CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AT M2
# 8.3 CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AT M3
# 8.4 CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AT M6
# 8.5 CHECK WHICH BASIC PHENOTYPES CONSISTENTLY AFFECTED INFANT FAECAL BETA DIVERSITY AT M1-M6
# 
# 9. ADONIS: CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY (TABLE S32A-D)
# 9.1 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M1 (TABLE S32A)
# 9.2 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M2 (TABLE S32B)
# 9.3 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M3 (TABLE S32C)
# 9.4 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M6 (TABLE S32D)
# 
# 10. ADONIS: CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY (NOT USED)
# 10.1 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M1 - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY
# 10.2 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M2 - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY
# 10.3 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M3 - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY
# 10.4 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M6 - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/")

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


##### =========================== 1. SELECT DATA OF INTEREST: INFANT FAECAL MICROBIOTA DATA AND MATERNAL MILK HMO DATA =========================== #####

## select infant faecal samples and select IDs, (potentially) relevant phenotype columns and bacterial data
inff <- data[data$sample_origin_type=="infant_faeces", c(1:15,504,507,            #IDs, incl. seq ID and sample type
                                                         510,526,                 #DNA isolation method, seq reads
                                                         25:26,473:475,           #maternal Le/Se status (genetic- and HMO-based)
                                                         427:429,                 #infant Se status (genetic)
                                                         394:395,430:434,442:445, #maternal/infant phenotypes (gestational age, infant sex, delivery place and mode, infant feeding)
                                                         476:503,                 #milk HMOs
                                                         531:1154)]               #milk relative bacterial abundances
#552x690


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

## check that data structure is correct
str(inff[,1:50])
str(inff[,c(51:67,689:690)])
# -> data structure is correct


##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

## Infant phenotypes: for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
inff$mother_birth_gestational_age_categories <- factor(inff$mother_birth_gestational_age_categories, levels = c("term", "preterm"))

inff$infant_birth_delivery_mode <- factor(inff$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
inff$infant_birth_delivery_mode_detailed <- factor(inff$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                        "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
inff$infant_ffq_feeding_type_delivery <- factor(inff$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))

# also change it for milk groups
inff$mother_milk_HMO_milk_group <- factor(inff$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


## No phenotypes were changed to numeric.


##### =========================== 4. INVERSE-RANK TRANSFORM NUMERIC MILK HMO DATA AND NUMERIC PHENOTYPES =========================== #####

### ===== 4.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 4.2 INVERSE-RANK TRANSFORM MILK HMO DATA AND NUMERIC PHENOTYPES ===== ###

inff_invr <- inff
inff_invr[,c(29,39:66)] <- as.data.frame(apply(inff_invr[,c(29,39:66)], 2, invrank))
colnames(inff_invr)[c(29,39:66)] <- paste0(colnames(inff_invr)[c(29,39:66)], "_invr")


##### =========================== 5. CALCULATE INFANT FAECAL RELATIVE ABUNDANCES =========================== #####

## save only columns with absolute bacterial abundances
rownames(inff_invr) <- inff_invr$seq_16S_sample_ID
inffbact <- inff_invr[,67:ncol(inff_invr)] #552x624

## calculate relative bacterial abundances
inffbact_relab <- (inffbact/rowSums(inffbact))
table(rowSums(inffbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
inffbact_relab_tmp <- inffbact_relab[,-624]

## exclude columns for bacteria that were not identified in any sample
inffbact_relab_tmp2 <- inffbact_relab_tmp[,colSums(inffbact_relab_tmp)>0]
dim(inffbact_relab_tmp2)

## merge back with metadata
inff_invr_RA <- merge(inff_invr[,1:66], inffbact_relab_tmp2, by="row.names")
rownames(inff_invr_RA) <- inff_invr_RA$Row.names
inff_invr_RA <- inff_invr_RA[,-1]


##### =========================== 6. CLR-TRANSFORM INFANT FAECAL RELATIVE ABUNDANCES =========================== #####

### ===== 6.1 FUNCTION FOR CLR-TRANSFORMATION ===== ###

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


### ===== 6.2 CLR-TRANSFORM INFANT FAECAL RELATIVE ABUNDANCES ===== ###

inff_invr_RAclr <- inff_invr_RA
inff_invr_RAclr[,67:ncol(inff_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(inff_invr_RAclr[,67:ncol(inff_invr_RAclr)], inff_invr_RAclr[,67:ncol(inff_invr_RAclr)]))

colnames(inff_invr_RAclr)[67:ncol(inff_invr_RAclr)] <- c(paste0(colnames(inff_invr_RAclr)[67:ncol(inff_invr_RAclr)], "_clr"))


##### =========================== 7. FUNCTION FOR ADONIS =========================== #####

library(vegan)
library(ggplot2)

## function to run Adonis with correction for DNA isolation batch and clean read counts
run.adonis2 <- function(mydata, timepoint, covardetails, covariates, ncovar, columns, outputfolder){
  column <- c()
  R2 <- c()
  p <- c()
  n <- c()
  tp <- c()
  
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
    tp <- c(tp, timepoint)
    
    print("Calculating variance explained")
    beta.cmd <- cmdscale(as.matrix(ait), k=2, eig=T)
    PCoA1 <- beta.cmd$eig[1]/sum(beta.cmd$eig)*100
    PCoA2 <- beta.cmd$eig[2]/sum(beta.cmd$eig)*100
    print(paste0("PCoA1 explains ", PCoA1, "% of the variance"))
    print(paste0("PCoA2 explains ", PCoA2, "% of the variance"))
    
    ## ADONIS
    print("Running adonis")
    a <- table(new_df[,i])/nrow(new_df)
    
    if (any(a>0.97)){ #set to NA if >97% of the data is the same (this includes if there is only 1 factor because that is 100% the same)
      R2 <- c(R2, NA)
      p <- c(p, NA)
    }else{ ##only run adonis if there is sufficient variation in the answers (at least 3% need to have different answer), this also means that >=2 levels are present for the variable of interest
      my.adonis.function <- paste0("ait ~ ", covariates, "new_df[,i]")
      print(paste0("Adonis function: ", my.adonis.function))
      a1 <- adonis2(as.formula(my.adonis.function), by="margin")
      print(a1)
      R2 <- c(R2, a1$R2[ncovar+1]) #we pick the value x+1 as the first x are for the factors we correct for
      p <- c(p, a1$'Pr(>F)'[ncovar+1]) #we pick the value x+1 as the first x are for the factors we correct for
    }
    
    # PCoA plot for X1 (PC1) and X2 (PC2)
    if (class(new_df[,i]) == "numeric" | class(new_df[,i]) == "integer"){
      print("Plotting PCoA without ellipses")
      pcoa1 <- ggplot(new_df, aes(x=X1, y=X2, color=new_df[,i]))+
        geom_point(alpha=0.5)+
        labs(x=paste0("PC1 (", signif(PCoA1,5), "%)"), y=paste0("PC2 (", signif(PCoA2,5), "%)"), color=paste0(colnames(new_df)[i]),
             title=paste0("Beta diversity by ", colnames(new_df)[i], "\n",
                          timepoint, ", n=", nrow(new_df), "\n", 
                          covardetails, "\n",
                          " R2=", signif(R2[length(R2)],3), " and p=", signif(p[length(p)],3)))+
        theme_bw()
    }else{
      print("Plotting PCoA with ellipses")
      pcoa1 <- ggplot(new_df, aes(x=X1, y=X2, color=new_df[,i]))+
        geom_point(alpha=0.5)+
        stat_ellipse(aes(group=new_df[,i]))+
        labs(x=paste0("PC1 (", signif(PCoA1,5), "%)"), y=paste0("PC2 (", signif(PCoA2,5), "%)"), color=paste0(colnames(new_df)[i]),
             title=paste0("Beta diversity by ", colnames(new_df)[i], "\n",
                          timepoint, ", n=", nrow(new_df), "\n", 
                          covardetails, "\n",
                          " R2=", signif(R2[length(R2)],3), " and p=", signif(p[length(p)],3)))+
        theme_bw()
    }
    
    if (p[length(p)]<0.05){
      ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/beta_diversity/", outputfolder, "/", timepoint, "/nominally_significant/240216_infant_faeces_PCoA_PC1_PC2_by_", colnames(mydata)[i], "_", covardetails, "_", timepoint, ".pdf"), pcoa1, dpi=500, width=20, height=15, units="cm")
    }else{
      ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/beta_diversity/", outputfolder, "/", timepoint, "/not_significant/240216_infant_faeces_PCoA_PC1_PC2_by_", colnames(mydata)[i], "_", covardetails, "_", timepoint, ".pdf"), pcoa1, dpi=500, width=20, height=15, units="cm")
    }
    
  }
  
  print("Preparing results data frame")
  b <- data.frame(column=column, n=n, R2=R2, pvalue=p, time_point=tp)
  
  #correct for multiple testing
  b$FDR <- p.adjust(b$p, method="BH")
  
  #resort results based on FDR and p
  b <- b[order(b$FDR, b$p),]
  
  #add info on model run
  b$model <- paste0(covardetails)
  
  ##resort columns
  b2 <- b[,c(5,1:4,6,7)]
  
  ##save output table
  write.table(b2, file=paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/beta_diversity/", outputfolder, "/", timepoint , "/240216_beta_diversity_", covardetails, "_results_", timepoint, ".txt"), row.names=F, col.names=T, sep="\t", quote=F)
  
  ##return output table
  return(b2)
}


##### =========================== 8. ADONIS: CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AND DECIDE ON CORRECTION FACTORS =========================== #####

### ===== 8.1 CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AT M1 ===== ###

adonis2_results_by_basic_phenos_M1 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="1_month",],
                                                  timepoint = "M1",
                                                  covardetails = "adonis2_without_correction",
                                                  covariates = "",
                                                  ncovar = 0,
                                                  columns = c(18,19,29,30,31,33,37),
                                                  outputfolder = "PCoA_basic_phenotypes")


### ===== 8.2 CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AT M2 ===== ###

adonis2_results_by_basic_phenos_M2 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="2_months",],
                                                  timepoint = "M2",
                                                  covardetails = "adonis2_without_correction",
                                                  covariates = "",
                                                  ncovar = 0,
                                                  columns = c(18,19,29,30,31,33,37),
                                                  outputfolder = "PCoA_basic_phenotypes")


### ===== 8.3 CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AT M3 ===== ###

adonis2_results_by_basic_phenos_M3 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="3_months",],
                                                  timepoint = "M3",
                                                  covardetails = "adonis2_without_correction",
                                                  covariates = "",
                                                  ncovar = 0,
                                                  columns = c(18,19,29,30,31,33,37),
                                                  outputfolder = "PCoA_basic_phenotypes")


### ===== 8.4 CHECK EFFECT OF BASIC PHENOTYPES ON INFANT FAECAL BETA DIVERSITY AT M6 ===== ###

adonis2_results_by_basic_phenos_M6 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="6_months",],
                                                  timepoint = "M6",
                                                  covardetails = "adonis2_without_correction",
                                                  covariates = "",
                                                  ncovar = 0,
                                                  columns = c(18,19,29,30,31,33,37),
                                                  outputfolder = "PCoA_basic_phenotypes")


### ===== 8.5 CHECK WHICH BASIC PHENOTYPES CONSISTENTLY AFFECTED INFANT FAECAL BETA DIVERSITY AT M1-M6 ===== ###

adonis2_results_by_basic_phenos_M1[adonis2_results_by_basic_phenos_M1$p<0.05,"column"]
# "infant_birth_delivery_place"            
# "infant_birth_delivery_mode"             
# "mother_birth_gestational_age_weeks_invr"
# "DNA_isolation_method"                   
# "seq_16S_n_reads_clean" 

adonis2_results_by_basic_phenos_M2[adonis2_results_by_basic_phenos_M2$p<0.05,"column"]
# seq_16S_n_reads_clean"

adonis2_results_by_basic_phenos_M3[adonis2_results_by_basic_phenos_M3$p<0.05,"column"]
# "DNA_isolation_method"
# "seq_16S_n_reads_clean"      
# "infant_birth_delivery_place"
# "infant_birth_delivery_mode"

adonis2_results_by_basic_phenos_M6[adonis2_results_by_basic_phenos_M6$p<0.05,"column"]
# "DNA_isolation_method"
# "infant_birth_delivery_place"

## -> for the associations at each time point, correct for DNA_isolation_method, seq_16S_n_reads_clean and infant_birth_delivery_mode


##### =========================== 9. ADONIS: CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY (TABLE S32A-D) =========================== #####

### ===== 9.1 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M1 (TABLE S32A) ===== ###

adonis2_results_by_milkgroup_milk_HMOs_M1 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="1_month" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode),],
                                                         timepoint = "M1",
                                                         covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode",
                                                         covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + ",
                                                         ncovar = 3,
                                                         columns = c(22,39:66),
                                                         outputfolder = "PCoA_milkgroup_HMOs")

## check results for M1 when additionally correcting for mother_birth_gestational_age_weeks_invr
adonis2_results_by_milkgroup_milk_HMOs_M1_corGA <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="1_month" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode) & !is.na(inff_invr_RAclr$mother_birth_gestational_age_weeks_invr),],
                                                               timepoint = "M1",
                                                               covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode_GA",
                                                               covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + new_df$mother_birth_gestational_age_weeks_invr + ",
                                                               ncovar = 4,
                                                               columns = c(22,39:66),
                                                               outputfolder = "PCoA_milkgroup_HMOs")


### ===== 9.2 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M2 (TABLE S32B) ===== ###

adonis2_results_by_milkgroup_milk_HMOs_M2 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="2_months" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode),],
                                                         timepoint = "M2",
                                                         covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode",
                                                         covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + ",
                                                         ncovar = 3,
                                                         columns = c(22,39:66),
                                                         outputfolder = "PCoA_milkgroup_HMOs")


### ===== 9.3 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M3 (TABLE S32C) ===== ###

adonis2_results_by_milkgroup_milk_HMOs_M3 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="3_months" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode),],
                                                         timepoint = "M3",
                                                         covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode",
                                                         covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + ",
                                                         ncovar = 3,
                                                         columns = c(22,39:66),
                                                         outputfolder = "PCoA_milkgroup_HMOs")


### ===== 9.4 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M6 (TABLE S32D) ===== ###

adonis2_results_by_milkgroup_milk_HMOs_M6 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="6_months" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode),],
                                                         timepoint = "M6",
                                                         covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode",
                                                         covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + ",
                                                         ncovar = 3,
                                                         columns = c(22,39:66),
                                                         outputfolder = "PCoA_milkgroup_HMOs")


##### =========================== 10. ADONIS: CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY (NOT USED) =========================== #####

### ===== 10.1 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M1 - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY ===== ###

adonis2_results_by_milkgroup_milk_HMOs_BF_M1 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="1_month" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode) & !is.na(inff_invr_RAclr$infant_ffq_feeding_mode) & inff_invr_RAclr$infant_ffq_feeding_mode=="breastfeeding",],
                                                         timepoint = "M1",
                                                         covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode_BF_only",
                                                         covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + ",
                                                         ncovar = 3,
                                                         columns = c(22,39:66),
                                                         outputfolder = "PCoA_milkgroup_HMOs_BF_only")

## check results for M1 when additionally correcting for mother_birth_gestational_age_weeks_invr
adonis2_results_by_milkgroup_milk_HMOs_BF_M1_corGA <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="1_month" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode) & !is.na(inff_invr_RAclr$mother_birth_gestational_age_weeks_invr) & !is.na(inff_invr_RAclr$infant_ffq_feeding_mode) & inff_invr_RAclr$infant_ffq_feeding_mode=="breastfeeding",],
                                                               timepoint = "M1",
                                                               covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode_GA_BF_only",
                                                               covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + new_df$mother_birth_gestational_age_weeks_invr + ",
                                                               ncovar = 4,
                                                               columns = c(22,39:66),
                                                               outputfolder = "PCoA_milkgroup_HMOs_BF_only")


### ===== 10.2 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M2 - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY ===== ###

adonis2_results_by_milkgroup_milk_HMOs_BF_M2 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="2_months" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode) & !is.na(inff_invr_RAclr$infant_ffq_feeding_mode) & inff_invr_RAclr$infant_ffq_feeding_mode=="breastfeeding",],
                                                         timepoint = "M2",
                                                         covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode_BF_only",
                                                         covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + ",
                                                         ncovar = 3,
                                                         columns = c(22,39:66),
                                                         outputfolder = "PCoA_milkgroup_HMOs_BF_only")


### ===== 10.3 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M3 - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY ===== ###

adonis2_results_by_milkgroup_milk_HMOs_BF_M3 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="3_months" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode) & !is.na(inff_invr_RAclr$infant_ffq_feeding_mode) & inff_invr_RAclr$infant_ffq_feeding_mode=="breastfeeding",],
                                                         timepoint = "M3",
                                                         covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode_BF_only",
                                                         covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + ",
                                                         ncovar = 3,
                                                         columns = c(22,39:66),
                                                         outputfolder = "PCoA_milkgroup_HMOs_BF_only")


### ===== 10.4 CHECK EFFECT OF MATERNAL MILK GROUP AND MILK HMOs ON INFANT FAECAL BETA DIVERSITY AT M6 - IN EXCLUSIVELY BREAST MILK FED INFANTS ONLY ===== ###

adonis2_results_by_milkgroup_milk_HMOs_BF_M6 <- run.adonis2(mydata = inff_invr_RAclr[inff_invr_RAclr$time_point=="6_months" & !is.na(inff_invr_RAclr$infant_birth_delivery_mode) & !is.na(inff_invr_RAclr$infant_ffq_feeding_mode) & inff_invr_RAclr$infant_ffq_feeding_mode=="breastfeeding",],
                                                         timepoint = "M6",
                                                         covardetails = "adonis2_with_correction_for_DNAmethod_reads_delmode_BF_only",
                                                         covariates = "new_df$DNA_isolation_method + new_df$seq_16S_n_reads_clean + new_df$infant_birth_delivery_mode + ",
                                                         ncovar = 3,
                                                         columns = c(22,39:66),
                                                         outputfolder = "PCoA_milkgroup_HMOs_BF_only")


