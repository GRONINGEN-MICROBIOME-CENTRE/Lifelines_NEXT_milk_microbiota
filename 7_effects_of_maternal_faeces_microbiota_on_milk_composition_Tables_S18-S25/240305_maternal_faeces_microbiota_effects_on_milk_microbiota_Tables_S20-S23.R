#######################################################################################################################
### ANALYSE EFFECT OF MATERNAL FAECAL MICROBIOTA ON MILK MICROBIOTA (ALPHA DIV AND REL. ABUNDANCES)  ##################
#######################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessions type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 0. DATA IMPORT
# 
# 1. SELECT DATA OF INTEREST
# 1.1 MATERNAL FAECAL MICROBIOTA DATA
# 1.2 MILK MICROBIOTA DATA
# 1.3 SELECT ONLY MOTHERS WITH BOTH FAECAL AND MILK MICROBIOTA DATA
# 
# 2. ENSURE CORRECT DATA STRUCTURE
# 
# 3. INVERSE-RANK TRANSFORM ALPHA DIVERSITY MEASURES
# 3.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 3.2 INVERSE-RANK TRANSFORM MATERNAL FAECAL ALPHA DIVERSITY DATA
# 3.3 INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY DATA
# 
# 4. CALCULATE RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE
# 4.1 CALCULATE MATERNAL FAECAL RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA
# 4.2 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA
# 
# 5. CLR-TRANSFORM RELATIVE ABUNDANCES
# 5.1 FUNCTION FOR CLR-TRANSFORMATION
# 5.2 CLR-TRANSFORM MATERNAL FAECAL RELATIVE ABUNDANCES
# 5.3 CLR-TRANSFORM MILK RELATIVE ABUNDANCES
# 
# 6. COMBINE MATERNAL FAECAL AND MILK MICROBIOTA DATA INTO 1 DATA FRAME
# 
# 7. MODEL FOR ASSOCIATION OF SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK MICROBIOTA WITH CORRECTION FOR MILK DNA ISOLATION BATCH, MILK SEQUENCING READS AND FAECAL SEQUENCING READS
# 7.1 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK MICROBIOTA WITH CORRECTION FOR MILK DNA ISOLATION BATCH, MILK SEQUENCING READS AND FAECAL SEQUENCING READS
# 7.2 RUN ASSOCIATION OF MATERNAL SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK MICROBIOTA DATA
# 
# 8. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING
# 8.1 EXCLUDE RESULTS FOR SIMPSON DIVERSITY
# 8.2 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES
# 8.3 CORRECT FOR MULTIPLE TESTING
# 8.4 SAVE RESULTS (TABLES S20, S21, S22, S23)


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/")

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

### ===== 1.1 MATERNAL FAECAL MICROBIOTA DATA ===== ###

## select (1) maternal faecal samples and (2) every mother only 1x (i.e. remove duplication of data generated for twin babies) and select IDs and maternal faecal microbiota data
matf_tmp <- data[data$sample_origin_type=="mother_faeces",]
matf_tmp2 <- matf_tmp[grep(";1", matf_tmp$mother_sample_ID), c(1:9,11:15,504,526,528:1154)] #170x643

## Maternal faecal data set (matf_tmp2)
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16   = seq_16S_n_reads_clean (maternal faeces reads)
## 17:19  = maternal faecal alpha diversity measures
## 20:643 = maternal faecal bacterial relative abundances


### ===== 1.2 MILK MICROBIOTA DATA ===== ###

## select (1) human milk samples from M3 and (2) every mother only 1x (i.e. remove duplication of data generated for twin babies) and select IDs and maternal milk microbiota data
milk_tmp <- data[data$sample_origin_type=="mother_human_milk" & data$time_point=="3_months",]
milk_tmp2 <- milk_tmp[grep(";1", milk_tmp$mother_sample_ID), c(1:9,11:15,504,511,526,528:1154)] #268x644

## Maternal milk data set (milk_tmp2)
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:17  = DNA_isolation_batch and seq_16S_n_reads_clean (for milk)
## 18:20  = maternal faecal alpha diversity measures
## 21:644 = maternal faecal bacterial relative abundances


### ===== 1.3 SELECT ONLY MOTHERS WITH BOTH FAECAL AND MILK MICROBIOTA DATA ===== ###

motherIDs <- intersect(matf_tmp2$mother_sample_ID, milk_tmp2$mother_sample_ID) #164 mothers

matf <- matf_tmp2[matf_tmp2$mother_sample_ID %in% motherIDs,] #164x643
milk <- milk_tmp2[milk_tmp2$mother_sample_ID %in% motherIDs,] #164x644


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

## check that data structure is correct
str(matf[,c(1:20,642:643)])
# -> data structure is correct

## check that data structure is correct
str(milk[,c(1:21,643:644)])
# -> data structure is correct


##### =========================== 3. INVERSE-RANK TRANSFORM ALPHA DIVERSITY MEASURES =========================== #####

### ===== 3.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 3.2 INVERSE-RANK TRANSFORM MATERNAL FAECAL ALPHA DIVERSITY DATA ===== ###

matf_invr <- matf
matf_invr[,grep("alpha_div",colnames(matf_invr))] <- as.data.frame(apply(matf_invr[,grep("alpha_div",colnames(matf_invr))], 2, invrank))

## fix colnames of all invr-transformed columns
colnames(matf_invr)[17:19] <- paste0(colnames(matf_invr)[17:19], "_invr")


### ===== 3.3 INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY DATA ===== ###

milk_invr <- milk
milk_invr[,grep("alpha_div",colnames(milk_invr))] <- as.data.frame(apply(milk_invr[,grep("alpha_div",colnames(milk_invr))], 2, invrank))

## fix colnames of all invr-transformed columns
colnames(milk_invr)[18:20] <- paste0(colnames(milk_invr)[18:20], "_invr")


##### =========================== 4. CALCULATE RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE =========================== #####

### ===== 4.1 CALCULATE MATERNAL FAECAL RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA ===== ###

## save only columns with absolute bacterial abundances
rownames(matf_invr) <- matf_invr$seq_16S_sample_ID
matfbact <- matf_invr[,20:ncol(matf_invr)] #164x624

## calculate relative bacterial abundances
matfbact_relab <- (matfbact/rowSums(matfbact))
table(rowSums(matfbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
matfbact_relab_tmp <- matfbact_relab[,-624]

## select only bacterial genera with ≥10% prevalence as for other analyses
matf_RAprev10_sum_M3 <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/relative_abundances/231215_maternal_faeces_RA_prev10_summary_statistics.txt", header=T, sep="\t")
matf_RAprev10_bact <- c(paste0("g__", matf_RAprev10_sum_M3$bacterium)) #83

matfbact_relab_tmp2 <- matfbact_relab_tmp[,colnames(matfbact_relab_tmp) %in% matf_RAprev10_bact]
length(matfbact_relab_tmp2) #83
# colnames(matfbact_relab_tmp2)

## merge back with metadata
matf_invr_RA <- merge(matf_invr[,1:19], matfbact_relab_tmp2, by="row.names")
rownames(matf_invr_RA) <- matf_invr_RA$Row.names
matf_invr_RA <- matf_invr_RA[,-1]
#164x102


### ===== 4.2 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA ===== ###

## save only columns with absolute bacterial abundances
rownames(milk_invr) <- milk_invr$seq_16S_sample_ID
milkbact <- milk_invr[,21:ncol(milk_invr)] #164x624

## calculate relative bacterial abundances
milkbact_relab <- (milkbact/rowSums(milkbact))
table(rowSums(milkbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
milkbact_relab_tmp <- milkbact_relab[,-624]

## select only bacterial genera with ≥10% prevalence as for other analyses
milk_RAprev10_sum <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/231215_milk_RA_prev10_summary_statistics.txt", header=T, sep="\t")
milk_RAprev10_sum_M3 <- milk_RAprev10_sum[milk_RAprev10_sum$time_point=="3_months",]
milk_RAprev10_bact <- c(paste0("g__", milk_RAprev10_sum_M3$bacterium)) #36

milkbact_relab_tmp2 <- milkbact_relab_tmp[,colnames(milkbact_relab_tmp) %in% milk_RAprev10_bact]
length(milkbact_relab_tmp2) #36
# colnames(milkbact_relab_tmp2)

## merge back with metadata
milk_invr_RA <- merge(milk_invr[,1:20], milkbact_relab_tmp2, by="row.names")
rownames(milk_invr_RA) <- milk_invr_RA$Row.names
milk_invr_RA <- milk_invr_RA[,-1]
#164x56


##### =========================== 5. CLR-TRANSFORM RELATIVE ABUNDANCES =========================== #####

### ===== 5.1 FUNCTION FOR CLR-TRANSFORMATION ===== ###

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


### ===== 5.2 CLR-TRANSFORM MATERNAL FAECAL RELATIVE ABUNDANCES ===== ###

matf_invr_RAclr <- matf_invr_RA
matf_invr_RAclr[,20:ncol(matf_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(matf_invr_RAclr[,20:ncol(matf_invr_RAclr)], matf_invr_RAclr[,20:ncol(matf_invr_RAclr)]))

## fix colnames of clr-transformed columns
colnames(matf_invr_RAclr)[20:ncol(matf_invr_RAclr)] <- c(paste0(colnames(matf_invr_RAclr)[20:ncol(matf_invr_RAclr)], "_clr"))

## fix rownames
rownames(matf_invr_RAclr) <- matf_invr_RAclr$mother_sample_ID

## Maternal data set (matf_invr_RAclr), 164x102
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16   = maternal faecal seq_16S_n_reads_clean
## 17:19  = invr-transformed maternal faecal alpha diversity measures
## 20:102 = clr-transformed maternal faecal relative abundances of bacterial genera with ≥10% prevalence


### ===== 5.3 CLR-TRANSFORM MILK RELATIVE ABUNDANCES ===== ###

milk_invr_RAclr <- milk_invr_RA
milk_invr_RAclr[,21:ncol(milk_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(milk_invr_RAclr[,21:ncol(milk_invr_RAclr)], milk_invr_RAclr[,21:ncol(milk_invr_RAclr)]))

## fix colnames of clr-transformed columns
colnames(milk_invr_RAclr)[21:ncol(milk_invr_RAclr)] <- c(paste0(colnames(milk_invr_RAclr)[21:ncol(milk_invr_RAclr)], "_clr"))

## fix rownames
rownames(milk_invr_RAclr) <- milk_invr_RAclr$mother_sample_ID

## Maternal data set (milk_invr_RAclr), 164x56
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:17  = milk DNA_isolation_batch and milk seq_16S_n_reads_clean
## 18:20  = invr-transformed milk alpha diversity measures
## 21:56 = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence


##### =========================== 6. COMBINE MATERNAL FAECAL AND MILK MICROBIOTA DATA INTO 1 DATA FRAME =========================== #####

## ensure all columns have unique names, indicating from which data set they derived
colnames(matf_invr_RAclr)[c(2:102)] <- c(paste0("matf_", colnames(matf_invr_RAclr)[c(2:102)]))
colnames(milk_invr_RAclr)[c(2:56)] <- c(paste0("milk_", colnames(milk_invr_RAclr)[c(2:56)]))

matfmilk <- merge(matf_invr_RAclr, milk_invr_RAclr, by="mother_sample_ID")


##### =========================== 7. MODEL FOR ASSOCIATION OF SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK MICROBIOTA WITH CORRECTION FOR MILK DNA ISOLATION BATCH, MILK SEQUENCING READS AND FAECAL SEQUENCING READS =========================== #####

### ===== 7.1 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK MICROBIOTA WITH CORRECTION FOR MILK DNA ISOLATION BATCH, MILK SEQUENCING READS AND FAECAL SEQUENCING READS ===== ###

run.lmer.cor.DNAbatch.mreads.freads.1tp <- function(datasetname, inputdata, xcolumn, ycolumn){
  
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
      
      #mixed model with correction for milk_group and seq_16S_n_reads_clean as fixed effects and WITH the phenotype of interest (i)
      m1 <- lm(my_df[,j] ~ milk_DNA_isolation_batch + milk_seq_16S_n_reads_clean + matf_seq_16S_n_reads_clean + my_df[,i], data=my_df)
      sum1 <- summary(m1)
      print(sum1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      e <- c(e, c(sum1$coefficients[grep("my_df",rownames(sum1$coefficients)), "Estimate"])) #first rows are for mother_milk_HMO_milk_group and seq_16S_n_reads_clean, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      p <- c(p, c(sum1$coefficients[grep("my_df",rownames(sum1$coefficients)), "Pr(>|t|)"])) #first rows are for mother_milk_HMO_milk_group and seq_16S_n_reads_clean, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, paste0("lm_1_time_point_correction_DNAbatch_milkreads_faecalreads"))
      d <- c(d, datasetname)
      x <- c(x, colnames(inputdata)[i])
      y <- c(y, colnames(inputdata)[j])
      n <- c(n, nrow(my_df))
    }
  }
  
  #save results in data frame
  res <- data.frame(dataset=d, statistic=stat, x=x, y=y, n_total=n, estimate=e, p=p)
  
  return(res)
}


### ===== 7.2 RUN ASSOCIATION OF MATERNAL SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK MICROBIOTA DATA ===== ###

## Maternal data set (matfmilk), 164x157
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point) from faecal data frame
## 16   = maternal faecal seq_16S_n_reads_clean
## 17:19  = invr-transformed maternal faecal alpha diversity measures
## 20:102 = clr-transformed maternal faecal relative abundances of bacterial genera with ≥10% prevalence
## 103:116 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point) from milk data frame
## 117:118  = milk DNA_isolation_batch and milk seq_16S_n_reads_clean
## 119:121  = invr-transformed milk alpha diversity measures
## 122:157 = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence

matfmilk_results <- run.lmer.cor.DNAbatch.mreads.freads.1tp(datasetname = "mother_maternal_microbiota_milk_microbiota_1_time_point",
                                                            inputdata   = matfmilk,
                                                            xcolumn     = c(17:19,20:102),     #invr-transformed maternal faecal alpha diversity and clr-transformed maternal faecal relative abundances of bacterial genera with ≥10% prevalence
                                                            ycolumn     = c(119:121,122:157))  #invr-transformed milk alpha diversity and clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence
# write.table(matfmilk_results, file="association_model_results/240305_milk_microbiota_by_maternal_faecal_microbiota_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

## re-import results
# setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/")
# matfmilk_results <- read.table(file="association_model_results/240305_milk_microbiota_by_maternal_faecal_microbiota_model_results.txt", header=T, sep="\t")


##### =========================== 8. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING =========================== #####

### ===== 8.1 EXCLUDE RESULTS FOR SIMPSON DIVERSITY ===== ###

matfmilk_results_no_simpson <- matfmilk_results[matfmilk_results$x!="matf_alpha_div_simpson_invr" & matfmilk_results$y!="milk_alpha_div_simpson_invr",] #3230x7


### ===== 8.2 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES ===== ###

## milk alpha diversity as outcome
milk_alpha_results <- matfmilk_results_no_simpson[grep("milk_alpha_div", matfmilk_results_no_simpson$y),] #170x7
milk_alpha_div_by_matf_alpha_div_results <- milk_alpha_results[grep("matf_alpha_div", milk_alpha_results$x),] #4x7
milk_alpha_div_by_matf_RA_results <- milk_alpha_results[grep("matf_g__", milk_alpha_results$x),] #166x7

## milk relative abundances as outcome
milk_RA_results <- matfmilk_results_no_simpson[grep("milk_g__", matfmilk_results_no_simpson$y),] #3060x7
milk_RA_div_by_matf_alpha_div_results <- milk_RA_results[grep("matf_alpha_div", milk_RA_results$x),] #72x7
milk_RA_div_by_matf_RA_results <- milk_RA_results[grep("matf_g__", milk_RA_results$x),] #2988x7


### ===== 8.3 CORRECT FOR MULTIPLE TESTING ===== ###

## FDR correction
milk_alpha_div_by_matf_alpha_div_results$FDR <- p.adjust(milk_alpha_div_by_matf_alpha_div_results$p, method="BH")
milk_alpha_div_by_matf_RA_results$FDR <- p.adjust(milk_alpha_div_by_matf_RA_results$p, method="BH")

milk_RA_div_by_matf_alpha_div_results$FDR <- p.adjust(milk_RA_div_by_matf_alpha_div_results$p, method="BH")
milk_RA_div_by_matf_RA_results$FDR <- p.adjust(milk_RA_div_by_matf_RA_results$p, method="BH")

## sort results by FDR and p
milk_alpha_div_by_matf_alpha_div_results_ord <- milk_alpha_div_by_matf_alpha_div_results[order(milk_alpha_div_by_matf_alpha_div_results$FDR, milk_alpha_div_by_matf_alpha_div_results$p),]
milk_alpha_div_by_matf_RA_results_ord <- milk_alpha_div_by_matf_RA_results[order(milk_alpha_div_by_matf_RA_results$FDR, milk_alpha_div_by_matf_RA_results$p),]

milk_RA_div_by_matf_alpha_div_results_ord <- milk_RA_div_by_matf_alpha_div_results[order(milk_RA_div_by_matf_alpha_div_results$FDR, milk_RA_div_by_matf_alpha_div_results$p),]
milk_RA_div_by_matf_RA_results_ord <- milk_RA_div_by_matf_RA_results[order(milk_RA_div_by_matf_RA_results$FDR, milk_RA_div_by_matf_RA_results$p),]

## -> No significant effects of maternal faecal alpha diversity or maternal faecal relative bacterial abundances on human milk alpha diversity and human milk relative bacterial abundances.


### ===== 8.4 SAVE RESULTS (TABLES S20, S21, S22, S23) ===== ###

write.table(milk_alpha_div_by_matf_alpha_div_results_ord, file="association_model_results/240305_milk_alpha_diversity_by_maternal_faecal_alpha_diversity_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(milk_alpha_div_by_matf_RA_results_ord, file="association_model_results/240305_milk_alpha_diversity_by_maternal_faecal_relative_abundances_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

write.table(milk_RA_div_by_matf_alpha_div_results_ord, file="association_model_results/240305_milk_relative_abundances_by_maternal_faecal_alpha_diversity_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(milk_RA_div_by_matf_RA_results_ord, file="association_model_results/240305_milk_relative_abundances_by_maternal_faecal_relative_abundances_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)




