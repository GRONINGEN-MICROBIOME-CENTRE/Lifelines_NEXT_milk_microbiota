#######################################################################################################################
### ANALYSE EFFECT OF MATERNAL FAECAL MICROBIOTA ON MILK HMOs  ########################################################
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
# 1.1 MATERNAL FAECAL MICROBIOTA DATA AND MILK HMO DATA
# 
# 2. ENSURE CORRECT DATA STRUCTURE
# 
# 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
# 
# 4. INVERSE-RANK TRANSFORM MILK HMO DATA AND MATERNAL FAECAL ALPHA DIVERSITY MEASURES
# 4.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 4.2 INVERSE-RANK TRANSFORM MILK HMO DATA
# 4.3 INVERSE-RANK TRANSFORM MATERNAL FAECAL ALPHA DIVERSITY DATA
# 
# 5. CALCULATE MATERNAL FAECAL RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE
# 5.1 CALCULATE MATERNAL FAECAL RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA
# 
# 6. CLR-TRANSFORM MATERNAL FAECAL RELATIVE ABUNDANCES
# 6.1 FUNCTION FOR CLR-TRANSFORMATION
# 6.2 CLR-TRANSFORM MATERNAL FAECAL RELATIVE ABUNDANCES
# 
# 7. MODEL FOR ASSOCIATION OF SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK HMOS WITH CORRECTION FOR MILK GROUP AND SEQUENCING READS
# 7.1 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK HMOS WITH CORRECTION FOR MILK GROUP AND SEQUENCING READS
# 7.2 RUN ASSOCIATION OF MATERNAL SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK HMOS
# 
# 8. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING
# 8.1 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES
# 8.2 CORRECT FOR MULTIPLE TESTING
# 8.3 SAVE RESULTS (TABLE S25)
# 8.4 EXCLUDE SIMPSON DIVERSITY, RECALCULATE FDR AND SAVE RESULTS AGAIN (TABLE S24)


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

### ===== 1.1 MATERNAL FAECAL MICROBIOTA DATA AND MILK HMO DATA ===== ###

## select (1) maternal faecal samples and (2) every mother only 1x (i.e. remove duplication of data generated for twin babies) and select IDs, maternal faecal microbiota and milk HMO data
matf_tmp <- data[data$sample_origin_type=="mother_faeces",]
matf <- matf_tmp[grep(";1", matf_tmp$mother_sample_ID), c(1:9,11:15,504,473,526,476:503,528:1154)] #170x672

## Maternal data set (matf)
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:17  = mother_milk_HMO_milk_group and seq_16S_n_reads_clean
## 18:45  = milk HMOs
## 46:48  = maternal faecal alpha diversity measures
## 49:672 = maternal faecal bacterial relative abundances


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

## check that data structure is correct
str(matf[,1:50])
str(matf[,671:672])
# -> data structure is correct


##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

# change reference group for milk groups
matf$mother_milk_HMO_milk_group <- factor(matf$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 4. INVERSE-RANK TRANSFORM MILK HMO DATA AND MATERNAL FAECAL ALPHA DIVERSITY MEASURES =========================== #####

### ===== 4.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 4.2 INVERSE-RANK TRANSFORM MILK HMO DATA ===== ###

matf_invr <- matf
matf_invr[,grep("_ugml", colnames(matf_invr))] <- as.data.frame(apply(matf_invr[,grep("_ugml", colnames(matf_invr))], 2, invrank))


### ===== 4.3 INVERSE-RANK TRANSFORM MATERNAL FAECAL ALPHA DIVERSITY DATA ===== ###

matf_invr[,grep("alpha_div",colnames(matf_invr))] <- as.data.frame(apply(matf_invr[,grep("alpha_div",colnames(matf_invr))], 2, invrank))

## fix colnames of all invr-transformed columns
colnames(matf_invr)[18:48] <- paste0(colnames(matf_invr)[18:48], "_invr")


##### =========================== 5. CALCULATE MATERNAL FAECAL RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE =========================== #####

### ===== 5.1 CALCULATE MATERNAL FAECAL RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA ===== ###

## save only columns with absolute bacterial abundances
rownames(matf_invr) <- matf_invr$seq_16S_sample_ID
matfbact <- matf_invr[,49:ncol(matf_invr)] #725x624

## calculate relative bacterial abundances
matfbact_relab <- (matfbact/rowSums(matfbact))
table(rowSums(matfbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
matfbact_relab_tmp <- matfbact_relab[,-624]

## select only bacterial genera with ≥10% prevalence
matfbact_relab_tmp2 <- matfbact_relab_tmp[,colSums(matfbact_relab_tmp>0)>=17]
length(matfbact_relab_tmp2) #36
# colnames(matfbact_relab_tmp2)

## merge back with metadata
matf_invr_RA <- merge(matf_invr[,1:48], matfbact_relab_tmp2, by="row.names")
rownames(matf_invr_RA) <- matf_invr_RA$Row.names
matf_invr_RA <- matf_invr_RA[,-1]


##### =========================== 6. CLR-TRANSFORM MATERNAL FAECAL RELATIVE ABUNDANCES =========================== #####

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


### ===== 6.2 CLR-TRANSFORM MATERNAL FAECAL RELATIVE ABUNDANCES ===== ###

matf_invr_RAclr <- matf_invr_RA
matf_invr_RAclr[,49:ncol(matf_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(matf_invr_RAclr[,49:ncol(matf_invr_RAclr)], matf_invr_RAclr[,49:ncol(matf_invr_RAclr)]))

## fix colnames of clr-transformed columns
colnames(matf_invr_RAclr)[49:ncol(matf_invr_RAclr)] <- c(paste0(colnames(matf_invr_RAclr)[49:ncol(matf_invr_RAclr)], "_clr"))

## fix rownames
rownames(matf_invr_RAclr) <- matf_invr_RAclr$mother_sample_ID

## Maternal data set (matf_invr_RAclr), 170x131
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:17  = mother_milk_HMO_milk_group and seq_16S_n_reads_clean
## 18:45  = invr-transformed milk HMOs
## 46:48  = invr-transformed maternal faecal alpha diversity measures
## 49:131 = clr-transformed maternal faecal relative abundances of bacterial genera with ≥10% prevalence


##### =========================== 7. MODEL FOR ASSOCIATION OF SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK HMOS WITH CORRECTION FOR MILK GROUP AND SEQUENCING READS =========================== #####

### ===== 7.1 FUNCTION FOR ASSOCIATION OF SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK HMOS WITH CORRECTION FOR MILK GROUP AND SEQUENCING READS ===== ###

run.lmer.cor.milkgroup.reads.1tp <- function(datasetname, inputdata, xcolumn, ycolumn){
  
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
      m1 <- lm(my_df[,j] ~ mother_milk_HMO_milk_group + seq_16S_n_reads_clean + my_df[,i], data=my_df)
      sum1 <- summary(m1)
      
      #retrieve estimates and p values from model with phenotype of interest (from sum1)
      e <- c(e, c(sum1$coefficients[grep("my_df",rownames(sum1$coefficients)), "Estimate"])) #first rows are for mother_milk_HMO_milk_group and seq_16S_n_reads_clean, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      p <- c(p, c(sum1$coefficients[grep("my_df",rownames(sum1$coefficients)), "Pr(>|t|)"])) #first rows are for mother_milk_HMO_milk_group and seq_16S_n_reads_clean, the number of rows varies between phenotypes -> ensure that this picks from the correct row on to pull out the results for i
      stat <- c(stat, paste0("lm_1_time_point_correction_milkgroup_reads"))
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


### ===== 7.2 RUN ASSOCIATION OF MATERNAL SINGLE-TIME-POINT MATERNAL MICROBIOTA DATA WITH MILK HMOS ===== ###

## Maternal data set (matf_invr_RAclr), 170x131
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:17  = mother_milk_HMO_milk_group and seq_16S_n_reads_clean
## 18:45  = invr-transformed milk HMOs
## 46:48  = invr-transformed maternal faecal alpha diversity measures
## 49:131 = clr-transformed maternal faecal relative abundances of bacterial genera with ≥10% prevalence

matf_1tp_results <- run.lmer.cor.milkgroup.reads.1tp(datasetname = "mother_maternal_microbiota_milk_HMOs_1_time_point",
                                                     inputdata   = matf_invr_RAclr,
                                                     xcolumn     = c(46:48,49:131),   #invr-transformed maternal faecal alpha diversity and clr-transformed maternal faecal relative abundances of bacterial genera with ≥10% prevalence
                                                     ycolumn     = c(18:45))          #invr-transformed milk HMOs
# write.table(matf_1tp_results, file="association_model_results/240123_milk_HMOs_by_maternal_faecal_microbiota_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)


##### =========================== 8. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING =========================== #####

### ===== 8.1 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES ===== ###

matf_alpha_div_results <- matf_1tp_results[grep("alpha_div", matf_1tp_results$x),] #84x7
matf_RA_results <- matf_1tp_results[grep("g__", matf_1tp_results$x),] #2324x7


### ===== 8.2 CORRECT FOR MULTIPLE TESTING ===== ###

## FDR correction
matf_alpha_div_results$FDR <- p.adjust(matf_alpha_div_results$p, method="BH")
matf_RA_results$FDR <- p.adjust(matf_RA_results$p, method="BH")

## sort results by FDR and p
matf_alpha_div_results_ord <- matf_alpha_div_results[order(matf_alpha_div_results$FDR, matf_alpha_div_results$p),]
matf_RA_results_ord <- matf_RA_results[order(matf_RA_results$FDR, matf_RA_results$p),]

## -> No significant effects of maternal faecal alpha diversity or maternal faecal relative bacterial abundances (of bacterial genera with ≥10% prevalence) on huma milk HMOs.


### ===== 8.3 SAVE RESULTS (TABLE S25) ===== ###

## results of maternal faecal alpha diversity, including simpson, on milk HMOs
# write.table(matf_alpha_div_results_ord, file="association_model_results/240123_milk_HMOs_invr_real_by_maternal_faecal_alpha_diversity_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

## results of maternal faecal relative abundances on milk HMOs (TABLE S25)
write.table(matf_RA_results_ord, file="association_model_results/240123_milk_HMOs_invr_real_by_maternal_faecal_relative_abundances_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)


### ===== 8.4 EXCLUDE SIMPSON DIVERSITY, RECALCULATE FDR AND SAVE RESULTS AGAIN (TABLE S24) ===== ###

## reimport results
matf_alpha_div_results_ord <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/association_model_results/240123_milk_HMOs_invr_real_by_maternal_faecal_alpha_diversity_model_results.txt", header=T, sep="\t", stringsAsFactors=T)

## exclude results for Simpson diversity
matf_alpha_div_results_ord_no_simpson <- matf_alpha_div_results_ord[matf_alpha_div_results_ord$x!="alpha_div_simpson_invr",]

## recalculate FDR
matf_alpha_div_results_ord_no_simpson$FDR <- p.adjust(matf_alpha_div_results_ord_no_simpson$p, method="BH")

## resort results
matf_alpha_div_results_ord_no_simpson <- matf_alpha_div_results_ord_no_simpson[order(matf_alpha_div_results_ord_no_simpson$FDR, matf_alpha_div_results_ord_no_simpson$p),]

## save results (TABLE S24)
write.table(matf_alpha_div_results_ord_no_simpson, file="association_model_results/240222_milk_HMOs_invr_real_by_maternal_faecal_alpha_diversity_richness_and_shannon_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)







