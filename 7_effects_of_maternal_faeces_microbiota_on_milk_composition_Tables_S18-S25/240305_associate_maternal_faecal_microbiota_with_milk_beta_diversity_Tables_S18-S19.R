#######################################################################################################################
### TEST EFFECTS OF MATERNAL FAECAL MICROBIOTA ON MILK MICROBIOTA DATA (BETA DIVERSITY) FOR MILK COMPOSITION PAPER  ###
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
# 2. INVERSE-RANK TRANSFORM MATERNAL FAECAL ALPHA DIVERSITY DATA
# 2.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 2.2 INVERSE-RANK TRANSFORM MATERNAL FAECAL ALPHA DIVERSITY DATA
# 
# 3. CALCULATE RELATIVE ABUNDANCES
# 3.1 CALCULATE MATERNAL FAECAL RELATIVE ABUNDANCES AND SELECT ONLY GENERA WITH ≥10% PREVALENCE
# 3.2 CALCULATE MILK RELATIVE ABUNDANCES
# 
# 4. CLR-TRANSFORM RELATIVE ABUNDANCES
# 4.1 FUNCTION FOR CLR-TRANSFORMATION
# 4.2 CLR-TRANSFORM MATERNAL FAECAL RELATIVE ABUNDANCES
# 4.3 CLR-TRANSFORM MILK RELATIVE ABUNDANCES
# 
# 5. COMBINE MATERNAL FAECAL AND MILK MICROBIOTA DATA INTO 1 DATA FRAME
# 
# 6. ASSOCIATE MATERNAL FAECAL MICROBIOTA DATA WITH MILK BETA DIVERSITY
# 6.1 FUNCTION FOR ADONIS WITH CORRECTION FOR MILK DNA ISOLATION BATCH, MILK SEQUENCING READS AND FAECAL SEQUENCING READS
# 6.2 ASSOCIATE MATERNAL FAECAL MICROBIOTA DATA WITH MILK BETA DIVERSITY
# 
# 7. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING
# 7.1 EXCLUDE RESULTS FOR SIMPSON DIVERSITY
# 7.2 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES
# 7.3 CORRECT FOR MULTIPLE TESTING
# 7.4 SAVE RESULTS (TABLES S18 AND S19)


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


##### =========================== 2. INVERSE-RANK TRANSFORM MATERNAL FAECAL ALPHA DIVERSITY DATA =========================== #####

### ===== 2.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 2.2 INVERSE-RANK TRANSFORM MATERNAL FAECAL ALPHA DIVERSITY DATA ===== ###

matf_invr <- matf
matf_invr[,17:19] <- as.data.frame(apply(matf_invr[,17:19], 2, invrank))
colnames(matf_invr)[17:19] <- c(paste0(colnames(matf_invr)[17:19], "_invr"))


##### =========================== 3. CALCULATE RELATIVE ABUNDANCES =========================== #####

### ===== 3.1 CALCULATE MATERNAL FAECAL RELATIVE ABUNDANCES AND SELECT ONLY GENERA WITH ≥10% PREVALENCE ===== ###

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


### ===== 3.2 CALCULATE MILK RELATIVE ABUNDANCES ===== ###

## save only columns with absolute bacterial abundances
rownames(milk) <- milk$seq_16S_sample_ID
milkbact <- milk[,21:ncol(milk)] #164x624

## calculate relative bacterial abundances
milkbact_relab <- (milkbact/rowSums(milkbact))
table(rowSums(milkbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
milkbact_relab_tmp <- milkbact_relab[,-624]

## exclude columns for bacteria that were not identified in any sample
milkbact_relab_tmp2 <- milkbact_relab_tmp[,colSums(milkbact_relab_tmp)>0]
dim(milkbact_relab_tmp2) #164x300

## merge back with metadata
milk_RA <- merge(milk[,1:20], milkbact_relab_tmp2, by="row.names")
rownames(milk_RA) <- milk_RA$Row.names
milk_RA <- milk_RA[,-1]
#164x320


##### =========================== 4. CLR-TRANSFORM RELATIVE ABUNDANCES =========================== #####

### ===== 4.1 FUNCTION FOR CLR-TRANSFORMATION ===== ###

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


### ===== 4.2 CLR-TRANSFORM MATERNAL FAECAL RELATIVE ABUNDANCES ===== ###

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


### ===== 4.3 CLR-TRANSFORM MILK RELATIVE ABUNDANCES ===== ###

milk_RAclr <- milk_RA
milk_RAclr[,21:ncol(milk_RAclr)] <- as.data.frame(do_clr_externalWeighting(milk_RAclr[,21:ncol(milk_RAclr)], milk_RAclr[,21:ncol(milk_RAclr)]))

colnames(milk_RAclr)[21:ncol(milk_RAclr)] <- c(paste0(colnames(milk_RAclr)[21:ncol(milk_RAclr)], "_clr"))


##### =========================== 5. COMBINE MATERNAL FAECAL AND MILK MICROBIOTA DATA INTO 1 DATA FRAME =========================== #####

## ensure all columns have unique names, indicating from which data set they derived
colnames(matf_invr_RAclr)[c(2:102)] <- c(paste0("matf_", colnames(matf_invr_RAclr)[c(2:102)]))
colnames(milk_RAclr)[c(2:320)] <- c(paste0("milk_", colnames(milk_RAclr)[c(2:320)]))

matfmilk <- merge(matf_invr_RAclr, milk_RAclr, by="mother_sample_ID")


##### =========================== 6. ASSOCIATE MATERNAL FAECAL MICROBIOTA DATA WITH MILK BETA DIVERSITY =========================== #####

library(vegan)


### ===== 6.1 FUNCTION FOR ADONIS WITH CORRECTION FOR MILK DNA ISOLATION BATCH, MILK SEQUENCING READS AND FAECAL SEQUENCING READS ===== ###

## function to run Adonis with correction for DNA isolation batch and clean read counts
run.adonis2.cor.DNAbatch.mreads.freads <- function(mydata, columns){
  column <- c()
  R2 <- c()
  p <- c()
  n <- c()

  note <- c()
  
  for (i in columns){
    print(paste0("Starting analysis for ", colnames(mydata)[i]))
    
    d <- mydata[!is.na(mydata[,i]),]

    ## AITCHISON DISTANCE MATRIX + METADATA
    print("Creating Aitchison distance matrix")
    ait <- vegdist(d[,grep("milk_g__", colnames(d))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
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
    
    if (any(a>0.97)){ #set to NA if >97% of the data is the same (this includes if there is only 1 factor because that is 100% the same)
      R2 <- c(R2, NA)
      p <- c(p, NA)
      note <- c(note, paste0("Adonis_could_not_be_run_as_there_was_too_little_variation_in_the_data"))
    }else{ ##only run adonis if there is sufficient variation in the answers (at least 3% need to have different answer), this also means that >=2 levels are present for the variable of interest
      a1 <- adonis2(ait ~ new_df$milk_DNA_isolation_batch + new_df$milk_seq_16S_n_reads_clean + new_df$matf_seq_16S_n_reads_clean + new_df[,i], by="margin")
      print(a1)
      R2 <- c(R2, a1$R2[4]) #we pick the 4th value as the first 3 are for the factors we correct for
      p <- c(p, a1$'Pr(>F)'[4]) #we pick the 4th value as the first 3 are for the factors we correct for
      note <- c(note, NA)
    }
    }
  
  print("Preparing results data frame")
  b <- data.frame(column=column, n=n, R2=R2, pvalue=p, note=note)
  return(b)
}


### ===== 6.2 ASSOCIATE MATERNAL FAECAL MICROBIOTA DATA WITH MILK BETA DIVERSITY ===== ###

adonis2_results <- run.adonis2.cor.DNAbatch.mreads.freads(mydata = matfmilk,
                                                          columns = c(17:19,20:102)) #invr-transformed maternal faecal alpha diversity and clr-transformed relative abundances of maternal faecal bacteria with ≥10% prevalence

adonis2_results <- adonis2_results[,-5] #exclude note column as it was only NA


##### =========================== 7. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING =========================== #####

### ===== 7.1 EXCLUDE RESULTS FOR SIMPSON DIVERSITY ===== ###

adonis2_results_no_simpson <- adonis2_results[adonis2_results$column!="matf_alpha_div_simpson_invr",] #85x4


### ===== 7.2 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES ===== ###

adonis2_results_by_alpha_div <- adonis2_results_no_simpson[grep("matf_alpha_div", adonis2_results_no_simpson$column),] #2x4
adonis2_results_by_RA <- adonis2_results_no_simpson[grep("matf_g__", adonis2_results_no_simpson$column),] #82x4


### ===== 7.3 CORRECT FOR MULTIPLE TESTING ===== ###

## FDR correction
adonis2_results_by_alpha_div$FDR <- p.adjust(adonis2_results_by_alpha_div$p, method="BH")
adonis2_results_by_RA$FDR <- p.adjust(adonis2_results_by_RA$p, method="BH")

## sort results by FDR and p
adonis2_results_by_alpha_div_ord <- adonis2_results_by_alpha_div[order(adonis2_results_by_alpha_div$FDR, adonis2_results_by_alpha_div$p),]
adonis2_results_by_RA_ord <- adonis2_results_by_RA[order(adonis2_results_by_RA$FDR, adonis2_results_by_RA$p),]

## -> No significant effects of maternal faecal alpha diversity or maternal faecal relative bacterial abundances on human milk beta diversity.


### ===== 7.4 SAVE RESULTS (TABLES S18 AND S19) ===== ###

write.table(adonis2_results_by_alpha_div_ord, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/association_model_results/240305_milk_beta_diversity_by_maternal_faecal_alpha_diversity_adonis_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis2_results_by_RA_ord, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/association_model_results/240305_milk_beta_diversity_by_maternal_faecal_relative_abundances_adonis_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

