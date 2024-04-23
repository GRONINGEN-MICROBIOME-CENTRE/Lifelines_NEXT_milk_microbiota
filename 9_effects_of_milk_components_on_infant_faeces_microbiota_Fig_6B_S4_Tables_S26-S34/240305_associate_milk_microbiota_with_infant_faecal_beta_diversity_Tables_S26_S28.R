###########################################################################################################################################
### TEST EFFECTS OF MILK MICROBIOTA (ALPHA DIV AND REL ABUNDANCES) ON INFANT FAECAL MICROBIOTA DATA (BETA DIVERSITY)  (TABLES S26, S28) ###
###########################################################################################################################################

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
# 1.1 MILK MICROBIOTA DATA
# 1.2 INFANT FAECAL MICROBIOTA DATA
# 1.3 SELECT ONLY PAIRS WITH BOTH MILK AND INFANT FAECAL MICROBIOTA DATA (PER TIME POINT)
# 
# 2. INVERSE-RANK TRANSFORM ALPHA DIVERSITY MEASURES
# 2.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 2.2 INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY DATA
# 
# 3. CALCULATE RELATIVE ABUNDANCES
# 3.1 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA
# 3.2 CALCULATE INFANT FAECAL RELATIVE ABUNDANCES AND SELECT MOST PRVEALENT BACTERIA
# 
# 4. CLR-TRANSFORM RELATIVE ABUNDANCES
# 4.1 FUNCTION FOR CLR-TRANSFORMATION
# 4.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES
# 4.3 CLR-TRANSFORM INFANT RELATIVE ABUNDANCES
# 
# 5. COMBINE MILK AND INFANT FAECAL MICROBIOTA DATA INTO 1 DATA FRAME
# 
# 6. ASSOCIATE MILK MICROBIOTA DATA WITH INFANT FAECAL BETA DIVERSITY
# 6.1 FUNCTION FOR ADONIS WITH CORRECTION FOR MILK DNA ISOLATION BATCH, MILK SEQUENCING READS AND FAECAL SEQUENCING READS
# 6.2 ASSOCIATE MILK MICROBIOTA DATA WITH INFANT FAECAL BETA DIVERSITY
# 
# 7. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING
# 7.1 EXCLUDE RESULTS FOR SIMPSON DIVERSITY
# 7.2 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES
# 7.3 CORRECT FOR MULTIPLE TESTING
# 7.4 CHECK RESULTS
# 7.5 SAVE RESULTS (TABLES S26, S28)


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

### ===== 1.1 MILK MICROBIOTA DATA ===== ###

## select human milk samples keep duplicated maternal data
milk_tmp <- data[data$sample_origin_type=="mother_human_milk",]
milk_tmp2 <- milk_tmp[, c(1:9,11:15,504,511,526,528:1154)] #737x644

## Maternal milk data set (milk_tmp2)
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:17  = DNA_isolation_batch and seq_16S_n_reads_clean (for milk)
## 18:20  = maternal faecal alpha diversity measures
## 21:644 = maternal faecal bacterial relative abundances


### ===== 1.2 INFANT FAECAL MICROBIOTA DATA ===== ###

## select infant faecal samples
inff_tmp <- data[data$sample_origin_type=="infant_faeces",]
inff_tmp2 <- inff_tmp[, c(1:9,11:15,504,510,526,433,528:1154)] #552x645

## Maternal faecal data set (inff_tmp2)
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:18  = DNA_isolation_method, seq_16S_n_reads_clean (infant faeces reads), infant_birth_delivery_mode
## 19:21  = infant faecal alpha diversity measures
## 22:645 = infant faecal bacterial relative abundances


### ===== 1.3 SELECT ONLY PAIRS WITH BOTH MILK AND INFANT FAECAL MICROBIOTA DATA (PER TIME POINT) ===== ###

pairIDs <- intersect(milk_tmp2$mother_sample_ID, inff_tmp2$mother_sample_ID) #509 pairs (these IDs include the time points :))

milk <- milk_tmp2[milk_tmp2$mother_sample_ID %in% pairIDs,] #509x644
inff <- inff_tmp2[inff_tmp2$mother_sample_ID %in% pairIDs,] #509x645


##### =========================== 2. INVERSE-RANK TRANSFORM ALPHA DIVERSITY MEASURES =========================== #####

### ===== 2.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 2.2 INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY DATA ===== ###

milk_invr <- milk
milk_invr[,grep("alpha_div",colnames(milk_invr))] <- as.data.frame(apply(milk_invr[,grep("alpha_div",colnames(milk_invr))], 2, invrank))

## fix colnames of all invr-transformed columns
colnames(milk_invr)[18:20] <- paste0(colnames(milk_invr)[18:20], "_invr")


##### =========================== 3. CALCULATE RELATIVE ABUNDANCES =========================== #####

### ===== 3.1 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA ===== ###

## save only columns with absolute bacterial abundances
rownames(milk_invr) <- milk_invr$seq_16S_sample_ID
milkbact <- milk_invr[,21:ncol(milk_invr)] #509x624

## calculate relative bacterial abundances
milkbact_relab <- (milkbact/rowSums(milkbact))
table(rowSums(milkbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
milkbact_relab_tmp <- milkbact_relab[,-624]

## select only bacterial genera with ≥10% prevalence as for other analyses
milk_RAprev10_sum <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/relative_abundances/231215_milk_RA_prev10_summary_statistics.txt", header=T, sep="\t")
milk_RAprev10_bact <- unique(milk_RAprev10_sum$bacterium)
milk_RAprev10_bact <- c(paste0("g__", milk_RAprev10_bact)) #36

milkbact_relab_tmp2 <- milkbact_relab_tmp[,colnames(milkbact_relab_tmp) %in% milk_RAprev10_bact]
length(milkbact_relab_tmp2) #36
# colnames(milkbact_relab_tmp2)

## merge back with metadata
milk_invr_RA <- merge(milk_invr[,1:20], milkbact_relab_tmp2, by="row.names")
rownames(milk_invr_RA) <- milk_invr_RA$Row.names
milk_invr_RA <- milk_invr_RA[,-1]
#509x56


### ===== 3.2 CALCULATE INFANT FAECAL RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA ===== ###

## save only columns with absolute bacterial abundances
rownames(inff) <- inff$seq_16S_sample_ID
inffbact <- inff[,22:ncol(inff)] #509x624

## calculate relative bacterial abundances
inffbact_relab <- (inffbact/rowSums(inffbact))
table(rowSums(inffbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
inffbact_relab_tmp <- inffbact_relab[,-624]

## exclude columns for bacteria that were not identified in any sample
inffbact_relab_tmp2 <- inffbact_relab_tmp[,colSums(inffbact_relab_tmp)>0]
dim(inffbact_relab_tmp2) #509x212

## merge back with metadata
inff_RA <- merge(inff[,1:21], inffbact_relab_tmp2, by="row.names")
rownames(inff_RA) <- inff_RA$Row.names
inff_RA <- inff_RA[,-1]
#509x233


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


### ===== 4.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES ===== ###

milk_invr_RAclr <- milk_invr_RA
milk_invr_RAclr[,21:ncol(milk_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(milk_invr_RAclr[,21:ncol(milk_invr_RAclr)], milk_invr_RAclr[,21:ncol(milk_invr_RAclr)]))

colnames(milk_invr_RAclr)[21:ncol(milk_invr_RAclr)] <- c(paste0(colnames(milk_invr_RAclr)[21:ncol(milk_invr_RAclr)], "_clr"))


### ===== 4.3 CLR-TRANSFORM INFANT RELATIVE ABUNDANCES ===== ###

inff_RAclr <- inff_RA
inff_RAclr[,22:ncol(inff_RAclr)] <- as.data.frame(do_clr_externalWeighting(inff_RAclr[,22:ncol(inff_RAclr)], inff_RAclr[,22:ncol(inff_RAclr)]))

colnames(inff_RAclr)[22:ncol(inff_RAclr)] <- c(paste0(colnames(inff_RAclr)[22:ncol(inff_RAclr)], "_clr"))


##### =========================== 5. COMBINE MILK AND INFANT FAECAL MICROBIOTA DATA INTO 1 DATA FRAME =========================== #####

## ensure all columns have unique names, indicating from which data set they derived
colnames(milk_invr_RAclr)[c(2:56)] <- c(paste0("milk_", colnames(milk_invr_RAclr)[c(2:56)]))
colnames(inff_RAclr)[c(2:233)] <- c(paste0("inff_", colnames(inff_RAclr)[c(2:233)]))

milkinff <- merge(milk_invr_RAclr, inff_RAclr, by="mother_sample_ID")


##### =========================== 6. ASSOCIATE MILK MICROBIOTA DATA WITH INFANT FAECAL BETA DIVERSITY =========================== #####

library(vegan)


### ===== 6.1 FUNCTION FOR ADONIS WITH CORRECTION FOR MILK DNA ISOLATION BATCH, MILK SEQUENCING READS AND FAECAL SEQUENCING READS ===== ###

## function to run Adonis with correction for DNA isolation batch and clean read counts
run.adonis2.cor.DNAbatch.mreads.DNAmethod.freads.delmode <- function(mydata, columns){
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
    ait <- vegdist(d[,grep("inff_g__", colnames(d))], method="euclidian") #only include columns with relative bacterial abundances on the genus level
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
      a1 <- adonis2(ait ~ new_df$milk_DNA_isolation_batch + new_df$milk_seq_16S_n_reads_clean + new_df$inff_DNA_isolation_method + new_df$inff_seq_16S_n_reads_clean + new_df$inff_infant_birth_delivery_mode + new_df[,i], by="margin")
      print(a1)
      R2 <- c(R2, a1$R2[6]) #we pick the 6th value as the first 5 are for the factors we correct for
      p <- c(p, a1$'Pr(>F)'[6]) #we pick the 6th value as the first 5 are for the factors we correct for
      note <- c(note, NA)
    }
    }
  
  print("Preparing results data frame")
  b <- data.frame(column=column, n=n, R2=R2, pvalue=p, note=note)
  return(b)
}


### ===== 6.2 ASSOCIATE MILK MICROBIOTA DATA WITH INFANT FAECAL BETA DIVERSITY ===== ###

## M1
adonis2_results_M1 <- run.adonis2.cor.DNAbatch.mreads.DNAmethod.freads.delmode(mydata = milkinff[milkinff$milk_time_point=="1_month" & !is.na(milkinff$inff_infant_birth_delivery_mode),],
                                                                            columns = c(18:20,21:56)) #invr-transformed milk alpha diversity and clr-transformed relative abundances of milk bacteria with ≥10% prevalence

adonis2_results_M1 <- adonis2_results_M1[,-5] #exclude note column as it was only NA


## M2
adonis2_results_M2 <- run.adonis2.cor.DNAbatch.mreads.DNAmethod.freads.delmode(mydata = milkinff[milkinff$milk_time_point=="2_months" & !is.na(milkinff$inff_infant_birth_delivery_mode),],
                                                                            columns = c(18:20,21:56)) #invr-transformed milk alpha diversity and clr-transformed relative abundances of milk bacteria with ≥10% prevalence

adonis2_results_M2 <- adonis2_results_M2[,-5] #exclude note column as it was only NA


## M3
adonis2_results_M3 <- run.adonis2.cor.DNAbatch.mreads.DNAmethod.freads.delmode(mydata = milkinff[milkinff$milk_time_point=="3_months" & !is.na(milkinff$inff_infant_birth_delivery_mode),],
                                                                            columns = c(18:20,21:56)) #invr-transformed milk alpha diversity and clr-transformed relative abundances of milk bacteria with ≥10% prevalence

adonis2_results_M3 <- adonis2_results_M3[,-5] #exclude note column as it was only NA


## M6
adonis2_results_M6 <- run.adonis2.cor.DNAbatch.mreads.DNAmethod.freads.delmode(mydata = milkinff[milkinff$milk_time_point=="6_months" & !is.na(milkinff$inff_infant_birth_delivery_mode),],
                                                                            columns = c(18:20,21:56)) #invr-transformed milk alpha diversity and clr-transformed relative abundances of milk bacteria with ≥10% prevalence

adonis2_results_M6 <- adonis2_results_M6[,-5] #exclude note column as it was only NA


##### =========================== 7. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING =========================== #####

### ===== 7.1 EXCLUDE RESULTS FOR SIMPSON DIVERSITY ===== ###

adonis2_results_M1_no_simpson <- adonis2_results_M1[adonis2_results_M1$column!="milk_alpha_div_simpson_invr",] #38x4
adonis2_results_M2_no_simpson <- adonis2_results_M2[adonis2_results_M2$column!="milk_alpha_div_simpson_invr",] #38x4
adonis2_results_M3_no_simpson <- adonis2_results_M3[adonis2_results_M3$column!="milk_alpha_div_simpson_invr",] #38x4
adonis2_results_M6_no_simpson <- adonis2_results_M6[adonis2_results_M6$column!="milk_alpha_div_simpson_invr",] #38x4


### ===== 7.2 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES ===== ###

adonis2_results_M1_alpha_div <- adonis2_results_M1_no_simpson[grep("milk_alpha_div", adonis2_results_M1_no_simpson$column),] #2x4
adonis2_results_M1_RA <- adonis2_results_M1_no_simpson[grep("milk_g__", adonis2_results_M1_no_simpson$column),] #36x4

adonis2_results_M2_alpha_div <- adonis2_results_M2_no_simpson[grep("milk_alpha_div", adonis2_results_M2_no_simpson$column),] #2x4
adonis2_results_M2_RA <- adonis2_results_M2_no_simpson[grep("milk_g__", adonis2_results_M2_no_simpson$column),] #36x4

adonis2_results_M3_alpha_div <- adonis2_results_M3_no_simpson[grep("milk_alpha_div", adonis2_results_M3_no_simpson$column),] #2x4
adonis2_results_M3_RA <- adonis2_results_M3_no_simpson[grep("milk_g__", adonis2_results_M3_no_simpson$column),] #36x4

adonis2_results_M6_alpha_div <- adonis2_results_M6_no_simpson[grep("milk_alpha_div", adonis2_results_M6_no_simpson$column),] #2x4
adonis2_results_M6_RA <- adonis2_results_M6_no_simpson[grep("milk_g__", adonis2_results_M6_no_simpson$column),] #36x4


### ===== 7.3 CORRECT FOR MULTIPLE TESTING ===== ###

## FDR correction
adonis2_results_M1_alpha_div$FDR <- p.adjust(adonis2_results_M1_alpha_div$p, method="BH")
adonis2_results_M1_RA$FDR <- p.adjust(adonis2_results_M1_RA$p, method="BH")

adonis2_results_M2_alpha_div$FDR <- p.adjust(adonis2_results_M2_alpha_div$p, method="BH")
adonis2_results_M2_RA$FDR <- p.adjust(adonis2_results_M2_RA$p, method="BH")

adonis2_results_M3_alpha_div$FDR <- p.adjust(adonis2_results_M3_alpha_div$p, method="BH")
adonis2_results_M3_RA$FDR <- p.adjust(adonis2_results_M3_RA$p, method="BH")

adonis2_results_M6_alpha_div$FDR <- p.adjust(adonis2_results_M6_alpha_div$p, method="BH")
adonis2_results_M6_RA$FDR <- p.adjust(adonis2_results_M6_RA$p, method="BH")


## sort results by FDR and p
adonis2_results_M1_alpha_div_ord <- adonis2_results_M1_alpha_div[order(adonis2_results_M1_alpha_div$FDR, adonis2_results_M1_alpha_div$p),]
adonis2_results_M1_RA_ord <- adonis2_results_M1_RA[order(adonis2_results_M1_RA$FDR, adonis2_results_M1_RA$p),]

adonis2_results_M2_alpha_div_ord <- adonis2_results_M2_alpha_div[order(adonis2_results_M2_alpha_div$FDR, adonis2_results_M2_alpha_div$p),]
adonis2_results_M2_RA_ord <- adonis2_results_M2_RA[order(adonis2_results_M2_RA$FDR, adonis2_results_M2_RA$p),]

adonis2_results_M3_alpha_div_ord <- adonis2_results_M3_alpha_div[order(adonis2_results_M3_alpha_div$FDR, adonis2_results_M3_alpha_div$p),]
adonis2_results_M3_RA_ord <- adonis2_results_M3_RA[order(adonis2_results_M3_RA$FDR, adonis2_results_M3_RA$p),]

adonis2_results_M6_alpha_div_ord <- adonis2_results_M6_alpha_div[order(adonis2_results_M6_alpha_div$FDR, adonis2_results_M6_alpha_div$p),]
adonis2_results_M6_RA_ord <- adonis2_results_M6_RA[order(adonis2_results_M6_RA$FDR, adonis2_results_M6_RA$p),]


### ===== 7.4 CHECK RESULTS ===== ###

## infant faecal beta diversity results by milk alpha diversity
nrow(adonis2_results_M1_alpha_div_ord[adonis2_results_M1_alpha_div_ord$FDR<0.05,]) #0 FDR-significant results
nrow(adonis2_results_M1_alpha_div_ord[adonis2_results_M1_alpha_div_ord$p<0.05,])   #1 nominally significant results

nrow(adonis2_results_M2_alpha_div_ord[adonis2_results_M2_alpha_div_ord$FDR<0.05,]) #0 FDR-significant results
nrow(adonis2_results_M2_alpha_div_ord[adonis2_results_M2_alpha_div_ord$p<0.05,])   #1 nominally significant results

nrow(adonis2_results_M3_alpha_div_ord[adonis2_results_M3_alpha_div_ord$FDR<0.05,]) #1 FDR-significant results
nrow(adonis2_results_M3_alpha_div_ord[adonis2_results_M3_alpha_div_ord$p<0.05,])   #1 nominally significant results

nrow(adonis2_results_M6_alpha_div_ord[adonis2_results_M6_alpha_div_ord$FDR<0.05,]) #0 FDR-significant results
nrow(adonis2_results_M6_alpha_div_ord[adonis2_results_M6_alpha_div_ord$p<0.05,])   #0 nominally significant results


## infant faecal beta diversity results by milk relative abundances
nrow(adonis2_results_M1_RA_ord[adonis2_results_M1_RA_ord$FDR<0.05,]) #0 FDR-significant results
nrow(adonis2_results_M1_RA_ord[adonis2_results_M1_RA_ord$p<0.05,])   #0 nominally significant results

nrow(adonis2_results_M2_RA_ord[adonis2_results_M2_RA_ord$FDR<0.05,]) #0 FDR-significant results
nrow(adonis2_results_M2_RA_ord[adonis2_results_M2_RA_ord$p<0.05,])   #0 nominally significant results

nrow(adonis2_results_M3_RA_ord[adonis2_results_M3_RA_ord$FDR<0.05,]) #0 FDR-significant results
nrow(adonis2_results_M3_RA_ord[adonis2_results_M3_RA_ord$p<0.05,])   #1 nominally significant results

nrow(adonis2_results_M6_RA_ord[adonis2_results_M6_RA_ord$FDR<0.05,]) #0 FDR-significant results
nrow(adonis2_results_M6_RA_ord[adonis2_results_M6_RA_ord$p<0.05,])   #0 nominally significant results


### ===== 7.5 SAVE RESULTS (TABLES S26, S28) ===== ###

setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/effects_of_milk_microbiota/")

## TABLES S26A-D
write.table(adonis2_results_M1_alpha_div_ord, file="240305_infant_faecal_beta_diversity_by_milk_alpha_diversity_adonis_results_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis2_results_M2_alpha_div_ord, file="240305_infant_faecal_beta_diversity_by_milk_alpha_diversity_adonis_results_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis2_results_M3_alpha_div_ord, file="240305_infant_faecal_beta_diversity_by_milk_alpha_diversity_adonis_results_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis2_results_M6_alpha_div_ord, file="240305_infant_faecal_beta_diversity_by_milk_alpha_diversity_adonis_results_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)

## TABLES S28A-D
write.table(adonis2_results_M1_RA_ord, file="240305_infant_faecal_beta_diversity_by_milk_relative_abundances_adonis_results_M1.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis2_results_M2_RA_ord, file="240305_infant_faecal_beta_diversity_by_milk_relative_abundances_adonis_results_M2.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis2_results_M3_RA_ord, file="240305_infant_faecal_beta_diversity_by_milk_relative_abundances_adonis_results_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)
write.table(adonis2_results_M6_RA_ord, file="240305_infant_faecal_beta_diversity_by_milk_relative_abundances_adonis_results_M6.txt", row.names=F, col.names=T, sep="\t", quote=F)

