#########################################################################################################################################################
### ANALYSE EFFECT OF MILK MICROBIOTA (ALPHA DIV AND REL ABUNDANCES) ON INFANT FAECAL MICROBIOTA (ALPHA DIV AND REL ABUNDANCES) (TABLES S27, S29-S31) ###
#########################################################################################################################################################

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
# 2. ENSURE CORRECT DATA STRUCTURE
# 
# 3. INVERSE-RANK TRANSFORM ALPHA DIVERSITY MEASURES
# 3.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 3.2 INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY DATA
# 3.3 INVERSE-RANK TRANSFORM INFANT FAECAL ALPHA DIVERSITY DATA
# 
# 4. CALCULATE RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE
# 4.1 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA
# 4.2 CALCULATE INFANT FAECAL RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA
# 
# 5. CLR-TRANSFORM RELATIVE ABUNDANCES
# 5.1 FUNCTION FOR CLR-TRANSFORMATION
# 5.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES
# 5.3 CLR-TRANSFORM INFANT FAECAL RELATIVE ABUNDANCES
# 
# 6. COMBINE MILK AND INFANT FAECAL MICROBIOTA DATA INTO 1 DATA FRAME
# 
# 7. MODEL FUNCTION (LMER/LM)
# 7.1 RUN MODELS FOR INFANT FAECAL ALPHA DIVERSITY
# 
# 8. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING
# 8.1 EXCLUDE RESULTS FOR SIMPSON DIVERSITY
# 8.2 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES
# 8.3 CORRECT FOR MULTIPLE TESTING
# 8.4 CHECK RESULTS
# 8.5 SAVE RESULTS (TABLES S27, S29, S30, S31)


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


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

## check that data structure is correct
str(milk[,c(1:21,643:644)])
# -> data structure is correct

## check that data structure is correct
str(inff[,c(1:22,644:645)])
# -> data structure is correct


##### =========================== 3. INVERSE-RANK TRANSFORM ALPHA DIVERSITY MEASURES =========================== #####

### ===== 3.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 3.2 INVERSE-RANK TRANSFORM MILK ALPHA DIVERSITY DATA ===== ###

milk_invr <- milk
milk_invr[,grep("alpha_div",colnames(milk_invr))] <- as.data.frame(apply(milk_invr[,grep("alpha_div",colnames(milk_invr))], 2, invrank))

## fix colnames of all invr-transformed columns
colnames(milk_invr)[18:20] <- paste0(colnames(milk_invr)[18:20], "_invr")


### ===== 3.3 INVERSE-RANK TRANSFORM INFANT FAECAL ALPHA DIVERSITY DATA ===== ###

inff_invr <- inff
inff_invr[,grep("alpha_div",colnames(inff_invr))] <- as.data.frame(apply(inff_invr[,grep("alpha_div",colnames(inff_invr))], 2, invrank))

## fix colnames of all invr-transformed columns
colnames(inff_invr)[19:21] <- paste0(colnames(inff_invr)[19:21], "_invr")


##### =========================== 4. CALCULATE RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE =========================== #####

### ===== 4.1 CALCULATE MILK RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA ===== ###

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


### ===== 4.2 CALCULATE INFANT FAECAL RELATIVE ABUNDANCES AND SELECT MOST PREVALENT BACTERIA ===== ###

## save only columns with absolute bacterial abundances
rownames(inff_invr) <- inff_invr$seq_16S_sample_ID
inffbact <- inff_invr[,22:ncol(inff_invr)] #509x624

## calculate relative bacterial abundances
inffbact_relab <- (inffbact/rowSums(inffbact))
table(rowSums(inffbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
inffbact_relab_tmp <- inffbact_relab[,-624]

## select only bacterial genera with ≥10% prevalence as for other analyses
inff_RAprev10_sum <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/231215_infant_faeces_RA_prev10_summary_statistics.txt", header=T, sep="\t")
inff_RAprev10_bact <- unique(inff_RAprev10_sum$bacterium)
inff_RAprev10_bact <- c(paste0("g__", inff_RAprev10_bact)) #30

inffbact_relab_tmp2 <- inffbact_relab_tmp[,colnames(inffbact_relab_tmp) %in% inff_RAprev10_bact]
length(inffbact_relab_tmp2) #30
# colnames(inffbact_relab_tmp2)

## merge back with metadata
inff_invr_RA <- merge(inff_invr[,1:21], inffbact_relab_tmp2, by="row.names")
rownames(inff_invr_RA) <- inff_invr_RA$Row.names
inff_invr_RA <- inff_invr_RA[,-1]
#509x51


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


### ===== 5.2 CLR-TRANSFORM MILK RELATIVE ABUNDANCES ===== ###

milk_invr_RAclr <- milk_invr_RA
milk_invr_RAclr[,21:ncol(milk_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(milk_invr_RAclr[,21:ncol(milk_invr_RAclr)], milk_invr_RAclr[,21:ncol(milk_invr_RAclr)]))

## fix colnames of clr-transformed columns
colnames(milk_invr_RAclr)[21:ncol(milk_invr_RAclr)] <- c(paste0(colnames(milk_invr_RAclr)[21:ncol(milk_invr_RAclr)], "_clr"))

## fix rownames
rownames(milk_invr_RAclr) <- milk_invr_RAclr$mother_sample_ID

## Maternal data set (milk_invr_RAclr), 509x56
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:17  = milk DNA_isolation_batch and milk seq_16S_n_reads_clean
## 18:20  = invr-transformed milk alpha diversity measures
## 21:56 = clr-transformed milk relative abundances of bacterial genera with ≥10% prevalence


### ===== 5.3 CLR-TRANSFORM INFANT FAECAL RELATIVE ABUNDANCES ===== ###

inff_invr_RAclr <- inff_invr_RA
inff_invr_RAclr[,22:ncol(inff_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(inff_invr_RAclr[,22:ncol(inff_invr_RAclr)], inff_invr_RAclr[,22:ncol(inff_invr_RAclr)]))

## fix colnames of clr-transformed columns
colnames(inff_invr_RAclr)[22:ncol(inff_invr_RAclr)] <- c(paste0(colnames(inff_invr_RAclr)[22:ncol(inff_invr_RAclr)], "_clr"))

## fix rownames
rownames(inff_invr_RAclr) <- inff_invr_RAclr$mother_sample_ID

## Maternal data set (inff_invr_RAclr), 164x102
## Columns:
## 1:15 = Information columns (mother_ID etc., incl. NEXT_participation_number and time_point)
## 16:18  = DNA_isolation_method, infant faecal seq_16S_n_reads_clean, infant_birth_delivery_mode
## 19:21  = invr-transformed infant faecal alpha diversity measures
## 22:51  = clr-transformed infant faecal relative abundances of bacterial genera with ≥10% prevalence


##### =========================== 6. COMBINE MILK AND INFANT FAECAL MICROBIOTA DATA INTO 1 DATA FRAME =========================== #####

## ensure all columns have unique names, indicating from which data set they derived
colnames(milk_invr_RAclr)[c(2:56)] <- c(paste0("milk_", colnames(milk_invr_RAclr)[c(2:56)]))
colnames(inff_invr_RAclr)[c(2:51)] <- c(paste0("inff_", colnames(inff_invr_RAclr)[c(2:51)]))

milkinff <- merge(milk_invr_RAclr, inff_invr_RAclr, by="mother_sample_ID")


##### =========================== 7. MODEL FUNCTION (LMER/LM) =========================== #####

## function to run lmer or lm models with/without correction for covariates
run.models <- function(datasetname, modeltype, covardetails, covariates, inputdata, xcolumn, ycolumn, outputfolder){
  
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
        
        my.lm.model <- paste0("my_df[,j] ~ ", covariates, " + my_df[,i]")
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
      jcolumntmp <- gsub("\\[", "", colnames(my_df)[j])
      jcolumn <- gsub("\\]", "", jcolumntmp)
      
      if (class(my_df[,i]) == "numeric" | class(my_df[,i]) == "integer"){
        print("Generating scatterplot")
        scatterplot1 <- ggplot(my_df, aes(x=my_df[,i], y=my_df[,j]))+
          geom_point(alpha=0.5)+
          theme_bw()+
          facet_grid(.~milk_time_point_numeric)+
          geom_smooth(method='lm')+
          labs(x=paste0(colnames(my_df)[i]), y=paste0(colnames(my_df)[j]),
               title=paste0(colnames(my_df)[j], " by ", colnames(my_df)[i], "\n",
                            "n=", nrow(my_df), "\n", 
                            covardetails, "\n",
                            "p=", signif(p[length(p)],3)))
        
        if (p[length(p)]<0.05){
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/nominally_significant/240305_infant_faeces_", jcolumn ,"_by_", colnames(my_df)[i], "_and_time.pdf"), scatterplot1, dpi=500, width=30, height=15, units="cm")
        }else{
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/not_significant/240305_infant_faeces_", jcolumn ,"_by_", colnames(my_df)[i], "_and_time.pdf"), scatterplot1, dpi=500, width=30, height=15, units="cm")
        }
        
      }else{
        print("Generating boxplots")
        boxplot1 <- ggplot(my_df, aes(x=my_df[,i], y=my_df[,j]))+
          geom_jitter(alpha=0.5)+
          geom_boxplot(alpha=0)+
          theme_bw()+
          facet_grid(.~milk_time_point_numeric)+
          labs(x=paste0(colnames(my_df)[i]), y=paste0(colnames(my_df)[j]),
               title=paste0(colnames(my_df)[j], " by ", colnames(my_df)[i], "\n",
                            "n=", nrow(my_df), "\n", 
                            covardetails, "\n",
                            "p=", signif(p[length(p)],3)))
        
        if (p[length(p)]<0.05){
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/nominally_significant/240305_infant_faeces_", jcolumn ,"_by_", colnames(my_df)[i], "_and_time.pdf"), boxplot1, dpi=500, width=30, height=15, units="cm")
        }else{
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/not_significant/240305_infant_faeces_", jcolumn ,"_by_", colnames(my_df)[i], "_and_time.pdf"), boxplot1, dpi=500, width=30, height=15, units="cm")
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
  write.table(res, file=paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/240305_", datasetname, "_", outputfolder, "_", covardetails, "_results.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  
  #return output table
  return(res)
}


### ===== 7.1 RUN MODELS FOR INFANT FAECAL ALPHA DIVERSITY ===== ###

inff_results <- run.models(datasetname = "infant_faeces_microbiota_by_maternal_microbiota",
                           modeltype = "lmer",
                           covardetails = "lmer_correction_for_DNAbatch_mreads_DNAmethod_freads_delmode_time_ID",
                           covariates = "milk_DNA_isolation_batch + milk_seq_16S_n_reads_clean + inff_DNA_isolation_method + inff_seq_16S_n_reads_clean + inff_infant_birth_delivery_mode + milk_time_point_numeric + (1|milk_infant_ID) + ",
                           inputdata   = milkinff[!is.na(milkinff$inff_infant_birth_delivery_mode),],
                           xcolumn     = c(18:20,21:56),  #invr-transformed milk alpha diversity measures and clr-transformed relative abundances of milk bacteria with ≥10% prevalence
                           ycolumn     = c(74:76,77:106),  #invr-transformed infant faecal alpha diversity measures and clr-transformed relative abundances of infant faecal bacteria with ≥10% prevalence
                           outputfolder = "effects_of_milk_microbiota")


##### =========================== 8. SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES AND CORRECT FOR MULTIPLE TESTING =========================== #####

### ===== 8.1 EXCLUDE RESULTS FOR SIMPSON DIVERSITY ===== ###

inff_results_no_simpson <- inff_results[inff_results$x!="milk_alpha_div_simpson_invr" & inff_results$y!="inff_alpha_div_simpson_invr",] #1216x9


### ===== 8.2 SEPARATE RESULTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES ===== ###

## infant faecal alpha diversity as outcome
inff_alpha_results <- inff_results_no_simpson[grep("inff_alpha_div", inff_results_no_simpson$y),] #76x9
inff_alpha_div_by_milk_alpha_div_results <- inff_alpha_results[grep("milk_alpha_div", inff_alpha_results$x),] #4x9
inff_alpha_div_by_milk_RA_results <- inff_alpha_results[grep("milk_g__", inff_alpha_results$x),] #72x9

## infant faecal relative abundances as outcome
inff_RA_results <- inff_results_no_simpson[grep("inff_g__", inff_results_no_simpson$y),] #1140x9
inff_RA_div_by_milk_alpha_div_results <- inff_RA_results[grep("milk_alpha_div", inff_RA_results$x),] #60x9
inff_RA_div_by_milk_RA_results <- inff_RA_results[grep("milk_g__", inff_RA_results$x),] #1080x9


### ===== 8.3 CORRECT FOR MULTIPLE TESTING ===== ###

## FDR correction
inff_alpha_div_by_milk_alpha_div_results$FDR <- p.adjust(inff_alpha_div_by_milk_alpha_div_results$p, method="BH")
inff_alpha_div_by_milk_RA_results$FDR <- p.adjust(inff_alpha_div_by_milk_RA_results$p, method="BH")

inff_RA_div_by_milk_alpha_div_results$FDR <- p.adjust(inff_RA_div_by_milk_alpha_div_results$p, method="BH")
inff_RA_div_by_milk_RA_results$FDR <- p.adjust(inff_RA_div_by_milk_RA_results$p, method="BH")

## sort results by FDR and p
inff_alpha_div_by_milk_alpha_div_results_ord <- inff_alpha_div_by_milk_alpha_div_results[order(inff_alpha_div_by_milk_alpha_div_results$FDR, inff_alpha_div_by_milk_alpha_div_results$p),]
inff_alpha_div_by_milk_RA_results_ord <- inff_alpha_div_by_milk_RA_results[order(inff_alpha_div_by_milk_RA_results$FDR, inff_alpha_div_by_milk_RA_results$p),]

inff_RA_div_by_milk_alpha_div_results_ord <- inff_RA_div_by_milk_alpha_div_results[order(inff_RA_div_by_milk_alpha_div_results$FDR, inff_RA_div_by_milk_alpha_div_results$p),]
inff_RA_div_by_milk_RA_results_ord <- inff_RA_div_by_milk_RA_results[order(inff_RA_div_by_milk_RA_results$FDR, inff_RA_div_by_milk_RA_results$p),]


### ===== 8.4 CHECK RESULTS ===== ###

## effect of milk alpha diversity on infant faecal alpha diversity
nrow(inff_alpha_div_by_milk_alpha_div_results_ord[inff_alpha_div_by_milk_alpha_div_results_ord$FDR<0.05,]) #0 FDR significant results
nrow(inff_alpha_div_by_milk_alpha_div_results_ord[inff_alpha_div_by_milk_alpha_div_results_ord$p<0.05,])   #1 nominally significant result

## effect of milk relative abundances on infant faecal alpha diversity
nrow(inff_alpha_div_by_milk_RA_results_ord[inff_alpha_div_by_milk_RA_results_ord$FDR<0.05,]) #0 FDR significant results
nrow(inff_alpha_div_by_milk_RA_results_ord[inff_alpha_div_by_milk_RA_results_ord$p<0.05,])   #3 nominally significant results

## effect of milk alpha diversity on infant faecal relative abundances
nrow(inff_RA_div_by_milk_alpha_div_results_ord[inff_RA_div_by_milk_alpha_div_results_ord$FDR<0.05,]) #0 FDR significant results
nrow(inff_RA_div_by_milk_alpha_div_results_ord[inff_RA_div_by_milk_alpha_div_results_ord$p<0.05,])   #3 nominally significant results

## effect of milk relative abundances on infant faecal relative abundances
nrow(inff_RA_div_by_milk_RA_results_ord[inff_RA_div_by_milk_RA_results_ord$FDR<0.05,]) #  8 FDR significant results
nrow(inff_RA_div_by_milk_RA_results_ord[inff_RA_div_by_milk_RA_results_ord$p<0.05,])   #100 nominally significant results

inff_RA_div_by_milk_RA_results_ord[inff_RA_div_by_milk_RA_results_ord$FDR<0.05,c(3:5,7:9)]
#                               x                            y n_total
# 853   milk_g__Lactobacillus_clr    inff_g__Lactobacillus_clr     500
# 815     milk_g__Haemophilus_clr      inff_g__Haemophilus_clr     500
# 807     milk_g__Haemophilus_clr inff_g__Conservatibacter_clr     500
# 748         milk_g__Gemella_clr          inff_g__Gemella_clr     500
# 749         milk_g__Gemella_clr      inff_g__Haemophilus_clr     500
# 922       milk_g__Neisseria_clr    inff_g__Streptococcus_clr     500
# 122   milk_g__Acinetobacter_clr      inff_g__Haemophilus_clr     500
# 426 milk_g__Bifidobacterium_clr   inff_g__Staphylococcus_clr     500
#        estimate            p          FDR
# 853  0.34093377 1.097757e-18 1.185578e-15
# 815  0.32970079 6.613132e-16 3.571092e-13
# 807  0.10287963 2.741145e-07 9.868121e-05
# 748  0.07907869 4.262309e-05 1.150823e-02
# 749  0.16015871 1.258062e-04 2.717414e-02
# 922 -0.19084733 1.677101e-04 2.737151e-02
# 122 -0.11564628 1.774080e-04 2.737151e-02
# 426  0.18868371 3.501406e-04 4.726898e-02


### ===== 8.5 SAVE RESULTS (TABLES S27, S29, S30, S31) ===== ###

setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/effects_of_milk_microbiota")

## TABLE S29
write.table(inff_alpha_div_by_milk_alpha_div_results_ord, file="240305_infant_faecal_alpha_diversity_by_milk_alpha_diversity_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

## TABLE S30
write.table(inff_alpha_div_by_milk_RA_results_ord, file="240305_infant_faecal_alpha_diversity_by_milk_relative_abundances_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

## TABLE S31
write.table(inff_RA_div_by_milk_alpha_div_results_ord, file="240305_infant_faecal_relative_abundances_by_milk_alpha_diversity_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

## TABLE S27
write.table(inff_RA_div_by_milk_RA_results_ord, file="240305_infant_faecal_relative_abundances_by_milk_relative_abundances_model_results.txt", row.names=F, col.names=T, sep="\t", quote=F)




