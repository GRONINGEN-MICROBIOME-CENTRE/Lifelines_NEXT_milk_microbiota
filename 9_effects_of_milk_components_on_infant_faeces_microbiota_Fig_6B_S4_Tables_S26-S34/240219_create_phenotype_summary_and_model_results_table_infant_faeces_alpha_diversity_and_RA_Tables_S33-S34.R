################################################################################################################################
### ASSOCIATION OF MILK GROUPS AND HMOS WITH INFANT FAECAL MICROBIOTA DATA (ALPHA DIV AND REL. ABUNDANCES) (TABLES S33, S34) ###
################################################################################################################################

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
# 4. CREATE SUMMARY STATISTICS FOR NUMERIC/INTEGER PHENOTYPES
# 4.1 FUNCTION TO CREATE SUMMARY STATISTICS
# 4.2 SUMMARY STATISTICS FOR MILK HMOs
# 
# 5. INVERSE-RANK TRANSFORM MILK HMOs AND INFANT FAECAL ALPHA DIVERSITY DATA
# 5.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION
# 5.2. INVERSE-RANK TRANSFORM MILK HMOs, INFANT FAECAL ALPHA DIVERSITY MEASURES AND NUMERIC PHENOTYPES
# 
# 6. CALCULATE INFANT FAECAL RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE
# 
# 7. CLR-TRANSFORM MILK RELATIVE ABUNDANCES
# 7.1 FUNCTION FOR CLR-TRANSFORMATION
# 7.2 CLR-TRANSFORM INFANT FAECAL RELATIVE ABUNDANCES
# 
# 8. MODEL FUNCTION (LMER/LM)
# 
# 9. DETAILS ABOUT inff_invr_RAclr DATA FRAME FOR ASSOCIATION OF BASIC PHENOTYPES AND MILK HMOs WITH INFANT FAECAL MICROBIOTA
# 
# 10. RUN MODELS FOR ASSOCIATION OF BASIC PHENOTYPES WITH INFANT FAECAL MICROBIOTA WITH CORRECTION FOR DNA ISOLATION METHOD, READS, TIME AND ID
# 10.1 RUN MODELS FOR INFANT FAECAL ALPHA DIVERSITY
# 10.2 RUN MODELS FOR INFANT FAECAL RELATIVE ABUNDANCES
# 
# 11. RUN MODELS FOR ASSOCIATION OF MILK GROUP AND MILK HMOs WITH INFANT FAECAL MICROBIOTA WITH CORRECTION FOR DNA ISOLATION METHOD, READS, TIME, ID AND DELIVERY MODE
# 11.1 RUN MODELS FOR INFANT FAECAL ALPHA DIVERSITY (TABLE S33)
# 11.2 RUN MODELS FOR INFANT FAECAL RELATIVE ABUNDANCES (TABLE S34)
# 
# [12. RUN MODELS FOR ASSOCIATION OF MILK GROUP AND MILK HMOs WITH INFANT FAECAL MICROBIOTA WITH CORRECTION FOR DNA ISOLATION METHOD, READS AND DELIVERY MODE - PER TIME POINT]
# [12.1 RUN MODELS FOR INFANT FAECAL ALPHA DIVERSITY]
# [12.2 RUN MODELS FOR INFANT FAECAL RELATIVE ABUNDANCES]


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
                                                         528:530,                 #milk alpha diversity
                                                         531:1154)]               #milk relative bacterial abundances
#552x693


##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

## check that data structure is correct
str(inff[,1:50])
str(inff[,c(51:70,692:693)])
# -> data structure is correct


##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

## Infant phenotypes: for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
inff$infant_birth_delivery_mode <- factor(inff$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
inff$infant_birth_delivery_mode_detailed <- factor(inff$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                        "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
inff$infant_ffq_feeding_type_delivery <- factor(inff$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))

# also change it for milk groups
inff$mother_milk_HMO_milk_group <- factor(inff$mother_milk_HMO_milk_group, levels=c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


## No phenotypes were changed to numeric.


##### =========================== 4. CREATE SUMMARY STATISTICS FOR NUMERIC/INTEGER PHENOTYPES  =========================== #####

### ===== 4.1 FUNCTION TO CREATE SUMMARY STATISTICS ===== ###

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


### ===== 4.2 SUMMARY STATISTICS FOR MILK HMOs ===== ###

## calculate summary statistics
HMO_sumstats <- create.summary.numeric(inff[,39:66])

## resort columns
HMO_sumstats <- HMO_sumstats[,c(1,9:13,2:8)]

## save phenotype summary statistics table
write.table(HMO_sumstats, file="phenotype_summary_statistics/240208_infant_faecal_microbiota_milk_HMOs_summary_statistics_n28.txt", row.names=F, col.names=T, sep="\t", quote=F)


##### =========================== 5. INVERSE-RANK TRANSFORM MILK HMOs AND INFANT FAECAL ALPHA DIVERSITY DATA =========================== #####

### ===== 5.1 FUNCTION FOR INVERSE-RANK TRANSFORMATION ===== ###

invrank <- function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


### ===== 5.2. INVERSE-RANK TRANSFORM MILK HMOs, INFANT FAECAL ALPHA DIVERSITY MEASURES AND NUMERIC PHENOTYPES ===== ###

inff_invr <- inff
inff_invr[,c(29,39:66,67:69)] <- as.data.frame(apply(inff_invr[,c(29,39:66,67:69)], 2, invrank))
colnames(inff_invr)[c(29,39:66,67:69)] <- paste0(colnames(inff_invr)[c(29,39:66,67:69)], "_invr")


##### =========================== 6. CALCULATE INFANT FAECAL RELATIVE ABUNDANCES AND SELECT BACTERIA WITH ≥10% PREVALENCE =========================== #####

## save only columns with absolute bacterial abundances
rownames(inff_invr) <- inff_invr$seq_16S_sample_ID
inffbact <- inff_invr[,70:ncol(inff_invr)] #552x624

## calculate relative bacterial abundances
inffbact_relab <- (inffbact/rowSums(inffbact))
table(rowSums(inffbact_relab), useNA="ifany") #all samples have 1

## exclude column 'Other'
inffbact_relab_tmp <- inffbact_relab[,-624]

## select only bacterial genera with ≥10% prevalence
inffbact_relab_tmp2 <- inffbact_relab_tmp[,colSums(inffbact_relab_tmp>0)>=56]
length(inffbact_relab_tmp2) #30
# colnames(inffbact_relab_tmp2)

## merge back with metadata
inff_invr_RA <- merge(inff_invr[,1:69], inffbact_relab_tmp2, by="row.names")
rownames(inff_invr_RA) <- inff_invr_RA$Row.names
inff_invr_RA <- inff_invr_RA[,-1]


##### =========================== 7. CLR-TRANSFORM MILK RELATIVE ABUNDANCES =========================== #####

### ===== 7.1 FUNCTION FOR CLR-TRANSFORMATION ===== ###

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


### ===== 7.2 CLR-TRANSFORM INFANT FAECAL RELATIVE ABUNDANCES ===== ###

inff_invr_RAclr <- inff_invr_RA
inff_invr_RAclr[,70:ncol(inff_invr_RAclr)] <- as.data.frame(do_clr_externalWeighting(inff_invr_RAclr[,70:ncol(inff_invr_RAclr)], inff_invr_RAclr[,70:ncol(inff_invr_RAclr)]))

colnames(inff_invr_RAclr)[70:ncol(inff_invr_RAclr)] <- c(paste0(colnames(inff_invr_RAclr)[70:ncol(inff_invr_RAclr)], "_clr"))


##### =========================== 8. MODEL FUNCTION (LMER/LM) =========================== #####

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
          facet_grid(.~time_point_numeric)+
          geom_smooth(method='lm')+
          labs(x=paste0(colnames(my_df)[i]), y=paste0(colnames(my_df)[j]),
               title=paste0(colnames(my_df)[j], " by ", colnames(my_df)[i], "\n",
                            "n=", nrow(my_df), "\n", 
                            covardetails, "\n",
                            "p=", signif(p[length(p)],3)))
        
        if (p[length(p)]<0.05){
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/nominally_significant/240306_infant_faeces_", jcolumn ,"_by_", colnames(my_df)[i], "_and_time.pdf"), scatterplot1, dpi=500, width=30, height=15, units="cm")
        }else{
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/not_significant/240306_infant_faeces_", jcolumn ,"_by_", colnames(my_df)[i], "_and_time.pdf"), scatterplot1, dpi=500, width=30, height=15, units="cm")
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
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/nominally_significant/240306_infant_faeces_", jcolumn ,"_by_", colnames(my_df)[i], "_and_time.pdf"), boxplot1, dpi=500, width=30, height=15, units="cm")
        }else{
          ggsave(paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/not_significant/240306_infant_faeces_", jcolumn ,"_by_", colnames(my_df)[i], "_and_time.pdf"), boxplot1, dpi=500, width=30, height=15, units="cm")
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
  write.table(res, file=paste0("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/", outputfolder, "/240306_", datasetname, "_", outputfolder, "_", covardetails, "_results.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  
  #return output table
  return(res)
}


##### =========================== 9. DETAILS ABOUT inff_invr_RAclr DATA FRAME FOR ASSOCIATION OF BASIC PHENOTYPES AND MILK HMOs WITH INFANT FAECAL MICROBIOTA =========================== #####

## Infant data set (inff_invr_RAclr)
## Columns:
## 1:17 = Information columns (infant_ID etc., incl. NEXT_participation_number and time_point)
## 18:19  = DNA_isolation_method and seq_16S_n_reads_clean
## 20:38  = maternal and infant phenotypes (check effect of these: c(22:24,26,28:31,33,37))
## 39:66  = invr-transformed milk HMOs
## 67:69  = invr-transformed infant faecal alpha diversity
## 70:99  = clr-transformed infant faecal relative abundances of bacterial genera with ≥10% prevalence


##### =========================== 10. RUN MODELS FOR ASSOCIATION OF BASIC PHENOTYPES WITH INFANT FAECAL MICROBIOTA WITH CORRECTION FOR DNA ISOLATION METHOD, READS, TIME AND ID =========================== #####

## -> Run models for association of basic maternal/infant phenotypes to decide if any phenotypes should be included as covariates in the other models.


### ===== 10.1 RUN MODELS FOR INFANT FAECAL ALPHA DIVERSITY ===== ###

inff_alpha_by_basic_phenotypes_results <- run.models(datasetname = "infant_faeces_microbiota_by_basic_phenotypes",
                                                     modeltype = "lmer",
                                                     covardetails = "lmer_correction_for_time_ID",
                                                     covariates = "time_point_numeric + (1|mother_ID) + ",
                                                     inputdata   = inff_invr_RAclr,
                                                     xcolumn     = c(18,19,29,30,31,33,37),  #basic phenotypes
                                                     ycolumn     = c(67:69),  #invr-transformed infant alpha diversity measures
                                                     outputfolder = "alpha_diversity")

## exclude results for Simpson, recalculate FDR and save
inff_alpha_by_basic_phenotypes_results_no_simpson <- inff_alpha_by_basic_phenotypes_results[inff_alpha_by_basic_phenotypes_results$y!="alpha_div_simpson_invr",]
inff_alpha_by_basic_phenotypes_results_no_simpson$FDR <- p.adjust(inff_alpha_by_basic_phenotypes_results_no_simpson$p, method="BH")
inff_alpha_by_basic_phenotypes_results_no_simpson <- inff_alpha_by_basic_phenotypes_results_no_simpson[order(inff_alpha_by_basic_phenotypes_results_no_simpson$FDR, inff_alpha_by_basic_phenotypes_results_no_simpson$p),]

write.table(inff_alpha_by_basic_phenotypes_results_no_simpson, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/240306_infant_faeces_microbiota_richness_shannon_by_basic_phenotypes_alpha_diversity_lmer_correction_for_time_ID_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

## check results
nrow(inff_alpha_by_basic_phenotypes_results_no_simpson[inff_alpha_by_basic_phenotypes_results_no_simpson$FDR<0.05,]) #2 FDR significant associations
nrow(inff_alpha_by_basic_phenotypes_results_no_simpson[inff_alpha_by_basic_phenotypes_results_no_simpson$p<0.05,])   #3 nominally significant associations


### ===== 10.2 RUN MODELS FOR INFANT FAECAL RELATIVE ABUNDANCES ===== ###

inff_RA_by_basic_phenotypes_results <- run.models(datasetname = "infant_faeces_microbiota_by_basic_phenotypes",
                                                  modeltype = "lmer",
                                                  covardetails = "lmer_correction_for_time_ID",
                                                  covariates = "time_point_numeric + (1|mother_ID) + ",
                                                  inputdata   = inff_invr_RAclr,
                                                  xcolumn     = c(18,19,29,30,31,33,37),  #basic phenotypes
                                                  ycolumn     = c(70:99),  #clr-transformed relative abundances of infant faecal bacteria with ≥10% prevalence
                                                  outputfolder = "relative_abundances")

## check results
nrow(inff_RA_by_basic_phenotypes_results[inff_RA_by_basic_phenotypes_results$FDR<0.05,]) #26 FDR significant associations
nrow(inff_RA_by_basic_phenotypes_results[inff_RA_by_basic_phenotypes_results$p<0.05,])   #52 nominally significant associations


##### =========================== 11. RUN MODELS FOR ASSOCIATION OF MILK GROUP AND MILK HMOs WITH INFANT FAECAL MICROBIOTA WITH CORRECTION FOR DNA ISOLATION METHOD, READS, TIME, ID AND DELIVERY MODE =========================== #####

### ===== 11.1 RUN MODELS FOR INFANT FAECAL ALPHA DIVERSITY (TABLE S33) ===== ###

inff_alpha_by_milkgroup_HMOs_results <- run.models(datasetname = "infant_faeces_microbiota_by_milkgroup_and_milk_HMOs",
                                                  modeltype = "lmer",
                                                  covardetails = "lmer_correction_for_DNAmethod_reads_time_ID_delmode",
                                                  covariates = "DNA_isolation_method + seq_16S_n_reads_clean + time_point_numeric + (1|mother_ID) + infant_birth_delivery_mode + ",
                                                  inputdata   = inff_invr_RAclr[!is.na(inff_invr_RAclr$infant_birth_delivery_mode),],
                                                  xcolumn     = c(22,39:66),  #maternal milk group and milk HMOs
                                                  ycolumn     = c(67:69),     #invr-transformed infant alpha diversity measures
                                                  outputfolder = "alpha_diversity") 

## exclude results for Simpson, recalculate FDR and save
inff_alpha_by_milkgroup_HMOs_results_no_simpson <- inff_alpha_by_milkgroup_HMOs_results[inff_alpha_by_milkgroup_HMOs_results$y!="alpha_div_simpson_invr",]
inff_alpha_by_milkgroup_HMOs_results_no_simpson$FDR <- p.adjust(inff_alpha_by_milkgroup_HMOs_results_no_simpson$p, method="BH")
inff_alpha_by_milkgroup_HMOs_results_no_simpson <- inff_alpha_by_milkgroup_HMOs_results_no_simpson[order(inff_alpha_by_milkgroup_HMOs_results_no_simpson$FDR, inff_alpha_by_milkgroup_HMOs_results_no_simpson$p),]

write.table(inff_alpha_by_milkgroup_HMOs_results_no_simpson, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/240306_infant_faeces_microbiota_richness_shannon_by_milkgroup_and_milk_HMOs_alpha_diversity_lmer_correction_for_DNAmethod_reads_time_ID_delmode_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# inff_alpha_by_milkgroup_HMOs_results_no_simpson <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/240306_infant_faeces_microbiota_richness_shannon_by_milkgroup_and_milk_HMOs_alpha_diversity_lmer_correction_for_DNAmethod_reads_time_ID_delmode_results.txt", header=T, sep="\t")

## check results
nrow(inff_alpha_by_milkgroup_HMOs_results_no_simpson[inff_alpha_by_milkgroup_HMOs_results_no_simpson$FDR<0.05,]) # 3 FDR significant associations
nrow(inff_alpha_by_milkgroup_HMOs_results_no_simpson[inff_alpha_by_milkgroup_HMOs_results_no_simpson$p<0.05,])   #11 nominally significant associations


## fix n for milk groups
df <- inff_invr_RAclr[!is.na(inff_invr_RAclr$infant_birth_delivery_mode),]
table(df$mother_milk_HMO_milk_group, useNA="ifany")
# Le+Se+ Le+Se- Le-Se+ Le-Se- 
#    376    121     36     10

inff_alpha_by_milkgroup_HMOs_results_no_simpson[!is.na(inff_alpha_by_milkgroup_HMOs_results_no_simpson$levels) & inff_alpha_by_milkgroup_HMOs_results_no_simpson$levels=="Le+Se-", "n_total"] <- 376+121
inff_alpha_by_milkgroup_HMOs_results_no_simpson[!is.na(inff_alpha_by_milkgroup_HMOs_results_no_simpson$levels) & inff_alpha_by_milkgroup_HMOs_results_no_simpson$levels=="Le-Se+", "n_total"] <- 376+36
inff_alpha_by_milkgroup_HMOs_results_no_simpson[!is.na(inff_alpha_by_milkgroup_HMOs_results_no_simpson$levels) & inff_alpha_by_milkgroup_HMOs_results_no_simpson$levels=="Le-Se-", "n_total"] <- 376+10
  
## fix levels for milk groups and add reference level
inff_alpha_by_milkgroup_HMOs_results_no_simpson[inff_alpha_by_milkgroup_HMOs_results_no_simpson$levels=="", "levels"] <- NA

inff_alpha_by_milkgroup_HMOs_results_no_simpson$reference_level <- NA
inff_alpha_by_milkgroup_HMOs_results_no_simpson[inff_alpha_by_milkgroup_HMOs_results_no_simpson$x=="mother_milk_HMO_milk_group", "reference_level"] <- "Le+Se+"

inff_alpha_by_milkgroup_HMOs_results_no_simpson <- inff_alpha_by_milkgroup_HMOs_results_no_simpson[,c(1:5,10,6:9)]

write.table(inff_alpha_by_milkgroup_HMOs_results_no_simpson, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/240306_infant_faeces_microbiota_richness_shannon_by_milkgroup_and_milk_HMOs_alpha_diversity_lmer_correction_for_DNAmethod_reads_time_ID_delmode_results.txt", row.names=F, col.names=T, sep="\t", quote=F)

  
### ===== 11.2 RUN MODELS FOR INFANT FAECAL RELATIVE ABUNDANCES (TABLE S34) ===== ###

inff_RA_by_milkgroup_HMOs_results <- run.models(datasetname = "infant_faeces_microbiota_by_milkgroup_and_milk_HMOs",
                                                   modeltype = "lmer",
                                                   covardetails = "lmer_correction_for_DNAmethod_reads_time_ID_delmode",
                                                   covariates = "DNA_isolation_method + seq_16S_n_reads_clean + time_point_numeric + (1|mother_ID) + infant_birth_delivery_mode + ",
                                                   inputdata   = inff_invr_RAclr[!is.na(inff_invr_RAclr$infant_birth_delivery_mode),],
                                                   xcolumn     = c(22,39:66),  #maternal milk group and milk HMOs
                                                   ycolumn     = c(70:99),     #clr-transformed relative abundances of infant faecal bacteria with ≥10% prevalence
                                                   outputfolder = "relative_abundances") 

# inff_RA_by_milkgroup_HMOs_results <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/240306_infant_faeces_microbiota_by_milkgroup_and_milk_HMOs_relative_abundances_lmer_correction_for_DNAmethod_reads_time_ID_delmode_results.txt", header=T, sep="\t")


## check results
nrow(inff_RA_by_milkgroup_HMOs_results[inff_RA_by_milkgroup_HMOs_results$FDR<0.05,]) # 4 FDR significant associations
nrow(inff_RA_by_milkgroup_HMOs_results[inff_RA_by_milkgroup_HMOs_results$p<0.05,])   #88 nominally significant associations


## fix n for milk groups
# df <- inff_invr_RAclr[!is.na(inff_invr_RAclr$infant_birth_delivery_mode),]
table(df$mother_milk_HMO_milk_group, useNA="ifany")
# Le+Se+ Le+Se- Le-Se+ Le-Se- 
#    376    121     36     10

inff_RA_by_milkgroup_HMOs_results[!is.na(inff_RA_by_milkgroup_HMOs_results$levels) & inff_RA_by_milkgroup_HMOs_results$levels=="Le+Se-", "n_total"] <- 376+121
inff_RA_by_milkgroup_HMOs_results[!is.na(inff_RA_by_milkgroup_HMOs_results$levels) & inff_RA_by_milkgroup_HMOs_results$levels=="Le-Se+", "n_total"] <- 376+36
inff_RA_by_milkgroup_HMOs_results[!is.na(inff_RA_by_milkgroup_HMOs_results$levels) & inff_RA_by_milkgroup_HMOs_results$levels=="Le-Se-", "n_total"] <- 376+10

## fix levels for milk groups and add reference level
inff_RA_by_milkgroup_HMOs_results[inff_RA_by_milkgroup_HMOs_results$levels=="", "levels"] <- NA

inff_RA_by_milkgroup_HMOs_results$reference_level <- NA
inff_RA_by_milkgroup_HMOs_results[inff_RA_by_milkgroup_HMOs_results$x=="mother_milk_HMO_milk_group", "reference_level"] <- "Le+Se+"

inff_RA_by_milkgroup_HMOs_results <- inff_RA_by_milkgroup_HMOs_results[,c(1:5,10,6:9)]


write.table(inff_RA_by_milkgroup_HMOs_results, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/240306_infant_faeces_microbiota_by_milkgroup_and_milk_HMOs_relative_abundances_lmer_correction_for_DNAmethod_reads_time_ID_delmode_results.txt", row.names=F, col.names=T, sep="\t", quote=F)


# ##### =========================== 12. RUN MODELS FOR ASSOCIATION OF MILK GROUP AND MILK HMOs WITH INFANT FAECAL MICROBIOTA WITH CORRECTION FOR DNA ISOLATION METHOD, READS AND DELIVERY MODE - PER TIME POINT =========================== #####
# 
# ### ===== 12.1 RUN MODELS FOR INFANT FAECAL ALPHA DIVERSITY ===== ###
# 
# inff_alpha_by_milkgroup_HMOs_results <- run.models(datasetname = "infant_faeces_microbiota_by_milkgroup_and_milk_HMOs",
#                                                    modeltype = "lmer",
#                                                    covardetails = "lmer_correction_for_DNAmethod_reads_time_ID_delmode",
#                                                    covariates = "DNA_isolation_method + seq_16S_n_reads_clean + time_point_numeric + (1|mother_ID) + infant_birth_delivery_mode + ",
#                                                    inputdata   = inff_invr_RAclr,
#                                                    xcolumn     = c(22,39:66),  #maternal milk group and milk HMOs
#                                                    ycolumn     = c(67:69),     #invr-transformed infant alpha diversity measures
#                                                    outputfolder = "alpha_diversity") 
# 
# ## exclude results for Simpson, recalculate FDR and save
# inff_alpha_by_milkgroup_HMOs_results_no_simpson <- inff_alpha_by_milkgroup_HMOs_results[inff_alpha_by_milkgroup_HMOs_results$y!="alpha_div_simpson_invr",]
# inff_alpha_by_milkgroup_HMOs_results_no_simpson$FDR <- p.adjust(inff_alpha_by_milkgroup_HMOs_results_no_simpson$p, method="BH")
# inff_alpha_by_milkgroup_HMOs_results_no_simpson <- inff_alpha_by_milkgroup_HMOs_results_no_simpson[order(inff_alpha_by_milkgroup_HMOs_results_no_simpson$FDR, inff_alpha_by_milkgroup_HMOs_results_no_simpson$p),]
# 
# write.table(inff_alpha_by_milkgroup_HMOs_results_no_simpson, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/240219_infant_faeces_microbiota_richness_shannon_by_milkgroup_and_milk_HMOs_alpha_diversity_lmer_correction_for_DNAmethod_reads_time_ID_delmode_results.txt", row.names=F, col.names=T, sep="\t", quote=F)
# 
# ## check results
# nrow(inff_alpha_by_milkgroup_HMOs_results_no_simpson[inff_alpha_by_milkgroup_HMOs_results_no_simpson$FDR<0.05,]) # 3 FDR significant associations
# nrow(inff_alpha_by_milkgroup_HMOs_results_no_simpson[inff_alpha_by_milkgroup_HMOs_results_no_simpson$p<0.05,])   #11 nominally significant associations
# 
# 
# ### ===== 12.2 RUN MODELS FOR INFANT FAECAL RELATIVE ABUNDANCES ===== ###
# 
# inff_RA_by_milkgroup_HMOs_results <- run.models(datasetname = "infant_faeces_microbiota_by_milkgroup_and_milk_HMOs",
#                                                 modeltype = "lmer",
#                                                 covardetails = "lmer_correction_for_DNAmethod_reads_time_ID_delmode",
#                                                 covariates = "DNA_isolation_method + seq_16S_n_reads_clean + time_point_numeric + (1|mother_ID) + infant_birth_delivery_mode + ",
#                                                 inputdata   = inff_invr_RAclr,
#                                                 xcolumn     = c(22,39:66),  #maternal milk group and milk HMOs
#                                                 ycolumn     = c(70:99),     #clr-transformed relative abundances of infant faecal bacteria with ≥10% prevalence
#                                                 outputfolder = "relative_abundances") 
# 
# ## check results
# nrow(inff_RA_by_milkgroup_HMOs_results[inff_RA_by_milkgroup_HMOs_results$FDR<0.05,]) # 4 FDR significant associations
# nrow(inff_RA_by_milkgroup_HMOs_results[inff_RA_by_milkgroup_HMOs_results$p<0.05,])   #88 nominally significant associations



## run per time point and check if trends/results remain
## check in breast-fed only infants
