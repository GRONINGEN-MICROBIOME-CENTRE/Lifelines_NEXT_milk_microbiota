#########################################################
### CHECK PHENOTYPES USED FOR MILK COMPOSITION PAPER  ###
#########################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE
# 0. DATA IMPORT
# 1. CHECK UNIQUE PHENOTYPES USED FOR ASSOCIATION STUDIES


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information")

## import files with included phenotypes

phenos_HMO <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results/240229_phenotypes_associated_with_milk_HMOs_real_n162.txt", header=T, sep="\t", stringsAsFactors=T)
#162x2

phenos_microbiota <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/association_model_results/240229_phenotypes_associated_with_milk_microbiota_real_n185.txt", header=T, sep="\t", stringsAsFactors=T)
#185x2


##### =========================== 1. CHECK UNIQUE PHENOTYPES USED FOR ASSOCIATION STUDIES =========================== #####

## combine data frames and only keep unique rows
phenos <- as.data.frame(rbind(phenos_HMO, phenos_microbiota))
phenos <- phenos[!duplicated(phenos),]

length(phenos$x) #191
length(unique(phenos$x)) #191
nrow(phenos[phenos$phenotype_group!="Human_milk_HMO_concentrations",]) #163
nrow(phenos[phenos$phenotype_group!="Human_milk_HMO_concentrations" & phenos$phenotype_group!="Maternal_genetics",]) #161


table(phenos$phenotype_group, useNA="ifany")
# Human_milk_collection               Infant_feeding 
# 3                                   4 
# Maternal_and_infant_anthropometrics Maternal_diet 
# 8                                   72 
# Maternal_genetics                   Maternal_health_and_diseases 
# 2                                   21 
# Maternal_lifestyle_and_exposures    Maternal_medication_use 
# 12                                  22 
# Pregnancy_and_birth_characteristics Human_milk_HMO_concentrations 
# 19                                  28

## Figure 1 includes 11 more phenotypes, thesee include infant growth parameters (n=4) and infant crying and gastrointestinal health (n=7).

