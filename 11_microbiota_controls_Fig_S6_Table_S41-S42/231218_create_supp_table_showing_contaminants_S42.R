### ========================== SUPPLEMENTARY TABLE: IDENTIFIED CONTAMINANTS PER SAMPLE TYPE (TABLE S42) ========================== ###

## SCRIPT:           CREATE SUPPLEMENTARY TABLE SHOOWING IDENTIFIED CONTAMINANTS PER SAMPLE TYPE / CONTROL TYPE
## DESCRIPTION:      Here, I create the supplementary table showing which bacteria SCRuB identified as contaminants.
##                   The data used here is from the dada2 pipeline WITH filtering to only include ASVs with lengths 400-431 bp.
## AUTHORS:          Johanne Spreckels
## NOTES:            
## DATE OF CREATION: December 2023

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23):     R

## CONTENTS OF THIS FILE

# 1. IMPORT DATA
# 1.1 IMPORT FILES WITH INFO ON CONTAMINANTS PER DATA SET
# 
# 2. POLISH DATA FILE ENTRIES
# 
# 3. COMBINE DATA FILES
# 
# 4. SAVE COMBINED TABLE (TABLE S42)


### ==================== 1. IMPORT DATA ==================== ###

### ===== 1.1 IMPORT FILES WITH INFO ON CONTAMINANTS PER DATA SET ===== ###

## load files with contaminants identified in milk isolation negative controls (c1) and library prep negative controls (c2)
c1milk <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination_milk_v2/seq_16S_milk_step1_contaminants.rds") #290x4
c2milk <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination_milk_v2/seq_16S_milk_step2_contaminants.rds") #221x4

## load file with contaminants identified in maternal faecal library prep negative controls
cmatf <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/seq_16S_mother_feces_contaminants.rds") #131x4

## load file with contaminants identified in infant faecal library prep negative controls
cinff <- readRDS("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/data_after_decontamination/seq_16S_infant_feces_contaminants.rds") #131x4


### ==================== 2. POLISH DATA FILE ENTRIES ==================== ###

## change colnames
colnames(c1milk) <- c("Bacterial_ASV", "Relative_abundance", "Contamination_source", "Batch")
colnames(c2milk) <- c("Bacterial_ASV", "Relative_abundance", "Contamination_source", "Batch")
colnames(cmatf) <- c("Bacterial_ASV", "Relative_abundance", "Contamination_source", "Batch")
colnames(cinff) <- c("Bacterial_ASV", "Relative_abundance", "Contamination_source", "Batch")

## change data structure
for (i in c(1,3:4)){c1milk[,i] <- as.factor(as.character(c1milk[,i]))}
for (i in c(1,3:4)){c2milk[,i] <- as.factor(as.character(c2milk[,i]))}
for (i in c(1,3:4)){cmatf[,i] <- as.factor(as.character(cmatf[,i]))}
for (i in c(1,3:4)){cinff[,i] <- as.factor(as.character(cinff[,i]))}

## rename data entries
levels(c1milk$Batch) <- c("Milk_DNA_isolation_batch_1", "Milk_DNA_isolation_batch_2")
levels(c2milk$Batch) <- c(paste0("Milk_sequencing_run_ID_", levels(c2milk$Batch)))
levels(cmatf$Batch) <- c(paste0("Maternal_faeces_sequencing_run_ID_", levels(cmatf$Batch)))
levels(cinff$Batch) <- c(paste0("Infant_faeces_sequencing_run_ID_", levels(cinff$Batch)))

## add column showing in from which sample type contaminants were removed
c1milk$Sample_type <- "Mother_human_milk"
c2milk$Sample_type <- "Mother_human_milk"
cmatf$Sample_type <- "Mother_faeces"
cinff$Sample_type <- "Infant_faeces"


### ==================== 3. COMBINE DATA FILES ==================== ###

## combine files
c <- as.data.frame(rbind(c1milk, c2milk, cmatf, cinff))

## fix data structure
c$Sample_type <- as.factor(as.character(c$Sample_type))

## resort columns
c2 <- c[,c(5,4,3,1,2)]


### ==================== 4. SAVE COMBINED TABLE (TABLE S42) ==================== ###

write.table(c2, "/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/microbiota_data/231218_contaminants_identified_in_all_sample_types.txt", sep="\t", row.names=F, quote=F)




