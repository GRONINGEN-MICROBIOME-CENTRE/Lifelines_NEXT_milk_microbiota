### ========================== SUPPLEMENTARY DATA FOR MICROBIOTA CONTROLS (FIGURE S6, TABLE S41) ========================== ###

## SCRIPT:           CREATE SUPPLEMENTARY FIGURES/TABLES FOR POSITIVE AND NEGATIVE MICROBIOTA CONTROLS FOR MILK COMPOSITION PAPER
## DESCRIPTION:      The data used here is from the dada2 pipeline WITH filtering to only include ASVs with lengths 400-431 bp.
## AUTHORS:          Johanne Spreckels
## NOTES:            Note that there are samples with 0 reads, i.e. they failed 16S sequencing / had no reads left after processing the 16S sequencing data.
## DATE OF CREATION: October 2023

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23):     R

## CONTENTS OF THIS FILE

# 1. DATA IMPORT AND DATA PREPARATION
# 1.1. IMPORT DATA
# 1.2. CHECK AND ENSURE CORRECT DATA STRUCTURE
# 1.3. SELECT DATA OF INTEREST
# 
# 2. PREPARE POSITIVE CONTROL DATA SET WITH RELATIVE BACTERIAL ABUNDANCES
# 2.1. POSITIVE CONTROLS
# 
# 3. PLOT RELATIVE BACTERIAL ABUNDANCES IN POSITIVE CONTROLS (FIGURE S6)
# 
# 4. SAVE TABLE SHOWING RELATIVE BACTERIAL ABUNDANCES IN POSITIVE CONTROLS
# 
# 5. PREPARE POSITIVE CONTROL DATA SET WITH RELATIVE BACTERIAL ABUNDANCES
# 
# 6. PLOT RELATIVE BACTERIAL ABUNDANCES IN NEGATIVE CONTROLS
# 
# 7. SAVE TABLE SHOWING RELATIVE BACTERIAL ABUNDANCES IN NEGATIVE CONTROLS (TABLE S41)


##### =========================== 1. DATA IMPORT AND DATA PREPARATION =========================== #####

### ===== 1.1. IMPORT DATA ===== ###

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


### ===== 1.2. CHECK AND ENSURE CORRECT DATA STRUCTURE ===== ###

## check data structure
str(data[,1:50])
str(data[,51:100])
str(data[,100:150])
str(data[,151:200])
str(data[,200:250])
str(data[,251:300])
str(data[,300:350])
str(data[,351:400])
str(data[,400:450])
str(data[,451:500])
str(data[,500:530])
str(data[,c(531,1153,1154)])
# -> looks good

## check sample numbers per sample type in data frame 'data' -> note that this data frame includes duplicated maternal data if a mother had twins!
dim(data)
table(data$sample_origin_type)
# 737 human milk
# 552 infant faeces
# 172 maternal faeces
# 7 isolation negative controls
# 25 lib prep negative controls
# 1 isolation positive control
# 21 lib prep positive controls


### ===== 1.3. SELECT DATA OF INTEREST ===== ###

## select only positive and negative controls
pos <- data[data$sample_origin=="positive_control",] #22x1154
neg <- data[data$sample_origin=="negative_control",] #32x1154

table(pos$sample_origin_type, useNA="ifany")
# 1 isolation mock
# 21 DNA mocks

table(neg$sample_origin_type, useNA="ifany")
# 7 isolation negative controls
# 25 library preparation negative controls


##### =========================== 2. PREPARE POSITIVE CONTROL DATA SET WITH RELATIVE BACTERIAL ABUNDANCES =========================== #####

### ===== 2.1. POSITIVE CONTROLS ===== ###

## save only columns with absolute bacterial abundances
rownames(pos) <- pos$seq_16S_sample_ID
posbact <- pos[,531:ncol(pos)] #22x624

## calculate relative bacterial abundances
posbact_relab <- (posbact/rowSums(posbact))
table(rowSums(posbact_relab), useNA="ifany") #21 have 1, DNA_mock_9 has NA (this sample failed seq)

## set relative abundances to 0 for DNA_mock_9
posbact_relab[rownames(posbact_relab)=="DNA_mock_9",] <- 0
table(rowSums(posbact_relab), useNA="ifany") #now 21 have 1 and DNA_mock_9 has 0

## remove bacterial genera not detected in any positive control
## exclude column 'Other'
posbact_relab_v2 <- posbact_relab[,-624]

## exclude columns for bacteria that were not identified in any sample
posbact_relab_v2 <- posbact_relab_v2[,colSums(posbact_relab_v2)>0]
dim(posbact_relab_v2) #22x11

## merge back with metadata
posRA <- merge(pos[,c(504:509,525,528:530)], posbact_relab_v2, by="row.names")
rownames(posRA) <- posRA$Row.names
posRA <- posRA[,-1]
posRA <- droplevels(posRA)

# All mocks contain 8 expected genera, except for DNA_mock_4, which contains 3 additional bacteria (Gemella, Streptococcus and Veillonella, each with <0.1% relative abundance)

## shorten colnames
colnames(posRA) <- gsub(".*g__", "", colnames(posRA))

## add theoretical mock community composition
# theoretical mock composition in 16S sequencing:
#    % Bacterial species         Note
#  4.2 Pseudomonas aeruginosa    
# 10.1 Escherichia coli          
# 10.4 Salmonella enterica       
# 18.4 Lactobacillus fermentum   now called Limosilactobacillus fermentum
#  9.9 Enterococcus faecalis     
# 15.5 Staphylococcus aureus     
# 14.1 Listeria monocytogenes    
# 17.4 Bacillus subtilis    

for (i in c(1,4,6)){posRA[,i] <- as.character(as.factor(posRA[,i]))}
posRA[23,] <- c("expected", rep(NA,2), "expected", NA, "expected",  rep(NA,4), 0.174, 0.099, 0.101, 0, 0.184, 0.141, 0.042, 0.104, 0.155, 0, 0)
for (i in c(1,4,6)){posRA[,i] <- as.factor(as.character(posRA[,i]))}
rownames(posRA)[23] <- "expected"
for (i in c(7:ncol(posRA))){posRA[,i] <- as.numeric(as.character(posRA[,i]))}

## convert data frame from wide to long format
library(reshape)
posRAlong <- melt(posRA, id.vars=colnames(posRA)[c(1:10)], variable_name="Genus")

## resort positive controls
posRAlong$seq_16S_sample_ID <- factor(posRAlong$seq_16S_sample_ID, levels = c("expected",
                                                                              "isolation_mock",
                                                                              "DNA_mock_1", "DNA_mock_2", "DNA_mock_3", "DNA_mock_4", "DNA_mock_5",
                                                                              "DNA_mock_6", "DNA_mock_7", "DNA_mock_8", "DNA_mock_9", "DNA_mock_10",
                                                                              "DNA_mock_11", "DNA_mock_12", "DNA_mock_13", "DNA_mock_14", "DNA_mock_15",
                                                                              "DNA_mock_16", "DNA_mock_17", "DNA_mock_18", "DNA_mock_19", "DNA_mock_20",
                                                                              "DNA_mock_21"))


##### =========================== 3. PLOT RELATIVE BACTERIAL ABUNDANCES IN POSITIVE CONTROLS (FIGURE S6) =========================== #####

library(ggplot2)

## create stacked barplot for isolation positive controls + theoretical mock
plotRAposctrl1 <- ggplot(posRAlong, aes(x=seq_16S_sample_ID, y=value, fill=Genus))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90),
        text = element_text(size=20),
        legend.position = "right")+
  labs(x="", y="Relative abundance", title="Relative abundances in positive controls")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  # facet_grid(.~sample_type, scales = "free_x", space = "free_x")+
  scale_fill_manual(values=c("#CC79A7","#009E73","#56B4E9","indianred3","#F0E442","#D55E00","#0072B2","darkred","#E69F00","palegreen","lightskyblue1"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/results_controls_16S_seq/231218_suppl_relative_abundances_plot_positive_controls.pdf", plotRAposctrl1, dpi=500, width=25, height=15, units="cm", useDingbats=F)


##### =========================== 4. SAVE TABLE SHOWING RELATIVE BACTERIAL ABUNDANCES IN POSITIVE CONTROLS =========================== #####

write.table(posRA, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/results_controls_16S_seq/231218_suppl_relative_abundances_table_positive_controls.txt", sep="\t", row.names=F, quote=F)


##### =========================== 5. PREPARE POSITIVE CONTROL DATA SET WITH RELATIVE BACTERIAL ABUNDANCES   =========================== #####

## save only columns with absolute bacterial abundances
rownames(neg) <- neg$seq_16S_sample_ID
negbact <- neg[,531:ncol(neg)] #32x624

## calculate relative bacterial abundances
negbact_relab <- (negbact/rowSums(negbact))
table(rowSums(negbact_relab), useNA="ifany") #24 have 1, 8 have NA (these samples failed seq)
neg_ctrls_0_reads <- neg[neg$seq_16S_dada2_n_reads_final==0,"seq_16S_sample_ID"]
# Batch2_NTC_DEPC_1 LibPrep_NTC_11    LibPrep_NTC_12    LibPrep_NTC_13   
# LibPrep_NTC_16    LibPrep_NTC_23    LibPrep_NTC_8     LibPrep_NTC_9 

## set relative abundances to 0 for DNA_mock_9
negbact_relab[rownames(negbact_relab) %in% neg_ctrls_0_reads,] <- 0
table(rowSums(negbact_relab), useNA="ifany") #now 24 have 1, 8 have 0

## remove bacterial genera not detected in any negative control
## exclude column 'Other'
negbact_relab_v2 <- negbact_relab[,-624]

## exclude columns for bacteria that were not identified in any sample
negbact_relab_v2 <- negbact_relab_v2[,colSums(negbact_relab_v2)>0]
dim(negbact_relab_v2) #32x82

## abundance filtering: filter on bacteria with at least 2% relative abundance in any sample
maxvector <- c()
for (i in 1:ncol(negbact_relab_v2)){maxvector[i] <- max(negbact_relab_v2[,i])}
negbact_relab_v2[33,] <- maxvector
negbact_relab_v2filt <- negbact_relab_v2[,negbact_relab_v2[33,]>0.02] #64 genera kept, 18 removed
negbact_relab_v2 <- negbact_relab_v2[-33,]
negbact_relab_v2filt <- negbact_relab_v2filt[-33,]

# prevalence filtering: filter on bacteria present in at least 2 samples
negbact_relab_v2filt2 <- negbact_relab_v2filt[,colSums(negbact_relab_v2filt>0)>=2] #24 genera kept, 40 removed

## sum up excluded genera in 'Other' (filtered genera)
negbact_relab_v2filt2$Other <- 1-rowSums(negbact_relab_v2filt2)

## set 'Other' to 0 for the 8 negative controls that had 0 reads
negbact_relab_v2filt2[rownames(negbact_relab_v2filt2) %in% neg_ctrls_0_reads, "Other"] <- 0

## merge back with metadata
negRA <- merge(neg[,c(504:509,525,528:530)], negbact_relab_v2filt2, by="row.names")
rownames(negRA) <- negRA$Row.names
negRA <- negRA[,-1]
negRA <- droplevels(negRA)

## shorten colnames
colnames(negRA) <- gsub(".*g__", "", colnames(negRA))

## combine seq_16S_sample_ID and seq_16S_dada2_n_reads_final in 1 column
negRA$sample_ID_reads <- as.factor(as.character(c(paste0(negRA$seq_16S_sample_ID, "__", negRA$seq_16S_dada2_n_reads_final, "_reads"))))
negRA <- negRA[,c(36,1:35)]

## convert data frame from wide to long format
# library(reshape)
negRAlong <- melt(negRA, id.vars=colnames(negRA)[c(1:11)], variable_name="Genus")

## resort positive controls
negRAlong$sample_ID_reads <- factor(negRAlong$sample_ID_reads, levels = c("Batch1_NTC_CD1_1__161_reads", "Batch1_NTC_CD2_1__217_reads", "Batch1_NTC_CD2_2__101_reads", "Batch1_NTC_DEPC_1__348_reads",
                                                                          "Batch2_NTC_CD1_1__6102_reads", "Batch2_NTC_CD1_2__397_reads", "Batch2_NTC_DEPC_1__0_reads",
                                                                          "LibPrep_NTC_1__44_reads", "LibPrep_NTC_2__528_reads", "LibPrep_NTC_3__25_reads", "LibPrep_NTC_4__8_reads", "LibPrep_NTC_5__28_reads",
                                                                          "LibPrep_NTC_6__2608_reads", "LibPrep_NTC_7__1864_reads", "LibPrep_NTC_8__0_reads", "LibPrep_NTC_9__0_reads", "LibPrep_NTC_10__1_reads",
                                                                          "LibPrep_NTC_11__0_reads", "LibPrep_NTC_12__0_reads", "LibPrep_NTC_13__0_reads", "LibPrep_NTC_14__1_reads", "LibPrep_NTC_15__1_reads",
                                                                          "LibPrep_NTC_16__0_reads", "LibPrep_NTC_17__1341_reads", "LibPrep_NTC_18__140_reads", "LibPrep_NTC_19__377_reads", "LibPrep_NTC_20__260_reads",
                                                                          "LibPrep_NTC_21__4602_reads", "LibPrep_NTC_22__28_reads", "LibPrep_NTC_23__0_reads", "LibPrep_NTC_24__31_reads", "LibPrep_NTC_25__17_reads"))


##### =========================== 6. PLOT RELATIVE BACTERIAL ABUNDANCES IN NEGATIVE CONTROLS =========================== #####

## create stacked barplot for isolation negative controls
plotRAnegctrl1 <- ggplot(negRAlong, aes(x=sample_ID_reads, y=value, fill=Genus))+
  geom_bar(stat="identity", position="fill")+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_text(angle=90),
        text = element_text(size=20),
        legend.position = "right")+
  labs(x="", y="Relative abundance", title="Relative abundances in negative controls")+
  scale_y_continuous(limits=c(0,1),
                     breaks=seq(from=0,to=1,by=0.2),
                     labels=c("0%","20%","40%","60%","80%","100%"))+
  facet_grid(.~sample_type, scales="free_x", space="free_x")+
  scale_fill_manual(values=c("honeydew4", "lightpink3", "lightgreen", "slategray3", "orange",
                             "blue", "goldenrod4", "red", "seagreen", "olivedrab1",
                             "#56B4E9", "rosybrown1", "turquoise3", "yellow", "lightblue1",
                             "seagreen1", "darkslateblue", "lightgoldenrod", "orchid2", "sienna",
                             "#E69F00", "navy", "hotpink3", "dodgerblue", "grey"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/results_controls_16S_seq/231218_suppl_relative_abundances_plot_negative_controls.pdf", plotRAnegctrl1, dpi=500, width=50, height=20, units="cm", useDingbats=F)


##### =========================== 7. SAVE TABLE SHOWING RELATIVE BACTERIAL ABUNDANCES IN NEGATIVE CONTROLS (TABLE S41) =========================== #####

## extract negative control RA incl. from Other
negbact_relab_with_other <- negbact_relab[,colSums(negbact_relab)>0]
dim(negbact_relab_with_other) #32x83
table(rowSums(negbact_relab_with_other)) #24 have 1, 8 have 0 (failed seq)

## merge back with metadata
negRA_with_other <- merge(neg[,c(507,504:505,525)], negbact_relab_with_other, by="row.names")
# rownames(negRA_with_other) <- negRA_with_other$Row.names
negRA_with_other <- negRA_with_other[,-1]
negRA_with_other <- droplevels(negRA_with_other)

## shorten colnames
colnames(negRA_with_other) <- gsub(".*g__", "", colnames(negRA_with_other))

## resort controls
negRA_with_other_v2 <- negRA_with_other[c(1:7,8,19,26:32,9:18,20:25),]

## add "run_" to run ID
negRA_with_other_v2$seq_16S_run_ID <- c(paste0("run_ID_", negRA_with_other_v2$seq_16S_run_ID))

## change colnames
colnames(negRA_with_other_v2)[1:4] <- c("Control_type", "sample_ID", "sequencing_run_ID", "n_reads")

write.table(negRA_with_other_v2, file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/results_controls_16S_seq/231218_suppl_relative_abundances_table_negative_controls.txt", sep="\t", row.names=F, quote=F)





