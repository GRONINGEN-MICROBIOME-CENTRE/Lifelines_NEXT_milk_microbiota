##################################################################################################
### CREATE HEATMAP PLOTS FOR MODEL RESULTS TABLE FOR REAL HMO DATA FOR MILK COMPOSITION PAPER  ###
##################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE
# 0. DATA IMPORT
# 1. CHECK AND PREPARE RESULTS FOR PLOTTING HEATMAPS
# 1.1 CHECK RESULTS
# 1.2 CHECK N PER PHENOTYPE GROUP
# 1.3 PREPARE DATA FOR PLOTTING
# 2. PLOT HEATMAPS
# 2.1 FUNCTION TO CREATE HEATMAPS
# 2.2 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED
# 2.3 CREATE HEATMAPS FOR ALL FDR-SIGNIFICANT ASSOCIATIONS
# 2.4 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED - ONLY INCLUDING PHENOTYPES WITH AT LEAST X NOMINALLY SIGNIFICANT ASSOCIATIONS (FIGURE 3A)
# 2.5 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED - ONLY INCLUDING FDR-SIGNIFICANT ASSOCIATIONS OR PHENOTYPES WITH AT LEAST X NOMINALLY SIGNIFICANT ASSOCIATIONS
# 2.6 CREATE HEATMAPS PER PHENOTYPE GROUP


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/")
# setwd("/Users/johanne/Desktop/PhD-Groningen/Thesis_chapters/Chapter3_milk_composition/Results/4_HMO_phenotype_associations/")

## import results table
# model_results_short <- readRDS(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results/240227_model_results_milk_HMOs_real_by_phenotypes_n162.rds")
# model_results_short <- readRDS(file="/Users/johanne/Desktop/PhD-Groningen/Thesis_chapters/Chapter3_milk_composition/Results/4_HMO_phenotype_associations/association_model_results/240227_model_results_milk_HMOs_real_by_phenotypes_n162.rds")
resHMO <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/4_HMO_phenotype_associations/association_model_results/240229_polished_model_results_milk_HMOs_real_by_phenotypes_n162.txt",  header=T, sep="\t", stringsAsFactors=T)
# resHMO <- read.table(file="association_model_results/240227_polished_model_results_milk_HMOs_real_by_phenotypes_n162.txt",  header=T, sep="\t", stringsAsFactors=T)
#4816x15


##### =========================== 1. CHECK AND PREPARE RESULTS FOR PLOTTING HEATMAPS =========================== #####

### ===== 1.1 CHECK RESULTS ===== ###

nrow(resHMO[resHMO$FDR<0.05,]) # 12
nrow(resHMO[resHMO$p<0.05,])   #419


### ===== 1.2 CHECK N PER PHENOTYPE GROUP ===== ###

## check unique phenotypes
# unique(resHMO$Phenotype)[order(unique(resHMO$Phenotype))] #n=162
# unique(resHMO$Phenotype_label_in_plots)[order(unique(resHMO$Phenotype_label_in_plots))] #n=172

## check n per phenotype group
phenos_info <- resHMO[,c(1,2)]
phenos_info <- phenos_info[!duplicated(phenos_info),]
table(phenos_info$Phenotype_group, useNA="ifany")
# Human_milk_collection               Infant_feeding   Maternal_and_infant_anthropometrics  Maternal_diet                   Maternal_genetics 
# 3                                   4                8                                    72                              1 
# Maternal_health_and_diseases    Maternal_lifestyle_and_exposures    Maternal_medication_use   Pregnancy_and_birth_characteristics 
# 21                              12                                  22                        19


### ===== 1.3 PREPARE DATA FOR PLOTTING ===== ###

## save data separately to adapt it for plotting
resHMO_for_plotting <- resHMO

## set all estimates of associations with p>0.05 to 0 so they show white in the plot
resHMO_for_plotting[resHMO_for_plotting$p>0.05, "Estimate"] <- 0

## determine labels for FDR and nominally significant results in plot
resHMO_for_plotting$my_label <- NA
resHMO_for_plotting[resHMO_for_plotting$p>0.05, "my_label"] <- ""
resHMO_for_plotting[resHMO_for_plotting$p<0.05 & resHMO_for_plotting$Estimate<0, "my_label"] <- "-"
resHMO_for_plotting[resHMO_for_plotting$p<0.05 & resHMO_for_plotting$Estimate>0, "my_label"] <- "+"
resHMO_for_plotting[resHMO_for_plotting$FDR<0.05 & resHMO_for_plotting$Estimate<0, "my_label"] <- "-*"
resHMO_for_plotting[resHMO_for_plotting$FDR<0.05 & resHMO_for_plotting$Estimate>0, "my_label"] <- "+*"
table(resHMO_for_plotting$my_label, useNA="ifany")

## fix labels from trait_name_in_plot: replace _ with space
resHMO_for_plotting$Phenotype_label_in_plots <- gsub("_", " ", resHMO_for_plotting$Phenotype_label_in_plots)

## add n to trait label
resHMO_for_plotting$Phenotype_label_in_plots <- c(paste0(resHMO_for_plotting$Phenotype_label_in_plots, " ", resHMO_for_plotting$n))

## reverse the order of the levels so they show them alphabetically from top to bottom in the plot
resHMO_for_plotting$Phenotype_label_in_plots <- as.factor(as.character(resHMO_for_plotting$Phenotype_label_in_plots))

# library(forcats)
levels(resHMO_for_plotting$Phenotype_label_in_plots)
resHMO_for_plotting$Phenotype_label_in_plots <- factor(resHMO_for_plotting$Phenotype_label_in_plots,
                                                       levels = rev(levels(resHMO_for_plotting$Phenotype_label_in_plots)))
levels(resHMO_for_plotting$Phenotype_label_in_plots)


## fix HMO labels and set the order of milk HMOs as in other figures
resHMO_for_plotting$plot_HMO <- resHMO_for_plotting$HMO
# resHMO_for_plotting$plot_HMO <- as.factor(as.character(gsub("mother_milk_HMO_", "", as.character(as.factor(resHMO_for_plotting$plot_HMO)))))
# resHMO_for_plotting$plot_HMO <- as.factor(as.character(gsub("_ugml_invr", "", as.character(as.factor(resHMO_for_plotting$plot_HMO)))))
resHMO_for_plotting$plot_HMO <- factor(resHMO_for_plotting$plot_HMO,
                                                 levels=c("Total", "Neut", "Fuc", "Sia", #total and groups
                                                          "3GL","6GL","LNH","LNnT","LNT", #neutral HMOs
                                                          "2FL","3FL","A_tetra","DFLNHa","LDFT", #neutral+fucosylated HMOs
                                                          "LNDH_I","LNFP_I","LNFP_II","LNFP_III","LNFP_V",
                                                          "LNnDFH","LNnFP_V","MFLNH_III",
                                                          "3F3SL",
                                                          "3SL","6SL","DSLNT","LSTb","LSTc"))
levels(resHMO_for_plotting$plot_HMO) <- c("Total HMOs", "Neutral HMOs", "Fucosylated HMOs", "Sialylated HMOs",
                                                           "3'GL", "6'GL", "LNH", "LNnT", "LNT",
                                                           "2'FL", "3'FL", "A-tetra", "DFLNHa", "LDFT", 
                                                           "LNDH I", "LNFP I", "LNFP II", "LNFP III",
                                                           "LNFP V", "LNnDFH", "LNnFP V", "MFLNH III", "3'F3'SL",
                                                           "3'SL", "6'SL", "DSLNT", "LSTb", "LSTc")


##### =========================== 2. PLOT HEATMAPS =========================== #####

library(ggplot2)

### ===== 2.1 FUNCTION TO CREATE HEATMAPS ===== ###

create.heatmap <- function(inputdata, outputname, w, h){
  
  ## create plot
  myplot <- ggplot(inputdata, aes(x=plot_HMO, y=Phenotype_label_in_plots, fill=Estimate))+
  geom_tile(color="lightgrey", linewidth=0.5, linetype=1)+
  scale_fill_gradient2(low="#0072B2", mid="white", high="maroon",
                       midpoint=0,
                       limits=c(-1.2,1.2))+
  scale_x_discrete(position = "top")+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0, vjust=1.5))+
  labs(x="", y="")+
    geom_text(label=inputdata$my_label)+
    facet_grid(Phenotype_group ~ ., scales = "free_y", space = "free")

  ## save plot
  ggsave(paste0("association_model_results_plots/", outputname), myplot, dpi=500, width=w, height=h, units="cm")
  
}

create.heatmap.ppt <- function(inputdata, outputname, w, h){

  ## create plot
  myplot <- ggplot(inputdata, aes(x=Phenotype_label_in_plots, y=plot_HMO, fill=Estimate))+
    geom_tile()+
    scale_fill_gradient2(low="#0072B2", mid="white", high="maroon",
                         midpoint=0,
                         limits=c(-1.2,1.2))+
    theme_bw()+
    theme(axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))+
    labs(x="", y="")+
    geom_text(label=inputdata$my_label)+
    facet_grid(. ~ Phenotype_group, scales = "free_x", space = "free")

  ## save plot
  ggsave(paste0("association_model_results_plots/", outputname), myplot, dpi=500, width=w, height=h, units="cm")

}


### ===== 2.2 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED ===== ###

create.heatmap(resHMO_for_plotting,
               outputname = "240229_milk_real_HMOs_by_shortlist_phenotypes_n172.pdf",
               w=40, h=125) #activate facet_grid in function for this


### ===== 2.3 CREATE HEATMAPS FOR ALL FDR-SIGNIFICANT ASSOCIATIONS ===== ###

## only including phenotype associations that are FDR significant
create.heatmap.ppt(resHMO_for_plotting[resHMO_for_plotting$FDR<0.05,],
                   outputname = "240229_milk_real_HMOs_by_shortlist_phenotypes_with_FDR_sign_associations_n12_turned.pdf",
                   w=40, h=20)


### ===== 2.4 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED - ONLY INCLUDING PHENOTYPES WITH AT LEAST X NOMINALLY SIGNIFICANT ASSOCIATIONS (FIGURE 3A) ===== ###

resHMO_for_plotting_nom_tmp <- resHMO_for_plotting[,c(10,6)]
resHMO_for_plotting_nom_tmp <- resHMO_for_plotting_nom_tmp[resHMO_for_plotting_nom_tmp$p<0.05,]
column_name <- c()
n_associations <- c()
for (i in unique(resHMO_for_plotting_nom_tmp$Phenotype_label_in_plots)){
  column_name <- c(column_name, i)
  n_associations <- c(n_associations, nrow(resHMO_for_plotting_nom_tmp[resHMO_for_plotting_nom_tmp$Phenotype_label_in_plots==i,]))
}
resHMO_for_plotting_nom_tmp2 <- data.frame(column_name = column_name, n_associations = n_associations)
resHMO_for_plotting_nom_tmp2$FDR_significant_association <- "no"
resHMO_for_plotting_nom_tmp2[resHMO_for_plotting_nom_tmp2$column_name %in% resHMO_for_plotting[resHMO_for_plotting$FDR<0.05,"Phenotype_label_in_plots"],"FDR_significant_association"] <- "yes"
resHMO_for_plotting_nom_tmp2 <- resHMO_for_plotting_nom_tmp2[order(-resHMO_for_plotting_nom_tmp2$n_associations),]
# to automatically capture all FDR-significant associations, include phenotypes with at least 4 p<0.05 associations


## only including phenotype that have at least 4 nominally/FDR significant associations
phenotypes_with_at_least_4_nom_associations <- resHMO_for_plotting_nom_tmp2[resHMO_for_plotting_nom_tmp2$n_associations>=4, "column_name"]
length(phenotypes_with_at_least_4_nom_associations) #40

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_label_in_plots %in% phenotypes_with_at_least_4_nom_associations,],
               outputname = "240229_milk_real_HMOs_by_shortlist_phenotypes_with_at_least_4_nom_sign_associations_n40.pdf",
               w=30, h=30) #activate facet_grid in function for this

create.heatmap.ppt(resHMO_for_plotting[resHMO_for_plotting$Phenotype_label_in_plots %in% phenotypes_with_at_least_4_nom_associations,],
                   outputname = "240229_milk_real_HMOs_by_shortlist_phenotypes_with_at_least_4_nom_sign_associations_n40_turned.pdf",
                   w=50, h=35)


### ===== 2.5 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED - ONLY INCLUDING FDR-SIGNIFICANT ASSOCIATIONS OR PHENOTYPES WITH AT LEAST X NOMINALLY SIGNIFICANT ASSOCIATIONS ===== ###

## only including phenotype that have at least 5 nominally/FDR significant associations
phenotypes_with_at_least_5_nom_associations <- resHMO_for_plotting_nom_tmp2[resHMO_for_plotting_nom_tmp2$n_associations>=5 | resHMO_for_plotting_nom_tmp2$FDR_significant_association=="yes", "column_name"]
length(phenotypes_with_at_least_5_nom_associations) #32

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_label_in_plots %in% phenotypes_with_at_least_5_nom_associations,],
               outputname = "240229_milk_real_HMOs_by_shortlist_phenotypes_with_FDR_or_at_least_5_nom_sign_associations_n32.pdf",
               w=30, h=30) #activate facet_grid in function for this

create.heatmap.ppt(resHMO_for_plotting[resHMO_for_plotting$Phenotype_label_in_plots %in% phenotypes_with_at_least_5_nom_associations,],
                   outputname = "240229_milk_real_HMOs_by_shortlist_phenotypes_with_FDR_or_at_least_5_nom_sign_associations_n32_turned.pdf",
                   w=50, h=35)


# setdiff(phenotypes_with_at_least_4_nom_associations, phenotypes_with_at_least_5_nom_associations)
# [1] "Prolonged rupture of membranes 1352"                         
# [2] "Avoiding gluten during pregnancy 963"                        
# [3] "Alpha-linolenic acid 512"                                    
# [4] "Collection season (autumn vs spring) 1540"                   
# [5] "Oxytocin during birth 408"                                   
# [6] "Gestational age (weeks) 1157"                                
# [7] "Infant birth delivery mode (C-section vs vaginal birth) 1411"
# [8] "Vitamin C 512"  


### ===== 2.6 CREATE HEATMAPS PER PHENOTYPE GROUP ===== ###

# Human_milk_collection               Infant_feeding   Maternal_and_infant_anthropometrics  Maternal_diet                   Maternal_genetics 
# 3                                   4                8                                    72                              1 
# Maternal_health_and_diseases    Maternal_lifestyle_and_exposures    Maternal_medication_use   Pregnancy_and_birth_characteristics 
# 21                              12                                  22                        19

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Human_milk_collection",],
               outputname = "240229_milk_real_HMOs_by_human_milk_collection_phenotypes_n3.pdf",
               w=30, h=15)

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Infant_feeding",],
               outputname = "240229_milk_real_HMOs_by_infant_feeding_phenotypes_n4.pdf",
               w=30, h=15)

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Maternal_and_infant_anthropometrics",],
               outputname = "240229_milk_real_HMOs_by_maternal_and_infant_anthropometrics_phenotypes_n8.pdf",
               w=30, h=15)

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Maternal_diet",],
               outputname = "240229_milk_real_HMOs_by_maternal_diet_phenotypes_n72.pdf",
               w=30, h=50)

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Maternal_genetics",],
               outputname = "240229_milk_real_HMOs_by_maternal_genetics_phenotypes_n1.pdf",
               w=30, h=15)

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Maternal_health_and_diseases",],
               outputname = "240229_milk_real_HMOs_by_maternal_health_and_diseases_phenotypes_n21.pdf",
               w=30, h=20)

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Maternal_lifestyle_and_exposures",],
               outputname = "240229_milk_real_HMOs_by_maternal_lifestyle_and_exposures_phenotypes_n12.pdf",
               w=30, h=15)

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Maternal_medication_use",],
               outputname = "240229_milk_real_HMOs_by_maternal_medication_use_phenotypes_n22.pdf",
               w=30, h=20)

create.heatmap(resHMO_for_plotting[resHMO_for_plotting$Phenotype_group=="Pregnancy_and_birth_characteristics",],
               outputname = "240229_milk_real_HMOs_by_pregnancy_and_birth_characteristics_phenotypes_n19.pdf",
               w=30, h=20)






