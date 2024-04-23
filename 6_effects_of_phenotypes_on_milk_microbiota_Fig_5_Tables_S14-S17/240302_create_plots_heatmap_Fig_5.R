#######################################################################################################
### PLOT RESULTS OF PHENOTYPES ON MILK MICROBIOTA COMPOSITION FOR MILK COMPOSITION PAPER (FIGURE 5) ###
#######################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 0. DATA IMPORT
# 
# 1. CHECK AND PREPARE RESULTS FOR PLOTTING HEATMAPS
# 1.1 CHECK RESULTS
# 1.2 CHECK N PER PHENOTYPE GROUP
# 1.3 PREPARE DATA FOR PLOTTING
# 
# 2. PLOT HEATMAPS
# 2.1 FUNCTION TO CREATE HEATMAPS
# 2.2 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED
# 2.3 CREATE HEATMAPS FOR FDR-SIGNIFICANT ASSOCIATIONS
# 2.4 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED - ONLY INCLUDING PHENOTYPES WITH AT LEAST X NOMINALLY SIGNIFICANT ASSOCIATIONS (PREPARE PANELS FOR FIGURE 5)
# 2.5 COMBINE PLOTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES (FIGURE 5)
# [2.6 CREATE HEATMAPS PER PHENOTYPE GROUP]


##### =========================== 0. DATA IMPORT =========================== #####

## cluster
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/association_model_results/")

model_results_alpha <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/association_model_results/240302_polished_model_results_milk_alpha_diversity_by_phenotypes_n185.txt", header=T)
model_results_RA <- read.table(file="/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/association_model_results/240301_polished_model_results_milk_relative_abundances_by_phenotypes_n185.txt", header=T)


# ## local
# setwd("/Users/johanne/Desktop/PhD-Groningen/Thesis_chapters/Chapter3_milk_composition/Results/6_milk_microbiota/")
# 
# model_results_alpha <- read.table(file="association_model_results/240302_polished_model_results_milk_alpha_diversity_by_phenotypes_n185.txt", header=T)
# model_results_RA <- read.table(file="association_model_results/240301_polished_model_results_milk_relative_abundances_by_phenotypes_n185.txt", header=T)


##### =========================== 1. CHECK AND PREPARE RESULTS FOR PLOTTING HEATMAPS =========================== #####

### ===== 1.1 CHECK RESULTS ===== ###

nrow(model_results_alpha[model_results_alpha$FDR<0.05,]) # 0
nrow(model_results_alpha[model_results_alpha$p<0.05,])   #30

nrow(model_results_RA[model_results_RA$FDR<0.05,]) #  8
nrow(model_results_RA[model_results_RA$p<0.05,])   #613


### ===== 1.2 CHECK N PER PHENOTYPE GROUP ===== ###

## check unique phenotypes
# unique(model_results_alpha$Phenotype)[order(unique(model_results_alpha$Phenotype))] #n=185
# unique(model_results_RA$Phenotype)[order(unique(model_results_RA$Phenotype))] #n=185
# length(intersect(unique(model_results_alpha$Phenotype), unique(model_results_RA$Phenotype))) #n=185

## check n per phenotype group
phenos_info_alpha <- model_results_alpha[,c(1:2)]
phenos_info_alpha <- phenos_info_alpha[!duplicated(phenos_info_alpha),]
table(phenos_info_alpha$Phenotype_group, useNA="ifany")

phenos_info_RA <- model_results_RA[,c(1:2)]
phenos_info_RA <- phenos_info_RA[!duplicated(phenos_info_RA),]
table(phenos_info_RA$Phenotype_group, useNA="ifany")

# Human_milk_collection           Human_milk_HMO_concentrations 
# 3                               28 
# Infant_feeding                  Maternal_and_infant_anthropometrics 
# 4                               8 
# Maternal_diet                   Maternal_genetics 
# 72                              2 
# Maternal_health_and_diseases    Maternal_lifestyle_and_exposures 
# 21                              11 
# Maternal_medication_use         Pregnancy_and_birth_characteristics 
# 17                              19


### ===== 1.3 PREPARE DATA FOR PLOTTING ===== ###

## ensure outcome column has the same colname in both data frames
colnames(model_results_alpha)[7] <- "y"
colnames(model_results_RA)[7] <- "y"

## combine results for alpha diversity and relative abundances
model_results_for_plotting <- as.data.frame(rbind(model_results_alpha, model_results_RA))

## fix outcome entries for alpha diversity
model_results_for_plotting$y <- as.factor(as.character(model_results_for_plotting$y))
levels(model_results_for_plotting$y)[21] <- "Genera richness"
levels(model_results_for_plotting$y)[33] <- "Shannon diversity index"

## set all estimates of associations with p>0.05 to 0 so they show white in the plot
model_results_for_plotting[model_results_for_plotting$p>0.05, "Estimate"] <- 0

## determine labels for FDR and nominally significant results in plot
model_results_for_plotting$my_label <- NA
model_results_for_plotting[model_results_for_plotting$p>0.05, "my_label"] <- ""
model_results_for_plotting[model_results_for_plotting$p<0.05 & model_results_for_plotting$Estimate<0, "my_label"] <- "-"
model_results_for_plotting[model_results_for_plotting$p<0.05 & model_results_for_plotting$Estimate>0, "my_label"] <- "+"
model_results_for_plotting[model_results_for_plotting$FDR<0.05 & model_results_for_plotting$Estimate<0, "my_label"] <- "-*"
model_results_for_plotting[model_results_for_plotting$FDR<0.05 & model_results_for_plotting$Estimate>0, "my_label"] <- "+*"
table(model_results_for_plotting$my_label, useNA="ifany")

## fix labels from trait_name_in_plot: replace _ with space
model_results_for_plotting$Phenotype_label_in_plots <- gsub("_", " ", model_results_for_plotting$Phenotype_label_in_plots)

## add n to trait label
model_results_for_plotting$Phenotype_label_in_plots <- c(paste0(model_results_for_plotting$Phenotype_label_in_plots, " ", model_results_for_plotting$n))

## reverse the order of the levels so they show them alphabetically from top to bottom in the plot
model_results_for_plotting$Phenotype_label_in_plots <- as.factor(as.character(model_results_for_plotting$Phenotype_label_in_plots))

# library(forcats)
levels(model_results_for_plotting$Phenotype_label_in_plots)
model_results_for_plotting$Phenotype_label_in_plots <- factor(model_results_for_plotting$Phenotype_label_in_plots,
                                                              levels = rev(levels(model_results_for_plotting$Phenotype_label_in_plots)))
levels(model_results_for_plotting$Phenotype_label_in_plots)

## fix label for DFLNHa
levels(model_results_for_plotting$Phenotype_label_in_plots)[151] <- "DFLNHa (µg/ml) 725"

## fix label for fucosylated HMOs
levels(model_results_for_plotting$Phenotype_label_in_plots)[129] <- "Fucosylated HMOs (µg/ml) 725"

## add column for outcome type
model_results_for_plotting$outcome_type <- NA
model_results_for_plotting[model_results_for_plotting$y=="Genera richness" | model_results_for_plotting$y=="Shannon diversity index", "outcome_type"] <- "Alpha diversity"
model_results_for_plotting[is.na(model_results_for_plotting$outcome_type), "outcome_type"] <- "Relative bacterial abundance"

## fix labels for y -> replace _ with space
levels(model_results_for_plotting$y)[c(4,17,30)] <- c("Allo-/Neo-/Para-/Rhizobium", "Escherichia-Shigella", "Prevotella 7")


##### =========================== 2. PLOT HEATMAPS =========================== #####

setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/6_milk_microbiota/")

library(ggplot2)


### ===== 2.1 FUNCTION TO CREATE HEATMAPS ===== ###

create.heatmap <- function(inputdata, mylimits, outputname, w, h){
  
  ## create plot
  myplot <- ggplot(inputdata, aes(x=y, y=Phenotype_label_in_plots, fill=Estimate))+
  geom_tile(color="lightgrey", linewidth=0.5, linetype=1)+
  scale_fill_gradient2(low="#0072B2", mid="white", high="maroon",
                       midpoint=0,
                       limits=mylimits)+ #c(-1,1) for alpha div, c(-5,5) for RA
  scale_x_discrete(position = "top")+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0, vjust=1.5))+
  labs(x="", y="")+
    geom_text(label=inputdata$my_label)+
    facet_grid(Phenotype_group ~ ., scales = "free_y", space = "free")

  ## save plot
  ggsave(paste0("association_model_results_plots/", outputname), myplot, dpi=500, width=w, height=h, units="cm")
  
  ## return plot
  return(myplot)
  
}


### ===== 2.2 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED ===== ###

## heatmap for alpha diversity measures
create.heatmap(model_results_for_plotting[model_results_for_plotting$outcome_type=="Alpha diversity",],
               mylimits = c(-1,1),
               outputname = "240302_milk_alpha_diversity_by_shortlist_phenotypes_n185.pdf",
               w=40, h=125) #activate facet_grid in function for this

## heatmap for relative abundances
create.heatmap(model_results_for_plotting[model_results_for_plotting$outcome_type=="Relative bacterial abundance",],
               mylimits = c(-5,5),
               outputname = "240302_milk_relative_abundances_by_shortlist_phenotypes_n185.pdf",
               w=40, h=125) #activate facet_grid in function for this


### ===== 2.3 CREATE HEATMAPS FOR FDR-SIGNIFICANT ASSOCIATIONS ===== ###

## heatmap for alpha diversity measures
# No FDR sign. results.

## heatmap for relative abundances
create.heatmap(model_results_for_plotting[model_results_for_plotting$outcome_type=="Relative bacterial abundance" & model_results_for_plotting$FDR<0.05,],
               mylimits = c(-5,5),
               outputname = "240302_milk_relative_abundances_by_shortlist_phenotypes_with_FDR_sign_associations_n8.pdf",
               w=40, h=15) #activate facet_grid in function for this


### ===== 2.4 CREATE HEATMAPS FOR ALL PHENOTYPE GROUPS COMBINED - ONLY INCLUDING PHENOTYPES WITH AT LEAST X NOMINALLY SIGNIFICANT ASSOCIATIONS (PREPARE PANELS FOR FIGURE 5) ===== ###

model_results_for_plotting_nom_tmp <- model_results_for_plotting[,c(10,6)]
model_results_for_plotting_nom_tmp <- model_results_for_plotting_nom_tmp[model_results_for_plotting_nom_tmp$p<0.05,]
column_name <- c()
n_associations <- c()
for (i in unique(model_results_for_plotting_nom_tmp$Phenotype_label_in_plots)){
  column_name <- c(column_name, i)
  n_associations <- c(n_associations, nrow(model_results_for_plotting_nom_tmp[model_results_for_plotting_nom_tmp$Phenotype_label_in_plots==i,]))
}
model_results_for_plotting_nom_tmp2 <- data.frame(column_name = column_name, n_associations = n_associations)
model_results_for_plotting_nom_tmp2$FDR_significant_association <- "no"
model_results_for_plotting_nom_tmp2[model_results_for_plotting_nom_tmp2$column_name %in% model_results_for_plotting[model_results_for_plotting$FDR<0.05,"Phenotype_label_in_plots"],"FDR_significant_association"] <- "yes"
model_results_for_plotting_nom_tmp2 <- model_results_for_plotting_nom_tmp2[order(-model_results_for_plotting_nom_tmp2$n_associations),]
# to automatically capture all FDR-significant associations, include phenotypes with at least 6 p<0.05 associations

table(model_results_for_plotting_nom_tmp2$n_associations)
# 1  2  3  4  5  6  7  8  9 10 11 
# 26 42 29 25 19 19  7  5  3  1  1 


## only including phenotype that have at least 6 nominally/FDR significant associations
phenotypes_with_at_least_6_nom_associations <- model_results_for_plotting_nom_tmp2[model_results_for_plotting_nom_tmp2$n_associations>=6, "column_name"]
length(phenotypes_with_at_least_6_nom_associations) #36

data_6_nom <- model_results_for_plotting[model_results_for_plotting$Phenotype_label_in_plots %in% phenotypes_with_at_least_6_nom_associations,]

plot_alpha_nom6 <- create.heatmap(data_6_nom[data_6_nom$outcome_type=="Alpha diversity",],
               mylimits = c(-0.3,0.3),
               outputname = "240310_milk_alpha_diversity_by_shortlist_phenotypes_with_at_least_6_nom_sign_associations_n36.pdf",
               w=10, h=30) #activate facet_grid in function for this

plot_RA_nom6 <- create.heatmap(data_6_nom[data_6_nom$outcome_type=="Relative bacterial abundance",],
               mylimits = c(-5,5),
               outputname = "240310_milk_relative_abundances_by_shortlist_phenotypes_with_at_least_6_nom_sign_associations_n36.pdf",
               w=30, h=30) #activate facet_grid in function for this


### ===== 2.5 COMBINE PLOTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES (FIGURE 5) ===== ###

library(cowplot)


### COMBINED PLOTS WITH X AXIS LABELS

## fix plots before combining
plot_alpha_nom6_v2 <- plot_alpha_nom6 + theme(strip.text.y = element_blank())

plot_RA_nom6_v2 <- plot_RA_nom6 + theme(axis.text.y = element_blank(),
                                        strip.text.y = element_blank())

## combine plots
plots_nom6 <- plot_grid(plot_alpha_nom6_v2, plot_RA_nom6_v2, rel_widths = c(1,2), nrow=1)

## save combined plots
ggsave("association_model_results_plots/240310_milk_alpha_div_and_relative_abundances_by_shortlist_phenotypes_with_at_least_6_nom_sign_associations_n36.pdf", plots_nom6, dpi=500, width=40, height=60, units="cm")


### COMBINED AND ALIGNED PLOTS WITHOUT X AXIS LABELS

## fix plots before combining
plot_alpha_nom6_v3 <- plot_alpha_nom6 +theme(axis.text.x = element_blank(),
                                              strip.text.y = element_blank())

plot_RA_nom6_v3 <- plot_RA_nom6 + theme(axis.text.x = element_blank(),
                                        axis.text.y = element_blank(),
                                        strip.text.y = element_blank())

## combine plots
plots_nom6_aligned <- plot_grid(plot_alpha_nom6_v3, plot_RA_nom6_v3, rel_widths = c(0.97,3), nrow=1)

## save combined plots
ggsave("association_model_results_plots/240310_milk_alpha_div_and_relative_abundances_by_shortlist_phenotypes_with_at_least_6_nom_sign_associations_n36_aligned.pdf", plots_nom6_aligned, dpi=500, width=56, height=40, units="cm")



# ### ===== 2.6 CREATE HEATMAPS PER PHENOTYPE GROUP ===== ###
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Human_milk_collection",],
#                outputname = "240220_milk_relative_abundances_by_human_milk_collection_phenotypes_n4.pdf",
#                w=30, h=15)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Human_milk_HMO_concentrations",],
#                outputname = "240220_milk_relative_abundances_by_human_milk_HMO_concentrations_phenotypes_n28.pdf",
#                w=30, h=30)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Infant_feeding",],
#                outputname = "240220_milk_relative_abundances_by_infant_feeding_phenotypes_n4.pdf",
#                w=30, h=15)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Maternal_and_infant_anthropometrics",],
#                outputname = "240220_milk_relative_abundances_by_maternal_and_infant_anthropometrics_phenotypes_n8.pdf",
#                w=30, h=15)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Maternal_diet",],
#                outputname = "240220_milk_relative_abundances_by_maternal_diet_phenotypes_n72.pdf",
#                w=30, h=50)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Maternal_genetics",],
#                outputname = "240220_milk_relative_abundances_by_maternal_genetics_phenotypes_n2.pdf",
#                w=30, h=15)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Maternal_health_and_diseases",],
#                outputname = "240220_milk_relative_abundances_by_maternal_health_and_diseases_phenotypes_n21.pdf",
#                w=30, h=20)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Maternal_lifestyle_and_exposures",],
#                outputname = "240220_milk_relative_abundances_by_maternal_lifestyle_and_exposures_phenotypes_n11.pdf",
#                w=30, h=15)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Maternal_medication_use",],
#                outputname = "240220_milk_relative_abundances_by_maternal_medication_use_phenotypes_n17.pdf",
#                w=30, h=20)
# 
# create.heatmap(model_results_RA_for_plotting[model_results_RA_for_plotting$phenotype_group=="Pregnancy_and_birth_characteristics",],
#                outputname = "240220_milk_relative_abundances_by_pregnancy_and_birth_characteristics_phenotypes_n19.pdf",
#                w=30, h=20)






