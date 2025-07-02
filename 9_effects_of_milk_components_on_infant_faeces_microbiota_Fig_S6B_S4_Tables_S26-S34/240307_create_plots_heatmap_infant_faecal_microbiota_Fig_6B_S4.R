####################################################################################
### PLOT RESULTS OF MILK COMPONENTS ON INFANT FAECAL MICROBIOTA (FIGURES 6B, S4) ###
####################################################################################

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
# 1.2 PREPARE DATA FOR PLOTTING
# 
# 2. PLOT HEATMAPS
# 2.1 FUNCTION TO CREATE HEATMAPS
# 2.2 CREATE HEATMAPS FOR ALL HMO ASSOCIATIONS (EXCL. MILK GROUP) (FIGURE 6B)
# 2.3 CREATE HEATMAPS FOR ALL MILK RELATIVE ABUNDANCES ASSOCIATIONS (EXCL. MILK ALPHA DIVERSITY) (FIGURE S4)
# [2.4 CREATE HEATMAPS FOR FDR-SIGNIFICANT ASSOCIATIONS]


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/")

## import results of maternal milk groups and milk HMOs on infant faecal alpha diversity and relative abundances
inff_alpha_by_HMOs <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/alpha_diversity/240306_infant_faeces_microbiota_richness_shannon_by_milkgroup_and_milk_HMOs_alpha_diversity_lmer_correction_for_DNAmethod_reads_time_ID_delmode_results.txt", header=T, sep="\t")
#62x9

inff_RA_by_HMOs <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/relative_abundances/240306_infant_faeces_microbiota_by_milkgroup_and_milk_HMOs_relative_abundances_lmer_correction_for_DNAmethod_reads_time_ID_delmode_results.txt", header=T, sep="\t")
#930x9


## import results of milk alpha diversity and relative abundances on infant faecal alpha diversity
inff_alpha_by_milk_alpha <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/effects_of_milk_microbiota/240305_infant_faecal_alpha_diversity_by_milk_alpha_diversity_model_results.txt", header=T, sep="\t")
#4x9

inff_alpha_by_milk_RA <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/effects_of_milk_microbiota/240305_infant_faecal_alpha_diversity_by_milk_relative_abundances_model_results.txt", header=T, sep="\t")
#72x9


## import results of milk alpha diversity and relative abundances on infant faecal relative abundances
inff_RA_by_milk_alpha <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/effects_of_milk_microbiota/240305_infant_faecal_relative_abundances_by_milk_alpha_diversity_model_results.txt", header=T, sep="\t")
#60x9

inff_RA_by_milk_RA <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/effects_of_milk_microbiota/240305_infant_faecal_relative_abundances_by_milk_relative_abundances_model_results.txt", header=T, sep="\t")
#1080x9


##### =========================== 1. CHECK AND PREPARE RESULTS FOR PLOTTING HEATMAPS =========================== #####

### ===== 1.1 CHECK RESULTS ===== ###

## results of maternal milk groups and milk HMOs on infant faecal alpha diversity and relative abundances
nrow(inff_alpha_by_HMOs[inff_alpha_by_HMOs$FDR<0.05,]) # 3
nrow(inff_alpha_by_HMOs[inff_alpha_by_HMOs$p<0.05,])   #11

nrow(inff_RA_by_HMOs[inff_RA_by_HMOs$FDR<0.05,]) # 4
nrow(inff_RA_by_HMOs[inff_RA_by_HMOs$p<0.05,])   #88


## results of milk alpha diversity and relative abundances on infant faecal alpha diversity
nrow(inff_alpha_by_milk_alpha[inff_alpha_by_milk_alpha$FDR<0.05,]) #0
nrow(inff_alpha_by_milk_alpha[inff_alpha_by_milk_alpha$p<0.05,])   #1

nrow(inff_alpha_by_milk_RA[inff_alpha_by_milk_RA$FDR<0.05,]) #0
nrow(inff_alpha_by_milk_RA[inff_alpha_by_milk_RA$p<0.05,])   #3


## results of milk alpha diversity and relative abundances on infant faecal relative abundances
nrow(inff_RA_by_milk_alpha[inff_RA_by_milk_alpha$FDR<0.05,]) #0
nrow(inff_RA_by_milk_alpha[inff_RA_by_milk_alpha$p<0.05,])   #3

nrow(inff_RA_by_milk_RA[inff_RA_by_milk_RA$FDR<0.05,]) #  8
nrow(inff_RA_by_milk_RA[inff_RA_by_milk_RA$p<0.05,])   #100


### ===== 1.2 PREPARE DATA FOR PLOTTING ===== ###

## colnames are already the same in each of the 6 data frames

## ensure x and y columns have the same words for the same levels

# inff_alpha_by_HMOs and inff_RA_by_HMOs -> add inff_ for y entries
inff_alpha_by_HMOs$y <- c(paste0("inff_", inff_alpha_by_HMOs$y))
inff_RA_by_HMOs$y <- c(paste0("inff_", inff_RA_by_HMOs$y))

# inff_alpha_by_milk_alpha - already has milk_ and inff_ in beginning for x and y
# inff_alpha_by_milk_RA - already has milk_ and inff_ in beginning for x and y

# inff_RA_by_milk_alpha - already has milk_ and inff_ in beginning for x and y
# inff_RA_by_milk_RA - already has milk_ and inff_ in beginning for x and y


## add column indicating x type (HMO conc, milk alpha div, milk RA)
inff_alpha_by_HMOs$x_type <- "HMO_concentrations"
inff_alpha_by_milk_alpha$x_type <- "Milk_alpha_diversity"
inff_alpha_by_milk_RA$x_type <- "Milk_relative_bacterial_abundance"

inff_RA_by_HMOs$x_type <- "HMO_concentrations"
inff_RA_by_milk_alpha$x_type <- "Milk_alpha_diversity"
inff_RA_by_milk_RA$x_type <- "Milk_relative_bacterial_abundance"


## add column indicating y type (inff alpha div, inff RA)
inff_alpha_by_HMOs$y_type <- "Infant_faecal_alpha_diversity"
inff_alpha_by_milk_alpha$y_type <- "Infant_faecal_alpha_diversity"
inff_alpha_by_milk_RA$y_type <- "Infant_faecal_alpha_diversity"

inff_RA_by_HMOs$y_type <- "Infant_faecal_relative_bacterial_abundance"
inff_RA_by_milk_alpha$y_type <- "Infant_faecal_relative_bacterial_abundance"
inff_RA_by_milk_RA$y_type <- "Infant_faecal_relative_bacterial_abundance"


## combine results for alpha diversity and relative abundances
model_results_for_plotting <- as.data.frame(rbind(inff_alpha_by_HMOs, inff_alpha_by_milk_alpha, inff_alpha_by_milk_RA,
                                                  inff_RA_by_HMOs, inff_RA_by_milk_alpha, inff_RA_by_milk_RA))
#2208x11


## fix x entries
model_results_for_plotting$x <- gsub("mother_milk_HMO_", "", model_results_for_plotting$x)
model_results_for_plotting$x <- gsub("_ugml_invr", "", model_results_for_plotting$x)
model_results_for_plotting$x <- gsub("milk_alpha_div_", "", model_results_for_plotting$x)
model_results_for_plotting$x <- gsub("_invr", "", model_results_for_plotting$x)
model_results_for_plotting$x <- gsub("milk_g__", "", model_results_for_plotting$x)
model_results_for_plotting$x <- gsub("_clr", "", model_results_for_plotting$x)

## fix y entries
model_results_for_plotting$y <- gsub("inff_alpha_div_", "", model_results_for_plotting$y)
model_results_for_plotting$y <- gsub("_invr", "", model_results_for_plotting$y)
model_results_for_plotting$y <- gsub("inff_g__", "", model_results_for_plotting$y)
model_results_for_plotting$y <- gsub("_clr", "", model_results_for_plotting$y)


## set all estimates of associations with p>0.05 to 0 so they show white in the plot
colnames(model_results_for_plotting)[7] <- "Estimate"
model_results_for_plotting[model_results_for_plotting$p>0.05, "Estimate"] <- 0

## determine labels for FDR and nominally significant results in plot
model_results_for_plotting$my_label <- NA
model_results_for_plotting[model_results_for_plotting$p>0.05, "my_label"] <- ""
model_results_for_plotting[model_results_for_plotting$p<0.05 & model_results_for_plotting$Estimate<0, "my_label"] <- "-"
model_results_for_plotting[model_results_for_plotting$p<0.05 & model_results_for_plotting$Estimate>0, "my_label"] <- "+"
model_results_for_plotting[model_results_for_plotting$FDR<0.05 & model_results_for_plotting$Estimate<0, "my_label"] <- "-*"
model_results_for_plotting[model_results_for_plotting$FDR<0.05 & model_results_for_plotting$Estimate>0, "my_label"] <- "+*"
table(model_results_for_plotting$my_label, useNA="ifany")


# reverse the order of the x_type levels so they show them alphabetically from top to bottom in the plot
model_results_for_plotting$x <- as.factor(as.character(model_results_for_plotting$x))

# library(forcats)
levels(model_results_for_plotting$x)
model_results_for_plotting$x <- factor(model_results_for_plotting$x,
                                       levels = rev(levels(model_results_for_plotting$x)))
levels(model_results_for_plotting$x)


##### =========================== 2. PLOT HEATMAPS =========================== #####

setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/7_infant_faeces_microbiota/")

library(ggplot2)


### ===== 2.1 FUNCTION TO CREATE HEATMAPS ===== ###

create.heatmap <- function(inputdata, mylimits, mylabels, outputname, w, h){
  
  ## create plot
  myplot <- ggplot(inputdata, aes(x=y, y=x, fill=Estimate))+
  geom_tile(color="lightgrey", linewidth=0.5, linetype=1)+
  scale_fill_gradient2(low="#0072B2", mid="white", high="maroon",
                       midpoint=0,
                       limits=mylimits, breaks=mylabels, labels=mylabels)+ #c(-0.5,0.5) for alpha div, c(-3,3) for RA
  scale_x_discrete(position = "top")+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=0, vjust=1.5))+
  labs(x="", y="")+
    geom_text(label=inputdata$my_label)+
    facet_grid(x_type ~ ., scales = "free_y", space = "free")

  ## save plot
  ggsave(paste0("heatmap_plots/", outputname), myplot, dpi=500, width=w, height=h, units="cm")
  
  ## return plot
  return(myplot)
  
}


### ===== 2.2 CREATE HEATMAPS FOR ALL HMO ASSOCIATIONS (EXCL. MILK GROUP) (FIGURE 6B) ===== ###

## select results for HMO concentrations only
HMO_data <- model_results_for_plotting[model_results_for_plotting$x_type=="HMO_concentrations" & model_results_for_plotting$x!="milk_group",]

## fix order of phenotypes in plot
HMO_data$x <- as.factor(as.character(HMO_data$x))
levels(HMO_data$x)
HMO_data$x <- factor(HMO_data$x,
                     levels = c("LSTc","LSTb","DSLNT","6SL","3SL",
                                "3F3SL",
                                "MFLNH_III","LNnFP_V","LNnDFH","LNFP_V","LNFP_III",
                                "LNFP_II","LNFP_I","LNDH_I","LDFT","DFLNHa",
                                "A_tetra","3FL","2FL",
                                "LNT","LNnT","LNH","6GL","3GL",
                                "Sia","Fuc","Neut","Total"))
levels(HMO_data$x) <- c("LSTc",     "LSTb",     "DSLNT",    "6'SL",      "3'SL",      "3'F3'SL",   
                        "MFLNH III","LNnFP V",  "LNnDFH",   "LNFP V",   "LNFP III", "LNFP II", 
                        "LNFP I",   "LNDH I",   "LDFT",     "DFLNHa",   "A-tetra",  "3'FL",     
                        "2'FL",      "LNT",      "LNnT",     "LNH",      "6'GL",      "3'GL",     
                        "Sialylated HMOs",      "Fucosylated HMOs",      "Neutral HMOs",     "Total HMOs")

plot_alpha_by_HMO <- create.heatmap(HMO_data[HMO_data$y_type=="Infant_faecal_alpha_diversity",],
                             mylimits = c(-0.2,0.2),
                             mylabels = c(-0.2,0,0.2),
                             outputname = "240308_infant_faecal_alpha_div_by_HMOs_n28.pdf",
                             w=10, h=30) #activate facet_grid in function for this

plot_RA_by_HMO <- create.heatmap(HMO_data[HMO_data$y_type=="Infant_faecal_relative_bacterial_abundance",],
                                 mylimits = c(-0.6,0.6),
                                 mylabels = c(-0.6,0,0.6),
                                 outputname = "240308_infant_faecal_relative_abundances_by_HMOs_n28.pdf",
                                 w=30, h=30) #activate facet_grid in function for this


### === COMBINE PLOTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES (FIGURE 6B) === ###

library(cowplot)


### COMBINED PLOTS WITH X AXIS LABELS

## fix plots before combining
plot_alpha_by_HMO_v2 <- plot_alpha_by_HMO + theme(strip.text.y = element_blank())

plot_RA_by_HMO_v2 <- plot_RA_by_HMO + theme(axis.text.y = element_blank(),
                              strip.text.y = element_blank())

## combine plots
comb_plots_by_HMO <- plot_grid(plot_alpha_by_HMO_v2, plot_RA_by_HMO_v2, rel_widths = c(1,2), nrow=1)

## save combined plots
ggsave("heatmap_plots/240308_infant_faecal_alpha_div_and_RA_by_HMOs_n28.pdf", comb_plots_by_HMO, dpi=500, width=40, height=60, units="cm")


### COMBINED AND ALIGNED PLOTS WITHOUT X AXIS LABELS

## fix plots before combining
plot_alpha_by_HMO_v3 <- plot_alpha_by_HMO + theme(axis.text.x = element_blank(),
                                   strip.text.y = element_blank())

plot_RA_by_HMO_v3 <- plot_RA_by_HMO + theme(axis.text.x = element_blank(),
                              axis.text.y = element_blank(),
                              strip.text.y = element_blank())

## combine plots
plot_RA_by_HMO_aligned <- plot_grid(plot_alpha_by_HMO_v3, plot_RA_by_HMO_v3, rel_widths = c(0.69,3), nrow=1)

## save combined plots
ggsave("heatmap_plots/240308_infant_faecal_alpha_div_and_RA_by_HMOs_n28_aligned.pdf", plot_RA_by_HMO_aligned, dpi=500, width=45, height=30, units="cm")


### ===== 2.3 CREATE HEATMAPS FOR ALL MILK RELATIVE ABUNDANCES ASSOCIATIONS (EXCL. MILK ALPHA DIVERSITY) (FIGURE S4) ===== ###

RA_data <- model_results_for_plotting[model_results_for_plotting$x_type=="Milk_relative_bacterial_abundance",]
# plot_alpha_by_RA <- create.heatmap(RA_data[RA_data$y_type=="Infant_faecal_alpha_diversity",],
#                                     mylimits = c(-0.1,0.1),
#                                     mylabels = c(-0.1,0,0.1),
#                                     outputname = "240308_infant_faecal_alpha_div_by_milk_relative_abundances_n36.pdf",
#                                     w=20, h=30) #activate facet_grid in function for this

## FIGURE S4
plot_RA_by_RA <- create.heatmap(RA_data[RA_data$y_type=="Infant_faecal_relative_bacterial_abundance",],
                                 mylimits = c(-0.4,0.4),
                                 mylabels = c(-0.4,0,0.4),
                                 outputname = "240308_infant_faecal_relative_abundances_by_milk_relative_abundances_n36.pdf",
                                 w=33, h=30) #activate facet_grid in function for this


# ### ===== 2.4 CREATE HEATMAPS FOR FDR-SIGNIFICANT ASSOCIATIONS ===== ###
# 
# my_x_phenos <- unique(model_results_for_plotting[model_results_for_plotting$FDR<0.05,"x"])
# FDR_results_for_plotting <- model_results_for_plotting[model_results_for_plotting$x %in% my_x_phenos,]
# 
# my_y_phenos <- unique(FDR_results_for_plotting[FDR_results_for_plotting$p<0.05,"y"])
# FDR_results_for_plotting2 <- FDR_results_for_plotting[FDR_results_for_plotting$y %in% my_y_phenos,]
# 
# ## fix order of phenotypes in plot
# FDR_results_for_plotting2$x <- as.factor(as.character(FDR_results_for_plotting2$x))
# FDR_results_for_plotting2$x <- factor(FDR_results_for_plotting2$x,
#                                      levels = c("DSLNT","MFLNH_III","LNnFP_V","LNH","3GL",
#                                                 "Neisseria","Lactobacillus","Haemophilus",
#                                                 "Gemella","Bifidobacterium","Acinetobacter"))
# 
# plot_alpha <- create.heatmap(FDR_results_for_plotting2[FDR_results_for_plotting2$y_type=="Infant_faecal_alpha_diversity",],
#                              mylimits = c(-0.2,0.2),
#                              mylabels = c(-0.2,0,0.2),
#                              outputname = "240307_infant_faecal_alpha_div_by_FDR_sign_phenos_n12.pdf",
#                              w=10, h=30) #activate facet_grid in function for this
# 
# plot_RA <- create.heatmap(FDR_results_for_plotting2[FDR_results_for_plotting2$y_type=="Infant_faecal_relative_bacterial_abundance",],
#                           mylimits = c(-0.6,0.6),
#                           mylabels = c(-0.6,0,0.6),
#                           outputname = "240307_infant_faecal_RA_by_FDR_sign_phenos_n12.pdf",
#                           w=30, h=30) #activate facet_grid in function for this
# 
# 
# ### === COMBINE PLOTS FOR ALPHA DIVERSITY AND RELATIVE ABUNDANCES === ###
# 
# library(cowplot)
# 
# 
# ### COMBINED PLOTS WITH X AXIS LABELS
# 
# ## fix plots before combining
# plot_alpha_v2 <- plot_alpha + theme(strip.text.y = element_blank())
# 
# plot_RA_v2 <- plot_RA + theme(axis.text.y = element_blank(),
#                               strip.text.y = element_blank())
# 
# ## combine plots
# comb_plots <- plot_grid(plot_alpha_v2, plot_RA_v2, rel_widths = c(1,2), nrow=1)
# 
# ## save combined plots
# ggsave("heatmap_plots/240307_infant_faecal_alpha_div_and_RA_by_FDR_sign_phenos_n12.pdf", comb_plots, dpi=500, width=40, height=60, units="cm")
# 
# 
# ### COMBINED AND ALIGNED PLOTS WITHOUT X AXIS LABELS
# 
# ## fix plots before combining
# plot_alpha_v3 <- plot_alpha +theme(axis.text.x = element_blank(),
#                                               strip.text.y = element_blank())
# 
# plot_RA_v3 <- plot_RA + theme(axis.text.x = element_blank(),
#                               axis.text.y = element_blank(),
#                               strip.text.y = element_blank())
# 
# ## combine plots
# plot_RA_aligned <- plot_grid(plot_alpha_v3, plot_RA_v3, rel_widths = c(0.65,3), nrow=1)
# 
# ## save combined plots
# ggsave("heatmap_plots/240307_infant_faecal_alpha_div_and_RA_by_FDR_sign_phenos_n12_aligned.pdf", plot_RA_aligned, dpi=500, width=56, height=24, units="cm")


