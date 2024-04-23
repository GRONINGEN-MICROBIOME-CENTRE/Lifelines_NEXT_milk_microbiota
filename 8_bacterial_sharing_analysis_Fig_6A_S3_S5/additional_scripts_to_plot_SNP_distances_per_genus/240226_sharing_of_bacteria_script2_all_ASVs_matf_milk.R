################################################################################################################
### PLOT SNP DISTANCES FOR ASVs OF EACH SELECTED GENUS FOR MATERNAL FAECES - MILK - SCRIPT 2 - all ##############
################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y

### CONTENT OF THIS FILE
# 1. LOAD DATA AND CREATE DISTANCE MATRICES
# 2. EXTRACT DISTANCE MATRICES FOR EACH GENUS
# 3. PER GENUS, SAVE ALL SNP DISTANCES IN A VECTOR
# 4. CREATE DATA FRAMES WITH GENUS AND SNP DISTANCES
# 5. COMBINE ALL GENERA DATA FRAMES INTO 1 DATA FRAME
# 6. PREPARE DATA FOR PLOTTING AND PLOT SNP DISTANCES BY GENUS
# 7. EXTRACT NUMBER OF UNIQUE ASVs PER GENUS
# 8. CREATE DATA FRAMES WITH GENUS AND NUMBER OF UNIQUE ASVs
# 9. PREPARE DATA FOR PLOTTING AND PLOT NUMBER OF UNIQUE ASVs BY GENUS
# 10. CREATE COMBINED PLOT
# 11. CREATE COMBINED PLOT - TURNED


### ===== 1. LOAD DATA AND CREATE DISTANCE MATRICES ===== ###

# Import: Distance matrices created by **snp-dists** using the MSAs generated before
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/allfa/out")

# Load and List the SNP distance matrices
filenames_all <- list.files(pattern="*.tab", recursive = TRUE)

# Create list of data frame names excluding extension
names_shared_genera_all <- sub("_MA.*", "", filenames_all)
names_all <- paste(names_shared_genera_all , "all", sep = "_")

# Load all files (add 1st column as rownames)
SNP_distances_all <- lapply(filenames_all, function(file) read.delim(file, check.names = FALSE))
names(SNP_distances_all) <- names_all
for (i in 1:length(SNP_distances_all)){
  rownames(SNP_distances_all[[i]]) <- SNP_distances_all[[i]][,1]
  SNP_distances_all[[i]][,1] <- NULL
}

# Get half of the distance matrices (avoid duplicated values)
for (i in 1:length(SNP_distances_all)) {
  SNP_distances_all[[i]][upper.tri(SNP_distances_all[[i]])] <- NA
}


### ===== 2. EXTRACT DISTANCE MATRICES FOR EACH GENUS ===== ###

## extract distance matrix as data frame per genus
Bacteroides_distances <- SNP_distances_all$Bacteroides_ASVs_all
Bifidobacterium_distances <- SNP_distances_all$Bifidobacterium_ASVs_all
Blautia_distances <- SNP_distances_all$Blautia_ASVs_all
Clostridium_sensu_stricto_distances <- SNP_distances_all$Clostridium_sensu_stricto_ASVs_all
Escherichia_Shigella_distances <- SNP_distances_all$Escherichia_Shigella_ASVs_all
Faecalibacterium_distances <- SNP_distances_all$Faecalibacterium_ASVs_all
Haemophilus_distances <- SNP_distances_all$Haemophilus_ASVs_all
Prevotella_distances <- SNP_distances_all$Prevotella_ASVs_all
Streptococcus_distances <- SNP_distances_all$Streptococcus_ASVs_all
Veillonella_distances <- SNP_distances_all$Veillonella_ASVs_all


## set comparison of identical ASVs from 0 to NA
remove.identical.comparisons <- function(dist_matrix){
  for (i in 1:nrow(dist_matrix)){
    for (j in 1:ncol(dist_matrix)){
      
      if(i == j){
        print("Set comparison of identical ASVs to 0")
        dist_matrix[i,j] <- NA
      }else{
        print("Nothing was changed")
      }
    }
  }
  return(dist_matrix)
}

Bacteroides_distances_2 <- remove.identical.comparisons(Bacteroides_distances)
Bifidobacterium_distances_2 <- remove.identical.comparisons(Bifidobacterium_distances)
Blautia_distances_2 <- remove.identical.comparisons(Blautia_distances)
Clostridium_sensu_stricto_distances_2 <- remove.identical.comparisons(Clostridium_sensu_stricto_distances)
Escherichia_Shigella_distances_2 <- remove.identical.comparisons(Escherichia_Shigella_distances)
Faecalibacterium_distances_2 <- remove.identical.comparisons(Faecalibacterium_distances)
Haemophilus_distances_2 <- remove.identical.comparisons(Haemophilus_distances)
Prevotella_distances_2 <- remove.identical.comparisons(Prevotella_distances)
Streptococcus_distances_2 <- remove.identical.comparisons(Streptococcus_distances)
Veillonella_distances_2 <- remove.identical.comparisons(Veillonella_distances)


### ===== 3. PER GENUS, SAVE ALL SNP DISTANCES IN A VECTOR ===== ###

## save all SNP distances in vectors per genus, exclude NAs
Bacteroides_distances_vector <- unlist(Bacteroides_distances_2, use.names=F)
Bacteroides_distances_vector <- Bacteroides_distances_vector[!is.na(Bacteroides_distances_vector)]

Bifidobacterium_distances_vector <- unlist(Bifidobacterium_distances_2, use.names=F)
Bifidobacterium_distances_vector <- Bifidobacterium_distances_vector[!is.na(Bifidobacterium_distances_vector)]

Blautia_distances_vector <- unlist(Blautia_distances_2, use.names=F)
Blautia_distances_vector <- Blautia_distances_vector[!is.na(Blautia_distances_vector)]

Clostridium_sensu_stricto_distances_vector <- unlist(Clostridium_sensu_stricto_distances_2, use.names=F)
Clostridium_sensu_stricto_distances_vector <- Clostridium_sensu_stricto_distances_vector[!is.na(Clostridium_sensu_stricto_distances_vector)]

Escherichia_Shigella_distances_vector <- unlist(Escherichia_Shigella_distances_2, use.names=F)
Escherichia_Shigella_distances_vector <- Escherichia_Shigella_distances_vector[!is.na(Escherichia_Shigella_distances_vector)]

Faecalibacterium_distances_vector <- unlist(Faecalibacterium_distances_2, use.names=F)
Faecalibacterium_distances_vector <- Faecalibacterium_distances_vector[!is.na(Faecalibacterium_distances_vector)]

Haemophilus_distances_vector <- unlist(Haemophilus_distances_2, use.names=F)
Haemophilus_distances_vector <- Haemophilus_distances_vector[!is.na(Haemophilus_distances_vector)]

Prevotella_distances_vector <- unlist(Prevotella_distances_2, use.names=F)
Prevotella_distances_vector <- Prevotella_distances_vector[!is.na(Prevotella_distances_vector)]

Streptococcus_distances_vector <- unlist(Streptococcus_distances_2, use.names=F)
Streptococcus_distances_vector <- Streptococcus_distances_vector[!is.na(Streptococcus_distances_vector)]

Veillonella_distances_vector <- unlist(Veillonella_distances_2, use.names=F)
Veillonella_distances_vector <- Veillonella_distances_vector[!is.na(Veillonella_distances_vector)]


### ===== 4. CREATE DATA FRAMES WITH GENUS AND SNP DISTANCES ===== ###

## prepare data frames from vectors, per genus
Bacteroides_rows <- data.frame(genus=c(rep("Bacteroides", length(Bacteroides_distances_vector))),
                                        SNP_distance=Bacteroides_distances_vector)
Bifidobacterium_rows <- data.frame(genus=c(rep("Bifidobacterium", length(Bifidobacterium_distances_vector))),
                                        SNP_distance=Bifidobacterium_distances_vector)
Blautia_rows <- data.frame(genus=c(rep("Blautia", length(Blautia_distances_vector))),
                                        SNP_distance=Blautia_distances_vector)
Clostridium_sensu_stricto_rows <- data.frame(genus=c(rep("Clostridium_sensu_stricto", length(Clostridium_sensu_stricto_distances_vector))),
                                        SNP_distance=Clostridium_sensu_stricto_distances_vector)
Escherichia_Shigella_rows <- data.frame(genus=c(rep("Escherichia_Shigella", length(Escherichia_Shigella_distances_vector))),
                                        SNP_distance=Escherichia_Shigella_distances_vector)
Faecalibacterium_rows <- data.frame(genus=c(rep("Faecalibacterium", length(Faecalibacterium_distances_vector))),
                                        SNP_distance=Faecalibacterium_distances_vector)
Haemophilus_rows <- data.frame(genus=c(rep("Haemophilus", length(Haemophilus_distances_vector))),
                                        SNP_distance=Haemophilus_distances_vector)
Prevotella_rows <- data.frame(genus=c(rep("Prevotella", length(Prevotella_distances_vector))),
                                        SNP_distance=Prevotella_distances_vector)
Streptococcus_rows <- data.frame(genus=c(rep("Streptococcus", length(Streptococcus_distances_vector))),
                                        SNP_distance=Streptococcus_distances_vector)
Veillonella_rows <- data.frame(genus=c(rep("Veillonella", length(Veillonella_distances_vector))),
                                        SNP_distance=Veillonella_distances_vector)


### ===== 5. COMBINE ALL GENERA DATA FRAMES INTO 1 DATA FRAME ===== ###

## combine all per genus data frames into 1 data frame
df_SNP_dist_all <- as.data.frame(rbind(Bacteroides_rows, Bifidobacterium_rows, Blautia_rows, Clostridium_sensu_stricto_rows, Escherichia_Shigella_rows,
                                     Faecalibacterium_rows, Haemophilus_rows, Prevotella_rows, Streptococcus_rows, Veillonella_rows))


### ===== 6. PREPARE DATA FOR PLOTTING AND PLOT SNP DISTANCES BY GENUS ===== ###

library(ggplot2)

## check medians, sort from large to small manually
median(Prevotella_distances_vector) #54
median(Veillonella_distances_vector) #43
median(Bacteroides_distances_vector) #36
median(Clostridium_sensu_stricto_distances_vector) #19
median(Haemophilus_distances_vector) #16
median(Bifidobacterium_distances_vector) #14
median(Blautia_distances_vector) #13
median(Streptococcus_distances_vector) #13
median(Faecalibacterium_distances_vector) #6
median(Escherichia_Shigella_distances_vector) #2

## resort levels of genera to show nicely in plot
df_SNP_dist_all$genus <- as.factor(as.character(df_SNP_dist_all$genus))
df_SNP_dist_all$genus <- factor(df_SNP_dist_all$genus, levels = c("Prevotella", "Veillonella", "Bacteroides", "Clostridium_sensu_stricto", "Haemophilus",
                                                              "Bifidobacterium", "Blautia", "Streptococcus", "Faecalibacterium", "Escherichia_Shigella"))


## plot SNP distances by genus
SNP_dist_plot1 <- ggplot(df_SNP_dist_all, aes(x=genus, y=SNP_distance))+
  geom_jitter(alpha=0.5, size=0.2, color="#0072B2")+
  geom_boxplot(outlier.shape = NA, alpha=0)+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=0),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  labs(x="", y="SNP distance")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/240228_matf_milk_ASV_SNP_distances_by_genus.pdf", SNP_dist_plot1, dpi=500, width=20, height=20, unit="cm")


### ===== 7. EXTRACT NUMBER OF UNIQUE ASVs PER GENUS ===== ###

## extract n(unique ASVs) detected per genus
Bacteroides_n_ASVs <- length(rownames(Bacteroides_distances_2))
Bifidobacterium_n_ASVs <- length(rownames(Bifidobacterium_distances_2))
Blautia_n_ASVs <- length(rownames(Blautia_distances_2))
Clostridium_sensu_stricto_n_ASVs <- length(rownames(Clostridium_sensu_stricto_distances_2))
Escherichia_Shigella_n_ASVs <- length(rownames(Escherichia_Shigella_distances_2))
Faecalibacterium_n_ASVs <- length(rownames(Faecalibacterium_distances_2))
Haemophilus_n_ASVs <- length(rownames(Haemophilus_distances_2))
Prevotella_n_ASVs <- length(rownames(Prevotella_distances_2))
Streptococcus_n_ASVs <- length(rownames(Streptococcus_distances_2))
Veillonella_n_ASVs <- length(rownames(Veillonella_distances_2))


### ===== 8. CREATE DATA FRAMES WITH GENUS AND NUMBER OF UNIQUE ASVs ===== ###

## combine n(unique ASVs) detected in 1 data frame
df_n_ASVs_all <- data.frame(genus=c("Escherichia_Shigella", "Faecalibacterium", "Streptococcus", "Blautia", "Bifidobacterium",
                                    "Haemophilus", "Clostridium_sensu_stricto", "Bacteroides", "Veillonella", "Prevotella"),
                            n_ASVs=c(Escherichia_Shigella_n_ASVs, Faecalibacterium_n_ASVs, Streptococcus_n_ASVs, Blautia_n_ASVs, Bifidobacterium_n_ASVs,
                                     Haemophilus_n_ASVs, Clostridium_sensu_stricto_n_ASVs, Bacteroides_n_ASVs, Veillonella_n_ASVs, Prevotella_n_ASVs))


### ===== 9. PREPARE DATA FOR PLOTTING AND PLOT NUMBER OF UNIQUE ASVs BY GENUS ===== ###

# library(ggplot2)

## sort genera as in SNP distances plot
df_n_ASVs_all$genus <- as.factor(as.character(df_n_ASVs_all$genus))
df_n_ASVs_all$genus <- factor(df_n_ASVs_all$genus, levels = c("Prevotella", "Veillonella", "Bacteroides", "Clostridium_sensu_stricto", "Haemophilus",
                                                              "Bifidobacterium", "Blautia", "Streptococcus", "Faecalibacterium", "Escherichia_Shigella"))


## plot SNP distances by genus
n_ASVs_plot1 <- ggplot(df_n_ASVs_all, aes(x=genus, y=n_ASVs))+
  geom_bar(stat="identity", fill="#0072B2")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.x = element_text(angle=90, hjust=1, vjust=1),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  labs(x="", y="Number of unique ASVs")+
  scale_y_continuous(limits = c(0,400))+
  geom_text(aes(label=n_ASVs), position=position_dodge(0.9), vjust=-1)
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/240228_matf_milk_number_of_unique_ASVs_per_genus.pdf", n_ASVs_plot1, dpi=500, width=30, height=20, unit="cm")


### ===== 10. CREATE COMBINED PLOT ===== ###

library(cowplot)

## create combined plot
combined_plot1 <- plot_grid(SNP_dist_plot1, n_ASVs_plot1, ncol=1)
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/240228_matf_milk_SNP_distances_and_number_of_unique_ASVs_per_genus.pdf", combined_plot1, dpi=500, width=30, height=20, unit="cm")


### ===== 11. CREATE COMBINED PLOT - TURNED ===== ###

## reverse levels in both data frames so they show in the correct order from top to bottom later
df_n_ASVs_all2 <- df_n_ASVs_all
df_n_ASVs_all2$genus <- factor(df_n_ASVs_all2$genus, levels = c("Escherichia_Shigella", "Faecalibacterium", "Streptococcus", "Blautia", "Bifidobacterium",
                                                                "Haemophilus", "Clostridium_sensu_stricto", "Bacteroides", "Veillonella", "Prevotella"))

df_SNP_dist_all2 <- df_SNP_dist_all
df_SNP_dist_all2$genus <- factor(df_SNP_dist_all2$genus, levels = c("Escherichia_Shigella", "Faecalibacterium", "Streptococcus", "Blautia", "Bifidobacterium",
                                                                    "Haemophilus", "Clostridium_sensu_stricto", "Bacteroides", "Veillonella", "Prevotella"))


## plot number of unique ASVs by genus, turned
n_ASVs_plot2 <- ggplot(df_n_ASVs_all2, aes(x=n_ASVs, y=genus))+
  geom_bar(stat="identity", fill="#0072B2")+
  theme_bw()+
  theme(legend.position="none",
        axis.text.y = element_text(hjust=1, vjust=1),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  labs(x="", y="", title="Number of unique ASVs")+
  scale_x_continuous(limits = c(0,500))+
  geom_text(aes(label=n_ASVs), position=position_dodge(0.9), hjust=-0.5)
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/240228_matf_milk_number_of_unique_ASVs_per_genus_turned.pdf", n_ASVs_plot2, dpi=500, width=20, height=30, unit="cm")

## plot SNP distances by genus, turned
SNP_dist_plot2 <- ggplot(df_SNP_dist_all2, aes(x=SNP_distance, y=genus))+
  geom_jitter(alpha=0.5, size=0.2, color="#0072B2")+
  geom_boxplot(outlier.shape = NA, alpha=0)+
  theme_bw()+
  theme(legend.position="none",
        # axis.text.y = element_text(angle=90, hjust=0, vjust=1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank())+
  labs(x="", y="", title="SNP distance")
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/240228_matf_milk_ASV_SNP_distances_by_genus_turned.pdf", SNP_dist_plot2, dpi=500, width=20, height=30, unit="cm")


# library(cowplot)
combined_plot2 <- plot_grid(n_ASVs_plot2, SNP_dist_plot2, nrow=1)
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/all_matfmilk_ASVs/240228_matf_milk_SNP_distances_and_number_of_unique_ASVs_per_genus_turned.pdf", combined_plot2, dpi=500, width=20, height=30, unit="cm")





