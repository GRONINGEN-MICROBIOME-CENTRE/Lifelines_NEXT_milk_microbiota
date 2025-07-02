################################################################################
##### ASVs sharing: statistical analysis and plots 
### Author(s): Asier Fernández / Johanne E. Spreckels
### Last updated: 16th July, 2024
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(ggplot2)

# Set working directory
setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/")

# Read results files
final_PSE_data_milk_feces_infant <- read.delim("Mother_Milk-Baby_Feces/PSE_table_milk_infant_feces.txt", header = T, check.names = F) 
final_PSE_data_milk_feces_mother <- read.delim("Mother_Milk-Mother_Feces/PSE_table_milk_maternal_feces.txt", header = T, check.names = F) 
cont_tables_milk_feces_infant <- readRDS("Mother_Milk-Baby_Feces/Contingency_tables_milk_infant_feces.rds")
cont_tables_milk_feces_mother <- readRDS("Mother_Milk-Mother_Feces/Contingency_tables_milk_maternal_feces.rds")

###############################################################
# Generate PSE results table and test statistical significance 
###############################################################

#################################################
#A. 16S - Mother Milk vs Baby Feces 
#################################################

#******************************
# Test PSE differences in related vs unrelated pairs
#******************************
# Reformat the data table 
final_PSE_data_milk_feces_infant$Pair_ID <- c(rownames(final_PSE_data_milk_feces_infant))
final_PSE_data_milk_feces_infant_long <- final_PSE_data_milk_feces_infant %>%
  pivot_longer(cols = -c(Pair_ID, Relatedness), names_to = "Genus_Comparison", values_to = "Value") %>%
  separate(Genus_Comparison, into = c("Genus", "Comparison"), sep = "_", extra = "merge", fill = "right")

final_PSE_data_milk_feces_infant_summary <- final_PSE_data_milk_feces_infant_long %>%
  dplyr::group_by(Genus, Value, Relatedness) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::pivot_wider(names_from = Value, values_from = count, values_fill = 0)

final_PSE_data_milk_feces_infant_summary$all_no_PSE <- rowSums(data.frame(final_PSE_data_milk_feces_infant_summary$no_PSE, final_PSE_data_milk_feces_infant_summary$`NA`))


# Generate contingency tables and estimate p-values for each genus
Fischer_tests <- data.frame(matrix(nrow=length(unique(final_PSE_data_milk_feces_infant_summary$Genus)), ncol=1))
colnames(Fischer_tests) <- "p_value"
rownames(Fischer_tests) <- unique(final_PSE_data_milk_feces_infant_summary$Genus)

for (i in 1:length(rownames(Fischer_tests))) {
  genus_data <- final_PSE_data_milk_feces_infant_summary[final_PSE_data_milk_feces_infant_summary$Genus==rownames(Fischer_tests)[i],]
  cont_table <- data.frame(rbind(c(genus_data[genus_data$Relatedness=="Related", "all_no_PSE"], 
                                   genus_data[genus_data$Relatedness=="Related", "PSE"]),
                                 c(genus_data[genus_data$Relatedness=="Unrelated","all_no_PSE"], 
                                   genus_data[genus_data$Relatedness=="Unrelated","PSE"])))
  cont_table <- matrix(unlist(cont_table), 2)
  colnames(cont_table) <- c("no_PSE", "PSE")
  rownames(cont_table) <- c("Related", "Unrelated")
  test <- fisher.test(cont_table)
  Fischer_tests[i,1] <- test$p.value  
}

#Estimate the FDR 
Fischer_tests$FDR <- p.adjust(Fischer_tests$p_value, method = "BH")


# Reformat results table in long format
final_PSE_data_milk_feces_infant_summary_long <- final_PSE_data_milk_feces_infant_summary %>%
  pivot_longer(cols = c("PSE", "no_PSE", "NA"), names_to = "Value", values_to = "count")

# Add the total counts, percentages and propotion for each group to the dataframe (Relatedness, Genus, Technology)
total_count <- final_PSE_data_milk_feces_infant_summary_long %>% 
  dplyr::group_by(Relatedness, Genus) %>% 
  dplyr::summarize(total = sum(count))

final_PSE_data_milk_feces_infant_summary_long <- dplyr::left_join(final_PSE_data_milk_feces_infant_summary_long, total_count)
final_PSE_data_milk_feces_infant_summary_long$percent <- paste0(round(final_PSE_data_milk_feces_infant_summary_long$count / final_PSE_data_milk_feces_infant_summary_long$total * 100), "%")
final_PSE_data_milk_feces_infant_summary_long$proportion <- final_PSE_data_milk_feces_infant_summary_long$count / final_PSE_data_milk_feces_infant_summary_long$total * 100
final_PSE_data_milk_feces_infant_summary_long$all_no_PSE <- NULL

#################################################
#B. 16S - Mother Milk vs Mother Feces 
#################################################

#******************************
# Test PSE differences in related vs unrelated pairs
#******************************
# Reformat the data table 
final_PSE_data_milk_feces_mother$Pair_ID <- c(rownames(final_PSE_data_milk_feces_mother))
final_PSE_data_milk_feces_mother_long <- final_PSE_data_milk_feces_mother %>%
  pivot_longer(cols = -c(Pair_ID, Relatedness), names_to = "Genus_Comparison", values_to = "Value") %>%
  separate(Genus_Comparison, into = c("Genus", "Comparison"), sep = "_", extra = "merge", fill = "right")

final_PSE_data_milk_feces_mother_summary <- final_PSE_data_milk_feces_mother_long %>%
  dplyr::group_by(Genus, Value, Relatedness) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::pivot_wider(names_from = Value, values_from = count, values_fill = 0)

final_PSE_data_milk_feces_mother_summary$all_no_PSE <- rowSums(data.frame(final_PSE_data_milk_feces_mother_summary$no_PSE, final_PSE_data_milk_feces_mother_summary$`NA`))


# Generate contingency tables and estimate p-values for each genus
Fischer_tests <- data.frame(matrix(nrow=length(unique(final_PSE_data_milk_feces_mother_summary$Genus)), ncol=1))
colnames(Fischer_tests) <- "p_value"
rownames(Fischer_tests) <- unique(final_PSE_data_milk_feces_mother_summary$Genus)

for (i in 1:length(rownames(Fischer_tests))) {
  genus_data <- final_PSE_data_milk_feces_mother_summary[final_PSE_data_milk_feces_mother_summary$Genus==rownames(Fischer_tests)[i],]
  cont_table <- data.frame(rbind(c(genus_data[genus_data$Relatedness=="Related", "all_no_PSE"], 
                                   genus_data[genus_data$Relatedness=="Related", "PSE"]),
                                 c(genus_data[genus_data$Relatedness=="Unrelated","all_no_PSE"], 
                                   genus_data[genus_data$Relatedness=="Unrelated","PSE"])))
  cont_table <- matrix(unlist(cont_table), 2)
  colnames(cont_table) <- c("no_PSE", "PSE")
  rownames(cont_table) <- c("Related", "Unrelated")
  test <- fisher.test(cont_table)
  Fischer_tests[i,1] <- test$p.value  
}

#Estimate the FDR 
Fischer_tests$FDR <- p.adjust(Fischer_tests$p_value, method = "BH")


# Reformat results table in long format
final_PSE_data_milk_feces_mother_summary_long <- final_PSE_data_milk_feces_mother_summary %>%
  pivot_longer(cols = c("PSE", "no_PSE", "NA"), names_to = "Value", values_to = "count")

# Add the total counts, percentages and propotion for each group to the dataframe (Relatedness, Genus, Technology)
total_count <- final_PSE_data_milk_feces_mother_summary_long %>% 
  dplyr::group_by(Relatedness, Genus) %>% 
  dplyr::summarize(total = sum(count))

final_PSE_data_milk_feces_mother_summary_long <- dplyr::left_join(final_PSE_data_milk_feces_mother_summary_long, total_count)
final_PSE_data_milk_feces_mother_summary_long$percent <- paste0(round(final_PSE_data_milk_feces_mother_summary_long$count / final_PSE_data_milk_feces_mother_summary_long$total * 100), "%")
final_PSE_data_milk_feces_mother_summary_long$proportion <- final_PSE_data_milk_feces_mother_summary_long$count / final_PSE_data_milk_feces_mother_summary_long$total * 100
final_PSE_data_milk_feces_mother_summary_long$all_no_PSE <- NULL


#******************************
# Generate the plots
#******************************
# Generate the stacked bar plots:

#A) 16S - Mother Milk vs Baby Feces 
final_PSE_data_milk_feces_infant_summary_long$Value <- factor(final_PSE_data_milk_feces_infant_summary_long$Value, levels = c("NA", "no_PSE", "PSE"))

pdf('Mother_Milk-Baby_Feces/PSE_relatedness_mother_milk_baby_feces.pdf', width=12, height=4)
ggplot(final_PSE_data_milk_feces_infant_summary_long, aes(x = Relatedness, y = proportion, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Genus, switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#E69F00","#0072B2"), labels = c("Genus not present","No", "Yes")) +
  scale_y_continuous(breaks=c(0,25, 50, 75, 100), expand=c(0.07,0),limits=c(0,107)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  # Add percentage labels for the PSE group
  geom_label(aes(label = percent),
             fill = "#FFFFFF",
             data = subset(final_PSE_data_milk_feces_infant_summary_long, Value == "PSE"),
             position = position_fill(vjust = 50), size = 4.7) +
  labs(y = "Percentage of pairs with PSEs", fill = "Sharing")
dev.off()

#B) 16S - Mother Milk vs Mother Feces 
final_PSE_data_milk_feces_mother_summary_long$Value <- factor(final_PSE_data_milk_feces_mother_summary_long$Value, levels = c("NA", "no_PSE", "PSE"))

pdf('Mother_Milk-Mother_Feces//PSE_relatedness_mother_milk_mother_feces.pdf', width=9, height=4)
ggplot(final_PSE_data_milk_feces_mother_summary_long, aes(x = Relatedness, y = proportion, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Genus, switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#CCCCCC","#E69F00","#0072B2"), labels = c("Genus not present","No", "Yes")) +
  scale_y_continuous(breaks=c(0,25, 50, 75, 100), expand=c(0.07,0),limits=c(0,107)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", linewidth = 1),
        strip.text.x = element_text(size = 13),
        strip.background = element_rect(color = "black", linewidth = 1),
        axis.text.x = element_text(size =13, angle =90, hjust = 1),
        axis.text.y = element_text(size=14), 
        axis.title.y = element_text(size=13), 
        axis.title.x = element_blank(),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  # Add percentage labels for the PSE group
  geom_label(aes(label = percent),
             fill = "#FFFFFF",
             data = subset(final_PSE_data_milk_feces_mother_summary_long, Value == "PSE"),
             position = position_fill(vjust = 50), size = 4.7) +
  labs(y = "Percentage of pairs with PSEs", fill = "Sharing")
dev.off()

########################
# Save results
########################
setwd("~/Desktop/PhD/Projects/Johanne Milk paper/RESULTS/")

write.table(final_PSE_data_milk_feces_infant_summary_long, file ="Mother_Milk-Baby_Feces/Mother_milk_Baby_Feces_final_PSE_results.txt", sep="\t")
write.table(final_PSE_data_milk_feces_mother_summary_long, file ="Mother_Milk-Mother_Feces/Mother_milk_Mother_Feces_final_PSE_results.txt", sep="\t")

