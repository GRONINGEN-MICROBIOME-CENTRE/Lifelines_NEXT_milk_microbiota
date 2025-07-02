################################################################################################################
### INVESTIGATE SHARING OF BACTERIA BETWEEN MATERNAL FAECES AND MILK - SCRIPT 3 - M3 (FIGURE S3) ###############
################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.2.1 (2022-06-23)):    R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE

# 1. LOAD DATA FROM SCRIPT 2 AND LIBRARIES
# 2. GENERATE PSE RESULTS TABLE AND TEST STATISTICAL SIGNIFICANCE
# 3. GENERATE PLOTS (FIGURE S3)
# 4. SAVE RESULTS


##### =========================== 1. LOAD DATA FROM SCRIPT 2 AND LIBRARIES =========================== #####

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(ggplot2)

# Set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/8_maternal_faeces_microbiota/sharing_with_milk/")

# Read results files
final_PSE_data_M3 <- read.delim("M3/240219_PSE_table_maternal_faeces_milk_M3.txt", header = T, check.names = F) 
cont_tables_M3 <- readRDS("M3/240219_maternal_faeces_milk_contingency_tables_M3.rds")


##### =========================== 2. GENERATE PSE RESULTS TABLE AND TEST STATISTICAL SIGNIFICANCE =========================== #####

#******************************
# Test PSE differences in related vs unrelated pairs
#******************************
# Reformat the data table 
final_PSE_data_M3$Pair_ID <- c(rownames(final_PSE_data_M3))
colnames(final_PSE_data_M3)[grep("Escherichia_Shigella", colnames(final_PSE_data_M3))] <- "EscherichiaShigella_M3"
colnames(final_PSE_data_M3)[grep("Clostridium_sensu_stricto", colnames(final_PSE_data_M3))] <- "Clostridiumsensustricto_M3"
final_PSE_data_long_M3 <- final_PSE_data_M3 %>%
  pivot_longer(cols = -c(Pair_ID, Relatedness), names_to = "Genus_Timepoint", values_to = "Value") %>%
  separate(Genus_Timepoint, into = c("Genus", "Timepoint"), sep = "_")

final_PSE_data_summary_M3 <- final_PSE_data_long_M3 %>%
  dplyr::group_by(Genus, Timepoint, Value, Relatedness) %>%
  dplyr::summarize(count = n()) %>%
  tidyr::pivot_wider(names_from = Value, values_from = count, values_fill = 0)

final_PSE_data_summary_M3$all_no_PSE <- rowSums(data.frame(final_PSE_data_summary_M3$no_PSE, final_PSE_data_summary_M3$`NA`))
final_PSE_data_summary_M3$Genus_Timepoint <- paste(final_PSE_data_summary_M3$Genus, final_PSE_data_summary_M3$Timepoint, sep = "_")


# Generate contingency tables and estimate p-values for each genus
Fischer_tests_M3 <- data.frame(matrix(nrow=length(unique(final_PSE_data_summary_M3$Genus_Timepoint)), ncol=1))
colnames(Fischer_tests_M3) <- "p_value"
rownames(Fischer_tests_M3) <- unique(final_PSE_data_summary_M3$Genus_Timepoint)

for (i in 1:length(rownames(Fischer_tests_M3))) {
  genus_data <- final_PSE_data_summary_M3[final_PSE_data_summary_M3$Genus_Timepoint==rownames(Fischer_tests_M3)[i],]
  print(genus_data)
  cont_table <- data.frame(rbind(c(genus_data[genus_data$Relatedness=="Related","all_no_PSE"], #7
                                   genus_data[genus_data$Relatedness=="Related","PSE"]),
                                 c(genus_data[genus_data$Relatedness=="Unrelated","all_no_PSE"], #7
                                   genus_data[genus_data$Relatedness=="Unrelated","PSE"])))
  cont_table <- matrix(unlist(cont_table), 2)
  colnames(cont_table) <- c("no_PSE", "PSE")
  rownames(cont_table) <- c("Related", "Unrelated")
  print(cont_table)
  test <- fisher.test(cont_table)
  Fischer_tests_M3[i,1] <- test$p.value  
}

#Estimate the FDR per time point
Fischer_tests_M3$FDR <- p.adjust(Fischer_tests_M3$p_value, method="BH")

Fischer_tests_M3$Genus_Timepoint <- rownames(Fischer_tests_M3)
final_PSE_data_summary_M3 <- merge(final_PSE_data_summary_M3, Fischer_tests_M3, by="Genus_Timepoint")

# Reformat results table in long format
final_PSE_data_summary_long_M3 <- final_PSE_data_summary_M3 %>%
  pivot_longer(cols = c("PSE", "no_PSE", "NA"), names_to = "Value", values_to = "count")

# Add the total counts, percentages and proportion for each group to the dataframe (Relatedness, Genus, Timepoint)
total_count_M3 <- final_PSE_data_summary_long_M3 %>% 
  dplyr::group_by(Relatedness, Genus, Timepoint) %>% 
  dplyr::summarize(total = sum(count))

final_PSE_data_summary_long_M3 <- dplyr::left_join(final_PSE_data_summary_long_M3, total_count_M3)
final_PSE_data_summary_long_M3$percent <- paste0(round(final_PSE_data_summary_long_M3$count / final_PSE_data_summary_long_M3$total * 100), "%")
final_PSE_data_summary_long_M3$proportion <- final_PSE_data_summary_long_M3$count / final_PSE_data_summary_long_M3$total * 100
final_PSE_data_summary_long_M3$all_no_PSE <- NULL


##### =========================== 3. GENERATE PLOTS (FIGURE S3) =========================== #####

#******************************
# Generate the plots
#******************************
# Generate the stacked bar plot for each time point
# M3
final_PSE_data_summary_long_M3$Value <- factor(final_PSE_data_summary_long_M3$Value, levels = c("NA", "no_PSE", "PSE"))

pdf('M3/240219_plots_PSE_relatedness_maternal_faeces_milk_M3.pdf', width=30, height=5)
ggplot(final_PSE_data_summary_long_M3, aes(x = Relatedness, y = proportion, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Genus, switch = "y") +
  labs(y = "Number of pairs", fill = "Sharing") +
  scale_fill_manual(values = c("#F4F4F4","#CCCCCC","#0072B2"), labels = c("Genus not present","No", "Yes")) +  #"#F0E7DD","#CCCCCC","#009E73"
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
             data = subset(final_PSE_data_summary_long_M3, Value == "PSE"),
             position = position_fill(vjust = 50), size = 4.7) +
  labs(y = "Percentage of pairs with PSEs", fill = "Sharing")
dev.off()


##### =========================== 4. SAVE RESULTS =========================== #####

## detailed table
write.table(final_PSE_data_summary_long_M3, file ="M3/240219_maternal_faeces_milk_final_PSE_results_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)


## create shortened table
# split data by relatedness
final_PSE_data_summary_M3_rel <- final_PSE_data_summary_M3[final_PSE_data_summary_M3$Relatedness=="Related",]
final_PSE_data_summary_M3_unrel <- final_PSE_data_summary_M3[final_PSE_data_summary_M3$Relatedness=="Unrelated",]

# indicate relaetdness in colnames
colnames(final_PSE_data_summary_M3_rel)[5:8] <- c(paste0(colnames(final_PSE_data_summary_M3_rel)[5:8], "_related"))
colnames(final_PSE_data_summary_M3_unrel)[5:8] <- c(paste0(colnames(final_PSE_data_summary_M3_unrel)[5:8], "_unrelated"))

# re-merge related and unrelated data
final_PSE_data_summary_M3_v2 <- merge(final_PSE_data_summary_M3_rel[,c(2:3,5:8)], final_PSE_data_summary_M3_unrel[,c(2,5:8,9:10)], by="Genus")

# create columns for n_total_realted and n_total_unrelated
final_PSE_data_summary_M3_v2$n_total_related <- c(final_PSE_data_summary_M3_v2[,3] + final_PSE_data_summary_M3_v2[,4] + final_PSE_data_summary_M3_v2[,5])
final_PSE_data_summary_M3_v2$n_total_unrelated <- c(final_PSE_data_summary_M3_v2[,7] + final_PSE_data_summary_M3_v2[,8] + final_PSE_data_summary_M3_v2[,9])

# calculate proportions
final_PSE_data_summary_M3_v2$perc_no_PSE_related <- final_PSE_data_summary_M3_v2$no_PSE_related / final_PSE_data_summary_M3_v2$n_total_related
final_PSE_data_summary_M3_v2$perc_PSE_related <- final_PSE_data_summary_M3_v2$PSE_related / final_PSE_data_summary_M3_v2$n_total_related
final_PSE_data_summary_M3_v2$perc_NA_related <- final_PSE_data_summary_M3_v2$NA_related / final_PSE_data_summary_M3_v2$n_total_related
final_PSE_data_summary_M3_v2$perc_all_no_PSE_related <- final_PSE_data_summary_M3_v2$all_no_PSE_related / final_PSE_data_summary_M3_v2$n_total_related

final_PSE_data_summary_M3_v2$perc_no_PSE_unrelated <- final_PSE_data_summary_M3_v2$no_PSE_unrelated / final_PSE_data_summary_M3_v2$n_total_unrelated
final_PSE_data_summary_M3_v2$perc_PSE_unrelated <- final_PSE_data_summary_M3_v2$PSE_unrelated / final_PSE_data_summary_M3_v2$n_total_unrelated
final_PSE_data_summary_M3_v2$perc_NA_unrelated <- final_PSE_data_summary_M3_v2$NA_unrelated / final_PSE_data_summary_M3_v2$n_total_unrelated
final_PSE_data_summary_M3_v2$perc_all_no_PSE_unrelated <- final_PSE_data_summary_M3_v2$all_no_PSE_unrelated / final_PSE_data_summary_M3_v2$n_total_unrelated

# resort columns
final_PSE_data_summary_M3_v3 <- final_PSE_data_summary_M3_v2[,c(2,1,
                                                                13,4,16,3,15,5,17,6,18,
                                                                14,8,20,7,19,9,21,10,22,
                                                                11:12)]
# change colnames
colnames(final_PSE_data_summary_M3_v3)[c(4,6,8,10,13,15,17,19)] <- c(paste0("n_", colnames(final_PSE_data_summary_M3_v3)[c(4,6,8,10,13,15,17,19)]))
colnames(final_PSE_data_summary_M3_v3)[c(8,9,17,18)] <- c("n_not_present_related", "perc_not_present_related", "n_not_present_unrelated", "perc_not_present_unrelated")

## save shortened table
write.table(final_PSE_data_summary_M3_v3, file ="M3/240219_maternal_faeces_milk_final_PSE_results_short_M3.txt", row.names=F, col.names=T, sep="\t", quote=F)






