################################################################################################################
### CREATE BASIC COHORT SUMMARY STATISTICS (FOR HMO DATA) ######################################################
################################################################################################################

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=20gb --nodes=1 --qos=priority --time=6-23:59:59 --pty bash -i
#                                           srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=04:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.0.3 (2020-10-10):     R
#for long sessoins type "screen" before requesting the job, when coming back online, re-open the session using "screen -rd"
#if you want to remove the screen session, after entering it, type ctrl a (nothing shows up) and then press k and then confirm with y


### CONTENT OF THIS FILE
# 0. DATA IMPORT
# 1. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES
# 2. BASIC INFORMATION / STATISTICS
# 3. CREATE BASIC COHORT PHENOTYPE SUMMARY STATISTICS AND PLOTS FOR FIGURE 1


##### =========================== 0. DATA IMPORT =========================== #####

## set working directory
setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/phenotypes/")

## import files
h <- read.table("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/phenotypes/240105_milk_composition_paper_sel_phenotypes_hmo_data_n1563.txt", header=T, sep="\t", stringsAsFactors=T) #1563x511


##### =========================== 1. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

## for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
h$mother_milk_collection_season <- factor(h$mother_milk_collection_season, levels = c("Spring", "Summer", "Autumn", "Winter"))
h$mother_milk_collection_month <- factor(h$mother_milk_collection_month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",  "Oct",  "Nov", "Dec"))
h$mother_milk_collection_notes <- factor(h$mother_milk_collection_notes, levels = c("No_records_of_incorrect_sample_handling", levels(h$mother_milk_collection_notes)[c(3,5,7,6,2,1)]))
h$mother_milk_collection_breasts_sampled <- factor(h$mother_milk_collection_breasts_sampled, levels=c("one_breast", "both_breasts"))
h$mother_genetics_blood_group_genotype <- factor(h$mother_genetics_blood_group_genotype, levels = c("O/O", "A/O", "A/A", "B/O", "B/B", "A/B"))
h$mother_blood_group <- factor(h$mother_blood_group, levels = c("O", "A", "B", "AB"))
h$mother_exp_living_situation <- factor(h$mother_exp_living_situation, levels = c("big_house", "small_house", "flat", "farm", "other"))
h$mother_birth_gestational_age_categories <- factor(h$mother_birth_gestational_age_categories, levels = c("term", "preterm"))
h$mother_breastpump_brand <- factor(h$mother_breastpump_brand, levels = c("Philips", "Medela", "Ardo", "Other"))
h$infant_birth_delivery_mode <- factor(h$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
h$infant_birth_delivery_mode_detailed <- factor(h$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                      "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
h$infant_ffq_feeding_type_delivery <- factor(h$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))
h$mother_milk_HMO_milk_group <- factor(h$mother_milk_HMO_milk_group, levels = c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 2. BASIC INFORMATION / STATISTICS =========================== #####

## select every mother only 1x (i.e. remove duplication of data generated for twin babies)
hmo <- h[grep(";1", h$mother_sample_ID),] #1542x511

## number of total measured samples 
length(unique(hmo$mother_sample_ID)) #1542

## number of samples by time point
table(hmo$time_point)
# 0.5_months    1_month   2_months   3_months   6_months 
#        214        419        340        433        136

## number of mothers
length(unique(hmo$mother_ID)) #524 mother-participations

hmo$mother_ID_simple <- gsub("_1_", "_", hmo$mother_ID)
hmo$mother_ID_simple <- as.factor(as.character(gsub("_2_", "_", hmo$mother_ID_simple)))
hmo <- hmo[,c(1:6,512,7:511)]
length(unique(hmo$mother_ID_simple)) #500 -> 524-500=24 -> 24 mothers have HMO measurements for first and second NEXT pregnancy

## number of infants #!!!Note: use the data frame before exclusion of repeated maternal data (h)!!!
length(unique(h$infant_ID)) #532

## number of twins
twins <- h[h$mother_type_pregnancy=="twin_pregnancy",]
twin1 <- twins[grep("Infant1", twins$infant_ID),]
twin2 <- twins[grep("Infant2", twins$infant_ID),]
length(unique(twin1$infant_ID)) #8
length(unique(twin2$infant_ID)) #8

fam <- h[,c(3,7,9)]
fam <- fam[!duplicated(fam),]
table(fam$mother_type_pregnancy, useNA="ifany") #516 single pregnancies + 16 twin pregnancies (this counts every twin pregnancy 2x)
length(unique(fam$infant_ID)) #532 unique infants
length(fam[grep("Infant1", fam$infant_ID),2]) #524 unique infant1
length(fam[grep("Infant2", fam$infant_ID),2]) #8 unique infant2


## check number of samples per mother
table(table(hmo$mother_ID)) #manually exclude ;2 mothers - they don't have samples in this data frame; total number of mothers -> 524
#number of samples (n)              1   2   3   4   5 
#number of mothers with n samples  57 131 146 165  25 
#%                                11% 25% 28% 31%  5%

## average number of samples per mother: 1542/524 = 2.94 --> 3 samples/mother


##### =========================== 3. CREATE BASIC COHORT PHENOTYPE SUMMARY STATISTICS AND PLOTS FOR FIGURE 1 =========================== #####

## save maternal cross-sectional phenotypes for basic summary statistics separately, exclude duplicated data for twins
mc <- h[-grep(";2", h$mother_ID),c(6,30,31,364,395:396)]
colnames(mc)

mc2 <- mc[!duplicated(mc),] #524 unique mother-participations
mc2$xcolumn <- "mother"

for (i in 1:ncol(mc2)){
  print(colnames(mc2)[i])
  print(table(mc2[,i], useNA="ifany"))
}


summary(mc2$mother_age_birth)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 23.00   30.00   32.00   32.57   35.00   50.00      26

sd(mc2$mother_age_birth, na.rm=T)
#3.935097

matagedata <- mc2[!is.na(mc2$mother_age_birth),1:2]

library(ggplot2)

histmatage <- ggplot(matagedata, aes(x=mother_age_birth))+
  geom_histogram(binwidth=1, color="black", fill="seagreen1")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_x_continuous(limits = c(18,52), breaks = c(20,30,40,50), labels = c(20,30,40,50))+
  scale_y_continuous(limits = c(0,60), breaks = c(0,20,40,60), labels = c(0,20,40,60))+
  labs(x="", y="Number of mothers", title="Maternal age at birth")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information/Histogram_maternal_age.pdf", histmatage, dpi=500, width=10, height=10, units="cm", useDingbats=F)

# matageplot <- ggplot(mc2[!is.na(mc2$mother_age_birth),], aes(x=xcolumn, y=mother_age_birth, color=xcolumn))+
#   geom_jitter(alpha=0.5)+
#   scale_color_manual(values="seagreen1")+
#   geom_errorbar(aes(ymin=mean(mc2$mother_age_birth, na.rm=T)-sd(mc2$mother_age_birth, na.rm=T),
#                     ymax=mean(mc2$mother_age_birth, na.rm=T)+sd(mc2$mother_age_birth, na.rm=T)),
#                 width=0.5)+
#   geom_segment(aes(x=0.7, xend=1.3, y=mean(mc2$mother_age_birth, na.rm=T), yend=mean(mc2$mother_age_birth, na.rm=T)))+
#   theme_bw()+
#   theme(panel.grid = element_blank(),
#         legend.position="none")+
#   scale_y_continuous(limits = c(18,52), breaks = c(20,30,40,50), labels = c(20,30,40,50))
# ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information/Dotplot_maternal_age.pdf", matageplot, dpi=500, width=15, height=10, units="cm", useDingbats=F)


summary(mc2$mother_birth_gestational_age_weeks)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 29.85   38.86   39.86   39.63   40.71   42.43     135

sd(mc2$mother_birth_gestational_age_weeks, na.rm=T)
#1.598759

matGAdata <- mc2[!is.na(mc2$mother_birth_gestational_age_weeks),c(1,6)]

histmatGA <- ggplot(matGAdata, aes(x=mother_birth_gestational_age_weeks))+
  geom_histogram(binwidth=1, color="black", fill="seagreen1")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_x_continuous(limits = c(28,44), breaks = c(30,32,34,36,38,40,42))+
  scale_y_continuous(limits = c(0,120))+
  labs(x="", y="Number of mothers", title="Gestational age at birth")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information/Histogram_gestational_age.pdf", histmatGA, dpi=500, width=10, height=10, units="cm", useDingbats=F)


table(mc2$mother_education_level, useNA="ifany")
# 1__elementary_to_lower_secondary_education 
# 6 (2%)
# 2__upper_secondary_education 
# 66 (18%)
# 3__tertiary_education 
# 300 (81%)
# <NA> 
# 152


## plot mother_education_level

# library(ggplot2)

# prepare data for plotting
# group=NA (dummy)
# subgroup=mother_education_level
# value=percentage

edudata <- data.frame(A=as.factor(as.character(c("A"))),
                           mother_education=as.factor(as.character(c("elementary_or_secondary_education","tertiary_education"))),
                           percentage=c(0.19,0.81))

eduplot <- ggplot(edudata, aes(x=A, y=percentage, fill=mother_education))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("elementary_or_secondary_education"="#F7F1EB",
                             "tertiary_education"="navy"))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1), labels = c("0%","25%","50%","75%","100%"))+
  labs(x="", y="Percentage", fill="Maternal education", title="Maternal education")+
  annotate(geom="text", x=1, y=0.84, label=c("81%"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information/Stacked_barplot_maternal_education.pdf", eduplot, dpi=500, width=15, height=10, units="cm")



## save infant cross-sectional phenotypes for basic summary statistics separately
ic <- h[,c(7,437,440,442)]
colnames(ic)

ic2 <- ic[!duplicated(ic),] #532 unique infant-participations

for (i in 1:ncol(ic2)){
  print(colnames(ic2)[i])
  print(table(ic2[,i], useNA="ifany"))
}

table(ic2$infant_birth_delivery_mode, useNA="ifany")
# vaginal_birth     C-section          <NA> 
#           415            67            50

table(ic2$infant_sex, useNA="ifany")
# female   male   <NA> 
#    247    258     27

summary(ic2$infant_birthweight_g)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  800    3240    3578    3566    3940    5055     146

sd(ic2$infant_birthweight_g, na.rm=T)
#557.6062


infBWdata <- ic2[!is.na(ic2$infant_birthweight_g),c(1,4)]

histinfBW <- ggplot(infBWdata, aes(x=infant_birthweight_g))+
  geom_histogram(binwidth=100, color="black", fill="seagreen1")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_x_continuous(limits = c(500,5500))+
  scale_y_continuous(limits = c(0,40))+
  labs(x="", y="Number of infants", title="Infant birth weight")
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information/Histogram_infant_birthweight.pdf", histinfBW, dpi=500, width=10, height=10, units="cm", useDingbats=F)


## plot infant_birth_delivery_mode

# library(ggplot2)

# prepare data for plotting
# group=NA (dummy)
# subgroup=delivery_mode
# value=percentage

deliverydata <- data.frame(A=as.factor(as.character(c("A"))),
                           delivery_mode=as.factor(as.character(c("C-section","Vaginal birth"))),
                           percentage=c(0.14,0.86))

deliveryplot <- ggplot(deliverydata, aes(x=A, y=percentage, fill=delivery_mode))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("C-section"="#F7F1EB",
                             "Vaginal birth"="#CC79A7"))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1), labels = c("0%","25%","50%","75%","100%"))+
  labs(x="", y="Percentage", fill="Delivery mode", title="Delivery mode")+
  annotate(geom="text", x=1, y=0.89, label=c("86%"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information/Stacked_barplot_delivery_mode.pdf", deliveryplot, dpi=500, width=10, height=10, units="cm")


## plot infant_sex

# library(ggplot2)

# prepare data for plotting
# group=NA (dummy)
# subgroup=infant_sex
# value=percentage

sexdata <- data.frame(A=as.factor(as.character(c("A"))),
                      sex=as.factor(as.character(c("Girl","Boy"))),
                      percentage=c(0.49,0.51))

sexdata$sex <- factor(sexdata$sex, levels = c("Girl", "Boy"))

sexplot <- ggplot(sexdata, aes(x=A, y=percentage, fill=sex))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("Girl"="#F7F1EB",
                             "Boy"="#56B4E9"))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1), labels = c("0%","25%","50%","75%","100%"))+
  labs(x="", y="Percentage", fill="Infant sex", title="Infant sex")+
  annotate(geom="text", x=1, y=0.54, label=c("51%"))
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information/Stacked_barplot_infant_sex.pdf", sexplot, dpi=500, width=10, height=10, units="cm")



## save dynamic phenotypes for basic summary statistics separately
id <- h[,c(7,8,451)]

id_W2 <- id[id$time_point=="0.5_months",] #217
id_M1 <- id[id$time_point=="1_month",]    #424
id_M2 <- id[id$time_point=="2_months",]   #345
id_M3 <- id[id$time_point=="3_months",]   #437
id_M6 <- id[id$time_point=="6_months",]   #140

table(id_W2$infant_ffq_feeding_mode, useNA="ifany")
# breastfeeding mixed_feeding          <NA> 
#     131 (93%)       10 (7%)           76

table(id_M1$infant_ffq_feeding_mode, useNA="ifany")
# breastfeeding mixed_feeding          <NA> 
#     223 (78%)      63 (22%)           138

table(id_M2$infant_ffq_feeding_mode, useNA="ifany")
# breastfeeding mixed_feeding          <NA> 
#     179 (80%)      44 (20%)           122

table(id_M3$infant_ffq_feeding_mode, useNA="ifany")
# breastfeeding mixed_feeding          <NA> 
#     203 (79%)      53 (21%)           181

table(id_M6$infant_ffq_feeding_mode, useNA="ifany")
# breastfeeding mixed_feeding          <NA> 
#      65 (64%)      37 (36%)            38


## plot infant feeding type by time point

# library(ggplot2)

# prepare data for plotting
# group=time_point
# subgroup=feeding_mode
# value=percentage

feedingdata <- data.frame(time_point=c(rep(c("0.5_months", "1_month", "2_months", "3_months", "6_months"),2)),
                          feeding_mode=c(rep(("Exclusive breastfeeding"),5), rep(("Mixed feeding"),5)),
                          percentage=c(0.93,0.78,0.80,0.79,0.64,0.07,0.22,0.20,0.21,0.36))

for (i in 1:2){feedingdata[,i] <- as.factor(as.character(feedingdata[,i]))}
feedingdata$feeding_mode <- factor(feedingdata$feeding_mode, levels = c("Mixed feeding", "Exclusive breastfeeding"))
mylabel <- c(paste0(feedingdata[feedingdata$feeding_mode=="Exclusive breastfeeding","percentage"]*100, "%"))

feedingplot <- ggplot(feedingdata, aes(x=time_point, y=percentage, fill=feeding_mode))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("Mixed feeding"="#F7F1EB",
                              "Exclusive breastfeeding"="#E69F00"))+
  scale_y_continuous(breaks = c(0,0.25,0.50,0.75,1), labels = c("0%","25%","50%","75%","100%"))+
  labs(x="", y="Percentage", fill="Feeding mode", title="Feeding mode by time point")+
  annotate(geom="text", x=1, y=0.96, label=mylabel[1])+
  annotate(geom="text", x=2, y=0.81, label=mylabel[2])+
  annotate(geom="text", x=3, y=0.83, label=mylabel[3])+
  annotate(geom="text", x=4, y=0.82, label=mylabel[4])+
  annotate(geom="text", x=5, y=0.67, label=mylabel[5])
ggsave("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/1_cohort_basic_information/Stacked_barplot_infant_feeding_mode_by_time_point.pdf", feedingplot, dpi=500, width=15, height=10, units="cm")



## check number of mothers with ABO, FUT2, FUT3 genotypes (unique mothers and number of mothers with 2 participations)
gen <- hmo[!is.na(hmo$mother_genetics_FUT2_G_to_A_rs601338),c(1,6:7,4,22:27,29)]

for (i in c(22:27,29)){
  print(colnames(hmo)[i])
  print(table(hmo[,i], useNA="ifany"))
}

for (i in 5:ncol(gen)){
  print(colnames(gen)[i])
  print(table(gen[,i], useNA="ifany"))
}
# -> there are no NAs in the genotype data in data frame gen and the number of samples with entries match the numbers from the hmo data frame

length(unique(gen$mother_ID)) #249
length(unique(gen$mother_ID_simple)) #229 -> 249-229=20 mothers with genotype data participated with more than 1 pregnancy


## check number of mothers with info on subclinical mastitis (mineral measurements)
table(hmo$mother_dis_subclinical_mastitis, useNA="ifany")
#  no  yes <NA> 
# 554   20  968 
## -> 554+20=574 samples had mineral measurements

scm <- hmo[!is.na(hmo$mother_dis_subclinical_mastitis),]
table(scm$time_point, useNA="ifany")
# 0.5_months    1_month   2_months   3_months   6_months 
#          0        229        171        174          


# Missings data n/N (%):
# Maternal age 26/524 (5%)
# Gestational age 135/524 (26%)
# Birth weight 146/532 (27%)
# Maternal education 152/524 (29%)
# Delivery mode 50/532 (9%)
# Infant sex 27/532 (5%)
# Feeding mode:
# W2 76/217 (35%)
# M1 138/424 (33%)
# M2 122/345 (35%)
# M3 181/437 (41%)
# M6 38/140 (27%)

