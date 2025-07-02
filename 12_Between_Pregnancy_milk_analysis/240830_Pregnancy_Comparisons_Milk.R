################################################################################
### MILK COMPARIONS ###
################################################################################

##### =========================== 0. DATA IMPORT =========================== #####

library(ggplot2)
library(vegan)
library(ggpubr)
## set working directory
#setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/2_HMO_conc_by_milk_group_and_time/")
setwd("~/UMCG/2024_LLnext_Johanne_rebuttal/")

## import file
h <- read.table("240108_milk_composition_paper_sel_phenotypes_hmo_data_n1563.txt", header=T, sep="\t", stringsAsFactors=T) #1563x512

# Notes:
# This file contains the following:
# - Phenotypes linked to measured HMO data.
#   !! Note that maternal phenotype and HMO data was duplicated for twins. This is indicated by ";1" (single baby/first baby) and ";2" (twin baby/duplicated data) in the mother_ID and mother_sample_ID.
#      -> For association with maternal phenotypes (i.e. non-duplicated data), grep only ";1"
#      -> For association with infant phenotypes, use both ";1" and ";2"
# - Real measured levels of 24 single HMOs and 4 grouped HMOs in µg/ml and HMO-based Le and Se status and milk groups
# - A quality control sample was included in all HMO batches and if there was less than <15% variation, the quality of the run was considered good. No need to correct for a batch effect from HMO measurements.


##### =========================== 1. SELECT DATA OF INTEREST =========================== #####
## select every mother only 1x (i.e. remove duplication of data generated for twin babies)
hmo <- h[grep(";1", h$mother_sample_ID),] #1542x512

##### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

## check that data structure for relevant columns is correct
str(hmo[,1:28])
str(hmo[,481:512])
# -> data structure in real HMO data file is correct.

##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

## for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
hmo$mother_milk_collection_season <- factor(hmo$mother_milk_collection_season, levels = c("Spring", "Summer", "Autumn", "Winter"))
hmo$mother_milk_collection_month <- factor(hmo$mother_milk_collection_month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",  "Oct",  "Nov", "Dec"))
hmo$mother_milk_collection_notes <- factor(hmo$mother_milk_collection_notes, levels = c("No_records_of_incorrect_sample_handling", levels(hmo$mother_milk_collection_notes)[c(3,5,7,6,2,1)]))
hmo$mother_milk_collection_breasts_sampled <- factor(hmo$mother_milk_collection_breasts_sampled, levels=c("one_breast", "both_breasts"))
hmo$mother_genetics_blood_group_genotype <- factor(hmo$mother_genetics_blood_group_genotype, levels = c("O/O", "A/O", "A/A", "B/O", "B/B", "A/B"))
hmo$mother_blood_group <- factor(hmo$mother_blood_group, levels = c("O", "A", "B", "AB"))
hmo$mother_exp_living_situation <- factor(hmo$mother_exp_living_situation, levels = c("big_house", "small_house", "flat", "farm", "other"))
hmo$mother_birth_gestational_age_categories <- factor(hmo$mother_birth_gestational_age_categories, levels = c("term", "preterm"))
hmo$mother_breastpump_brand <- factor(hmo$mother_breastpump_brand, levels = c("Philips", "Medela", "Ardo", "Other"))
hmo$infant_birth_delivery_mode <- factor(hmo$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
hmo$infant_birth_delivery_mode_detailed <- factor(hmo$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                      "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
hmo$infant_ffq_feeding_type_delivery <- factor(hmo$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))
hmo$mother_milk_HMO_milk_group <- factor(hmo$mother_milk_HMO_milk_group, levels = c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 4. CHECK BASIC INFORMATION / STATISTICS =========================== #####
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
hmo <- hmo[,c(1:6,513,7:512)]
length(unique(hmo$mother_ID_simple)) #500 -> 524-500=24 -> 24 mothers have HMO measurements for first and second NEXT pregnancy

### ============================================================================
##### ======== COMPARE HMOs BETWEEN PREGNANCIES ======= ###
#### BRAY CURTIS !
### Same mother, different pregnancies, same time point
### ===========================================================================
## note: they will always be same milk group, so no need for stratification here
# > find mothers with multiple pregnancies
mmp1 <- hmo$family_ID[hmo$NEXT_participation_number %in% c('first_participation')]
mmp2 <- hmo$family_ID[hmo$NEXT_participation_number %in% c('second_participation')]
mmp12 <- intersect(mmp1,mmp2)
length(mmp12) # we should have 24 mothers with measurements in two pregnancies

# >> grab hmo names
cnames <- colnames(hmo)[grep('mother_milk_.*_ugml',colnames(hmo))]
# > drop summaries (total / Fuc / Netu)
cnames <- cnames[1:25]

# now calculate Bray-Curtis distance
# > loop over time points
results_hmo_p1p2_bc <- NULL
hmo_mp <- hmo[hmo$family_ID %in% mmp12,]
for (t in levels(hmo$time_point)) {
  print(t)
  # > get first pregnancy
  hmo_p1 <- hmo_mp[hmo_mp$NEXT_participation_number %in% c('first_participation') & 
                     hmo_mp$time_point %in% t,]
  # > get second pregnancy
  hmo_p2 <- hmo_mp[hmo_mp$NEXT_participation_number %in% c('second_participation') & 
                     hmo_mp$time_point %in% t,]
  # > further subset to keep mothers who have measurement @ timepoint = t
  hmo_p12_t <-intersect(hmo_p1$family_ID,hmo_p2$family_ID)
  # >> only do comparison if we have at least 3 samples
  if (length(hmo_p12_t) < 3) {
    print(paste0(' > W: NOT ENOUGH (',length(hmo_p12_t),') mothers with p1 & p2 measurements @ ',t))
  } else {
    # >> do comparisons; prepare data
    print(paste0(' > ',length(hmo_p12_t),' mothers with p1 & p2 measurements @ ',t))
    hmo_p1 <- hmo_p1[hmo_p1$family_ID %in% hmo_p12_t,]
    hmo_p2 <- hmo_p2[hmo_p2$family_ID %in% hmo_p12_t,]
    # >> get all HMOs
    cnames <- colnames(hmo_p1)[grep('mother_milk_.*_ugml',colnames(hmo_p1))]
    # >> iterate over mothers
    for (mtr in hmo_p12_t) {
      hmo_v_p1 <- hmo_p1[hmo_p1$family_ID == mtr,colnames(hmo_p1) %in% cnames]
      hmo_v_p2 <- hmo_p2[hmo_p1$family_ID == mtr,colnames(hmo_p2) %in% cnames]
      # >> distance
      hmo_bray <- as.numeric(vegan::vegdist(rbind.data.frame(hmo_v_p1,hmo_v_p2),method = 'bray'))
      # >> distance 2
      hmo_aitch <- as.numeric(vegan::vegdist(rbind.data.frame(hmo_v_p1,hmo_v_p2),method = 'robust.aitchison'))
      # >> add results
      one_row <- data.frame(Mother=mtr,
                            HMO_BC_Dist=hmo_bray,
                            HMO_Dist_aitch = hmo_aitch,
                            Milk_Group=hmo_p1$mother_milk_HMO_milk_group[hmo_p1$family_ID==mtr],
                            Timepoint=t)
      results_hmo_p1p2_bc <- rbind.data.frame(results_hmo_p1p2_bc,one_row)
    }
  }
}
results_hmo_p1p2_bc$Group <- "P2-vs-P1"
ggplot(results_hmo_p1p2_bc,aes(y=HMO_BC_Dist,col=Timepoint,x=Timepoint)) + 
  theme_classic() + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group) 
ggplot(results_hmo_p1p2_bc,aes(y=HMO_Dist_aitch,col=Timepoint,x=Timepoint)) + 
  theme_classic() + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group) 

### Between different mothers, 1st pregnancy only
### ===========================================================================
## note: we stratify per milk group (not Milk Group X - vs - Milk Group Y comparisons)
# > find mothers with multiple pregnancies
mmp1 <- hmo$mother_sample_ID[hmo$NEXT_participation_number %in% c('first_participation')]
# now calculate distances in each HMO
# > loop over time points
results_hmo_p1pw_bc <- NULL
results_hmo_p1pw_bc_l <- list()
hmo_sp <- hmo[hmo$mother_sample_ID %in% mmp1,]
ccc <- 0
# >> time points
for (t in levels(hmo_sp$time_point)) {
  print(t)
  # > loop over milk groups
  for (mg in levels(hmo_sp$mother_milk_HMO_milk_group)) {
    print(mg)
    # > get data (milk group = mg & time point = t)
    hmo_mg <- hmo_sp[hmo_sp$mother_milk_HMO_milk_group %in% mg &
                       hmo_sp$time_point %in% t,]
    hmo_mg$family_ID <- as.factor(as.character( hmo_mg$family_ID))
    all_mother_ids <- hmo_mg$family_ID
    # > now do pair-wise comparison between DIFFERENT mothers
    # > iterate over first mother IDs
    for (pn1 in seq(1,length(all_mother_ids)-1) ) {
      # > iterate over second mother IDs
      #print(paste0(m_id1))
      for (pn2 in seq(pn1+1,length(all_mother_ids)) ) {
        m_id1 <- all_mother_ids[pn1]
        m_id2 <- all_mother_ids[pn2]
        # > get mothers and HMO values
        hmo_p1 <- hmo_mg[hmo_mg$family_ID == m_id1,colnames(hmo_mg) %in% cnames]
        hmo_p2 <- hmo_mg[hmo_mg$family_ID == m_id2,colnames(hmo_mg) %in% cnames]
        # >> bray
        hmo_bray <- as.numeric(vegan::vegdist(rbind.data.frame(hmo_p1,hmo_p2),method = 'bray'))
        hmo_aitch <- as.numeric(vegan::vegdist(rbind.data.frame(hmo_p1,hmo_p2),method = 'robust.aitchison'))
        # >> add data
        one_row <- data.frame(Mother=m_id1,
                              Mother2=m_id2,
                              HMO_BC_Dist=hmo_bray,
                              HMO_Dist_aitch = hmo_aitch,
                              Milk_Group=mg,
                              Timepoint=t)
        ccc <- ccc + 1
        results_hmo_p1pw_bc_l[[ccc]] <- one_row
      }
    }
  }
}
results_hmo_p1pw_bc <- do.call(rbind, results_hmo_p1pw_bc_l)  
results_hmo_p1pw_bc$Group <- "RND.Pairs"
# plot it
ggplot(results_hmo_p1pw_bc,aes(y=HMO_BC_Dist,col=Timepoint)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group)
ggplot(results_hmo_p1pw_bc,aes(y=HMO_Dist_aitch,col=Timepoint)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group)

# merge with p1 - p2 distances
results_hmo_p1p2_bc$Mother2 <- results_hmo_p1p2_bc$Mother
results_hmo_mrg_bc <- rbind.data.frame(results_hmo_p1pw_bc,results_hmo_p1p2_bc)
# subset to stuff we have in both
results_hmo_mrg_bc_s <- results_hmo_mrg_bc[results_hmo_mrg_bc$Timepoint %in% c('1_month','2_months','3_months','6_months') & 
                                             results_hmo_mrg_bc$Milk_Group %in% c('Le+Se+','Le+Se-'),]
# plot it
g <- ggplot(results_hmo_mrg_bc_s,aes(y=HMO_BC_Dist,col=Group,x=Timepoint)) + 
  #geom_violin(position = position_dodge(1),width=1.75) + 
  geom_boxplot(width = 0.9,position = position_dodge(1)) + 
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group) + 
  theme(legend.position = 'bottom')
print(g)
ggsave(plot=g,filename = 'Milk_bray_curtis_comparison_detailed.svg',
       width=6,height=4)

g <- ggplot(results_hmo_mrg_bc_s,aes(y=HMO_Dist_aitch,col=Group,x=Timepoint)) + 
  #geom_violin(position = position_dodge(1),width=1.75) + 
  geom_boxplot(width = 0.9,position = position_dodge(1)) + 
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group) + 
  theme(legend.position = 'bottom')
print(g)
ggsave(plot=g,filename = 'Milk_Aitchinson_comparison_detailed.svg',
       width=6,height=4)

# stats [ per group & time point]
p1p2all <- results_hmo_mrg_bc_s[results_hmo_mrg_bc_s$Group == 'P2-vs-P1',]
for (g in c('Le+Se+','Le+Se-')) {
  for (t in c('1_month','2_months','3_months','6_months')) {
    p1p2 <- results_hmo_mrg_bc_s[results_hmo_mrg_bc_s$Milk_Group == g & 
                                   results_hmo_mrg_bc_s$Timepoint == t & 
                                   results_hmo_mrg_bc_s$Group == 'P2-vs-P1',]
    rndp <- results_hmo_mrg_bc_s[results_hmo_mrg_bc_s$Milk_Group == g & 
                                   results_hmo_mrg_bc_s$Timepoint == t & 
                                   results_hmo_mrg_bc_s$Group == 'RND.Pairs',]
    rndps <- rndp[sample(c(1:nrow(rndp)),100),]#nrow(p1p2all)*2),]
    wc <- wilcox.test(p1p2$HMO_BC_Dist,rndps$HMO_BC_Dist)
    print(paste0(' > time-point = ',t,' & group = ',g,'; N[p1-p2] = ',nrow(p1p2),'; p-value = ',wc$p.value))
    print(paste0('   >> mean[p1-p2] = ',mean(p1p2$HMO_BC_Dist),'; sd[p1-p2] = ',sd(p1p2$HMO_BC_Dist)))
    print(paste0('   >> mean[rnd.p] = ',mean(rndps$HMO_BC_Dist),'; sd[rnd.p] = ',sd(rndps$HMO_BC_Dist)))
  }
}

# stats [total]
p1p2 <- results_hmo_mrg_bc_s[results_hmo_mrg_bc_s$Group == 'P2-vs-P1',]
print(paste0('P1-P2: mean = ',mean(p1p2$HMO_BC_Dist),'; sd = ',sd(p1p2$HMO_BC_Dist)))
rndp <- results_hmo_mrg_bc_s[results_hmo_mrg_bc_s$Group == 'RND.Pairs',]
print(paste0('RND.P: mean = ',mean(rndp$HMO_BC_Dist),'; sd = ',sd(rndp$HMO_BC_Dist)))
# summaries
summary(p1p2$HMO_BC_Dist); sd(p1p2$HMO_BC_Dist)
summary(rndp$HMO_BC_Dist); sd(rndp$HMO_BC_Dist)
# bootstrapping p-value ranges for selecting only 55 pairs
bootp <- c()
for (boot in c(1:1000)) {
  rndps <- rndp[sample(c(1:nrow(rndp)),nrow(p1p2)),]
  wc <- wilcox.test(p1p2$HMO_Dist_aitch,rndps$HMO_Dist_aitch)
  bootp <- c(bootp,wc$p.value)
}
summary(bootp)
sd(bootp)
mean(bootp)
# p-value: 1.57 * 10^-6 +/- 1.17*10-5
# "P1-P2: mean = 0.103035562342203; sd = 0.0698541451849081"
# "RND.P: mean = 0.172936499201888; sd = 0.0797357467364295"
# plot 2 [simple variant]
g <- ggplot(results_hmo_mrg_bc_s,aes(y=HMO_BC_Dist,col=Group,x=Group)) + 
  geom_violin(position = position_dodge(1),width=1.25) + 
  geom_boxplot(width = 0.35,position = position_dodge(1),alpha=0.1) + 
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  ylim(c(0,0.9)) + 
  theme(legend.position = 'none')

g <- ggplot(results_hmo_mrg_bc_s,aes(y=HMO_Dist_aitch,col=Group,x=Group)) + 
  geom_violin(position = position_dodge(1),width=1.25) + 
  geom_boxplot(width = 0.35,position = position_dodge(1),alpha=0.1) + 
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  #ylim(c(0,0.9)) + 
  theme(legend.position = 'none')

tst_forplot <- compare_means(HMO_BC_Dist ~ Group, 
                             data = results_hmo_mrg_bc_s,
                             method = "wilcox.test")
tst_forplot$p <- mean(bootp); tst_forplot$p.adj <- formatC(mean(bootp), format = "e", digits = 2)
tst_forplot$y.position <- 0.80

g <- g + stat_pvalue_manual(tst_forplot)
print(g)
ggsave(plot=g,filename = 'Milk_bray_curtis_comparison.svg',width=4,height=6)

### ===========================================================================
##### ======== COMPARE HMOs BETWEEN PREGNANCIES ======= ###
#### PER HMO:
### Same mother, different pregnancies, same time point
### ===========================================================================
## note: they will always be same milk group, so no need for stratification here
# > find mothers with multiple pregnancies
mmp1 <- hmo$family_ID[hmo$NEXT_participation_number %in% c('first_participation')]
mmp2 <- hmo$family_ID[hmo$NEXT_participation_number %in% c('second_participation')]
mmp12 <- intersect(mmp1,mmp2)
length(mmp12) # we should have 24 mothers with measurements in two pregnancies

# what to test:
cnames <- colnames(hmo)[grep('mother_milk_.*_ugml',colnames(hmo))]
# > drop summaries (total / Fuc / Netu)
cnames <- cnames[1:25]

# now calculate distances in each HMO
# > loop over time points
results_hmo_p1p2 <- NULL
hmo_mp <- hmo[hmo$family_ID %in% mmp12,]
for (t in levels(hmo$time_point)) {
  print(t)
  # > get first pregnancy
  hmo_p1 <- hmo_mp[hmo_mp$NEXT_participation_number %in% c('first_participation') & 
                     hmo_mp$time_point %in% t,]
  # > get second pregnancy
  hmo_p2 <- hmo_mp[hmo_mp$NEXT_participation_number %in% c('second_participation') & 
                     hmo_mp$time_point %in% t,]
  # > further subset to keep mothers who have measurement @ timepoint = t
  hmo_p12_t <-intersect(hmo_p1$family_ID,hmo_p2$family_ID)
  # >> only do comparison if we have at least 3 samples
  if (length(hmo_p12_t) < 3) {
    print(paste0(' > W: NOT ENOUGH (',length(hmo_p12_t),') mothers with p1 & p2 measurements @ ',t))
  } else {
    # >> do comparisons; prepare data
    print(paste0(' > ',length(hmo_p12_t),' mothers with p1 & p2 measurements @ ',t))
    hmo_p1 <- hmo_p1[hmo_p1$family_ID %in% hmo_p12_t,]
    hmo_p2 <- hmo_p2[hmo_p2$family_ID %in% hmo_p12_t,]
    # >> iterate over HMOs
    for (hmo_cname in cnames) {
      print(hmo_cname)
      # >> iterate over mothers
      for (mtr in hmo_p12_t) {
        hmo_v_p1 <- hmo_p1[[hmo_cname]][hmo_p1$family_ID==mtr]
        hmo_v_p2 <- hmo_p2[[hmo_cname]][hmo_p2$family_ID==mtr]
        hmo_delta <- abs(hmo_v_p1 - hmo_v_p2)
        hmo_ratio <- hmo_v_p1 / hmo_v_p2
        one_row <- data.frame(Mother=mtr,HMO=hmo_cname,HMO_P1=hmo_v_p1,HMO_P2=hmo_v_p2,
                              HMO_delta=hmo_delta,
                              MILK_group=hmo_p1$mother_milk_HMO_milk_group[hmo_p1$family_ID==mtr],
                              HMO_ratio=hmo_ratio,SETUP='P1_vs_P2',Timepoint=t)
        results_hmo_p1p2 <- rbind.data.frame(results_hmo_p1p2,one_row)
      }
    }
  }
}
results_hmo_p1p2$HMO <- gsub('mother_milk_HMO_','',results_hmo_p1p2$HMO)
results_hmo_p1p2$HMO <- gsub('_ugml','',results_hmo_p1p2$HMO)
results_hmo_p1p2$HMO_ratio[is.nan(results_hmo_p1p2$HMO_ratio)] <- NA
results_hmo_p1p2$HMO_ratio[is.infinite(results_hmo_p1p2$HMO_ratio)] <- NA

ggplot(results_hmo_p1p2,aes(x=HMO,y=log(HMO_ratio),col=Timepoint)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

g <- ggplot(results_hmo_p1p2,aes(x=HMO,y=log(HMO_ratio),col=HMO)) +
  theme_classic() +
  geom_boxplot() +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave(plot = g,filename = 'Milk_HMOs_change_P2vsP1.svg',width = 6,height=4)

# test all VS null [no change]
tstHMO <- NULL
for (hmoN in unique(results_hmo_p1p2$HMO)) {
  totest <- log(results_hmo_p1p2$HMO_ratio[results_hmo_p1p2$HMO == hmoN])
  wc <- wilcox.test(totest,mu = 0,conf.level = 0.95)
  print(paste0(hmoN,' == 0 test; p-value: ',wc$p.value))
  tstHMO <- rbind.data.frame(tstHMO,data.frame(HMO=hmoN,pV=wc$p.value))
}
tstHMO$FDR <- p.adjust(tstHMO$pV)

### alternative test: paired test per HMO
tstHMOpaired <- NULL
tstHMOcor <- NULL
for (hmoN in unique(results_hmo_p1p2$HMO)) {
  totest <- results_hmo_p1p2[results_hmo_p1p2$HMO == hmoN,]
  wc <- wilcox.test((totest$HMO_P1),(totest$HMO_P2),paired = T)
  print(paste0(hmoN,' paired test VS Ho: mean1=mean2; p-value: ',wc$p.value))
  tstHMOpaired <- rbind.data.frame(tstHMOpaired,data.frame(HMO=hmoN,pV=wc$p.value))
  ct <- cor.test((totest$HMO_P1),(totest$HMO_P2))
  tstHMOcor <- rbind.data.frame(tstHMOcor,data.frame(HMO=hmoN,pV=ct$p.value,
                                                     cor=ct$estimate,
                                                     cor.ci.low=ct$conf.int[1],
                                                     cor.ci.high=ct$conf.int[2],
                                                     group='Le+Se+'))
#  ct$p.value
}
tstHMOpaired$FDR <- p.adjust(tstHMOpaired$pV)
tstHMOcor$FDR <- p.adjust(tstHMOcor$pV)
tstHMOcor$setup='p1_vs_p2'
write.table(tstHMOcor,file = 'Milk_correlations_P1P2.csv',sep=',',row.names = F)

### Correlation of HMOs between random mothers and P1P2 mothers
# ================================================================
## note: we stratify per milk group (not Milk Group X - vs - Milk Group Y comparisons)
# > find mothers with multiple pregnancies
tstHMOcorbm <- NULL
mmp1 <- hmo$mother_sample_ID[hmo$NEXT_participation_number %in% c('first_participation')]
# now calculate distances in each HMO
# > loop over time points
hmo_sp <- hmo[hmo$mother_sample_ID %in% mmp1,]
ccc <- 0
# >> iterate over HMOs
# >> grab hmo names
cnames <- colnames(hmo)[grep('mother_milk_.*_ugml',colnames(hmo))]
# > drop summaries (total / Fuc / Netu)
cnames <- cnames[1:25]
for (hmo_cname in cnames) {
  #for (t in levels(hmo_sp$time_point)) {
  #print(t)
  # > loop over milk groups
  #for (mg in levels(hmo_sp$mother_milk_HMO_milk_group)) {
  for (mg in c('Le+Se+','Le+Se-')) {
    print(mg)
    # > get data (milk group = mg & time point = t)
    hmo_mg <- hmo_sp[hmo_sp$mother_milk_HMO_milk_group %in% mg,]
    #& hmo_sp$time_point %in% t,]
    hmo_mg$family_ID <- as.factor(as.character( hmo_mg$family_ID))
    all_mother_ids <- hmo_mg$family_ID
    # > get data [p1p2 mothers]
    hmo_mp1p2 <- hmo_mp[hmo_mp$mother_milk_HMO_milk_group %in% mg,]
    #& hmo_mp$time_point %in% t,]
    print(dim(hmo_mp1p2))
    
    hmo_p1p2 <- hmo_mp1p2[[hmo_cname]]
    hmo_other <- sample(hmo_mg[[hmo_cname]],size = length(hmo_p1p2))
    # > test
    ct <- cor.test(hmo_p1p2,hmo_other)
    tstHMOcorbm <- rbind.data.frame(tstHMOcorbm,data.frame(HMO=hmo_cname,pV=ct$p.value,
                                                           cor=ct$estimate,
                                                           group=mg,
                                                           setup='p1p2_vs_rnd_mothers', 
                                                           cor.ci.low=ct$conf.int[1],
                                                           cor.ci.high=ct$conf.int[2]))
  }
}
#}
# take average in case of two groups
library(tidyverse)
tstHMOcorbmr <- tstHMOcorbm %>% group_by(HMO) %>% summarise(cor = mean(cor),
                                                            )
# merge
tstHMOcorbm$FDR <- p.adjust(tstHMOcorbm$pV)
tstHMOcorBoth <- rbind(tstHMOcorbm[tstHMOcorbm$group=='Le+Se+',],tstHMOcor)
tstHMOcorBoth$HMO <- gsub('mother_milk_HMO_','',tstHMOcorBoth$HMO)
tstHMOcorBoth$HMO <- gsub('_ugml','',tstHMOcorBoth$HMO)
tstHMOcorBoth <- tstHMOcorBoth[tstHMOcorBoth$HMO != 'Total',]
tstHMOcorBoth$HMO <- factor(tstHMOcorBoth$HMO)

tsdf <- tstHMOcorBoth[tstHMOcorBoth$setup=='p1_vs_p2',]
tsdf <- tsdf[order(tsdf$cor),]
tsdf$HMO <- factor(tsdf$HMO,levels=tsdf$HMO)

tstHMOcorBoth$HMO <- factor(tstHMOcorBoth$HMO,levels=levels(tsdf$HMO))

g <- ggplot(tstHMOcorBoth,aes(x=setup,y=cor,col=setup)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  theme_classic() + theme(legend.position = 'top')

ggsave(plot = g,filename = 'Milk_HMO_correlation_p1p2_vs_rnd.svg',width = 4,height = 4)

g <- ggplot(tstHMOcorBoth,aes(x=HMO,y=cor,col=setup)) + 
  geom_point(position = position_dodge(0.5)) + 
  geom_errorbar(ymax=tstHMOcorBoth$cor.ci.high,
                ymin=tstHMOcorBoth$cor.ci.low,
                position = position_dodge(0.5)) + 
  theme_classic() + theme(legend.position = 'top') + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(plot = g,filename = 'Milk_HMO_correlation_p1p2_vs_rnd_perHMO.svg',width = 6,height = 4)

tstHMOcorBoth$FDR <- p.adjust(tstHMOcorBoth$pV)

### Between different mothers, 1st pregnancy only
### ===========================================================================
## note: we stratify per milk group (not Milk Group X - vs - Milk Group Y comparisons)
# > find mothers with multiple pregnancies
mmp1 <- hmo$mother_sample_ID[hmo$NEXT_participation_number %in% c('first_participation')]
# now calculate distances in each HMO
# > loop over time points
results_hmo_p1pw <- NULL
results_hmo_p1pw_l <- list()
hmo_sp <- hmo[hmo$mother_sample_ID %in% mmp1,]
ccc <- 0
# >> iterate over HMOs
cnames <- colnames(hmo_p1)[grep('mother_milk_.*_ugml',colnames(hmo_p1))]
for (hmo_cname in cnames[1:2]) {
  for (t in levels(hmo_sp$time_point)[1:2]) {
    print(t)
    # > loop over milk groups
    for (mg in levels(hmo_sp$mother_milk_HMO_milk_group)) {
      print(mg)
      # > get data (milk group = mg & time point = t)
      hmo_mg <- hmo_sp[hmo_sp$mother_milk_HMO_milk_group %in% mg &
                         hmo_sp$time_point %in% t,]
      hmo_mg$family_ID <- as.factor(as.character( hmo_mg$family_ID))
      all_mother_ids <- hmo_mg$family_ID
      # > now do pair-wise comparison between DIFFERENT mothers
      # > iterate over first mother IDs
      for (pn1 in seq(1,length(all_mother_ids)-1) ) {
        # > iterate over second mother IDs
        print(paste0(m_id1))
        for (pn2 in seq(pn1+1,length(all_mother_ids)) ) {
          m_id1 <- all_mother_ids[pn1]
          m_id2 <- all_mother_ids[pn2]
          # > get mothers
          hmo_p1 <- hmo_mg[hmo_mg$family_ID == m_id1,]
          hmo_p2 <- hmo_mg[hmo_mg$family_ID == m_id2,]
          
          #print(hmo_cname)
          # >> iterate over mothers
          hmo_v_p1 <- max(hmo_p1[[hmo_cname]],hmo_p2[[hmo_cname]])
          hmo_v_p2 <- min(hmo_p1[[hmo_cname]],hmo_p2[[hmo_cname]])
          hmo_delta <- sqrt((hmo_v_p1 - hmo_v_p2)^2)
          hmo_ratio <- hmo_v_p1 / hmo_v_p2
          # >>> reality check: milk groups
          if (hmo_p1$mother_milk_HMO_milk_group != hmo_p2$mother_milk_HMO_milk_group) {
            print('ERROR: milk groups got mixed up!')
            stop()
          }
          # add data
          one_row <- data.frame(Mother1=m_id1,
                                Mother2=m_id2,
                                HMO=hmo_cname,
                                HMO_P1=hmo_v_p1,
                                HMO_P2=hmo_v_p2,
                                HMO_delta=hmo_delta,
                                MILK_group=hmo_p1$mother_milk_HMO_milk_group,
                                HMO_ratio=hmo_ratio,
                                SETUP='P1_vs_P2',
                                Timepoint=t)
          ccc <- ccc + 1
          results_hmo_p1pw_l[[ccc]] <- one_row
          #results_hmo_p1pw <- rbind.data.frame(results_hmo_p1pw,one_row)
        }
      }
    }
  }
}  

results_hmo_p1pw <- do.call(rbind, results_hmo_p1pw_l)  
results_hmo_p1pw$HMO <- gsub('mother_milk_HMO_','',results_hmo_p1pw$HMO)
results_hmo_p1pw$HMO <- gsub('_ugml','',results_hmo_p1pw$HMO)
results_hmo_p1pw$HMO_ratio[is.nan(results_hmo_p1pw$HMO_ratio)] <- NA
results_hmo_p1pw$HMO_ratio[is.infinite(results_hmo_p1pw$HMO_ratio)] <- NA

ggplot(results_hmo_p1pw,aes(x=HMO,y=log(HMO_ratio),col=Timepoint)) + geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))


