################################################################################
### MILK COMPARIONS ###
################################################################################

##### =========================== 0. DATA IMPORT =========================== #####

library(ggplot2)
library(vegan)
library(ggpubr)
## set working directory
#setwd("/groups/umcg-llnext/tmp01/umcg-jspreckels/Milk_composition_paper/2_HMO_conc_by_milk_group_and_time/")
setwd("~/UMCG/2025_LLnext_Johanne_rebuttal_2/")

## microbiome and phenotypes
# > load file
in16S <- readRDS('milk_16_decontaminated/231129_milk_composition_paper_sel_phenotypes_sel_microbiota_data_genera_n1515.rds')

##### =========================== 1. SELECT DATA OF INTEREST =========================== #####
## select every mother only 1x (i.e. remove duplication of data generated for twin babies)
in16Sf <- in16S[grep(";1", in16S$mother_sample_ID),] #1542x512

#### =========================== 2. ENSURE CORRECT DATA STRUCTURE =========================== #####

##### =========================== 3. ENSURE REFERENCE CATEGORY SHOWS FIRST FOR CATEGORICAL PHENOTYPES  =========================== #####

## for categorical phenotypes, ensure that the order of the levels is as desired (reference category first)
in16Sf$mother_milk_collection_season <- factor(in16Sf$mother_milk_collection_season, levels = c("Spring", "Summer", "Autumn", "Winter"))
in16Sf$mother_milk_collection_month <- factor(in16Sf$mother_milk_collection_month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep",  "Oct",  "Nov", "Dec"))
in16Sf$mother_milk_collection_breasts_sampled <- factor(in16Sf$mother_milk_collection_breasts_sampled, levels=c("one_breast", "both_breasts"))
in16Sf$mother_genetics_blood_group_genotype <- factor(in16Sf$mother_genetics_blood_group_genotype, levels = c("O/O", "A/O", "A/A", "B/O", "B/B", "A/B"))
in16Sf$mother_blood_group <- factor(in16Sf$mother_blood_group, levels = c("O", "A", "B", "AB"))
in16Sf$mother_exp_living_situation <- factor(in16Sf$mother_exp_living_situation, levels = c("big_house", "small_house", "flat", "farm", "other"))
in16Sf$mother_birth_gestational_age_categories <- factor(in16Sf$mother_birth_gestational_age_categories, levels = c("term", "preterm"))
in16Sf$mother_breastpump_brand <- factor(in16Sf$mother_breastpump_brand, levels = c("Philips", "Medela", "Ardo", "Other"))
in16Sf$infant_birth_delivery_mode <- factor(in16Sf$infant_birth_delivery_mode, levels=c("vaginal_birth", "C-section"))
in16Sf$infant_birth_delivery_mode_detailed <- factor(in16Sf$infant_birth_delivery_mode_detailed, levels = c("vaginal_birth_spon_contrac", "vaginal_birth_spon_contrac_vaccum", "vaginal_birth_induction", "vaginal_birth_induction_vaccum",
                                                                                                      "C-section_planned", "C-section_pre_labor_emergency", "C-section_spon_contrac_emergency", "C-section_induction_emergency"))
in16Sf$infant_ffq_feeding_type_delivery <- factor(in16Sf$infant_ffq_feeding_type_delivery, levels = c("breastfeeding", "breast_milk_via_bottle", "mixed_feeding"))
in16Sf$mother_milk_HMO_milk_group <- factor(in16Sf$mother_milk_HMO_milk_group, levels = c("Le+Se+","Le+Se-","Le-Se+","Le-Se-"))


##### =========================== 4. CHECK BASIC INFORMATION / STATISTICS =========================== #####
## number of total measured samples 
length(unique(in16Sf$mother_sample_ID)) #768

## number of samples by time point
table(in16Sf$time_point)
# 1_month 2_months 3_months 6_months 
# 527       88      631      196

## number of mothers
length(unique(in16Sf$mother_ID)) #338 mother-participations

### ============================================================================
##### ======== COMPARE MILK MICROBIOME BETWEEN PREGNANCIES ======= ###
#### BRAY CURTIS !
### Same mother, different pregnancies, same time point
### ===========================================================================
## note: they will always be same milk group, so no need for stratification her
# > get milk only
in16Sf2 <- in16Sf[in16Sf$sample_type %in% c('human_milk'),]
# > find mothers with multiple pregnancies
mmp1 <- in16Sf2$family_ID[in16Sf2$NEXT_participation_number %in% c('first_participation')]
mmp2 <- in16Sf2$family_ID[in16Sf2$NEXT_participation_number %in% c('second_participation')]
mmp12 <- intersect(mmp1,mmp2)
length(mmp12) # we have 12 mothers with measurements in two pregnancies
in16Spp <- in16Sf2[in16Sf2$family_ID %in% mmp12,]
#  >> results DF
results_milk_p1p2_bc <- NULL
# > now go over time-points
for (tp in c('1_month','3_months','6_months')) {
  in16Sppt <- in16Spp[in16Spp$time_point %in% tp,]
  mmp1 <- in16Sppt$family_ID[in16Sppt$NEXT_participation_number %in% c('first_participation')]
  mmp2 <- in16Sppt$family_ID[in16Sppt$NEXT_participation_number %in% c('second_participation')]
  mmp12 <- intersect(mmp1,mmp2)
  in16Sp1p2 <- in16Sppt[in16Sppt$family_ID %in% mmp12]
  toKeep <- as.data.frame(table(in16Sp1p2$family_ID))
  toKeep <- toKeep[toKeep$Freq == 2,]
  print(paste0(tp,' -> ',length(toKeep$Var1),' mothers with 2 pregnancies!'))
  if (length(mmp12) >= 5 ) {
    print(' >> processing ...')
    # now actually process the data
    # > keep samples with Preg 1 and Preg 2 for this Time-point
    in16Sp1p3 <- in16Sp1p2[in16Sp1p2$family_ID %in% toKeep$Var1,]
    # > clean up 16S genera-level microbiome (drop all-0 and low prevalence genera)
    mbGenCols <- grep('^g__',colnames(in16Sp1p3))
    mbGens <- in16Sp1p3[,mbGenCols]
    mbGensSums <- colSums(mbGens)
    mbGensPrev <- apply(mbGens,MARGIN = 2,FUN = function(x) {sum(x>0)/length(x)})
    mbGensF <- mbGens[,mbGensSums > 0 & mbGensPrev > 0.05]
    generaToTest <- colnames(mbGensF)
    print(paste0('   >>> ',ncol(mbGensF),' genera found after filtering'))
    # > separate into first and second pregnancy  
    milk_p1 <- in16Sp1p3[in16Sp1p3$NEXT_participation_number %in% c('first_participation'),]
    milk_p2 <- in16Sp1p3[in16Sp1p3$NEXT_participation_number %in% c('second_participation'),]
    # > sanity check
    if (sum(!(milk_p1$family_ID %in% milk_p2$family_ID)) > 0 | 
        sum(!(milk_p2$family_ID %in% milk_p1$family_ID)) > 0)  {
      print('WARNING: something is wrong - mothers not matching!')
    }
    # >> iterate over mothers
    for (mtr in milk_p1$family_ID) {
      # >> split microbiome of motherX into P1 and P2 for distance calculation
      milk_p1_mtr <- milk_p1[milk_p1$family_ID == mtr,colnames(milk_p1) %in% generaToTest]
      milk_p2_mtr <- milk_p2[milk_p2$family_ID == mtr,colnames(milk_p2) %in% generaToTest]
      # >> distance
      milk_bray <- as.numeric(vegan::vegdist(rbind.data.frame(milk_p1_mtr,milk_p2_mtr),method = 'bray'))
      # >> distance 2
      mlik_aitch <- as.numeric(vegan::vegdist(rbind.data.frame(milk_p1_mtr,milk_p2_mtr),method = 'robust.aitchison'))
      # >> add results
      one_row <- data.frame(Mother=mtr,
                            Milk_BC_Dist=milk_bray,
                            Milk_Dist_aitch = mlik_aitch,
                            Milk_Group=milk_p1$mother_milk_HMO_milk_group[milk_p1$family_ID==mtr],
                            Timepoint=tp)
      results_milk_p1p2_bc <- rbind.data.frame(results_milk_p1p2_bc,one_row)
    }
  } else {
    print(' >> N < 5, dropping this time point!')
  }
}

results_milk_p1p2_bc$Group <- "P2-vs-P1"
ggplot(results_milk_p1p2_bc,aes(y=Milk_BC_Dist,col=Timepoint,x=Timepoint)) + 
  theme_classic() + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group) 
ggplot(results_milk_p1p2_bc,aes(y=Milk_BC_Dist,col=Timepoint,x=Timepoint)) + 
  theme_classic() + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggplot(results_milk_p1p2_bc,aes(y=Milk_Dist_aitch,col=Timepoint,x=Timepoint)) + 
  theme_classic() + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group) 
ggplot(results_milk_p1p2_bc,aes(y=Milk_Dist_aitch,col=Timepoint,x=Timepoint)) + 
  theme_classic() + 
  geom_boxplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

### ===========================================================================
### Between different mothers, 1st pregnancy only
### ===========================================================================
## note: we stratify per milk group (not Milk Group X - vs - Milk Group Y comparisons)
# > find mothers with multiple pregnancies
#mmp1 <- in16Sf2$family_ID[in16Sf2$NEXT_participation_number %in% c('first_participation')]

# now calculate distances between microbiomes
# > loop over time points
results_milk_p1pw_bc <- NULL
results_milk_p1pw_bc_l <- list()
#milk_sp <- in16Sf2[in16Sf2$mother_sample_ID %in% mmp1,]
#milk_sp <- 
mb16S <- in16SppSel[grep('g__',colnames(in16SppSel))]
mb16S <- mb16S[,colSums(mb16S) > 0]
cnames <- colnames(mb16S)
ccc <- 0
# >> time points
for (t in unique(results_milk_p1p2_bc$Timepoint)) {
  print(t)
  # > loop over milk groups
  #for (mg in unique(results_milk_p1p2_bc$Milk_Group)) {
  #  print(mg)
    # > get data (milk group = mg & time point = t)
#    in16SppSel <- in16Spp[in16Spp$mother_milk_HMO_milk_group %in% mg &
#                         in16Spp$time_point %in% t & in16Spp$NEXT_participation_number %in% c('first_participation'),]
    in16SppSel <- in16Sf2[#in16Sf2$mother_milk_HMO_milk_group %in% mg &
                            in16Sf2$time_point %in% t & in16Sf2$NEXT_participation_number %in% c('first_participation'),]
    in16SppSel$family_ID <- as.factor(as.character( in16SppSel$family_ID))
    all_mother_ids <- in16SppSel$family_ID
    # > now do pair-wise comparison between DIFFERENT mothers
    # > iterate over first mother IDs
    for (pn1 in seq(1,length(all_mother_ids)-1) ) {
      # > iterate over second mother IDs
      #print(paste0(m_id1))
      for (pn2 in seq(pn1+1,length(all_mother_ids)) ) {
        m_id1 <- all_mother_ids[pn1]
        m_id2 <- all_mother_ids[pn2]
        # > get mothers and HMO values
        milk_p1 <- in16SppSel[in16SppSel$family_ID == m_id1,colnames(in16SppSel) %in% cnames]
        milk_p2 <- in16SppSel[in16SppSel$family_ID == m_id2,colnames(in16SppSel) %in% cnames]
        # >> bray
        milk_bray <- as.numeric(vegan::vegdist(rbind.data.frame(milk_p1,milk_p2),method = 'bray'))
        milk_aitch <- as.numeric(vegan::vegdist(rbind.data.frame(milk_p1,milk_p2),method = 'robust.aitchison'))
        # >> add data
        one_row <- data.frame(Mother=m_id1,
                              Mother2=m_id2,
                              Milk_BC_Dist=milk_bray,
                              Milk_Dist_aitch = milk_aitch,
                              #Milk_Group=mg,
                              Timepoint=t)
        ccc <- ccc + 1
        results_milk_p1pw_bc_l[[ccc]] <- one_row
        if (ccc %% 1000 == 0) {print(ccc)}
      }
    }
  #}
}
results_milk_p1pw_bc <- do.call(rbind, results_milk_p1pw_bc_l)  
results_milk_p1pw_bc$Group <- "RND.Pairs"
# plot it
ggplot(results_milk_p1pw_bc,aes(y=Milk_BC_Dist,x=Timepoint,col=Timepoint)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + geom_jitter(width = 0.5,height = 0.01,alpha=0.1) + 
  geom_boxplot() +
  facet_wrap(~Milk_Group)
ggplot(results_milk_p1pw_bc,aes(y=Milk_Dist_aitch,x=Timepoint,col=Timepoint)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + geom_jitter(width = 0.5,height = 0.01,alpha=0.1) +
  geom_boxplot() + 
  facet_wrap(~Milk_Group)

# merge with p1 - p2 distances
results_milk_p1p2_bc$Mother2 <- results_milk_p1p2_bc$Mother
results_milk_mrg_bc <- rbind.data.frame(results_milk_p1pw_bc,results_milk_p1p2_bc)
# subset to stuff we have in both
results_milk_mrg_bc_s <- results_milk_mrg_bc[results_milk_mrg_bc$Timepoint %in% c('1_month','3_months') & 
                                                 results_milk_mrg_bc$Milk_Group %in% c('Le+Se+','Le+Se-'),]
# save
saveRDS(results_milk_p1pw_bc,'results_pairwise_bray.RDS')
# plot it
g <- ggplot(results_milk_mrg_bc_s,aes(y=Milk_BC_Dist,col=Group,x=Timepoint)) + 
  #geom_violin(position = position_dodge(1),width=1.75) + 
  geom_boxplot(width = 0.9,position = position_dodge(1)) + 
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group) + 
  theme(legend.position = 'bottom')
print(g)
ggsave(plot=g,filename = 'Milk_microbiome_bray_curtis_comparison_detailed.svg',
       width=6,height=4)

g <- ggplot(results_milk_mrg_bc_s,aes(y=Milk_Dist_aitch,col=Group,x=Timepoint)) + 
  #geom_violin(position = position_dodge(1),width=1.75) + 
  geom_boxplot(width = 0.9,position = position_dodge(1)) + 
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  facet_wrap(~Milk_Group) + 
  theme(legend.position = 'bottom')
print(g)
ggsave(plot=g,filename = 'Milk_microbiome_Aitchinson_comparison_detailed.svg',
       width=6,height=4)

# stats [ per group & time point]
# p1p2all <- results_milk_mrg_bc_s[results_milk_mrg_bc_s$Group == 'P2-vs-P1',]
# for (g in c('Le+Se+','Le+Se-')) {
#   for (t in c('1_month','3_months')) {
#     p1p2 <- results_milk_mrg_bc_s[results_milk_mrg_bc_s$Milk_Group == g & 
#                                    results_milk_mrg_bc_s$Timepoint == t & 
#                                    results_milk_mrg_bc_s$Group == 'P2-vs-P1',]
#     rndp <- results_milk_mrg_bc_s[results_milk_mrg_bc_s$Milk_Group == g & 
#                                    results_milk_mrg_bc_s$Timepoint == t & 
#                                    results_milk_mrg_bc_s$Group == 'RND.Pairs',]
#     rndps <- rndp[sample(c(1:nrow(rndp)),100),] #nrow(p1p2all)*2),]
#     #rndps <- rndp#[sample(c(1:nrow(rndp)),100),]#nrow(p1p2all)*2),]
#     wc <- wilcox.test(p1p2$Milk_Dist_aitch,rndps$Milk_Dist_aitch)
#     print(paste0(' > time-point = ',t,' & group = ',g,'; N[p1-p2] = ',nrow(p1p2),'; p-value = ',wc$p.value))
#     print(paste0('   >> mean[p1-p2] = ',mean(p1p2$Milk_BC_Dist),'; sd[p1-p2] = ',sd(p1p2$Milk_BC_Dist)))
#     print(paste0('   >> mean[rnd.p] = ',mean(rndps$Milk_BC_Dist),'; sd[rnd.p] = ',sd(rndps$Milk_BC_Dist)))
#   }
# }

p1p2all <- results_milk_mrg_bc_s[results_milk_mrg_bc_s$Group == 'P2-vs-P1',]
#for (g in c('Le+Se+','Le+Se-')) {
  for (t in c('1_month','3_months')) {
    p1p2 <- results_milk_mrg_bc_s[
                                    results_milk_mrg_bc_s$Timepoint == t & 
                                    results_milk_mrg_bc_s$Group == 'P2-vs-P1',]
    rndp <- results_milk_mrg_bc_s[
                                    results_milk_mrg_bc_s$Timepoint == t & 
                                    results_milk_mrg_bc_s$Group == 'RND.Pairs',]
    rndps <- rndp[sample(c(1:nrow(rndp)),100),] #nrow(p1p2all)*2),]
    #rndps <- rndp#[sample(c(1:nrow(rndp)),100),]#nrow(p1p2all)*2),]
    wc <- wilcox.test(p1p2$Milk_Dist_aitch,rndps$Milk_Dist_aitch)
    print(paste0(' > time-point = ',t,' & group = ',g,'; N[p1-p2] = ',nrow(p1p2),'; p-value = ',wc$p.value))
    print(paste0('   >> mean[p1-p2] = ',mean(p1p2$Milk_BC_Dist),'; sd[p1-p2] = ',sd(p1p2$Milk_BC_Dist)))
    print(paste0('   >> mean[rnd.p] = ',mean(rndps$Milk_BC_Dist),'; sd[rnd.p] = ',sd(rndps$Milk_BC_Dist)))
  }
#}

# stats [total]
p1p2 <- results_milk_mrg_bc_s[results_milk_mrg_bc_s$Group == 'P2-vs-P1',]
print(paste0('P1-P2: mean = ',mean(p1p2$Milk_Dist_aitch),'; sd = ',sd(p1p2$Milk_Dist_aitch)))
rndp <- results_milk_mrg_bc_s[results_milk_mrg_bc_s$Group == 'RND.Pairs',]
print(paste0('RND.P: mean = ',mean(rndp$Milk_Dist_aitch),'; sd = ',sd(rndp$Milk_Dist_aitch)))
# summaries
summary(p1p2$Milk_Dist_aitch); sd(p1p2$Milk_Dist_aitch)
summary(rndp$Milk_Dist_aitch); sd(rndp$Milk_Dist_aitch)
# bootstrapping p-value ranges for selecting only 12 pairs
bootp <- c()
for (boot in c(1:1000)) {
  rndps <- rndp[sample(c(1:nrow(rndp)),nrow(p1p2)),]
  wc <- wilcox.test(p1p2$Milk_Dist_aitch,rndps$Milk_Dist_aitch)
  bootp <- c(bootp,wc$p.value)
}
summary(bootp)
sd(bootp)
mean(bootp)

# plot 2 [simple variant]
g <- ggplot(results_milk_mrg_bc_s,aes(y=Milk_Dist_aitch,col=Group,x=Group)) + 
  geom_violin(position = position_dodge(1),width=1.25) + 
  geom_boxplot(width = 0.35,position = position_dodge(1),alpha=0.1) + 
  theme_classic2() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
  #ylim(c(0,0.9)) + 
  theme(legend.position = 'none')

tst_forplot <- compare_means(Milk_Dist_aitch ~ Group, 
                             data = results_milk_mrg_bc_s,
                             method = "wilcox.test")
tst_forplot$p <- mean(bootp); tst_forplot$p.adj <- formatC(mean(bootp), format = "e", digits = 2)
tst_forplot$y.position <- 15

g <- g + stat_pvalue_manual(tst_forplot)
print(g)
ggsave(plot=g,filename = 'Milk_bray_curtis_comparison_premutation_test.svg',width=4,height=6)


