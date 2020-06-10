##################################
# Alpha diversity - calculated through qiime2

## calculate alpha diversity through qiime2
observed_otus_vector_noFB_noFC_noRmBa003<-read_qza("observed_otus_vector_noFB_noFC_noRmBa003.qza")
write.csv(observed_otus_vector_noFB_noFC_noRmBa003$data, file = "observed_otus_vector_noFB_noFC_noRmBa003.csv")

## do not use Chao1 because it is significantly affected by singletons.
#Chao1_vector_noFB_noFC_noRmBa003<-read_qza("Chao1_vector_noFB-noFC_noRmBa003.qza")
#write.csv(Chao1_vector_noFB_noFC_noRmBa003$data, file = "Chao1_vector_noFB_noFC_noRmBa003.csv")

evenness_vector_noFB_noFC_noRmBa003<-read_qza("evenness_vector_noFB_noFC_noRmBa003.qza")
write.csv(evenness_vector_noFB_noFC_noRmBa003$data, file = "evenness_vector_noFB_noFC_noRmBa003.csv")

shannon_vector_noFB_noFC_noRmBa003<-read_qza("shannon_vector_noFB_noFC_noRmBa003.qza")
write.csv(shannon_vector_noFB_noFC_noRmBa003$data, file = "shannon_vector_noFB_noFC_noRmBa003.csv")

faith_pd_vector_noFB_noFC_noRmBa003<-read_qza("faith_pd_vector_noFB_noFC_noRmBa003.qza")
write.csv(faith_pd_vector_noFB_noFC_noRmBa003$data, file = "faith_pd_vector_noFB_noFC_noRmBa003.csv")

## add alpha diversity results to metadata file. 
alpha_wtPosC_wtNegC_wtBa<-read.csv("AlphaDiversity_withMeta_wtPosC_wtNegC_wtBa.csv", sep =',', header=T)
alpha_noPosC_noNegC_noBa<-read.csv("AlphaDiversity_withMeta_noPosC_noNegC_noBa.csv", sep =',', header=T)

##check data normality 
library(RVAideMemoire)
mshapiro.test(alpha_noPosC_noNegC_noBa$Shannon_index)
#W = 0.96403, p-value = 0.03743
mshapiro.test(alpha_noPosC_noNegC_noBa$pielou_evenness)
#W = 0.96906, p-value = 0.07314
mshapiro.test(alpha_noPosC_noNegC_noBa$observed_otus)
# W = 0.93, p-value = 0.000616
mshapiro.test(alpha_noPosC_noNegC_noBa$faith_pd)
# W = 0.60432, p-value = 1.244e-12
mshapiro.test(cbind(alpha_noPosC_noNegC_noBa$Shannon_index,alpha_noPosC_noNegC_noBa$pielou_evenness,alpha_noPosC_noNegC_noBa$observed_otus,alpha_noPosC_noNegC_noBa$faith_pd))
# W = 0.63204, p-value = 3.803e-12

##Should we try both parametric test ("t.test") to compare means 
##and nonparametric test ("wilcox.test") to compare medians for alpha diversity.
#The choosing of parametric test and nonparametric test based on the distribution of our data. 
#Parametric test requires data to be normally distributed, I need to check the histogram of my data distribution of the alpha index of the samples.
#check data normality through box-plots. 
#If your data comes from a normal distribution, the box will be symmetrical with the median in the center.
#Shapiro-Wilk’s method is widely recommended for normality test and it provides better power than K-S. 
#Shapiro-Wilk is based on the correlation between the data and the corresponding normal scores.
#The middle line of alpha box plots were medians. 
# make a pairwise list that we want to compare. (1) compare to large intestine,(2) compare to large intestine
my_comparisons <- list(c("large_intestine","cloacal"), c("large_intestine", "feces"), c("large_intestine", "small_intestine"))
my_comparisons_b <- list(c("small_intestine","cloacal"), c("small_intestine", "feces"), c("small_intestine", "large_intestine"))

a1 <- ggboxplot(alpha_noPosC_noNegC_noBa, x = "SampleTypes", y = "observed_otus",
                add = "jitter", fill = "SampleTypes")
# Pairwise comparison against reference and using non-parametric test (Wilcoxon test).
#a1 <- a1 + stat_compare_means(method = "wilcox.test", ref.group = "large_intestine", label = "p.format")

a1 <- a1 + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format")
print(a1)
a1_b <- a1 + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons_b, label = "p.format")
print(a1_b)


#because our data was not normally distributed. So we use nonparametric test. 
#a1 <- a1 + stat_compare_means(method = "t.test", ref.group = "large_intestine") #Pairwise comparision using non-parametric test (Wilcoxon test).
a2 <- ggboxplot(alpha_noPosC_noNegC_noBa, x = "SampleTypes", y = "pielou_evenness",
                add = "jitter", fill = "SampleTypes")
a2 <- a2 + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format") 
print(a2)
a2_b <- a2 + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons_b, label = "p.format") 
print(a2_b)

a3 <- ggboxplot(alpha_noPosC_noNegC_noBa, x = "SampleTypes", y = "Shannon_index",
                add = "jitter", fill = "SampleTypes")
a3 <- a3 + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format") 
print(a3)
a3_b <- a3 + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons_b, label = "p.format") 
print(a3_b)

a4 <- ggboxplot(alpha_noPosC_noNegC_noBa, x = "SampleTypes", y = "faith_pd",
                add = "jitter", fill = "SampleTypes")
a4 <- a4 + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons, label = "p.format") 
print(a4)
a4_b <- a4 + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons_b, label = "p.format") 
print(a4_b)


