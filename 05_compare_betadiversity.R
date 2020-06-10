#####################################################################
# ORDINATION & ADONIS / PERMANOVA

Toad5.hell.pco.bray <- ordinate(Toad5.hell, "MDS", "bray")

#distinguish between SampleTypes.
p1<-plot_ordination(Toad5.hell, Toad5.hell.pco.bray, shape="SampleTypes", color ="SampleTypes") +
  geom_point(size=5, alpha=0.75) +
  scale_color_hue(labels = c("Cloaca", "Feces", "Large Intestine","Small Intestine")) +  #use scale_color_manual to map color with value
  scale_shape_manual(labels = c("Cloaca", "Feces", "Large Intestine","Small Intestine"), values=c(19,17, 18, 16)) + # values to shape are 19=circle, 17 = triangle
  ggtitle("Bray Curtis PCoA")
p1 + theme_bw() # change background to white

#distinguish between ToadID and SampleTypes.
p2<-plot_ordination(Toad5.hell, Toad5.hell.pco.bray, color="ToadID", shape="SampleTypes") +
  geom_point(size=5, alpha=0.75) +
  ggtitle("Bray Curtis PCoA")
p2 + theme_bw() # change background to white

#distinguish between ToadID.
p3<-plot_ordination(Toad5.hell, Toad5.hell.pco.bray, color="ToadID", shape="ToadID") +
  geom_point(size=5, alpha=0.75) +
  scale_color_manual(labels = c("Toad1", "Toad2", "Toad3","Toad4", "Toad5", "Toad6", "Toad7", "Toad8", "Toad9", "Toad10", "Toad11", "Toad12", "Toad13", "Toad14", "Toad15", "Toad16", "Toad17", "Toad18"),values=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,19,17, 18, 16)) + 
  scale_shape_manual(labels = c("Toad1", "Toad2", "Toad3","Toad4", "Toad5", "Toad6", "Toad7", "Toad8", "Toad9", "Toad10", "Toad11", "Toad12", "Toad13", "Toad14", "Toad15", "Toad16", "Toad17", "Toad18"), values=c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,19,17, 18, 16)) + # values to shape are 19=circle, 17 = triangle
  ggtitle("Bray Curtis PCoA")
p3 + theme_bw() + theme(legend.position = "bottom") # change background to white
#p3 + guides(col=guide_legend(nrow = 2,byrow = TRUE))


# the sample data frame
Toad5.df <- data.frame(sample_data(Toad5))

# let's compare with bray curtis 
Toad5.hell.bc <- phyloseq::distance(Toad5.hell, method = "bray")

ToadSample.BC<-adonis2(Toad5.hell.bc~SampleTypes,data=Toad5.df,permutations=9999)
ToadSample.BC
#Df SumOfSqs     R2     F Pr(>F)    
#SampleTypes  3   2.2217 0.1164 2.986  1e-04 ***
#  Residual    68  16.8653 0.8836                 
#Total       71  19.0870 1.0000  
##use strata to set ToadID as fixed factor. 
#The results were not very different from the results without using strata. 
#adonis2(Toad5.hell.bc ~ SampleTypes, strata = Toad5.df$ToadID, data = Toad5.df, permutations=9999)
# Df SumOfSqs      R2      F Pr(>F)    
# SampleTypes  3   2.2856 0.11729 3.0117  1e-04 ***
#   Residual    68  17.2018 0.88271                  
# Total       71  19.4874 1.00000  

# with a significant result, now look to check the dispersions 
# (these are the cloud of samples in the plots above for factor groups) 
# are NON significant
HOV1.BC <- betadisper(Toad5.hell.bc, Toad5.df$SampleTypes)

## Display pvalue of HOV
anova(HOV1.BC)
#Response: Distances
#Df  Sum Sq  Mean Sq F value Pr(>F)
#Groups     3 0.04173 0.013910  1.3607 0.2622
#Residuals 68 0.69516 0.010223    # pvalue is nonsiginificant, therefore we have confidence the difference is true

# pairwise testing within groups
##pairwise comparison after adonis significant. 
library(RVAideMemoire)
# test between groups in permanova
# p adjust multiple testing using Wilks Test with Benjamini Hochberg adjustment for pval
pairwise.perm.manova(Toad5.hell.bc, Toad5.df$SampleTypes, test="Wilks",p.method="hochberg",nperm=9999)
#                  cloacal feces  large_intestine
#  feces           0.0108  -      -              
#  large_intestine 0.3022  0.0005 -              
#  small_intestine 0.0291  0.0005 0.3022   

## p value correction: Benjamini & Yekutieli (2001) ("BY") for dependent tests
pairwise.perm.manova(Toad5.hell.bc, Toad5.df$SampleTypes, test="Wilks",p.method="BY",nperm=9999)
#                  cloacal feces   large_intestine
#  feces           0.00980 -       -              
#   large_intestine 0.50009 0.00074 -              
#   small_intestine 0.03197 0.00074 0.56766   

## compare between individuals using all sample types (needs replicates).
# the sample data frame
Toad5.df <- data.frame(sample_data(Toad5))

# let's compare with bray curtis 
Toad5.hell  = microbiome::transform(Toad5, transform = 'hellinger')
Toad5.hell.bc <- phyloseq::distance(Toad5.hell, method = "bray")

ToadID.BC<-adonis2(Toad5.hell.bc~ToadID,data=Toad5.df,permutations=9999)
ToadID.BC
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Toad5.hell.bc ~ ToadID, data = Toad5.df, permutations = 9999)
# Df SumOfSqs      R2      F Pr(>F)    
# ToadID   17  11.4502 0.58757 4.5254  1e-04 ***
#   Residual 54   8.0372 0.41243                  
# Total    71  19.4874 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# with a significant result, now look to check the dispersions 
# (these are the cloud of samples in the plots above for factor groups) 
# are NON significant
HOV1.BC <- betadisper(Toad5.hell.bc, Toad5.df$ToadID)
## Display pvalue of HOV
# anova(HOV1.BC)
# Analysis of Variance Table
# 
# Response: Distances
# Df  Sum Sq  Mean Sq F value Pr(>F)
# Groups    17 0.13396 0.007880  0.4103 0.9771
# Residuals 54 1.03709 0.019205 

#group gut microbiota variation by sex in estern toads
Toad5_lge_clo<-subset_samples(Toad5, SampleTypes%in%c("large_intestine","cloacal"))
Toad5_lge_clo_east<-subset_samples(Toad5_lge_clo, Location%in%c("Port_Douglas"))
##distinguish between sex
library(microbiome)
# use hellinger - transforms data to rel abundance -- for bray curtis
# and then performs square root transform
Toad5_lge_clo_east.hell  = microbiome::transform(Toad5_lge_clo_east, transform = 'hellinger')

#install.packages('ggrepel')
library("ggrepel")
Toad5_lge_clo_east.hell.pco.bray <- ordinate(Toad5_lge_clo_east.hell, "MDS", "bray")
sex_east<-plot_ordination(Toad5_lge_clo_east.hell, Toad5_lge_clo_east.hell.pco.bray, shape="sex", color ="SampleTypes") +
  geom_point(size=5, alpha=0.75) +
  #geom_text (mapping = aes(label = ToadID), size = 3, colour="black", nudge_y = -0.02)+
  # use nudge_y to change your distance from data point
  scale_color_hue(labels = c("cloacal","large_intestine")) +  #use scale_color_manual to map color with value
  scale_shape_manual(labels = c("F", "M"), values=c(19,17, 18, 16)) + # values to shape are 19=circle, 17 = triangle
  ggtitle("Bray Curtis PCoA for Male and Female Toads from Eastern Australia")
sex_east + theme_bw() # change background to white
