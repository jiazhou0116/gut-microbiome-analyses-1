#############################################
# DIFFERENTIAL EXPRESSION WITH DESEQ to compare between samples
library(DESeq2)
# Deseq for differential differences of OTUs between groups
# set reference / control level 
# we idnetify 'control' as our common comparative element

# Create table for Genus
taxtab<-table(tax_table(Toad5)[,"Phylum"], exclude = NULL)
taxtab
# p__[Thermi]    p__Acidobacteria   p__Actinobacteria    p__Bacteroidetes       p__Chlamydiae 
# 2                   6                 123                 883                   4 
# p__Chloroflexi    p__Cyanobacteria  p__Deferribacteres    p__Elusimicrobia       p__Firmicutes 
# 18                  23                   5                   1                1256 
# p__Fusobacteria p__Gemmatimonadetes             p__GN02    p__Lentisphaerae              p__OD1 
# 55                   1                   9                  24                  54 
# p__Planctomycetes   p__Proteobacteria     p__Spirochaetes    p__Synergistetes      p__Tenericutes 
# 8                1072                  10                   1                  54 
# p__TM6              p__TM7  p__Verrucomicrobia                <NA> 
#   1                  23                  24                1641 

# You can see that there are NAs present so let's remove NAs.Removes 44.87% taxa that could not assigned to phylum. 
Toad6 <-subset_taxa(Toad5, ! is.na(Phylum) & ! Phylum %in%  c("","uncharacterized"))
taxtab2<-table(tax_table(Toad6)[,"Phylum"], exclude = NULL)
# p__[Thermi]    p__Acidobacteria   p__Actinobacteria    p__Bacteroidetes       p__Chlamydiae 
# 2                   6                 123                 883                   4 
# p__Chloroflexi    p__Cyanobacteria  p__Deferribacteres    p__Elusimicrobia       p__Firmicutes 
# 18                  23                   5                   1                1256 
# p__Fusobacteria p__Gemmatimonadetes             p__GN02    p__Lentisphaerae              p__OD1 
# 55                   1                   9                  24                  54 
# p__Planctomycetes   p__Proteobacteria     p__Spirochaetes    p__Synergistetes      p__Tenericutes 
# 8                1072                  10                   1                  54 
# p__TM6              p__TM7  p__Verrucomicrobia 
# 1                  23                  24 

###set "large_intestine" as common comparative element
sample_data(Toad6)$SampleTypes <- relevel(sample_data(Toad6)$SampleTypes, "large_intestine")
diagdds = phyloseq_to_deseq2(Toad6, ~ SampleTypes)
diagdds
diagdds_lg = DESeq(diagdds, test="Wald", fitType="parametric")
# lists the effects present in the result
# set contrast to extract results as needed
resultsNames(diagdds_lg)
# [1] "Intercept"                                     
# [2] "SampleTypes_cloacal_vs_large_intestine"        
# [3] "SampleTypes_feces_vs_large_intestine"          
# [4] "SampleTypes_small_intestine_vs_large_intestine"
# If you want, you can ake a look at the DESeq2 table.
DESeq2Table <- estimateSizeFactors(diagdds)
sizeFactors(DESeq2Table)


# individually look at the results for each pair
# Feces v Lge Intest
res_fec = results(diagdds_lg, cooksCutoff = FALSE, contrast=c("SampleTypes", "feces", "large_intestine"))
alpha = 0.05
sigtab_fec_lge = res_fec[which(res_fec$padj < alpha), ]
sigtab_fec_lge = cbind(as(sigtab_fec_lge, "data.frame"), as(tax_table(Toad5)[rownames(sigtab_fec_lge), ], "matrix"))

write.table(sigtab_fec_lge, "sigtab_fec_lge.txt") 

# Cloaca v Lge Intest
res_clo = results(diagdds_lg, cooksCutoff = FALSE, contrast=c("SampleTypes", "cloacal", "large_intestine"))
alpha = 0.05
sigtab_clo_lge = res_clo[which(res_clo$padj < alpha), ]
sigtab_clo_lge = cbind(as(sigtab_clo_lge, "data.frame"), as(tax_table(Toad5)[rownames(sigtab_clo_lge), ], "matrix"))

write.table(sigtab_clo_lge, "sigtab_clo_lge.txt") 

# Small Intest v Lge Intest
res_sml = results(diagdds_lg, cooksCutoff = FALSE, contrast=c("SampleTypes", "small_intestine", "large_intestine"))
alpha = 0.05
sigtab_sml_lge = res_sml[which(res_sml$padj < alpha), ]
sigtab_sml_lge = cbind(as(sigtab_sml_lge, "data.frame"), as(tax_table(Toad5)[rownames(sigtab_sml_lge), ], "matrix"))

write.table(sigtab_sml_lge, "sigtab_sml_lge.txt") 


#An example plot with class level
cols<-c("p__Actinobacteria" = "red", "p__Bacteroidetes" = "blue", "p__Chloroflexi" = "green", 
        "p__Firmicutes" = "yellow", "p__Proteobacteria" = "pink", "p__TM7" = "purple", "NA" = "grey")

# Phylum order
x = tapply(sigtab_fec_lge$log2FoldChange, sigtab_fec_lge$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_fec_lge$Phylum = factor(as.character(sigtab_fec_lge$Phylum), levels=names(x))
# Class order
x = tapply(sigtab_fec_lge$log2FoldChange, sigtab_fec_lge$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab_fec_lge$Class = factor(as.character(sigtab_fec_lge$Class), levels=c("NA", "c__Actinobacteria", "c__Alphaproteobacteria", 
                                                                           "c__Bacilli", "c__Bacteroidia", "c__Betaproteobacteria", 
                                                                           "c__Clostridia", "c__Deltaproteobacteria", "c__Epsilonproteobacteria", "c__Erysipelotrichi", "c__Flavobacteriia", "c__Fusobacteriia",
                                                                           "c__Gammaproteobacteria", "c__Sphingobacteriia", "c__Thermomicrobia",  "c__TM7-3", "c__[Saprospirae]"))

fec_lge_cla <- ggplot(sigtab_fec_lge, aes(y=Class, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  scale_color_manual(values = cols) +
  geom_point(size=2) + 
  xlim(-40, 30) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  ggtitle("Feces versus Large Intestine, p.value=0.05") ## \n gives new line

# Phylum order
x = tapply(sigtab_clo_lge$log2FoldChange, sigtab_clo_lge$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_clo_lge$Phylum = factor(as.character(sigtab_clo_lge$Phylum), levels=names(x))
# Class order
x = tapply(sigtab_clo_lge$log2FoldChange, sigtab_clo_lge$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab_clo_lge$Class = factor(as.character(sigtab_clo_lge$Class), levels=c("NA", "c__Actinobacteria", "c__Alphaproteobacteria", 
                                                                           "c__Bacilli", "c__Bacteroidia", "c__Betaproteobacteria", 
                                                                           "c__Clostridia", "c__Deltaproteobacteria", "c__Epsilonproteobacteria", "c__Erysipelotrichi", "c__Flavobacteriia", "c__Fusobacteriia",
                                                                           "c__Gammaproteobacteria", "c__Sphingobacteriia", "c__Thermomicrobia",  "c__TM7-3", "c__[Saprospirae]"))

clo_lge_cla <- ggplot(sigtab_clo_lge, aes(y=Class, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  scale_color_manual(values = cols) +
  geom_point(size=2) + 
  xlim(-40, 30) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  ggtitle("Cloaca versus Large Intestine, p.value=0.05") ## \n gives new line
# Phylum order
x = tapply(sigtab_sml_lge$log2FoldChange, sigtab_sml_lge$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_sml_lge$Phylum = factor(as.character(sigtab_sml_lge$Phylum), levels=names(x))
# Class order
x = tapply(sigtab_sml_lge$log2FoldChange, sigtab_sml_lge$Class, function(x) max(x))
x = sort(x, TRUE)

sml_lge_cla <- ggplot(sigtab_sml_lge, aes(y=Class, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  scale_color_manual(values = cols) +
  geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  ggtitle("Small Intestine versus Large Intestine, p.value=0.05") ## \n gives new line
# Error: Insufficient values in manual scale. 8 needed but only 7 provided.

library(gridExtra)
library(grid)
grid.arrange(fec_lge_cla, clo_lge_cla, ncol = 1, nrow = 2)


# then repeat comparisons with a different control 
#set "small_intestine" as common comparative element
sample_data(Toad6)$SampleTypes <- relevel(sample_data(Toad6)$SampleTypes, "small_intestine")
diagdds = phyloseq_to_deseq2(Toad6, ~ SampleTypes)
diagdds
diagdds_sml = DESeq(diagdds, test="Wald", fitType="parametric")
# lists the effects present in the result
# set contrast to extract results as needed
resultsNames(diagdds_sml)
# [1] "Intercept"                                     
# [2] "SampleTypes_cloacal_vs_large_intestine"        
# [3] "SampleTypes_feces_vs_large_intestine"          
# [4] "SampleTypes_small_intestine_vs_large_intestine"
# If you want, you can ake a look at the DESeq2 table.
DESeq2Table <- estimateSizeFactors(diagdds)
sizeFactors(DESeq2Table)

# individually look at the results for each pair
# Feces v Sml Intest
res_fec = results(diagdds_sml, cooksCutoff = FALSE, contrast=c("SampleTypes", "feces", "small_intestine"))
alpha = 0.05
sigtab_fec_sml = res_fec[which(res_fec$padj < alpha), ]
sigtab_fec_sml = cbind(as(sigtab_fec_sml, "data.frame"), as(tax_table(Toad5)[rownames(sigtab_fec_sml), ], "matrix"))

write.table(sigtab_fec_sml, "sigtab_fec_sml.txt") 

# Cloaca v Sml Intest
res_clo = results(diagdds_sml, cooksCutoff = FALSE, contrast=c("SampleTypes", "cloacal", "small_intestine"))
alpha = 0.05
sigtab_clo_sml = res_clo[which(res_clo$padj < alpha), ]
sigtab_clo_sml = cbind(as(sigtab_clo_sml, "data.frame"), as(tax_table(Toad5)[rownames(sigtab_clo_sml), ], "matrix"))

write.table(sigtab_clo_sml, "sigtab_clo_sml.txt") 

# Lge Intest V Small Intest
res_lge = results(diagdds_lg, cooksCutoff = FALSE, contrast=c("SampleTypes", "large_intestine","small_intestine"))
alpha = 0.05
sigtab_lge_sml = res_lge[which(res_lge$padj < alpha), ]
sigtab_lge_sml = cbind(as(sigtab_lge_sml, "data.frame"), as(tax_table(Toad5)[rownames(sigtab_lge_sml), ], "matrix"))

write.table(sigtab_sml_lge, "sigtab_lge_sml.txt") 


#An example plot with class level
cols<-c("p__Actinobacteria" = "red", "p__Bacteroidetes" = "blue", "p__Chloroflexi" = "green", 
        "p__Firmicutes" = "yellow", "p__Proteobacteria" = "pink", "p__TM7" = "purple", "NA" = "grey", "p__Fusobacteria" = "orange")
# Phylum order
x = tapply(sigtab_fec_sml$log2FoldChange, sigtab_fec_sml$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_fec_sml$Phylum = factor(as.character(sigtab_fec_sml$Phylum), levels=names(x))
# Class order
x = tapply(sigtab_fec_sml$log2FoldChange, sigtab_fec_sml$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab_fec_sml$Class = factor(as.character(sigtab_fec_sml$Class), levels=c("NA", "c__Actinobacteria", "c__Alphaproteobacteria", 
                                                                           "c__Bacilli", "c__Bacteroidia", "c__Betaproteobacteria", 
                                                                           "c__Clostridia", "c__Deltaproteobacteria", "c__Epsilonproteobacteria", "c__Erysipelotrichi", "c__Flavobacteriia", "c__Fusobacteriia",
                                                                           "c__Gammaproteobacteria", "c__Sphingobacteriia", "c__Thermomicrobia",  "c__TM7-3", "c__[Saprospirae]"))

fec_sml_cla <- ggplot(sigtab_fec_sml, aes(y=Class, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  scale_color_manual(values = cols) +
  geom_point(size=2) + 
  xlim(-40, 30) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  ggtitle("Feces versus Small Intestine, p.value=0.05") ## \n gives new line

# Phylum order
x = tapply(sigtab_clo_sml$log2FoldChange, sigtab_clo_sml$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab_clo_sml$Phylum = factor(as.character(sigtab_clo_sml$Phylum), levels=names(x))
# Class order
x = tapply(sigtab_clo_sml$log2FoldChange, sigtab_clo_sml$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab_clo_sml$Class = factor(as.character(sigtab_clo_sml$Class), levels=c("NA", "c__Actinobacteria", "c__Alphaproteobacteria", 
                                                                           "c__Bacilli", "c__Bacteroidia", "c__Betaproteobacteria", 
                                                                           "c__Clostridia", "c__Deltaproteobacteria", "c__Epsilonproteobacteria", "c__Erysipelotrichi", "c__Flavobacteriia", "c__Fusobacteriia",
                                                                           "c__Gammaproteobacteria", "c__Sphingobacteriia", "c__Thermomicrobia",  "c__TM7-3", "c__[Saprospirae]"))

clo_sml_cla <- ggplot(sigtab_clo_sml, aes(y=Class, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  scale_color_manual(values = cols) +
  geom_point(size=2) + 
  xlim(-40, 30) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  ggtitle("Cloaca versus Small Intestine, p.value=0.05") ## \n gives new line
sml_lge_cla <- ggplot(sigtab_sml_lge, aes(y=Class, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  scale_color_manual(values = cols) +
  geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  ggtitle("Small Intestine versus Large Intestine, p.value=0.05") ## \n gives new line

library(gridExtra)
library(grid)
grid.arrange(fec_lge_cla, clo_lge_cla, ncol = 1, nrow = 2)

###################################################################
#DIFFERENTIAL EXPRESSION WITH DESEQ to compare between male and female toads
library(DESeq2)
# Deseq for differential differences of OTUs between groups
# set reference / control level 
# in our case we will use feces as our 'control'
# we idnetify 'control' as our common comparative element

###set "large_intestine" as common comparative element
sample_data(Toad5_lge)$SampleTypes <- relevel(sample_data(Toad5_lge)$sex, "F")
diagdds = phyloseq_to_deseq2(Toad5_lge, ~ sex)
diagdds
diagdds_lge = DESeq(diagdds, test="Wald", fitType="parametric")
# lists the effects present in the result
# set contrast to extract results as needed
resultsNames(diagdds_lge)
#[1] "Intercept"  "sex_M_vs_F"

# male vs female
res = results(diagdds_lge, cooksCutoff = FALSE, contrast=c("sex", "M", "F"))
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Toad5_lge)[rownames(sigtab), ], "matrix"))

#An example plot 
#shows the OTUs that significantly differe in the male samples when compared with the female samples
sex_lge_class <- ggplot(sigtab, aes(y=Class, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  ggtitle("Male versus Female in Large Intestinal samples, p.value=0.05") ## \n gives new line
sex_lge_class

###set "cloaca" as common comparative element
sample_data(Toad5_clo)$SampleTypes <- relevel(sample_data(Toad5_clo)$sex, "F")
diagdds = phyloseq_to_deseq2(Toad5_clo, ~ sex)
diagdds
diagdds_clo = DESeq(diagdds, test="Wald", fitType="parametric")
# lists the effects present in the result
# set contrast to extract results as needed
resultsNames(diagdds_clo)
#[1] "Intercept"  "sex_M_vs_F"

# male vs female
res = results(diagdds_clo, cooksCutoff = FALSE, contrast=c("sex", "M", "F"))
alpha = 0.05
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(Toad5_clo)[rownames(sigtab), ], "matrix"))

#An example plot 
#shows the OTUs that significantly differe in the male samples when compared with the female samples
sex_clo_class <- ggplot(sigtab, aes(y=Class, x=log2FoldChange, color=Phylum)) + 
  geom_vline(xintercept = 0.0, color = "gray", size = 0.5) +
  geom_point(size=5) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=12)) +
  theme(axis.text.y = element_text(size=12)) +
  theme(axis.title=element_text(size=12)) +
  theme(legend.text=element_text(size=12), legend.title=element_text(size=12)) +
  ggtitle("Male versus Female in Cloacal Samples, p.value=0.05") ## \n gives new line
sex_clo_class