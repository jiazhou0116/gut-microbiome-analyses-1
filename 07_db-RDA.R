##############################################################################################
## CONSTRAINED ANALYSIS
## Using "CAP" in Phyloseq : vegan::capscale = 
## distance-based redundancy analysis (db-RDA)
head(sample_data(Toad5))
## because replicates are pseudo-replicates regarding the measurements
## this analysis needs replicates, so we kept both large intestine and cloaca data
## only compare within large intestine and cloaca
Toad5_lge_clo<-subset_samples(Toad5, SampleTypes%in%c("large_intestine","cloacal"))

Toad5.bc <- phyloseq::distance(physeq = Toad5, method = "bray")
Toad5_lge_clo.bc <- phyloseq::distance(physeq = Toad5_lge_clo, method = "bray")
?cor
## for variable selection
## check the correlation of host characteristics
metadata<-read.csv("host_characteristics_correlation_test.csv")
metadata = as.data.frame(metadata)
metadata.num<-lapply(metadata, as.numeric)
# a character string indicating which correlation coefficient (or covariance) is to be computed. 
#One of "pearson" (default), "kendall", or "spearman": can be abbreviated.
cor_metadata<-cor(metadata[,1:4])
write.table(cor_metadata, "cor_host_characteristics_chap1.txt", sep="\t")
#            SUL	        SVL	        BodyWeight	sex
# SUL	       1	          0.942306688	0.854387127	0.129717767
# SVL	       0.942306688	1	          0.828766593	0.083711432
# BodyWeight 0.854387127	0.828766593	1	          0.242487452
# sex	       0.129717767	0.083711432	0.242487452	1

### we do not need step function if the variables are not too many. 
cap_mod <- ordinate(
  physeq = Toad5, 
  method = "CAP",
  distance = Toad5.bc,
  formula = ~ BodyWeight+sex+SUL)

anova(cap_mod, by="axis", perm=500)
head(summary(cap_mod))

cap_all1Toad$anova
anova(cap_all1Toad, by="axis", perm=500)
anova(cap_all1Toad, by="term", perm=500)
anova(cap_all1Toad, by="mar", perm=500)
# Permutation test for capscale under reduced model
# Marginal effects of terms
# Permutation: free
# Number of permutations: 999
# 
# Model: capscale(formula = distance ~ sex + BodyWeight, data = data)
# Df SumOfSqs      F Pr(>F)    
# sex         1   2.5701 8.7014  0.001 ***
#   BodyWeight  1   0.5981 2.0250  0.011 *  
#   Residual   69  20.3802                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


## remove step function because it's beneficial when the factors are many. 
## In my case, with or without step function, the plot result is the same. 
#cap_mod
{
  title = "Variation in Cane Toad Intestinal Microbiome Explained by Host Characteristics"
  cap_plot <- plot_ordination(
    physeq = Toad5, 
    ordination = cap_mod, 
    title = title,
    color = "SampleTypes", 
    axes = c(1,2)
  ) +
    #aes(shape = P) + # cannot be a continuous var
    geom_point(aes(colour = SampleTypes), alpha = 0.4, size = 4)  +
    geom_point(colour = "grey90", size = 1.5)  
  #scale_color_manual(values = c("#a65628", "red", "#ffae19", "#4daf4a", 
  #                              "#1919ff", "darkorchid3", "magenta"))
  
  
  # Now add the environmental variables as arrows
  arrowmat <- vegan::scores(cap_mod, display = "bp")
  arrowmat
  ## rename
  
  
  # Add labels, make a data.frame
  arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
  # Define the arrow aesthetic mapping
  arrow_map <- aes(xend = CAP1, 
                   yend = CAP2, 
                   x = 0, 
                   y = 0, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  label_map <- aes(x = 1.3 * CAP1, 
                   y = 1.3 * CAP2, 
                   shape = NULL, 
                   color = NULL, 
                   label = labels)
  
  arrowhead = arrow(length = unit(0.02, "npc"))
  
  # Make a new graphic
  cap_plot + 
    geom_segment(
      mapping = arrow_map, 
      size = 0.5, 
      data = arrowdf, 
      color = "gray", 
      arrow = arrowhead
    ) + 
    geom_text(
      mapping = label_map, 
      size = 4.3, 
      #colour = "black", 
      fontface = "bold",
      data = arrowdf, 
      show.legend = FALSE
    )
}