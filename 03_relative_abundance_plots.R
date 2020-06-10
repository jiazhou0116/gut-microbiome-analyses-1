#######################################################################################
### BARPLOT ABUNDANCES for characterization of microbial community
## merge samples for each SampleTypes
Toad5m = merge_samples(Toad5, "SampleTypes")
## phylum abund
Toad5m_phy <- Toad5m %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770"
)
# Plot 
plotphy<-ggplot(Toad5m_phy, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance \n (Bacterial Phyla >2%) \n") +
  ggtitle("Phylum Composition of Different Sample Types for Cane Toad's Gut Microbiota") 
plotphy

## Notice that each sample doesn’t fully add up to 1. 
## This reflects the rare phyla that were removed. 
## If you want your plot to look like everything adds up, 
## you can add position = “fill” to the geom_bar() command.


## Now with plot abund with 'genus' or lowest classification
Toad5m_cla<- Toad5m %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Class)                                      # Sort data frame alphabetically by phylum
# Plot 
plotcla<-ggplot(Toad5m_cla, aes(x = Sample, y = Abundance, fill = Class)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1, title="Class")) +
  ylab("Relative Abundance \n (Bacterial Genera >2%) \n") +
  ggtitle("Class Composition of Different Sample Types for Cane Toad's Gut Microbiota")
plotcla

## Now with plot abund with 'genus' or lowest classification
Toad5m_gen<- Toad5m %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by phylum
# Plot 
gen_cols<-c("g__Achromobacter" = "dark salmon", "g__Acidiphilium" = "orange", "g__Acinetobacter" = "red", "g__Aeromicrobium" = "turquoise", "g__Aeromonas" = "blue", "g__Akkermansia" = "cadet blue", "g__Arcobacter" = "gold", "g__Anaerovorax" = "green", "g__Anaerorhabdus" = "blue violet",
            "g__Bacillus" = "yellow", "g__Bacteroides" = "pink", "g__Brachyspira" = "purple", "g__Cetobacterium" = "maroon", "g__Clostridium" = "olive drab", "g__Comamonas" = "spring green", "g__Elizabethkingia" = "corn flower blue", "g__Epulopiscium" = "light sea green", "g__Flavobacterium" = "violet", "g__Fluviicola" = "sky blue", "g__Microvirgula" = "powder blue","g__Niabella" = "aqua marine", "g__Odoribacter" = "lime green", "g__Oscillospira" = "khaki", "g__Parabacteroides" = "coral", "g__Plesiomonas" = "wheat", "g__Providencia" = "yellow green", "g__Pseudomonas" = "misty rose", "g__Shewanella" = "alice blue", "g__Sphingobacterium" = "tan", "g__Treponema" = "orchid", "g__Turicibacter" = "dark gray", "g__Vibrio" = "lavender")
plotgen<-ggplot(Toad5m_gen, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = gen_cols) +
  scale_x_discrete(
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1, title="Genus")) +
  ylab("Relative Abundance \n (Bacterial Genera >2%) \n") +
  ggtitle("Genus Composition of Different Sample Types for Cane Toad's Gut Microbiota")
plotgen

## Combine each plot into one graphic.
library(gridExtra)
grid.arrange(nrow=2, plotphy, plotcla)    

## merge samples for each sex
Toad5_sex = merge_samples(Toad5, "sex")
## phylum abund
Toad5_sex_phy <- Toad5_sex %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

# Set colors for plotting
phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770"
)
# Plot 
plotphy_sex<-ggplot(Toad5_sex_phy, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance \n (Bacterial Phyla >2%) \n") +
  ggtitle("Phylum Composition of Different Sex for Cane Toad's Gut Microbiota") 
plotphy_sex

## Now with plot abund with 'genus' or lowest classification
Toad5_sex_gen<- Toad5_sex %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Genus)                                      # Sort data frame alphabetically by phylum
# Plot 
gen_cols<-c("g__Achromobacter" = "dark salmon", "g__Acidiphilium" = "orange", "g__Acinetobacter" = "red", "g__Aeromicrobium" = "turquoise", "g__Aeromonas" = "blue", "g__Akkermansia" = "cadet blue", "g__Arcobacter" = "gold", "g__Anaerovorax" = "green", "g__Anaerorhabdus" = "blue violet",
            "g__Bacillus" = "yellow", "g__Bacteroides" = "pink", "g__Brachyspira" = "purple", "g__Cetobacterium" = "maroon", "g__Clostridium" = "olive drab", "g__Comamonas" = "spring green", "g__Elizabethkingia" = "corn flower blue", "g__Epulopiscium" = "light sea green", "g__Flavobacterium" = "violet", "g__Fluviicola" = "sky blue", "g__Microvirgula" = "powder blue","g__Niabella" = "aqua marine", "g__Odoribacter" = "lime green", "g__Oscillospira" = "khaki", "g__Parabacteroides" = "coral", "g__Plesiomonas" = "wheat", "g__Providencia" = "yellow green", "g__Pseudomonas" = "misty rose", "g__Shewanella" = "alice blue", "g__Sphingobacterium" = "tan", "g__Treponema" = "orchid", "g__Turicibacter" = "dark gray", "g__Vibrio" = "lavender")
plotgen_sex<-ggplot(Toad5_sex_gen, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = gen_cols) +
  scale_x_discrete(
    drop = FALSE
  ) +
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  #
  guides(fill = guide_legend(reverse = F, keywidth = 1, keyheight = 1, title="Genus")) +
  ylab("Relative Abundance \n (Bacterial Genera >2%) \n") +
  ggtitle("Genus Composition of Different Sex for Cane Toad's Gut Microbiota")
plotgen_sex

## Combine each plot into one graphic.
library(gridExtra)
grid.arrange(nrow=2, plotphy, plotcla) 
