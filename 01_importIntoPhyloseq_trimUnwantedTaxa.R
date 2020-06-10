#Import qiime2 output into phyloseq and trim the taxa that are unwanted
# set your working directory to where the files are located
setwd("")
# to check where your working directory is use
getwd()

# load your libraries (install packages if needed)
library(tidyverse)
#devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(vegan)
library(plyr)
library(phyloseq)
library(gridExtra)
library(ape)
library(data.table)
library(lme4)

physeq<-qza_to_phyloseq(features="table-dada2.qza", tree="rooted_tree.qza", taxonomy="taxonomy_gg.qza", metadata="MetaData_CaneToadData1_2.txt")
#physeq<-qza_to_phyloseq(features="table-dada2.qza", tree="unrooted-tree.qza", taxonomy="taxonomy_gg.qza", metadata="MetaData_CaneToadData1_2.txt")

physeq
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 8672 taxa and 89 samples ]
# sample_data() Sample Data:       [ 89 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 8672 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 8672 tips and 8546 internal nodes ]

rank_names(physeq)
#[1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"


# As we have previously tested and checked our controls lets remove them
sample_data(physeq)$SampleTypes

Toad1<-subset_samples(physeq, SampleTypes%in%c("cloacal","large_intestine", "feces", "small_intestine"))
Toad1
taxtab1<-table(tax_table(Toad1)[,"Phylum"], exclude = NULL)

# now let's trim thse taxa that are unwanted
# Create table, number of features for Kingdom
taxtab<-table(tax_table(Toad1)[,"Kingdom"], exclude = NULL)
taxtab
#k__Archaea k__Bacteria  Unassigned 
#2        8617          53 

# So let's remove Archaea and those unassigned
Toad2 <-subset_taxa(Toad1, ! is.na(Kingdom) & ! Kingdom %in%  c("","Unassigned", "k__Archaea"))

# Create table for phylum
taxtab2<-table(tax_table(Toad2)[,"Phylum"], exclude = NULL)
taxtab2

# You can see that there are NAs present so let's remove NAs.
#Toad3 <-subset_taxa(Toad2, ! is.na(Phylum) & ! Phylum %in%  c("","uncharacterized"))

taxtab3<-table(tax_table(Toad2)[,"Class"], exclude = NULL)
taxtab3

####Question2 Chloroplast should not be included in bacteria...Eukaryotic
## remove chloroplasts
filterClass = c("c__Chloroplast")
Toad3 = subset_taxa(Toad2, ! Class %in% filterClass)

# because we have trimmed many let's remove those OTUs that 
# may not have any counts in samples anymore
# are there any with 0 count?
any(taxa_sums(Toad3) == 0)
# [1] TRUE
## How many?
sum(taxa_sums(Toad3) == 0)
# [1] 898
# remove those “unobserved” OTUs, i.e. 
# prune OTUs that are not present in at least one sample
Toad4 = prune_taxa(taxa_sums(Toad3) > 0, Toad3)
# we want to look at read depth across all samples
# Take a look at spread of read depth across samples
sdt = data.table(as(sample_data(Toad4), "data.frame"),
                 TotalReads = sample_sums(Toad4), keep.rownames = TRUE)
pSeqDepth = ggplot(sdt, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
pSeqDepth
# lowest is around 70000 and largest is 180000 giving you an idea of spread.

# What about the total reads per sample
# and what does the distribution look like?
q1<-qplot(log10(rowSums(otu_table(Toad4))),binwidth=0.2, main="kept singletons data") +
  xlab("Logged counts-per-sample") 
q1

# now let's get a count of the number of single and double OTU counts
tdt = data.table(tax_table(Toad4),
                 TotalCounts = taxa_sums(Toad4),
                 OTU = taxa_names(Toad4))

# How many singletons?
tdt[(TotalCounts <= 1), .N]
# [1] 0

####Question. why remove the OTUs with less than 4 prevalence. 
#in order to reduce the tail of data distribution: 
#doubletons: more than 2 counts existing in more than two samples.
# How many singletons and doubletons?
tdt[(TotalCounts <= 4), .N]
# [1] 2402
# remove all OTUs with less than 4 prevalence
Toad5 <- prune_taxa(taxa_sums(Toad4) > 4, Toad4)
# compare those results in the table, note number removed
Toad5
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5298 taxa and 72 samples ]
# sample_data() Sample Data:       [ 72 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 5298 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 5298 tips and 5196 internal nodes ]
Toad4
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 7700 taxa and 72 samples ]
# sample_data() Sample Data:       [ 72 samples by 18 sample variables ]
# tax_table()   Taxonomy Table:    [ 7700 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 7700 tips and 7597 internal nodes ]

taxtab5<-table(tax_table(Toad5)[,"Phylum"], exclude = NULL)
taxtab5
##
# p__[Thermi]    p__Acidobacteria   p__Actinobacteria 
# 2                   6                 123 
# p__Bacteroidetes       p__Chlamydiae      p__Chloroflexi 
# 883                   4                  18 
# p__Cyanobacteria  p__Deferribacteres    p__Elusimicrobia 
# 23                   5                   1 
# p__Firmicutes     p__Fusobacteria p__Gemmatimonadetes 
# 1256                  55                   1 
# p__GN02    p__Lentisphaerae              p__OD1 
# 9                  24                  54 
# p__Planctomycetes   p__Proteobacteria     p__Spirochaetes 
# 8                1072                  10 
# p__Synergistetes      p__Tenericutes              p__TM6 
# 1                  54                   1 
# p__TM7  p__Verrucomicrobia                <NA> 
#   23                  24                1641 

# notice the skew of the data is reduced
q2<-qplot(log10(rowSums(otu_table(Toad5))),binwidth=0.2, main="remove doubletons from data") +
  xlab("Logged counts-per-sample")
q2

grid.arrange(nrow=2, q1, q2)

# to perform beta diversity testing we need to 
# compare either relative abundances or rarefied data
# first let's look at rarefaction
###question. the rarefaction based on the lowest sequencing depth of all the samples.
# we might lose a lot of information by remove sequencing based on the lowest sequence number. 
##???other option for rarefaction. or use data without rarefaction.
## hellinger has been used to transform data before calculating bray curtis

library(microbiome)

# use hellinger - transforms data to rel abundance -- for bray curtis
# and then performs square root transform
Toad5.hell  = microbiome::transform(Toad5, transform = 'hellinger')
q3<-qplot(log10(rowSums(otu_table(Toad5.hell))),binwidth=0.2, main = "Hellinger data") +
  xlab("Logged counts-per-sample") 
q3

grid.arrange(nrow=3, q1,q2,q3)


# check again do we need to remove any that have no counts
any(taxa_sums(Toad5) == 0)
# false so no need to take any other steps
