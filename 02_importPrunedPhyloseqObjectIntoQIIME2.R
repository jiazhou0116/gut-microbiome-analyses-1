# Regenerate new tree using pruned phyloseq object.
library(biomformat);packageVersion("biomformat")
#otu<-t(as(otu_table(Toad6),"matrix")) # 't' to transform if taxa_are_rows=FALSE
#if taxa_are_rows=TRUE
otu<-as(otu_table(Toad5),"matrix")
otu_biom<-make_biom(data=otu)
write_biom(otu_biom,"otu_biom_Toad5.biom")
##RUN CODE IN QIIME2
# qiime tools import \
# --input-path /Users/jiazhou/Box/CaneToadData1/Physeq/otu_biom_Toad5.biom \
# --type 'FeatureTable[Frequency]' \
# --input-format BIOMV100Format \
# --output-path feature-table_Toad5.qza

# qiime feature-table filter-seqs \
# --i-data rep-seqs-dada2.qza \
# --i-table feature-table_Toad5.qza \
# --o-filtered-data rep-seqs-Toad5.qza

# qiime phylogeny align-to-tree-mafft-fasttree \
# --i-sequences rep-seqs-Toad5.qza \
# --o-alignment aligned-rep-seqs-Toad5.qza \
# --o-masked-alignment masked-aligned-rep-seqs-Toad5.qza \
# --o-tree unrooted-tree-Toad5.qza \
# --o-rooted-tree rooted-tree-Toad5.qza

Toad6<-qza_to_phyloseq(features="feature-table_Toad5.qza", tree="rooted-tree-Toad5.qza", taxonomy="taxonomy_gg.qza", metadata="MetaData_CaneToadData1_2.txt")
Toad6