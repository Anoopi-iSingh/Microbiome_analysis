---
title: "16S_Phyloseq"
author: "Anoop Singh"
date: "08/03/2021"
output: "html_document"
---
#Phyloseq

###Load packages
```{r}
library(tidyverse); packageVersion("tidyverse")
library(phyloseq); packageVersion("phyloseq")
library(DESeq2); packageVersion("DESeq2")  
library(microbiome); packageVersion("microbiome")
library(vegan); packageVersion("vegan")  
library(picante); packageVersion("picante")
library(ALDEx2); packageVersion("ALDEx2")
library(metagenomeSeq); packageVersion("metagenomeSeq") 
library(HMP); packageVersion("HMP")  
library(dendextend); packageVersion("dendextend") 
library(selbal); packageVersion("selbal") 
library(rms); packageVersion("rms")
library(breakaway); packageVersion("breakaway")
library(ape)
```
###Construct phyloseq object
```{r}
abund_table<-read.delim("Data/table.tsv",row.names=1,check.names=FALSE)
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)
smpl<- read.delim("Data/sample-metadata.tsv",row.names=1,check.names=FALSE)
tax <- read.delim("Data/taxonomy.tsv",row.names=1,check.names=FALSE)
OTU_tree <- read.tree("Data/tree.nwk")
#phyloseq object
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX = tax_table(as.matrix(tax))
SAM = sample_data(smpl)
OTU_tree<-compute.brlen(OTU_tree,method="Grafen")
physeq<-merge_phyloseq(phyloseq(OTU,TAX),SAM,OTU_tree)
```
###Analysis
Relative abundance
```{r}
#relative abundance at phylum level
table(phyloseq::tax_table(physeq)[, "Phylum"])
ps_rel_abund = phyloseq::transform_sample_counts(physeq, function(x){x / sum(x)})
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  labs(x = "", y = "Relative Abundance\n") +
   facet_wrap(~ type, scales = "free") +
  theme(panel.background = element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
```
Boxplot
```{r}
ps_phylum <- phyloseq::tax_glom(ps_rel_abund, "Phylum")
phyloseq::taxa_names(ps_phylum) <- phyloseq::tax_table(ps_phylum)[, "Phylum"]
phyloseq::otu_table(ps_phylum)[1:5, 1:5]
phyloseq::psmelt(ps_phylum) %>%
  ggplot(data = ., aes(x = type, y = Abundance)) +
  geom_boxplot(outlier.shape  = NA) +
  geom_jitter(aes(color = OTU), height = 0, width = .2) +
  labs(x = "", y = "Abundance\n") +
  facet_wrap(~ OTU, scales = "free")
```


