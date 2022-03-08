##########################################
#
#       A.altamirnai skin micrbiome          
#       Martínez-Ugalde et al. 2022
#           emartug@gmail.com
#
##########################################

#Load libraries
library(qiime2R)
library(ggplot2)
library(phyloseq)
library(tidyverse)
library(ggforce)
library(fantaxtic)
library(readxl)
library(textshape)
library(tidyverse)
library(ggpubr)
library(UpSetR)
library(gridExtra)
library(grid)
library(vegan)


###Aplha and beta diversity stast and ploting between sample types##

#Data import to R
#Set work directory
setwd("~/A. altamirani microbiome/Imported qza files for R")

#Create a phyloseq object
ambysphy <- qza_to_phyloseq(
  features="filtrar_table.qza",  
  tree="Altamirani_rooted-tree.qza", 
  taxonomy= "Altamirani_taxonomy.qza", 
  metadata = "metadata_altamirani.txt", 
)

#Aplha diversity stats and figures

#Add alpha diveristy metrics as columns to metadata
metadata=read_q2metadata("metadata_altamirani_ultimaversion.txt")
shannon=read_qza("shannon_vector.qza")
faithpd=read_qza("faith_pd_vector.qza")
obsots=read_qza("observed_otus_vector.qza")

shannon=shannon$data %>% rownames_to_column("SampleID")
faithpd=faithpd$data %>% rownames_to_column("SampleID")
obsotu=obsotu$data %>% rownames_to_column("SampleID")

alfa=cbind(sample_data(ambysphy),shannon,faithpd,obsotu)

alfa2=alfa[,-28]
alfa2=alfa2[,-29]
alfa2=alfa2[,-30]
alfa2

#Add levels to metadata columns 
alfa2$Season<-factor(alfa2$Season,levels = c("Summer","Autumn","Winter","Spring"))
alfa2$Locality<-factor(alfa2$Locality,levels = c("LagunaSeca","Organillos","Sehuayan","Tecpan"))
alfa2$BdPresence<-factor(alfa2$BdPresence,levels = c("Positive","Negative"))
alfa2$SampleType<-factor(alfa2$SampleType,levels = c("Metamorphic","Non-metamorphic","Sediment","Water"))

#Stats for alpha diversity
kruskal.test(observed_otus ~ SampleType, data = alfa2)
pairwise.wilcox.test(x = alfa2$observed_otus, g = alfa2$SampleType, p.adjust.method = "BH")

kruskal.test(faith_pd ~ Season,data = alfa2)
pairwise.wilcox.test(x = alfa2$faith_pd, g = alfa2$Season, p.adjust.method = "BH")

kruskal.test(faith_pd ~ Locality,data = alfa2)
pairwise.wilcox.test(x = alfa2$faith_pd, g = alfa2$Season, p.adjust.method = "BH")

#Ploting alpha diversity between the four sample types
#Generate a comparission list acording to KW results in order to plot brackets for boxplots
alfacomps=list(c("Metamorphic","Non-metamorphic"),c("Metamorphic","Water"),c("Metamorphic","Sediment"),c("Non-metamorphic","Sediment"),
                c("Non-metamorphic","Water"),c("Sediment","Water"))

alfacomps1=list(c("Metamorphic","Non-metamorphic"),c("Metamorphic","Sediment"),c("Non-metamorphic","Sediment"),
                 c("Non-metamorphic","Water"),c("Sediment","Water"))

faisig<-ggplot(data=alfa2, aes(x=SampleType, y=faith_pd, fill=SampleType))+
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = alfacomps,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  labs(y="Faith's\nPhylogenetic Diversity")

shannsig<-ggplot(data=alfa2, aes(x=SampleType, y=shannon, fill=SampleType))+
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = alfacomps1,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  labs(y="Shannon \n Diversity Index",color="Sample Type")

otusig<-ggplot(data=alfa2, aes(x=SampleType, y=observed_otus, fill=SampleType))+
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = alfacomps1,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  labs(y="Observed ASV's",color="Sample Type")


#PCoA plots based on the weighted UniFrac matrix

#PCoA colred by sample type
betawuni<-wunifrac$data$Vectors %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=SampleType))+
  scale_color_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  geom_point(size=3) + 
  guides(color=guide_legend(ncol = 2))+
  labs(x="PCoA1 (47.86 %)",y="PCoA2 (23.13 %)")

#For beta diversity stats check A.-altamirani-16S-skin-microbiota/Scripts/A_altamirani_QIIME2_Workflow.txt 









