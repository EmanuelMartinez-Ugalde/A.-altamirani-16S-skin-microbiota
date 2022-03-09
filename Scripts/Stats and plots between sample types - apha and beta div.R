##########################################
#
#       A.altamirani skin microbiome          
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
library(ggpubr)

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

faisig=ggplot(data=alfa2, aes(x=SampleType, y=faith_pd, fill=SampleType))+
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = alfacomps,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  labs(y="Faith's\nPhylogenetic Diversity")

shannsig=ggplot(data=alfa2, aes(x=SampleType, y=shannon, fill=SampleType))+
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = alfacomps1,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  labs(y="Shannon \n Diversity Index",color="Sample Type")

otusig=ggplot(data=alfa2, aes(x=SampleType, y=observed_otus, fill=SampleType))+
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = alfacomps1,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  labs(y="Observed ASV's",color="Sample Type")


#PCoA plots based on the weighted UniFrac matrix

wunifrac=read_qza("weighted_unifrac_pcoa_results.qza")

#PCoA colred by sample type
betawuni=wunifrac$data$Vectors %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=SampleType))+
  scale_color_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  geom_point(size=3) + 
  guides(color=guide_legend(ncol = 2))+
  labs(x="PCoA1 (47.86 %)",y="PCoA2 (23.13 %)")

#For beta diversity stats check A.-altamirani-16S-skin-microbiota/Scripts/A_altamirani_QIIME2_Workflow.txt 

###Betadisper analysis 

ambysbetad=betadisper(ambyswuni, sample_data(ambysphy)$SampleType)

ki=permutest(ambysbetad,pairwise=TRUE) 

bdispplot=ambysbetad$distances
bdispplot1=as.data.frame(bdispplot)
alfadisp=cbind(alfa2,bdispplot1)

#Comparission list for brackets in boxplots
complist=list(c("Adult","Juvenile"),c("Adult","Sediment"),c("Juvenile","Sediment"),c("Juvenile","Water"),
             c("Sediment","Water"))

dispbox=ggplot(data=alfadisp,aes(x=SampleType, y=bdispplot,fill=SampleType)) +
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = complist,method = "t.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#E69F00", "#009E73", "#972D15", "#0072B2"))+
  labs(y="Distance to centrois", fill="Sample Type")

###UpSet plot

taxa_names(ambysphy) = paste0("ASV", seq(ntaxa(ambysphy)))#Change ASV id for a numericID 
ambysmerged <- merge_samples(ambysphy, "SampleType",fun = mean)#Merge ASV by sample type
write_phyloseq(ambysmerged, type = "all", path = "A. altamirani microbiome")

lvl7=read.csv(file="otu_table_by_ASV_presence_in_each_Sample.csv",header = T,sep = ";",
              stringsAsFactors = T, row.names = 1)#Counts were converted to 1 if present in any sample

j1=upset(lvl7,sets.bar.color = c("#972D15","#0072B2","#E69F00", "#009E73"),
          point.size = 2, line.size = 1, mainbar.y.label = "Shared ASV's",
          sets.x.label = "ASV's per Sample Type", text.scale = c(1.6, 1.4, 1.6, 1.4, 1.4, 2),
          queries = list(list(query = intersects,params = list("Metamorphic"), color = "#E69F00", active = T),
                         list(query = intersects,params = list("Non.metamorphic"), color = "#009E73", active = T),
                         list(query = intersects,params = list("Sediment"), color = "#972D15", active = T),
                         list(query = intersects,params = list("Water"), color ="#0072B2", active = T)))

###Alpha and beta stats and figures testing seasonal
###sampling location and bd influence in metamorphic and non-metamorphic skin bacterial communities

#Subset data and give levels to variables

Meta=subset(alfa2,SampleType=="Metamorphic")
Nometa=subset(alfa2,SampleType=="Non-metamorphic")

Meta$Season=factor(Meta$Season,levels = c("Summer","Autumn","Winter","Spring"))
Nometa$Season=factor(Nometa$Season,levels = c("Summer","Autumn","Winter","Spring"))

Meta$Locality=factor(Meta$Locality,levels = c("LagunaSeca","Organillos","Sehuayan","Tecpan"))
Nometa$BdPresence=factor(Nometa$BdPresence,levels = c("Positive","Negative"))

#Stats for metamorphic and non-metamorphic samples 
kruskal.test(faith_pd ~ Season,data = Meta)
pairwise.wilcox.test(x = Meta$faith_pd, g = Meta$Season, p.adjust.method = "BH" )

kruskal.test(faith_pd ~ Season,data = Nometa)
pairwise.wilcox.test(x = Nometa$faith_pd, g = Nometa$Season, p.adjust.method = "BH" )

kruskal.test(faith_pd ~ Locality,data = Meta)
pairwise.wilcox.test(x = Meta$faith_pd, g = Meta$Locality, p.adjust.method = "BH" )

kruskal.test(faith_pd ~ Locality,data = Nometa)
pairwise.wilcox.test(x = Nometa$faith_pd, g = Nometa$Locality, p.adjust.method = "BH" )

kruskal.test(faith_pd ~ BdPresence,data = Meta)
kruskal.test(faith_pd ~ BdPresence,data = Nometa)

#Alpha diversity plots

#Season
Seascomp1=list(c("Winter","Spring"))

faithsamet=Meta %>%
  filter(!is.na(faith_pd)) %>%
ggplot(aes(x=Season, y=faith_pd, fill=Season)) +
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = Seascomp1,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#855138ff", "#c58a57ff", "#8fc1d4ff", "#536c53ff"))+
  labs(y="Metamorphic Faith's\nPhylogenetic Diversity")

faithsnm=Nometa %>%
  filter(!is.na(faith_pd)) %>%
ggplot(aes(x=Season, y=faith_pd, fill=Season)) +
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  guides(fill=guide_legend(ncol = 2))+
  scale_color_manual(values=c("#855138ff", "#c58a57ff", "#8fc1d4ff", "#536c53ff"))+
  labs(y="Non-metamorphic Faith's\nPhylogenetic Diversity")

#Bd presence
bdcomps=list(c("Positive","Negative"))

faithbdm=Meta %>%
  filter(!is.na(faith_pd)) %>%
  filter(BdPresence=="Positive"|BdPresence=="Negative")%>%
  ggplot(aes(x=BdPresence, y=faith_pd, fill=BdPresence)) +
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  scale_fill_manual(values=c("#c8c6c6ff","#4b6587ff"))+
  labs(y="Metamorphic Faith's\nPhylogenetic Diversity")

faithbdnm=Nometa %>%
  filter(!is.na(faith_pd)) %>%
  filter(BdPresence=="Positive"|BdPresence=="Negative")%>%
  ggplot(aes(x=BdPresence, y=faith_pd, fill=BdPresence)) +
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  guides(fill=guide_legend(ncol = 2))+
  scale_fill_manual(values=c("#c8c6c6ff","#4b6587ff"))+
  labs(y="Non-metamorphic Faith's\nPhylogenetic Diversity")


#Locality

locomp6=list(c("LagunaSeca","Organillos"),c("LagunaSeca","Tecpan")
              ,c("Organillos","Tecpan"),c("Sehuayan","Tecpan"))

locomp1=list(c("Organillos","Sehuayan"))

faithlocm=Meta %>%
  filter(!is.na(faith_pd)) %>%
  ggplot(aes(x=Locality, y=faith_pd, fill=Locality)) +
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  stat_compare_means(comparisons = locomp1,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#ee9c39", "#186a62", "#942900", "#b46acd"))+
  labs(y="Metamorphic Faith's\nPhylogenetic Diversity")

faithlocnm=Nometa %>%
  filter(!is.na(faith_pd)) %>%
  ggplot(aes(x=Locality, y=faith_pd, fill=Locality)) +
  geom_boxplot(lwd=0.3,outlier.size = 0.7)+
  guides(fill=guide_legend(ncol = 2))+
  stat_compare_means(comparisons = locomp6,method = "wilcox.test",label = "p.signif",
                     size=3,tip.length = 0.01,step.increase = 0.1,vjust = 0.3,hide.ns = TRUE)+
  scale_fill_manual(values=c("#ee9c39", "#186a62", "#942900", "#b46acd"))+
  labs(y="Non-metamorphic Faith's\nPhylogenetic Diversity")


#Beta diversity plots
Metwunifrac=read_qza("Metamorphic-core-metrics/weighted_unifrac_pcoa_results.qza")
NMetwunifrac=read_qza("Nonmetamorphic-core-metrics/weighted_unifrac_pcoa_results.qza")

metbd=Metwunifrac$data$Vectors %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=BdPresence))+
  guides(color=guide_legend(ncol=2),shape=guide_legend(ncol=1))+
  scale_color_manual(values=c("#c8c6c6ff","#4b6587ff"))+
  geom_point(size=3) +
  labs(x="PCoA1 (49.12 %)",y="PCoA2 (14.29 %)")

nometbd=NMetwunifrac$data$Vectors %>%
  left_join(metadata) %>%
  filter(SampleType=="Non-metamorphic")%>%
  ggplot(aes(x=PC1, y=PC2, color=BdPresence))+
  guides(color=guide_legend(ncol=2),shape=guide_legend(ncol=1))+
  scale_color_manual(values=c("#c8c6c6ff","#4b6587ff"))+
  geom_point(size=3) + 
  labs(x="PCoA1 (64.56 %)",y="PCoA2 (8.009 %)")

metseas=Metwunifrac$data$Vectors %>%
  left_join(metadata) %>%
  filter(SampleType=="Metamorphic")%>%
  ggplot(aes(x=PC1, y=PC2, color=Season))+
  guides(color=guide_legend(ncol=2),shape=guide_legend(ncol=1))+
  scale_color_manual(values=c("#855138ff", "#c58a57ff", "#8fc1d4ff", "#536c53ff"))+
  geom_point(size=3) + 
  labs(x="PCoA1 (49.12 %)",y="PCoA2 (14.29 %)")

nometseas=NMetwunifrac$data$Vectors %>%
  left_join(metadata) %>%
  filter(SampleType=="Non-metamorphic")%>%
  ggplot(aes(x=PC1, y=PC2, color=Season))+
  guides(color=guide_legend(ncol=2),shape=guide_legend(ncol=1))+
  scale_color_manual(values=c("#855138ff", "#c58a57ff", "#8fc1d4ff", "#536c53ff"))+
  geom_point(size=3) + 
  labs(x="PCoA1 (64.56 %)",y="PCoA2 (8.009 %)")

metloc=Metwunifrac$data$Vectors %>%
  left_join(metadata) %>%
  filter(SampleType=="Metamorphic")%>%
  ggplot(aes(x=PC1, y=PC2, color=Locality))+
  guides(color=guide_legend(ncol=2),shape=guide_legend(ncol=1))+
  scale_color_manual(values=c("#ef7c8eff", "#638c80ff", "#85b7d8ff", "#887bb0ff"))+
  geom_point(size=3) + 
  labs(x="PCoA1 (49.12 %)",y="PCoA2 (14.29 %)")

nometloc=NMetwunifrac$data$Vectors %>%
  left_join(metadata) %>%
  filter(SampleType=="Non-metamorphic")%>%
  ggplot(aes(x=PC1, y=PC2, color=Locality))+
  guides(color=guide_legend(ncol=2),shape=guide_legend(ncol=1))+
  scale_color_manual(values=c("#ef7c8eff", "#638c80ff", "#85b7d8ff", "#887bb0ff"))+
  geom_point(size=3) + 
  labs(x="PCoA1 (64.56 %)",y="PCoA2 (8.009 %)")

bd=ggarrange(metbd,nometbd,ncol = 2,nrow = 1)
seas=ggarrange(metseas,nometseas,ncol = 2,nrow = 1)
loc=ggarrange(metloc,nometloc,ncol = 2,nrow = 1)









