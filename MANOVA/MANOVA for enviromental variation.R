##########################################
#
#       A.altamirani skin microbiome          
#       Martínez-Ugalde et al. 2022
#           emartug@gmail.com
#
##########################################

#Load libraries
library(dplyr)
library(readr)
library(tidyverse)
library(tidyr)
library(data.table)

#Loading data
dataMANOVA=read.csv(file="MANOVA_Table.csv",row.names = 1)

dataMANOVA$Season=factor(dataMANOVA$Season,levels = c("Summer","Autumn","Winter","Spring"))
dataMANOVA$Locality=factor(dataMANOVA$Locality,levels = c("Organillos","LagunaSeca","Sehuayan","Tecpan"))


am=manova(cbind(MeanTemp,SeasonMinTemp,SeasonMaxTemp,pH,
                DisolvedOxigen,Conductivity,DeltaTemp)~Season+Locality,data=dataMANOVA)


summary(am,test="Wilks")#Test statistic

summary.aov(am)#Summary of the variance 

