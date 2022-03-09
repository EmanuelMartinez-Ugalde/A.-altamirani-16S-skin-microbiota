##########################################
#
#       A.altamirani skin microbiome          
#       Martínez-Ugalde et al. 2022
#           emartug@gmail.com
#
##########################################

#RDA

UNIFRAC=read.csv(file = "wunifrac-matrix.csv",row.names = 1)

ENV1=read.csv(file="ENV-table1.csv",row.names = 1)
names(ENV1)

#z-score variables to get them at the same range
ENV1[,c(2:length(ENV1))]=scale(ENV1[,c(2:length(ENV1))]) 

ENV1$DevStage=factor(ENV1$DevStage,levels=c("Metamorphic","Non-metamorphic"))

#Run RDA
wuRDA=capscale(UNIFRAC ~ ., ENV1,sqrt.dist = TRUE)
summary(wuRDA)

#Select the best variables by forward selection using ordistep
step.forward=ordistep(wuRDA, psteps = 10000,direction = "forward",permutations = how(nperm = 999))

#Run RDA based on the model selection
wuRDA1=capscale( UNIFRAC ~ DevStage + pH + Weight + MeanTemp + DeltaTemp +      
                    DissolvedOxygen + Conductivity + Elevation , 
                  ENV1,sqrt.dist = TRUE)

#RDA stats
summary(wuRDA1)
variable=anova.cca(wuRDA1,by = "terms")
variable
axis=anova.cca(wuRDA1,by= "axis")
axis
coeficiente=coef(wuRDA1)
R2adj=RsquareAdj(wuRDA1)
R2adj
coeficiente

#Extract data form the RDA object in order to get a plot using ggplot
smry=summary(wuRDA1)
scrs=scores(wuRDA1)
df1=data.frame(smry$sites[,1:2]) # site scores for RDA1 and RDA2
df1$site=rownames(df1)  #add site names
df2=data.frame(smry$biplot[,1:2])  # mapping environmental variables
DevS=as.data.frame(ENV1$DevStage)
colnames(DevS)
names(DevS)[names(DevS) == "ENV1$DevStage"] = "DevStage"
df1=cbind(df1,DevS)


scor = scores(wuRDA1, display=c("lc", "cn", "bp"), scaling=2)
numeric_centroids=data.frame(scor$centroids)
numeric_centroids
str(numeric_centroids)
numeric_centroids$DevStage=rownames(numeric_centroids)
numeric_centroids


continuous_arrows=data.frame(scor$biplot)
continuous_arrows
continuous_arrows$class=rownames(continuous_arrows) #turn rownames into a variable
continuous_arrows=continuous_arrows[-c(1),]
continuous_arrows
site_centroids=data.frame(scor$constraints)
site_centroids
site_centroids$site_names=rownames(site_centroids)
site_centroids


basplot=plot(wuRDA1)
mult=attributes(basplot$biplot)$arrow.mul
mult=attributes(scores(wuRDA1))$const 

names(numeric_centroids)[names(numeric_centroids) == "DevStage"] <- "Centroids"
DEVcentroids=numeric_centroids[,-4]
DEVcentroids=numeric_centroids[1:2,]
DEVcentroids[1,3]="Metamorphic"
DEVcentroids[2,3]="Non-metamorphic"



RDAFig<-ggplot(numeric_centroids, aes(x = CAP1, y = CAP2))+
  geom_point(data=DEVcentroids,aes(shape=Centroids,color="darkred"),
             size=5,fill="darkred",color="darkred")+
  scale_shape_manual(values = c(17,18,19,22,25))+
  scale_size(guide = 'none')+
  geom_point(data = df1,aes(x=CAP1, y=CAP2,colour=DevStage),alpha=0.6,size=3)+
  scale_color_manual(values=c("#E69F00", "#009E73"))+
  geom_segment(data = continuous_arrows,
               arrow = arrow(length = unit(0.35, "cm")), colour = "black")+ 
  geom_text(data = continuous_arrows,
            aes(x= (mult + mult/60) * CAP1, y = (mult + mult/60) * CAP2, 
                label = class), 
            size = 3,parse = TRUE,
            hjust = 0.5)+
  guides(shape=FALSE,color=guide_legend(ncol=1))+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  xlim(-2.5, 3.5) +
  ylim(-2.5, 3.5) +
  xlab("CAP1 (17.33%)") + # Procentaje de CAP1
  ylab("CAP2 (5.04%)") + # Procentaje de CAP2
  labs(colour="Developmental Stage")+
  theme(axis.line.x = element_line(colour = 'black'),
        axis.line.y = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'black'),
        axis.title = element_text(size = 12,face = "bold"),
        axis.text = element_text(angle = 0,face = "bold",
                                 size = 11, hjust = 0.5,vjust = 0.5,colour = "black"),
        legend.background = element_rect(fill = 'transparent', colour = NA),
        legend.key = element_rect(fill = "transparent"),
        legend.key.size = unit(0.4, "cm"),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12,face = "bold"),
        panel.background = element_rect(fill = 'transparent', colour = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_line(colour = NA),
        plot.background = element_rect(fill = 'transparent', colour = NA),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        text = element_text(size = 10,face = "bold"))

