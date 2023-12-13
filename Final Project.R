#Cmd+Shift+R on the Mac

library(stringr)
library(ggplot2)
library(vegan)
library(tidyverse)
library(ggpubr)
library(vegan)
library(ggplot2)
library(MASS)
library(viridis)
library(patchwork)

# Data --------------------------------------------------------------------
setwd("/Users/logan/Fatty Acid Paper")
Merged <- read_csv("fatty_acid_isotope_merged.csv")
#separate table based on column
Saturated <- c("ID","Season", "Site", "C12.0", "C14.0", "C15.0", "C16.0", "C17.0", "C18.0", "C20.0")
Monounsaturated <- c("ID", "Season", "Site", "C16.1n.9", "C16.1n.7", "C17.1", "C18.1n.9", "C18.1n.7", "C20.1.1")
Polyunsaturated <- c("ID", "Season", 'Site', "C18.2n.6", "C18.3n.3", "C18.4n.3", "C20.4n.6", "C20.5n.3", "C22.2n.9", "C22.4n.6", "C22.5n.6", "C22.5n.3", "C22.6n.3")
Isotope <- c("ID", "Season", "Site", "N2.Amp", "perN", "del15N", "CO2.Amp", "perC", "del13C", "Weight..mg..1", "SO2.Area", "perS", "del34S", "weight..mg.", "H2.Amp", "perH", "del2H", "CO.Amp", "perO", "del18O")

#create new tables
Saturated <- Merged[, Saturated, with = FALSE]
Monounsaturated <- Merged[, Monounsaturated, with = FALSE]
Polyunsaturated <- Merged[, Polyunsaturated, with = FALSE]
Isotope <- Merged[, Isotope, with = FALSE]



# Saturated LDA (with color blind friendly plot)---------------------------------------------------------------------

Saturated_Filter <- Saturated[,c(4:10)]
Saturated_LDA <- Saturated_Filter[,colMeans(Saturated[,c(4:10)])>1]
Saturated_LDA$Site_Season <- paste(Saturated$Site, Saturated$Season, sep = "_")

#run discriminant analysis
SaturatedFA_LDA <- lda(Site_Season~.,Saturated_LDA)

#extract data on fatty acid loading data
Saturated_LDA_Loading <- as.data.frame(SaturatedFA_LDA$scaling)
Saturated_LDA_Loading$FA <- row.names(Saturated_LDA_Loading)
Saturated_LDA_Loading$Importance <- sqrt(Saturated_LDA_Loading$LD1^2 + Saturated_LDA_Loading$LD2^2)

#sort in order of importance
Saturated_LDA_Loading <- Saturated_LDA_Loading[order(-Saturated_LDA_Loading$Importance), ]

#Extract ordinations, tells how samples are related to each fatty acid and to each other
Saturated_Sample_Ordination <- as.data.frame(predict(SaturatedFA_LDA)$x)
Saturated_Sample_Ordination$Site <- Saturated$Site
Saturated_Sample_Ordination$Season <- Saturated$Season
#run lda to get eigenvalue, eigenvalue is proportion of trace

#LDA Plot
Saturated_LDA_Plot <- ggplot()+
  geom_text(data=Saturated_LDA_Loading,aes(LD1,LD2,label=FA),size=2)+
  geom_point(data=Saturated_Sample_Ordination,aes(LD1,LD2,color=Site,shape=Season),size=2)+
  xlab("Proportion Eigenvalue 0.587")+
  ylab("Proportion Eigenvalue 0.344")+
  ggtitle("Saturated FA LDA") +
  theme_bw()+
  scale_color_brewer(palette = "Dark2") +
  theme(aspect.ratio = 1)

Saturated_LDA_Plot

# Polyunsaturated LDA (with color blind friendly plot) -----------------------------------------------------

PUFA_Filter <- Polyunsaturated[,c(4:13)]
PUFA_LDA <- PUFA_Filter[,colMeans(Polyunsaturated[,c(4:13)])>1]
PUFA_LDA$Site_Season <- paste(Polyunsaturated$Site, Polyunsaturated$Season, sep = "_")

#run discriminant analysis
PolyunsatFA_LDA <- lda(Site_Season~.,PUFA_LDA)

#extract data on fatty acid loading data
PUFALDA_Loading <- as.data.frame(PolyunsatFA_LDA$scaling)
PUFALDA_Loading$FA <- row.names(PUFALDA_Loading)
PUFALDA_Loading$Importance <- sqrt(PUFALDA_Loading$LD1^2 + PUFALDA_Loading$LD2^2)

#sort in order of importance
PUFALDA_Loading <- PUFALDA_Loading[order(-PUFALDA_Loading$Importance), ]

#Extract ordinations, tells how samples are related to each fatty acid and to each other
PUFASample_Ordination <- as.data.frame(predict(PolyunsatFA_LDA)$x)
PUFASample_Ordination$Site <- Polyunsaturated$Site
PUFASample_Ordination$Season <- Polyunsaturated$Season
#run lda to get eigenvalue, eigenvalue is proportion of trace

#LDA Plot
PUFA_LDA_plot <- ggplot()+
  geom_text(data=PUFALDA_Loading,aes(LD1,LD2,label=FA),size=2)+
  geom_point(data=PUFASample_Ordination,aes(LD1,LD2,color=Site,shape=Season),size=2)+
  #scale_color_viridis(discrete=TRUE)+
  xlab("Proportion Eigenvalue 0.456")+
  ylab("Proportion Eigenvalue 0.323")+
  ggtitle("PUFA LDA") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  theme(aspect.ratio = 1)

PUFA_LDA_plot
# Monounsaturated LDA (with colorblind friendly plot) -----------------------------------------------------

MUFA_Filter <- Monounsaturated[,c(4:9)]
MUFA_LDA <- MUFA_Filter[,colMeans(Monounsaturated[,c(4:9)])>1]
MUFA_LDA$Site_Season <- paste(Monounsaturated$Site, Monounsaturated$Season, sep = "_")

#run discriminant analysis
MonounsatFA_LDA <- lda(Site_Season~.,MUFA_LDA)

#extract data on fatty acid loading data
MUFALDA_Loading <- as.data.frame(MonounsatFA_LDA$scaling)
MUFALDA_Loading$FA <- row.names(MUFALDA_Loading)
MUFALDA_Loading$Importance <- sqrt(MUFALDA_Loading$LD1^2 + MUFALDA_Loading$LD2^2)

#sort in order of importance
MUFALDA_Loading <- MUFALDA_Loading[order(-MUFALDA_Loading$Importance), ]

#Extract ordinations, tells how samples are related to each fatty acid and to each other
MUFASample_Ordination <- as.data.frame(predict(MonounsatFA_LDA)$x)
MUFASample_Ordination$Site <- Monounsaturated$Site
MUFASample_Ordination$Season <- Monounsaturated$Season
#run lda to get eigenvalue, eigenvalue is proportion of trace

#LDA Plot
MUFA_LDA_plot <- ggplot()+
  geom_text(data=MUFALDA_Loading,aes(LD1,LD2,label=FA),size=2)+
  geom_point(data=MUFASample_Ordination,aes(LD1,LD2,color=Site,shape=Season),size=2)+
  #scale_color_viridis(discrete=TRUE)+
  xlab("Proportion Eigenvalue 0.753")+
  ylab("Proportion Eigenvalue 0.179")+
  ggtitle("MUFA LDA") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  theme(aspect.ratio = 1)
MUFA_LDA_plot

# All Fatty Acid LDA (with color blind friendly plot)------------------------------------------------------


FA_Filter <- Merged[,c(7:29)]
FA_LDA <- FA_Filter[,colMeans(Merged[,c(7:29)])>1]
FA_LDA$Site_Season <- paste(Merged$Site, Merged$Season, sep = "_")

#run discriminant analysis
AllFA_LDA <- lda(Site_Season~.,FA_LDA)

#extract data on fatty acid loading data
FALDA_Loading <- as.data.frame(AllFA_LDA$scaling)
FALDA_Loading$FA <- row.names(FALDA_Loading)
FALDA_Loading$Importance <- sqrt(FALDA_Loading$LD1^2 + FALDA_Loading$LD2^2)

#sort in order of importance
FALDA_Loading <- FALDA_Loading[order(-FALDA_Loading$Importance), ]

#Extract ordinations, tells how samples are related to each fatty acid and to each other
AllSample_Ordination <- as.data.frame(predict(AllFA_LDA)$x)
AllSample_Ordination$Site <- Merged$Site
AllSample_Ordination$Season <- Merged$Season
#run lda to get eigenvalue, eigenvalue is proportion of trace

#LDA Plot
All_FA_LDA_plot <- ggplot()+
  geom_text(data=FALDA_Loading,aes(LD1,LD2,label=FA),size=2)+
  geom_point(data=AllSample_Ordination,aes(LD1,LD2,color=Site,shape=Season),size=2)+
  #scale_color_viridis(discrete=TRUE)+
  xlab("Proportion Eigenvalue 0.358")+
  ylab("Proportion Eigenvalue 0.276")+
  ggtitle("All Fatty Acids LDA") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  theme(aspect.ratio = 1)


# Isotope LDA -------------------------------------------------------------

Isotope_Filter <- Isotope[,c(6, 9, 13, 17, 20)]
Iso_LDA <- Isotope_Filter[,colMeans(Isotope[,c(6, 9, 13, 17, 20)])>-200]
Iso_LDA$Site_Season <- paste(Isotope$Site, Isotope$Season, sep = "_")

#run discriminant analysis
Isotope_LDA <- lda(Site_Season~.,Iso_LDA)

#extract data on stable isotope loading data
ISOLDA_Loading <- as.data.frame(Isotope_LDA$scaling)
ISOLDA_Loading$Isotope <- row.names(ISOLDA_Loading)
ISOLDA_Loading$Importance <- sqrt(ISOLDA_Loading$LD1^2 + ISOLDA_Loading$LD2^2)

#sort in order of importance
ISOLDA_Loading <- ISOLDA_Loading[order(-ISOLDA_Loading$Importance), ]

#Extract ordinations, tells how samples are related to each fatty acid and to each other
ISOSample_Ordination <- as.data.frame(predict(Isotope_LDA)$x)
ISOSample_Ordination$Site <- Isotope$Site
ISOSample_Ordination$Season <- Isotope$Season
#run lda to get eigenvalue, eigenvalue is proportion of trace

#LDA Plot
Isotope_LDA_Plot <- ggplot()+
  geom_text(data=ISOLDA_Loading,aes(LD1,LD2,label=Isotope),size=5)+
  geom_point(data=ISOSample_Ordination,aes(LD1,LD2,color=Site,shape=Season),size=5)+
  #scale_color_viridis(discrete=TRUE)+
  xlab("Proportion Eigenvalue 0.447")+
  ylab("Proportion Eigenvalue 0.399")+
  ggtitle("Isotope LDA") +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  theme(aspect.ratio = 1)

Isotope_LDA_Plot

#combine fatty acid lda plots
Saturated_LDA_Plot + MUFA_LDA_plot + PUFA_LDA_plot + All_FA_LDA_plot + plot_layout(nrow = 2, byrow = F)

#all plots have color blind friendly palettes

# Stable Isotope Boxplot (color blind friendly palettes) --------------------------------------------------
library(ggplot2)
library(viridisLite)
library(RColorBrewer)


#create oxygen stable isotope boxplot
Oboxplot <- ggplot(Isotope, aes(x = Site, y = del18O, color = Season, label = ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
  theme_bw()+
  scale_color_brewer(palette = "Dark2")

Oboxplot

#create hydrogen stable isotope boxplot
Hboxplot <- ggplot(Isotope, aes(x = Site, y = del2H, color = Season, label = ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  theme_bw()+
  scale_color_brewer(palette = "Dark2")


Hboxplot

#create nitrogen stable isotope boxplot
Nboxplot <- ggplot(Isotope, aes(x = Site, y = del15N, color = Season, label = ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  theme_bw()+
  scale_color_brewer(palette = "Dark2")
  
  
Nboxplot

#create carbon stable isotope boxplot
Cboxplot <- ggplot(Isotope, aes(x = Site, y = del13C, color = Season, label = ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  theme_bw()+
  scale_color_brewer(palette = "Dark2")
  
  
Cboxplot

#create sulfur stable isotope boxplot
Sboxplot <- ggplot(Isotope, aes(x = Site, y = del34S, color = Season, label = ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  theme_bw()+
  scale_color_brewer(palette = "Dark2")



Sboxplot


Cboxplot + Nboxplot + Hboxplot + Oboxplot + Sboxplot



# Fatty Acid Boxplot ------------------------------------------------------
DSAboxplot <- ggplot(Merged, aes(x = Site, y = C22.5n.3, color = Season, label = ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  ggtitle("DSA")+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_color_brewer(palette = "Dark2")

DSAboxplot

EPAboxplot <- ggplot(Merged, aes(x = Site, y = C20.5n.3, color = Season, label = ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  ggtitle("EPA")+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_color_brewer(palette = "Dark2")

EPAboxplot

DHAboxplot <- ggplot(Merged, aes(x = Site, y = C22.6n.3, color = Season, label = ID)) +
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + 
  ggtitle("DHA")+
  theme_bw()+
  theme(aspect.ratio = 1)+
  scale_color_brewer(palette = "Dark2")

DHAboxplot


DSAboxplot + EPAboxplot + DHAboxplot 

# Fatty Acid and Isotope RDA ----------------------------------------------


#run these once:
#install.packages("remotes")
#remotes::install_github("gavinsimpson/ggvegan")
setwd("/Users/logan/Stable Isotope Analysis")

library(stringr)
library(ggplot2)
library(vegan)
library(tidyverse)
library(ggpubr)


Merged$Site <- str_extract(Merged$ID, "[A-Z]+")


#reorder columns
#pivot longer delta values
delta_col_order <- c("ID", "Lipid.ID", "C12.0", "C14.0", "C15.0", "C16.0", "C16.1n.9", "C16.1n.7", "C17.0", "C17.1", "C18.0", "C18.1n.9", "C18.1n.7", "C18.2n.6", "C18.3n.3", "C18.4n.3", "C20.0", "C20.1.1", "C20.4n.6", "C20.5n.3", "C22.2n.9", "C22.4n.6", "C22.5n.6", "C22.5n.3", "C22.6n.3", "Season", "Site", "del15N", "del13C", "del34S", "del2H", "del18O", "perN", "perC", "perS", "perH", "perO", "sites")

CN_delta_col_order <- c("ID", "Lipid.ID", "C12.0", "C14.0", "C15.0", "C16.0", "C16.1n.9", "C16.1n.7", "C17.0", "C17.1", "C18.0", "C18.1n.9", "C18.1n.7", "C18.2n.6", "C18.3n.3", "C18.4n.3", "C20.0", "C20.1.1", "C20.4n.6", "C20.5n.3", "C22.2n.9", "C22.4n.6", "C22.5n.6", "C22.5n.3", "C22.6n.3", "Season", "Site", "del15N", "del13C", "perN", "perC", "perS", "perH", "perO", "sites")

Merged_edited_wide <- Merged[, delta_col_order]

Merged_CN_edited <- Merged[, CN_delta_col_order]


fa_rda<-Merged[,c("C12.0","C14.0","C15.0","C16.0","C16.1n.9","C16.1n.7","C17.0","C17.1","C18.0","C18.1n.9","C18.1n.7","C18.2n.6","C18.3n.3","C18.4n.3","C20.0","C20.1.1","C20.4n.6","C20.5n.3","C22.2n.9","C22.4n.6","C22.5n.6","C22.5n.3","C22.6n.3")]
pufa_fa_rda<-Merged[,c("C18.2n.6", "C18.3n.3", "C18.4n.3", "C20.4n.6", "C20.5n.3", "C22.2n.9", "C22.4n.6", "C22.5n.6", "C22.5n.3", "C22.6n.3")]

isotope_rda<-Merged[,c("del15N","del13C","del34S","del2H","del18O")]

fa_iso_rda<-rda(fa_rda~.,data=isotope_rda, scale = TRUE)
pufa_fa_iso_rda<-rda(pufa_fa_rda~.,data=isotope_rda)

pufa_scores <- as.data.frame(pufa_fa_iso_rda$CA$v)
pufa_scores$fattyacids <- row.names(pufa_scores)

isotope_scores <- as.data.frame(fa_iso_rda$CCA$biplot)
isotope_scores$isotopes <- row.names(isotope_scores)

fattyacid_scores <- as.data.frame(fa_iso_rda$CA$v)
fattyacid_scores$fattyacids <- row.names(fattyacid_scores)

site_scores <- as.data.frame(fa_iso_rda$CA$u)
site_scores$Site <- Merged$Site
site_scores$Season <- Merged$Season


#All fatty acid and isotope RDA plot
All_FA_RDA_plot <- ggplot()+
  geom_text(data = isotope_scores, aes(x = RDA1, y = RDA2, label = isotopes, size = 8), color = "PURPLE") +
  geom_text(data = fattyacid_scores, aes(x = PC1, y = PC2, label = fattyacids, size = 8), color = "darkGREEN") +
  #geom_text(size = 2.5, data = site_scores, aes(x = PC1, y = PC2, label = Site))+
  geom_point(size = 2.5, data = site_scores, aes(x = PC1, y = PC2, color = Site, shape = Season)) +
  theme_classic()+
  scale_color_brewer(palette = "Dark2")+
  theme(aspect.ratio = 1)

All_FA_RDA_plot

#PUFA and Isotope RDA plot
PUFA_RDA_plot <- ggplot()+
  geom_text(data = isotope_scores, aes(x = RDA1, y = RDA2, label = isotopes), color = "PURPLE") +
  geom_text(data = pufa_scores, aes(x = PC1, y = PC2, label = fattyacids), color = "darkGREEN") +
  #geom_text(size = 2.5, data = site_scores, aes(x = PC1, y = PC2, label = Site))
  geom_point(size = 2.5, data = site_scores, aes(x = PC1, y = PC2, color = Site, shape = Season))+
  theme_classic()+
  scale_color_brewer(palette = "Dark2")+
  theme(aspect.ratio = 1)

All_FA_RDA_plot + PUFA_RDA_plot

# LMM ---------------------------------------------------------------------

library(ggplot2)
library(GGally)
setwd("/Users/logan/Stable Isotope Analysis")

library(stringr)
library(ggplot2)
library(vegan)
library(tidyverse)
library(ggpubr)



delta_col_order <- c("ID", "Lipid.ID", "C12.0", "C14.0", "C15.0", "C16.0", "C16.1n.9", "C16.1n.7", "C17.0", "C17.1", "C18.0", "C18.1n.9", "C18.1n.7", "C18.2n.6", "C18.3n.3", "C18.4n.3", "C20.0", "C20.1.1", "C20.4n.6", "C20.5n.3", "C22.2n.9", "C22.4n.6", "C22.5n.6", "C22.5n.3", "C22.6n.3", "Season", "Site", "del15N", "del13C", "del34S", "del2H", "del18O", "perN", "perC", "perS", "perH", "perO")

#CN_delta_col_order <- c("ID", "Lipid.ID", "C12.0", "C14.0", "C15.0", "C16.0", "C16.1n.9", "C16.1n.7", "C17.0", "C17.1", "C18.0", "C18.1n.9", "C18.1n.7", "C18.2n.6", "C18.3n.3", "C18.4n.3", "C20.0", "C20.1.1", "C20.4n.6", "C20.5n.3", "C22.2n.9", "C22.4n.6", "C22.5n.6", "C22.5n.3", "C22.6n.3", "Season", "Site", "del15N", "del13C", "perN", "perC", "perS", "perH", "perO", "sites")

Merged_edited_wide <- Merged[, delta_col_order]

updown<-data.frame(Site=c("MDS","TDS","PDS","MUR","TUS","PUS"),DamPosition=c("Downstream","Downstream","Downstream","Upstream","Upstream","Upstream"), River = c("Miss", "Tipp", "Prairie_Creek", "Miss", "Tipp", "Prairie_Creek"))

Merged_edited_2 <- merge(Merged_edited_wide, updown, by = "Site", all.x = TRUE)

library(lmerTest)
library(emmeans)
library(broom)
library(stargazer)

#repeated measures mixed effects model.
#your sampling unit (each site) is the repeated measure and should be identified as the variable in the parenthesis, which sets the variable as random.

#DPA linear mixed model
lmerDPA <- lmer(C22.5n.3~DamPosition + Season + (1|Site), data = Merged_edited_2)
summary(lmerDPA)
plot(lmerDPA)



#EPA linear mixed model
lmerEPA <- lmer(C20.5n.3~DamPosition + Season + (1|Site), data = Merged_edited_2)
summary(lmerEPA)
plot(lmerEPA)

#DHA linear mixed model
lmerDHA <- lmer(C22.6n.3~DamPosition + Season + (1|Site), data = Merged_edited_2)
summary(lmerDHA)

#nitrogen linear mixed model
lmerdel15N <- lmer(del15N~DamPosition+Season+(1|Site), data=Merged_edited_2)
summary(lmerdel15N)
plot(lmerdel15N)


emmeans_season <- emmeans(lmerdel15N, "Season")
emmeans_damposition <- emmeans(lmerdel15N, "DamPosition")
contrast(emmeans_damposition, 'tukey')
contrast(emmeans_season, 'tukey')
plot(emmeans_season)
plot(emmeans_damposition)

#carbon linear mixed model
lmerdel13C <- lmer(del13C~DamPosition + Season + (1|Site), data = Merged_edited_2)
summary(lmerdel13C)
plot(lmerdel13C)

emmeans_season <- emmeans(lmerdel13C, "Season")
emmeans_damposition <- emmeans(lmerdel13C, "DamPosition")
contrast(emmeans_damposition, 'tukey')
contrast(emmeans_season, 'tukey')
plot(emmeans_season)

#oxygen linear mixed model 
lmerdel18O <- lmer(del18O~DamPosition + Season + (1|Site), data = Merged_edited_2)
summary(lmerdel18O)
plot(lmerdel18O)

emmeans_season <- emmeans(lmerdel18O, "Season")
emmeans_damposition <- emmeans(lmerdel18O, "DamPosition")
contrast(emmeans_damposition, 'tukey')
contrast(emmeans_season, 'tukey')

#hydrogen linear mixed model
lmerdel2H <- lmer(del2H~DamPosition + Season + (1|Site), data = Merged_edited_2)
summary(lmerdel2H)
plot(lmerdel2H)

emmeans_season <- emmeans(lmerdel2H, "Season")
emmeans_damposition <- emmeans(lmerdel2H, "DamPosition")
contrast(emmeans_damposition, 'tukey')
contrast(emmeans_season, 'tukey')

#sulfur linear mixed model
lmerdel34S <- lmer(del34S~DamPosition + Season + (1|Site), data = Merged_edited_2)
summary(lmerdel34S)
plot(lmerdel34S)

emmeans_season <- emmeans(lmerdel34S, "Season")
emmeans_damposition <- emmeans(lmerdel34S, "DamPosition")
contrast(emmeans_damposition, 'tukey')
contrast(emmeans_season, 'tukey')






