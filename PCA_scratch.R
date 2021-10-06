#dependcies
library(tidyverse)
library(ggpubr) 
library("FactoMineR")
library("factoextra")
library(plot3D)

#read in data
dat <- read_csv("MouseMasterDate_GXE.csv")
# get default color hex codes from ggpubr in case you need to share with someone.
get_palette("default", 4)
sum_df <- PCA_dat %>%   
  group_by(Strain, Diet, Week) %>%  
  summarise(glucose_mean = mean(Glucose),   
            glucose_se = sd(Glucose)/sqrt(n()), PercBWincrease_mean = mean(PercBWincrease), PercBWincrease_se = sd(PercBWincrease)/sqrt(n()),
            Ingest_bymouse_mean = mean(Ingest_bymouse), Ingest_bymouse_se = sd(Ingest_bymouse)/sqrt(n()), Trig_mean = mean(Trig), Trig_se = sd(Trig)/sqrt(n()),
            NEFA_mean = mean(NEFA), NEFA_se = sd(NEFA)/sqrt(n()), Insulin_mean = mean(Insulin), Insulin_se = sd(Insulin)/sqrt(n()))

#Figure 2, PCA
#get rid of average data we do not need
dat <- dat %>% select(-BW_Average, -Food_Avg, -Food_Mass)
# select columns to start setting up percent difference in weight
weight <- dat %>% select(`Mouse Cage`, `Mouse Number`, Strain, Diet, Week, PercBWincrease)
View(weight)
# take week 0 for weight to make a baseline
week0 <- weight %>% filter(Week == 0) %>% mutate(W0_mean = PercBWincrease) %>% select(W0_mean, `Mouse Number`, Strain, Diet)
View(week0)
#calculate percent difference in weight gain
weight1 <- weight %>% left_join(week0, by = c("Mouse Number", "Strain", "Diet")) %>% mutate(percdiff = abs((W0_mean - PercBWincrease)/W0_mean)*100)
View(weight1)
#check length of PercBWincrease to be sure the vector is the same length to replace raw BW data with percent increase data
length(dat$PercBWincrease)
length(weight1$percdiff)
dat$PercBWincrease <- weight1$percdiff
#need only week 4 and 8 as week 0 is a baseline, use week 0 only for line plots (see figure 4)
week0 <- dat %>% filter(Week == 0)
week4 <- dat %>% filter(Week == 4)
 week8 <- dat %>% filter(Week == 8)
# create PCA data by binding week 4 and 8 data 
# PCA for 2D plotting
PCA_dat_na <- week4 %>% bind_rows(week8)
# PCA for 3D plotting
PCA_dat <- week4 %>% bind_rows(week8)
# select all diets except for Chow for PCA
PCA_dat_na <- PCA_dat_na %>% filter(Diet != "Chow")
# create active numeric columns for PCA
# 2D plotting
active_na <- PCA_dat_na %>% select(PercBWincrease, Ingest_bymouse, Glucose, Trig, NEFA, Insulin)
# 3D plotting
active <- PCA_dat %>% select(PercBWincrease, Ingest_bymouse, Glucose, Trig, NEFA, Insulin)
#perform scaled PCA
res.pca <- PCA(active_na, scale.unit = TRUE, graph = FALSE)
# 2D plotting
# scree plot
fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))
# contributions of variables to PC1 and PC2
# extract contributions of variables to respective PCs in a table
var <- get_pca_var(res.pca)
var <- var$contrib
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
# begin to build biplot
fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca, geom.ind = "point", col.ind = PCA_dat$Diet)
# PCA biplot by diet, and strain
# also order diet appropriately 
PCA_dat_na$Diet <- factor(PCA_dat_na$Diet, levels = c("American", "Mediterranean", "Vegetarian", "Vegan"))
PCA_dat_na$Strain <- factor(PCA_dat_na$Strain, levels = c("C57BL/6J", "A/J", "DBA/2J", "S/J"))
fviz_pca_ind(res.pca, geom.ind = "point", col.ind = PCA_dat_na$Diet) + 
  theme(text = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20)) #, addEllipses = TRUE
fviz_pca_ind(res.pca, geom.ind = "point", col.ind = PCA_dat_na$Strain) + 
  theme(text = element_text(size = 20), axis.title = element_text(size = 20), axis.text = element_text(size = 20)) #, addEllipses = TRUE
fviz_pca_ind(res.pca, axes = c(2, 3), geom.ind = "point", col.ind = PCA_dat_na$Diet) +
  theme(text = element_text(size = 15), axis.title = element_text(size = 15), axis.text = element_text(size = 15)) #, addEllipses = TRUE #, addEllipses = TRUE
fviz_pca_ind(res.pca, axes = c(2, 3), geom.ind = "point", col.ind = PCA_dat_na$Strain) + 
  theme(text = element_text(size = 15), axis.title = element_text(size = 15), axis.text = element_text(size = 15)) #, addEllipses = TRUE#, addEllipses = TRUE

# 3D plotting
# also must create data frame with no na's invovling character to string to parse by diet and strain
# after PCA.
PCA_dat <- PCA_dat %>% filter(!is.na(PCA_dat$PercBWincrease))
PCA_dat <- PCA_dat %>% filter(!is.na(PCA_dat$Ingest_bymouse))
PCA_dat <- PCA_dat %>% filter(!is.na(PCA_dat$Glucose))
PCA_dat <- PCA_dat %>% filter(!is.na(PCA_dat$Trig))
PCA_dat <- PCA_dat %>% filter(!is.na(PCA_dat$NEFA))
PCA_dat <- PCA_dat %>% filter(!is.na(PCA_dat$Insulin))

active_prcomp <- as.data.frame(active)
active_prcomp <- na.omit(active_prcomp)
PCA_3d <- prcomp(active_prcomp[, 1:6], scale. = TRUE)
# this works
# read in numerica data so diet and strain are a numeric vector, this is the only way scatter3D
# can read classes for color plotting (what year is it?)
PCA_dat_numeric <- read_csv("PCA_dat_numeric.csv")
PCs <- as.data.frame(PCA_3d$x)
scatter3D(x = PCs$PC1, y = PCs$PC2, z = PCs$PC3, main = "3D Biplot of GXE Mice by Strain", 
          xlab = "PC1, 37.5%, %WeightGain", 
          ylab = "PC2,22.5%,Tri", 
          zlab = "PC3,14.2%,Glucose", 
          colvar = PCA_dat_numeric$Strain_number) # try this once you have PCA_dat again

# Barplots (Figure 3)
# make diets in order requested

fin <- PCA_dat %>% filter(Week == 8)
fin$Diet <- factor(fin$Diet, levels = c("American", "Mediterranean", "Vegetarian", "Vegan"))

glucose <- ggbarplot(fin, x = "Strain", y = "Glucose", color = "Diet", fill = "Diet", title = "Glucose (mg/dL)", xlab = "Strain", ylab = "glucose (mg/dL)", position = position_dodge(0.9), add = "mean_se")
glucose <- glucose + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
trig <- ggbarplot(fin, x = "Strain", y = "Trig", color = "Diet", fill = "Diet", title = "Triglycerides (mg/dL)", xlab = "Strain", ylab = "triglyceride (mg/dL)", position = position_dodge(0.9), add = "mean_se")
trig <- trig + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
percdiff <- ggbarplot(fin, x = "Strain", y = "PercBWincrease", color = "Diet", fill = "Diet", title = "Body Weight Percent Increase", xlab = "Strain", ylab = "Percent BW From Week 0", position = position_dodge(0.9), add = "mean_se")
percdiff <- percdiff + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
insulin <- ggbarplot(fin, x = "Strain", y = "Insulin", color = "Diet", fill = "Diet", title = "Insulin", xlab = "Strain", ylab = "Insulin (need units)", position = position_dodge(0.9), add = "mean_se")
insulin <- insulin + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
nefa <- ggbarplot(fin, x = "Strain", y = "NEFA", color = "Diet", fill = "Diet", title = "NEFA (mEq/L)", xlab = "Strain", ylab = "NEFA (mEq/L)", position = position_dodge(0.9), add = "mean_se")
nefa <- nefa + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
topPCs <- ggarrange(glucose, trig, percdiff, nefa,ncol = 2, nrow =2, common.legend = TRUE)
annotate_figure(topPCs, top = text_grob("Metabolic Phenotypes at Week 8"))
# look at ingestion of food by mouse.
ingest <- ggbarplot(dat, x = "Strain", y = "Ingest_bymouse", color = "Diet", fill = "Diet", title = "Ingestion of Food By Mouse (g/day) at Week 4", xlab = "Strain", ylab = "g/day", position = position_dodge(0.9), add = "mean_se")
ingest <- ingest + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")

# line plot of weight gain (Figure 4)
dat$Diet <- factor(dat$Diet, levels = c("American", "Mediterranean", "Vegetarian", "Vegan"))
dat$Strain <- factor(dat$Strain, levels = c("C57BL/6J", "A/J", "DBA/2J", "S/J"))
line <- ggline(dat, x = "Week", y = "PercBWincrease", color = "Diet", title = "Percent Body Weight Gain over Week 0", xlab = "Week", ylab = "Percent Weight Gain", add = "mean_se", facet.by = "Strain")
line + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")

#line plot of top variables by time (Supplementary Figure 4)
PCA_dat$Diet <- factor(PCA_dat$Diet, levels = c("American", "Mediterranean", "Vegetarian", "Vegan"))
PCA_dat$Strain <- factor(PCA_dat$Strain, levels = c("C57BL/6J", "A/J", "DBA/2J", "S/J"))
glucose <- ggline(PCA_dat, x = "Week", y = "Glucose", color = "Diet", title = "Glucose (mg/dL)", xlab = "Week", ylab = "glucose (mg/dL)", position = position_dodge(0.9), add = "mean_se", facet.by = "Strain")
glucose <- glucose + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
trig <- ggline(PCA_dat, x = "Week", y = "Trig", color = "Diet", title = "Triglycerides (mg/dL)", xlab = "Week", ylab = "triglycerides (mg/dL)", position = position_dodge(0.9), add = "mean_se", facet.by = "Strain")
trig <- trig + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
percdiff <- ggline(PCA_dat, x = "Week", y = "PercBWincrease", color = "Diet", title = "Body Weight Percent Increase", xlab = "Week", ylab = "Percent BW Gain", position = position_dodge(0.9), add = "mean_se", facet.by = "Strain")
percdiff <- percdiff + font("title", size = 15, face = "bold") + font("subtitle", size = 15) + font("xlab", size = 15, face = "bold") + font("ylab", size = 15, face = "bold")+ stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
insulin <- ggline(PCA_dat, x = "Week", y = "Insulin", color = "Diet", title = "Insulin (ng/mL)", xlab = "Week", ylab = "Insulin (ng/mL)", position = position_dodge(0.9), add = "mean_se", facet.by = "Strain", scales = "free_y")
insulin <- insulin + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
nefa <- ggline(PCA_dat, x = "Week", y = "NEFA", color = "Diet", title = "NEFA (mg/dL)", xlab = "Week", ylab = "NEFA (mg/dL)", position = position_dodge(0.9), add = "mean_se", facet.by = "Strain")
nefa <- nefa + stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif")
ingest <- ggline(PCA_dat, x = "Week", y = "Ingest_bymouse", color = "Diet",  title = "Ingestion of Food By Mouse", xlab = "Week", ylab = "g/g BW/day", position = position_dodge(0.9), add = "mean_se", facet.by = "Strain")
ingest <- ingest + font("title", size = 15, face = "bold") + font("subtitle", size = 15) + font("xlab", size = 15, face = "bold") + font("ylab", size = 15, face = "bold")+ stat_compare_means(aes(group = Diet), method = "anova", label = "p.signif", label.y = 0.175)
topPCs <- ggarrange(glucose, trig, percdiff, nefa,ncol = 2, nrow =2, common.legend = TRUE)
annotate_figure(topPCs, top = text_grob("Metabolic Phenotypes through Time"))
#health outcome plot
diffeffect <- ggarrange(line, insulin, widths = c(2, 1), common.legend = TRUE)
annotate_figure(diffeffect, top = text_grob("Weight Gain Does Not Always Negatively Effect Health"))

# two-way anova for all phenotypes
# body weight percent gain
require(emmeans)
lm <- lm(dat$PercBWincrease ~ dat$Diet * dat$Strain)
ml <- emmeans(lm, ~ Diet | Strain)
ml
pairs(ml)
lm <- lm(dat$Glucose ~ dat$Diet * dat$Strain)
ml <- emmeans(lm, ~ Diet | Strain)
ml
pairs(ml)
lm <- lm(dat$Trig ~ dat$Diet * dat$Strain)
ml <- emmeans(lm, ~ Diet | Strain)
ml
pairs(ml)
lm <- lm(dat$Insulin ~ dat$Diet * dat$Strain)
ml <- emmeans(lm, ~ Diet | Strain)
ml
pairs(ml)
lm <- lm(dat$NEFA ~ dat$Diet * dat$Strain)
ml <- emmeans(lm, ~ Diet | Strain)
ml
pairs(ml)
lm <- lm(dat$Ingest_bymouse ~ dat$Diet * dat$Strain)
ml <- emmeans(lm, ~ Diet | Strain)
ml
pairs(ml)

# Figure 5, "Health Score"
# read in data with strain and diet as numeric 
PCA_dat_numeric <- read_csv("PCA_dat_numeric.csv")
# calculate multi-factor ANOVA
lm <- lm(PCA_dat_numeric$Diet_numbe ~ PCA_dat_numeric$Strain_number + 
           PCA_dat_numeric$PercBWincrease + 
           PCA_dat_numeric$Glucose + 
           PCA_dat_numeric$Trig + PCA_dat_numeric$NEFA)
aov <- aov(lm)
# dependency for aov effect sizes
require(lsr)
# calculate effect sizes for aov above
etaSquared(aov, type = 2, anova = FALSE)

# plot ingestion through time with all strains present.

PCA_dat <- PCA_dat %>% select(Strain, Diet, Week, Ingest_bymouse)
PCA_dat <- PCA_dat %>% group_by(Strain, Diet, Week)
PCA_dat <- na.omit(PCA_dat)
PCA_dat$Week <- as.factor(PCA_dat$Week)

PCA_dat <- PCA_dat %>% 
  summarise(ingest_mean = mean(Ingest_bymouse), Ingest.min = ingest_mean - sd(Ingest_bymouse)/sqrt(length(Ingest_bymouse)), Ingest.max =  ingest_mean + sd(Ingest_bymouse)/sqrt(length(Ingest_bymouse)))

PCA_dat$Diet <- factor(PCA_dat$Diet, levels = c("American", "Mediterranean", "Vegetarian", "Vegan"))
PCA_dat$Strain <- factor(PCA_dat$Strain, levels = c("C57BL/6J", "A/J", "DBA/2J", "S/J"))

p <- ggplot(PCA_dat, aes(x = Week, y= ingest_mean, color = Strain)) + 
  labs(title = "Ingestion of Diet by Strain through Time", x = "Week", y = "g/g BW/day") 

p + geom_point(aes(shape = Strain)) + 
  geom_line(aes(linetype = Diet)) + theme_bw()  + 
  geom_errorbar( aes(ymin= Ingest.min, ymax = Ingest.max)) + theme(title = element_text(size = 15, face = "bold"), axis.title.x = element_text(size = 15, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"))



