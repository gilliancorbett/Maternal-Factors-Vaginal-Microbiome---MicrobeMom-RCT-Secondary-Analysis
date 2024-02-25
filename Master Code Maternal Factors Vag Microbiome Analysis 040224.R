library(vegan)
library(dplyr)
library(tidyverse)
library(vegan)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggnewscale)
library (tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RVAideMemoire)
library(ggpubr)
library(Hmisc)
library(scales)
library(matrixStats)
library (tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RVAideMemoire)
library(ggpubr)
library(Hmisc)
library(scales)
library(matrixStats)
library(car)
library(quantmod)
library(patchwork)
library(cowplot)
library(BiodiversityR) 
library(ggplot2)
library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)
library(BiodiversityR)
library(viridis)
library(RColorBrewer)
library(ggpubr)
library (tidyverse)
library(vegan)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RVAideMemoire)
library(ggpubr)
library(Hmisc)
library(scales)
library(matrixStats)
library(car)
library(quantmod)
library(patchwork)
library(cowplot)
library(BiodiversityR) 
library(ggplot2)
library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)
library(BiodiversityR)
# Figure 2 - Species NMDS Ordination, envfit and RDA Analysis ------------------------
#Import species and metadata datasets 

MM_Selected_Metadata_030224 <- read_excel("Desktop/Documents/PhD/Low Risk Vaginal Microbial Profile - MicroMOM Analysis/MM Datasets/MM_Early_Vaginal_Jan24_Analysis/MM_Selected_Metadata_030224.xlsx")
View(MM_Selected_Metadata_030224)

MM_Vaginal_Species_Jan24 <- read_excel("Desktop/Documents/PhD/Low Risk Vaginal Microbial Profile - MicroMOM Analysis/MM Datasets/MM_Early_Vaginal_Jan24_Analysis/MM_Vaginal_Species_Jan24.xlsx")
View(MM_Vaginal_Species_Jan24)

#nmds of species dataframe 
df = MM_Selected_Metadata_030224
#create metadata set excluding CST 
metadataSpec <- df[, 3:24]
View (metadataSpec)
df2 = MM_Vaginal_Species_Jan24
com = df2[,3:63] #Remove Sample ID and CST to leave only species columns

#Remove any columns where all samples have zero of that species 
com <- com %>% 
  select_if(negate(function(col) is.numeric(col) && sum(col) < 0.000001))

#Create matrix, set seed and run nmds of species 
m_com = as.matrix(com)
set.seed(123)
nmds_spec = metaMDS(m_com, distance = "bray", autotransform = FALSE)
nmds_spec
plot(nmds) #Record stress no. 0.0672

#Extract nmds coordinates into a dataframe and add CST for colour coding on plots 
data.scoresSpec = as.data.frame(scores(nmds_spec)$sites) #Extract sites values 
data.scoresSpec$CST = df2$CST #Add CST column for ggplot labelling 
#Order the CST column to allow legend order
data.scoresSpec$CST <- factor(data.scoresSpec$CST, levels = c("CST I", "CST II", 
                                                              "CST III", "CST IV", "CST V", "CST VIII", 
                                                              "CST IX"))
View(data.scoresSpec) #ensure no NAs in CST column 

#Stats for NMDS_spec Model 
anoSpec = anosim(m_com, df2$CST, distance = "bray", permutations = 9999)
anoSpec
CST.anoSpec <- with(df, anosim(m_com, df2$CST, distance = "bray", permutations = 9999, strata = NULL, parallel = getOption("mc.cores")))
summary(CST.anoSpec) #ANOSIM R0.9967, p<0.001
plot(CST.anoSpec)  #Supplementatry File - shows distance ranking of each CST Group for NMDS Model 

#Plot for nmds ordination 
species_ordination = ggplot(data.scoresSpec, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = df2$CST)) + 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  annotate("text", size = 3.0, x=7.0, y=2.5, label= " Stress = 0.0672") +
  annotate("text", size = 3.0, x=7.0, y=1.5, label= " ANOSIM R = 0.997 
           p<0.001") + ##Add NMDS Stress ANOSIM Stat R 0.9967, p <0.0001, stress no. 0.06722932
  labs(x = "NMDS1", colour = "CST", y = "NMDS2")  +
  stat_ellipse(type = "t", level = 0.95, aes(color = df2$CST))
species_ordination 

#Envfit of NMDS Species against env covariates

envfit_species <- envfit(nmds_spec, metadataSpec, perm=1000, na.rm=TRUE)
envfit_species
str(envfit_species)

#Plotting the envfit model on species NMDS Ordination plot
env <- df
env2 <- env[,3:24] #Remove Sample ID and CST
View(env)
View(env2) #No sample ID or CST 
en = envfit(nmds_spec, env2, permutations = 999, na.rm = TRUE)
en
plot(nmds_spec)
plot(en)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scoresEnv = as.data.frame(scores(nmds_spec)$sites)
#add 'season' column as before (from example) - just change df$CST_Group to df$class
data.scoresEnv$CST = env$CST

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#Filter covariates to just plot those with p<0.05. 
significant_vectors <- names(envfit_species$vectors$pvals)[envfit_species$vectors$pvals < 0.05]

#Create ggplot of envfit values 
ggSpecEnvfit <- ggplot(data = data.scoresEnv, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = CST), size = 3, alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = as.data.frame(scores(envfit_species, "vectors")),
               size = 2, alpha = 0.5, colour = "grey30") +
  geom_point(data = as.data.frame(scores(envfit_species, "vectors")),
             aes(x = NMDS1, y = NMDS2), shape = "diamond", size = 2, alpha = 0.6,
             colour = "grey30") +
  geom_text(data = as.data.frame(scores(envfit_species, "vectors"))[significant_vectors, ],
            aes(x = NMDS1, y = NMDS2 + 0.04), size = 3.5,
            label = significant_vectors,
            colour = "grey30", fontface = "bold") +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 12),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black", face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  annotate("text", x=5.0, y=1.0, label= " Leucine (Adj-R 0.086)
           Lysine (Adj-R 0.114) 
           Phenylalanine (Adj-R 0.085)
           Valine (Adj-R 0.004)
           Previous Breastfeeding (Adj-R 0.048)
           Previous Vaginal Delivery (Adj-R 0.048)") +
  labs(x = "NMDS1", colour = "CST", y = "NMDS2")
ggSpecEnvfit


#to Create horizontal barchart 

envfit_species ## enShort  #Generated above 

str(envfit_species)

#Make tibble of R and p values for Vectors then Factors 
envfit_species_Vectors <- tibble((envfit_species$vectors)$r, (envfit_species$vectors)$pvals, )
envfit_species_Factors <- tibble((envfit_species$factors)$r, (envfit_species$factors)$pvals)

##How to keep covariate names for these

#Rename columns as R-Squared and p value
names(envfit_species_Vectors)[1] <- "R2"
names(envfit_species_Vectors)[2] <- "pvalue"
names(envfit_species_Factors)[1] <- "R2"
names(envfit_species_Factors)[2] <- "pvalue"

View(envfit_species_Vectors)
View(envfit_species_Factors)

#Create this table as an object of all R2 and p Values for covariates in species envfit model 
EnvFitSpeciesMetadata <- (rbind(envfit_species_Vectors, envfit_species_Factors)) #Covariate Row names lost before this 
View(EnvFitSpeciesMetadata)

#compare both envfit results and your R2/pval dataframe to ensure correct
envfit_species
View(EnvFitSpeciesMetadataShort)

#get list from running enShort to display results 
covariate_categories <- data.frame("Covariate" =c('Maternal Age',
                                                  'Education', 
                                                  'HP Index', 
                                                  'BMI', 
                                                  'MET per week', 
                                                  'WHO Wellbeing Score', 
                                                  'Fat',
                                                  'Starch',
                                                  'Sugars',
                                                  'NSP',
                                                  'Vitamin K1',
                                                  'Niacin',
                                                  'Selenium',
                                                  'Leucine',
                                                  'Lysine',
                                                  'Methionine',
                                                  'Phenylalanine',
                                                  'Valine', 
                                                  'Ethnicity', 
                                                  'Prev BF', 
                                                  'Prev VD', 
                                                  'Smoking'),
                                   "Covariate type" = c('Demography',
                                                        'Demography',
                                                        'Demography',
                                                        'Demography',
                                                        'Lifestyle',
                                                        'Lifestyle',
                                                        'Macronutrient',
                                                        'Macronutrient',
                                                        'Macronutrient',
                                                        'Macronutrient',
                                                        'Micronutrient',
                                                        'Micronutrient',
                                                        'Micronutrient',
                                                        'Amino Acid',
                                                        'Amino Acid',
                                                        'Amino Acid',
                                                        'Amino Acid',
                                                        'Amino Acid',
                                                        'Demography',
                                                        'Obstetric', 
                                                        'Obstetric', 
                                                        'Lifestyle'))

envfit_spec_results <- cbind(covariate_categories,EnvFitSpeciesMetadata) 
View(envfit_spec_results)
#Confirm these are correct corresponding covariate names and r/p by looking at enShort result 
envfit_species
#Relabel R2 as Vaginal beta diversity for labelling 
colnames(envfit_spec_results)<- c('Covariate',
                        'Covariate_type',
                        'Vaginal beta diversity', 
                        'p-value')
names(envfit_spec_results)

envfitR2_long <- gather(envfit_spec_results, Sample_type, R2, 'Vaginal beta diversity', factor_key=TRUE)
envfitR2_long$Covariate_type = factor(envfitR2_long$Covariate_type, levels=c('Demography', 'Lifestyle', 'Obstetric', 'Macronutrient','Micronutrient', 'Amino Acid'))
envfitR2_long$Sample_type = factor(envfitR2_long$Sample_type, levels=c('Vaginal beta diversity'))
head(envfitR2_long)

#Horizontal Barchart ggplot 
ggEnvfitSpecChart <- ggplot(envfitR2_long, aes(x = fct_reorder(Covariate, +R2, max), y = R2, fill = Covariate_type, label = ifelse(envfitR2_long$`p-value` < 0.05, "*", "")))+
  facet_wrap(~Sample_type, nrow = 2) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0) +
  scale_color_brewer(palette = "Set2")+
  theme (axis.text.x = element_text(size = 7, colour = "black", angle = 45)) +
  theme (axis.text.y = element_text(size = 7, colour = "black"))+
  theme(axis.title = element_text(size = 15)) +
  theme(legend.position = "bottom",
        legend.text = element_text(color = "black", size = 12)) +
  theme(strip.text = element_text(face="bold", size=10),
        strip.background = element_rect(linewidth=2)) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        linewidth = 0.25, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white")) +
  xlab('')+ #Removed Axis label 
  ylab(bquote(''~R^2)) +
  labs(color=NULL) + 
  coord_flip()+ 
  theme(plot.title = element_textbox(size =15, hjust = 0.5, margin = margin(t = 5, b = 5)))
ggEnvfitSpecChart

### RDA of species ordination ###

#examine species nmds first nmds_spec
#species data object with Sample ID / CST Removed - just species data, as created above 
com
#Create Environmental MetadataRDA object where all covariates are numeric
MM_Selected_Metadata_RDA_030224 <- read_excel("Desktop/Documents/PhD/Low Risk Vaginal Microbial Profile - MicroMOM Analysis/MM Datasets/MM_Early_Vaginal_Jan24_Analysis/MM_Selected_Metadata_RDA_030224.xlsx")
View(MM_Selected_Metadata_RDA_030224)

metadataRDA <- MM_Selected_Metadata_RDA_030224
metadataRDA <- metadataRDA[, 3:24] #Remove Sample ID and CST Group 
#Explore species dataframe
names(com)
dim(com)  # dataset dimensions
head(com)  # look at first 5 rows
str(com)  # structure of objects in dataset
summary(com)  # summary statistics for all objects (min, mean, max, etc.)
# Count number of species frequencies in each abundance class
ab <- table(unlist(com))
# Plot distribution of species frequencies
barplot(ab, las = 1, # make axis labels perpendicular to axis
        xlab = "Abundance class", ylab = "Frequency", # label axes
        col = grey(5:0/5)) # 5-colour gradient for the bars
# Count the number of zeros in the dataset
sum(com == 0)
# Calculate proportion of zeros in the dataset = 96.5%
sum(com == 0)/(nrow(com) * ncol(com))

# Apply Hellinger transformation to correct for the double zero problem (The Hellinger transformation expresses abundances as the square-root of their relative abundance at each site (Borcard, Gillet, and Legendre 2011), solving the issue with double zeros)
com.hel <- decostand(com, method = "hellinger")

#Explore environmental dataset 
names(metadataRDA)
dim(metadataRDA)
head(metadataRDA)
str(metadataRDA)
summary(metadataRDA)

# good idea to check for correlations between variables, as the constrained ordination methods we will be using are highly sensitive to collinearities in the explanatory matrix.
# We can visually look for correlations between variables:

#Species covariate heatmap for supplementary File 
heatmap(abs(cor(metadataRDA)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topleft", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))

#From this you can observe colinearity - alot with amino acids and energy/macronutrients, protein&Amino acids, B vitamin group 
#Did this step prior and removed colinear covariates from raw data set 

#Standardise environmental variable 
# Scale and center variables
metadataRDA.z <- decostand(metadataRDA, method = "standardize")
# Variables are now centered around a mean of 0
round(apply(metadataRDA.z, 2, mean), 1)
# and scaled to have a standard deviation of 1
apply(metadataRDA.z, 2, sd)
#Now we have centred standardised our quantitative variables. Should be able to add back in binary/qualitative variables 

# Model the effect of all environmental variables on fish
# community composition
com.rda <- rda(com.hel ~ ., data = metadataRDA.z)
summary(com.rda)
sqrt(vif.cca(com.rda)) #Check colinearity 
#Interpreting results - Constrained Proportion: variance of  Yexplained by  X (16.7%) - he included environmental variables explain 40.71% of the variation in vaginal community composition across site
#Examine the RDA Model 
global_r2 <- RsquareAdj(com.rda)$adj.r.squared
global_r2 #r -0.0254 - Model is not significant - no forward selection performed 

### Plotting the Species RDA 
spec_rda <- rda(com.hel ~ metadataRDA$Lysine + metadataRDA$Valine + metadataRDA$Leucine + metadata$Phenylalanine + metadataRDA$Prev_VD + metadataRDA$Prev_BF)
summary(spec_rda)
plot(spec_rda)
RsquareAdj(spec_rda) #-0.027
anova.cca(spec_rda) #p 0.977 
anova.cca(spec_rda, by = "axis") # Nil significant
test2 <- anova.cca(spec_rda, by = "terms") # Nil significant 
sqrt(vif.cca(spec_rda)) #No colinearity all <2 

#### Creating Spec RDA Plot 
# To plot data via ggplot2, information on the locations of sites (circles in the ordiplot) is obtained via function sites.long. Information on species is extracted by function species.long.
spec_rda_plot1 <- ordiplot(spec_rda, scaling = 1, type = "text")
sites.longSpec <- sites.long(spec_rda_plot1, env.data=metadataRDA)
head(sites.longSpec)
species.longSpec <- species.long(spec_rda_plot1)
species.longSpec
#Information on the labeling of the axes is obtained with function axis.long. This information is extracted from the ordination model and not the ordiplot, hence it is important to select the same axes via argument ‘choices’. The extracted information includes information on the percentage of explained variation. I suggest that you cross-check with the summary of the redundancy analysis, where information on proportion explained is given.
Ordination.modelSpec <- rda(com.hel ~ ., data = metadataRDA, scaling = "species")
axis.longSpec <- axis.long(Ordination.modelSpec, choices=c(1, 2))
axis.longSpec ### Gives x axis RDA 1 8.5% and y axis RDA 2 2.6% 

library(BiodiversityR)
library(ggThemeAssist) #Find code for BioR Theme at: Page 7 of https://cran.r-project.org/web/packages/BiodiversityR/BiodiversityR.pdf
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

species.longSpec$p<-as.numeric(test2$`Pr(>F)`)
species.longSpec<-subset(species.longSpec, p < 0.05)

##Species RDA Plot 
plotRDA_Spec <- ggplot() + 
  xlab(axis.longSpec[1, "label"]) +
  ylab(axis.longSpec[2, "label"]) +  
  labs(colour = "CST")  +
  annotate("text", x=-1.0, y=0.4, label= " Species RDA p = 0.807") +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  geom_point(data=sites.longSpec, 
             aes(x=axis1, y=axis2, colour=df2$CST), 
             size=5) +
  geom_point(data=species.longSpec, 
             aes(x=axis1, y=axis2)) + 
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") + 
  geom_segment(data=species.longSpec, aes(x=0, xend=axis1, y=0, yend=axis2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) 

plotRDA_Spec

#Summary of Plots 
species_ordination
ggSpecEnvfit
ggEnvfitSpecChart
plotRDA_Spec

#Assimilating all species ordination, envfit and RDA plots 
Figure2 <- ggarrange(species_ordination, ggSpecEnvfit, plotRDA_Spec, ggEnvfitSpecChart, 
                      labels = c("A", "B", "C", "D"),
                      common.legend = TRUE, legend = "left",
                      ncol = 2, nrow = 2)
Figure2 


# Biological Processes Functional Analysis  ----------------------------------------------

#Create bp functional dataframe
BP_GO_Data <- read_excel("Desktop/Documents/PhD/Low Risk Vaginal Microbial Profile - MicroMOM Analysis/MM Datasets/MM_Early_Vaginal_Jan24_Analysis/BP_GO_Data.xlsx")
View(BP_GO_Data)
bp = BP_GO_Data 
bp = bp[,2:1090] #Remove Sample ID from BP data

#Create bp metadata dataframe
bp_metadata1 <- MM_Selected_Metadata_030224 #All 22 covariates with Sample ID and CST
bp_metadata = bp_metadata1[, 2:24] #Remove Sample ID, Keep CST

#Perform bp nmds ordination 
#Remove any columns where all samples have zero of that function
bp <- bp %>% 
  select_if(negate(function(col) is.numeric(col) && sum(col) < 0.000001))

m_bp = as.matrix(bp) #Create bp Matrix to run NMDS
set.seed(123)
nmds_bp = metaMDS(m_bp, distance = "bray", autotransform = FALSE)

nmds_bp
plot(nmds_bp) #Record stress no. 0.1098758 

data.scoresBP = as.data.frame(scores(nmds_bp)$sites) #Creates object of BP-nmds coordination points 
data.scoresBP$CST = bp_metadata1$CST #Adds CST labels 

#If you want the legend ordered
data.scoresBP$CST <- factor(data.scoresBP$CST, levels = c("CST I", "CST II", 
                                                          "CST III", "CST IV", "CST V", "CST VIII", 
                                                          "CST IX"))
View(data.scoresBP) #Ensure no NA Values 

#Stats for NMDS Biological Processes vs Metadata ordination 
anoBP = anosim(m_bp, bp_metadata1$CST, distance = "bray", permutations = 9999)
anoBP
#Plot for NMDS Biological Processes vs Metadata ordination 
CST.anoBP <- with(bp_metadata, anosim(m_bp, bp_metadata1$CST, distance = "bray", permutations = 9999, strata = NULL, parallel = getOption("mc.cores")))
summary(CST.anoBP) #ANOSIM 0.885 p <0.0001
plot(CST.anoBP) #Shows distance ranking per CST Group 

#How to add text to a plot - https://www.statology.org/add-text-to-ggplot/
#Creating plot of envfit for large meta-data set against BP nmds plot 
bp_ordination = ggplot(data.scoresBP, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = data.scoresBP$CST)) + 
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  annotate("text", size = 3.0, x=-1.0, y=0.0, label= " Stress = 0.110") +
  annotate("text", size = 3.0, x=-1.0, y=-0.1, label= " ANOSIM R = 0.885, 
           p <0.0001 ") + ##Add NMDS Stress ANOSIM Stat R 0.9967, p <0.0001, stress no. 0.06722932
  labs(x = "NMDS1", colour = "CST", y = "NMDS2")  +
  ggtitle( "Biological Process Functions", size = 12) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_ellipse(type = "t", level = 0.95, aes(color = data.scoresBP$CST))
bp_ordination


#Envfit of BP Functional Data against metadata covariates  
nmds_bp
envfit_bp <- envfit(nmds_bp, bp_metadata, perm=1000, na.rm=TRUE)
envfit_bp   
str(envfit_bp)
plot(envfit_bp)

#to annotate envfit_BP plot
envfit_BP_coord_cont = as.data.frame(scores(envfit_bp, "vectors")) * ordiArrowMul(envfit_bp)
envfit_BP_coord_cat = as.data.frame(scores(envfit_bp, "factors")) * ordiArrowMul(envfit_bp)

significant_bp_vectors <- names(envfit_bp$vectors$pvals)[envfit_bp$vectors$pvals < 0.05]

ggBPEnvfit <- ggplot(data = data.scoresBP, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = CST), size = 3, alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = as.data.frame(scores(envfit_bp, "vectors")),
               size = 2, alpha = 0.5, colour = "grey30") +
  geom_point(data = as.data.frame(scores(envfit_bp, "vectors")),
             aes(x = NMDS1, y = NMDS2), shape = "diamond", size = 2, alpha = 0.6,
             colour = "grey30") +
  geom_text(data = as.data.frame(scores(envfit_bp, "vectors"))[significant_bp_vectors, ],
            aes(x = NMDS1, y = NMDS2 + 0.04), size = 3.5,
            label = significant_bp_vectors,
            colour = "grey30", fontface = "bold") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 10),
        legend.text = element_text(size = 10, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  annotate("text", size = 3.0, x=-1.2, y=0.2, label= " CST (Adj-R 0.502)
           Smoking (Adj-R 0.042) 
           Prev BF (Adj-R 0.030)
           Prev VD (Adj-R 0.034)
           Vit K1 (Adj-R 0.047)") +
  labs(x = "NMDS1", colour = "CST", y = "NMDS2")
ggBPEnvfit


#to Create horizontal barchart - START HERE AGAIN!! went down to my_funBP but the code is wrong

envfit_bp  #Generated above 
str(envBPShort)

#Make tibble of R and p values for Vectors then Factors 
envfit_bp_Vectors <- tibble((envfit_bp$vectors)$r, (envfit_bp$vectors)$pvals, )
envfit_bp_Factors <- tibble((envfit_bp$factors)$r, (envfit_bp$factors)$pvals)

##How to keep covariate names for these - answer: matched below 

#Rename columns as R-Squared and p value
names(envfit_bp_Vectors)[1] <- "R2"
names(envfit_bp_Vectors)[2] <- "pvalue"
View(envfit_bp_Vectors)

#Repeat for factors 
names(envfit_bp_Factors)[1] <- "R2"
names(envfit_bp_Factors)[2] <- "pvalue"
View(envfit_bp_Factors)

#Create this table as an object of all R2 and p Values for covariates in species envfit model 
EnvFitBPMetadata <- (rbind(envfit_bp_Vectors, envfit_bp_Factors)) #Covariate Row names lost before this 
View(EnvFitBPMetadata)

#compare both envfit results and your R2/pval dataframe to ensure correct
envfit_bp
View(EnvFitBPMetadata)

#get list from running enShort to display results 
covariate_bp_categories <- data.frame("Covariate" =c('Maternal Age',
                                                  'Education', 
                                                  'HP Index', 
                                                  'BMI', 
                                                  'MET per week', 
                                                  'WHO Wellbeing Score', 
                                                  'Fat',
                                                  'Starch',
                                                  'Sugars',
                                                  'NSP',
                                                  'Vitamin K1',
                                                  'Niacin',
                                                  'Selenium',
                                                  'Leucine',
                                                  'Lysine',
                                                  'Methionine',
                                                  'Phenylalanine',
                                                  'Valine', 
                                                  'CST',
                                                  'Ethnicity', 
                                                  'Prev BF', 
                                                  'Prev VD', 
                                                  'Smoking'),
                                   "Covariate type" = c('Demography',
                                                        'Demography',
                                                        'Demography',
                                                        'Demography',
                                                        'Lifestyle',
                                                        'Lifestyle',
                                                        'Macronutrient',
                                                        'Macronutrient',
                                                        'Macronutrient',
                                                        'Macronutrient',
                                                        'Micronutrient',
                                                        'Micronutrient',
                                                        'Micronutrient',
                                                        'Amino Acid',
                                                        'Amino Acid',
                                                        'Amino Acid',
                                                        'Amino Acid',
                                                        'Amino Acid',
                                                        'Species',
                                                        'Demography',
                                                        'Obstetric', 
                                                        'Obstetric', 
                                                        'Lifestyle'))

envfit_bp_results <- cbind(covariate_bp_categories,EnvFitBPMetadata) 
View(envfit_bp_results)
#Confirm these are correct corresponding covariate names and r/p by looking at enShort result 
envfit_bp
#Relabel R2 as Vaginal beta diversity for labelling 
colnames(envfit_bp_results)<- c('Covariate',
                                  'Covariate_type',
                                  'Biological Process', 
                                  'p-value')
names(envfit_bp_results)
envfitR2_longBP <- gather(envfit_bp_results, Sample_type, R2, 'Biological Process', factor_key=TRUE)
envfitR2_longBP$Covariate_type = factor(envfitR2_longBP$Covariate_type, levels=c('Demography', 'Lifestyle', 'Obstetric', 'Macronutrient','Micronutrient', 'Amino Acid', 'Species'))
envfitR2_longBP$Sample_type = factor(envfitR2_longBP$Sample_type, levels=c('Biological Process'))
head(envfitR2_longBP)

#Horizontal Barchart ggplot for BP Function and metadata 
ggEnvfitBPChart <- ggplot(envfitR2_longBP, aes(x = fct_reorder(Covariate, +R2, max), y = R2, fill = Covariate_type, label = ifelse(envfitR2_longBP$`p-value` < 0.05, "*", "")))+
  facet_wrap(~Sample_type, nrow = 2) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0) +
  scale_color_brewer(palette = "Set2")+
  theme (axis.text.x = element_text(size = 7, colour = "black", angle = 45)) +
  theme (axis.text.y = element_text(size = 7, colour = "black"))+
  theme(axis.title = element_text(size = 15)) +
  theme(legend.position = "bottom",
        legend.text = element_text(color = "black", size = 8)) +
  theme(strip.text = element_text(face="bold", size=10),
        strip.background = element_rect(linewidth=2)) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        linewidth = 0.25, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white")) +
  xlab('')+ #Removed Axis label 
  ylab(bquote(''~R^2)) +
  labs(color=NULL) + 
  coord_flip()+ 
  theme(plot.title = element_textbox(size =15, hjust = 0.5, margin = margin(t = 5, b = 5)))
ggEnvfitBPChart

#RDA of BP Function

#Functional metadata set to include CST 
fun_metadataRDA <- MM_Selected_Metadata_RDA_030224
fun_metadataRDA <- fun_metadataRDA[, 2:24] #Remove Sample ID, keep CST

specBP <- BP_GO_Data
specBP <- specBP[, -1] #Remove non-species data columns (Sample ID)

#Explore species dataframe
names(specBP)
dim(specBP)  # dataset dimensions
head(specBP)  # look at first 5 rows
str(specBP)  # structure of objects in dataset
summary(specBP)  # summary statistics for all objects (min, mean, max, etc.)

# Count number of species frequencies in each abundance class
abBP <- table(unlist(specBP))
# Plot distribution of species frequencies
barplot(abBP, las = 1, # make axis labels perpendicular to axis
        xlab = "Abundance class", ylab = "Frequency", # label axes
        col = grey(5:0/5)) # 5-colour gradient for the bars
# Count the number of zeros in the dataset
sum(specBP == 0)
# Calculate proportion of zeros in the dataset = 68.0%
sum(specBP == 0)/(nrow(specBP) * ncol(specBP))
# Apply Hellinger transformation to correct for the double zero problem (The Hellinger transformation expresses abundances as the square-root of their relative abundance at each site (Borcard, Gillet, and Legendre 2011), solving the issue with double zeros)
specBP.hel <- decostand(specBP, method = "hellinger")

#Explore environmental dataset 
names(fun_metadataRDA)
dim(fun_metadataRDA)
head(fun_metadataRDA)
str(fun_metadataRDA)
summary(fun_metadataRDA)

# good idea to check for correlations between variables, as the constrained ordination methods we will be using are highly sensitive to collinearities in the explanatory matrix.

# We can visually look for colinearity between variables - HEAT MAP 
heatmap(abs(cor(fun_metadataRDA)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topleft", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
#From this you can observe colinearity  

#Standardise environmental variable 
# Scale and center variables
fun_metadataRDA.z <- decostand(fun_metadataRDA, method = "standardize")
# Variables are now centered around a mean of 0
round(apply(fun_metadataRDA.z, 2, mean), 1)
# and scaled to have a standard deviation of 1
apply(fun_metadataRDA.z, 2, sd)
#Now we have centred standardised our  variables. 

# Model the effect of all environmental variables on bp community composition
bp.rda <- rda(specBP.hel ~ ., data = fun_metadataRDA.z)
sqrt(vif.cca(bp.rda)) #Examine degree of colinearity - aim <2 for each 
summary(bp.rda)

# global R2 of the Model 
global_r2bp <- RsquareAdj(bp.rda)$adj.r.squared
global_r2bp #r2 0.00347

###Examining just covariates identified in envfit - CST, BF, VD, Smoking, Vitamin K 
bp_rda <- rda(specBP.hel ~ fun_metadataRDA$CST + fun_metadataRDA$Prev_BF + fun_metadataRDA$Prev_VD + fun_metadataRDA$Smoking + fun_metadataRDA$`Vitamin K1 (Phylloquinone)`)
summary(bp_rda)
plot(bp_rda)
RsquareAdj(bp_rda) ###-0.002
anova.cca(bp_rda) #p 0.491 - do not perform forward selection
anova.cca(bp_rda, by = "axis") # Nil significant
anova.cca(bp_rda, by = "terms") # Nil significant 
sqrt(vif.cca(bp_rda)) #No colinearity all <2 

#### Creating bp RDA Plot 
# To plot data via ggplot2, information on the locations of sites (circles in the ordiplot) is obtained via function sites.long. Information on species is extracted by function species.long.
bp_rda_plot1 <- ordiplot(bp_rda, scaling = 1, type = "text")
sites.longBP <- sites.long(bp_rda_plot1, env.data=fun_metadataRDA)
head(sites.longBP)

species.longBP <- species.long(bp_rda_plot1)
species.longBP

#Information on the labeling of the axes is obtained with function axis.long. This information is extracted from the ordination model and not the ordiplot, hence it is important to select the same axes via argument ‘choices’. The extracted information includes information on the percentage of explained variation. I suggest that you cross-check with the summary of the redundancy analysis, where information on proportion explained is given.
Ordination.modelBP <- rda(specBP.hel ~ ., data = fun_metadataRDA, scaling = "species")
axis.longBP <- axis.long(Ordination.modelBP, choices=c(1, 2))
axis.longBP ### Gives x axis RDA 1 11% and y axis RDA 2 3.3% 

library(BiodiversityR)
library(ggThemeAssist) #Find code for BioR Theme at: Page 7 of https://cran.r-project.org/web/packages/BiodiversityR/BiodiversityR.pdf
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 10, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.key = element_blank())

## bp RDA Plot - many black spots - can this be fixed by using different colour dataframe/CST for geom_point? 
plotRDA_BP <- ggplot() + 
  xlab(axis.longBP[1, "label"]) +
  ylab(axis.longBP[2, "label"]) +  
  labs(colour = "CST")  +
  annotate("text", x=0.11, y=-0.7, label= "p = 0.807") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  geom_point(data=sites.longBP, 
             aes(x=axis1, y=axis2, colour=df2$CST), 
             size=5) +
  geom_point(data=species.longBP,color = "#3300CC", 
             aes(x=axis1, y=axis2), size=5)+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") + 
  geom_segment(data=species.longBP, aes(x=0, xend=axis1, y=0, yend=axis2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) 


plotRDA_BP


#Summary of bp plots 
bp_ordination
ggBPenvfit 
ggEnvfitBPChart
plotRDA_BP

#Assimilating all species ordination, envfit and RDA plots 
Figure3a <- ggarrange(bp_ordination, ggBPenvfit, plotRDA_BP, ggEnvfitBPChart, 
                     labels = c("A", "B", "C", "D"),
                     common.legend = TRUE, legend = "left",
                     ncol = 4, nrow = 1)
Figure3a


# Cellular Component Functional Analysis ----------------------------------

#Create cc functional dataframe

CC_GO_Data <- read_excel("Desktop/Documents/PhD/Low Risk Vaginal Microbial Profile - MicroMOM Analysis/MM Datasets/MM_Early_Vaginal_Jan24_Analysis/CC_GO_Data.xlsx")
View(CC_GO_Data)
View(CC_GO_Data)
cc = CC_GO_Data 
cc = cc[,2:224] #Remove Sample ID from BP data

#Create bp metadata dataframe
cc_metadata1 <- MM_Selected_Metadata_030224 #All 22 covariates with Sample ID and CST
cc_metadata = cc_metadata1[, 2:24] #Remove Sample ID, Keep CST

#Perform bp nmds ordination 
#Remove any columns where all samples have zero of that function
cc <- cc %>% 
  select_if(negate(function(col) is.numeric(col) && sum(col) < 0.000001))

m_cc = as.matrix(cc) #Create bp Matrix to run NMDS
set.seed(123)
nmds_cc = metaMDS(m_cc, distance = "bray", autotransform = FALSE)
nmds_cc
plot(nmds_cc) #Record stress no. 0.05789

data.scoresCC = as.data.frame(scores(nmds_cc)$sites) #Creates object of BP-nmds coordination points 
data.scoresCC$CST = cc_metadata1$CST #Adds CST labels 

#If you want the legend ordered
data.scoresCC$CST <- factor(data.scoresCC$CST, levels = c("CST I", "CST II", 
                                                          "CST III", "CST IV", "CST V", "CST VIII", 
                                                          "CST IX"))
View(data.scoresCC) #Ensure no NA Values 

#Stats for NMDS Biological Processes vs Metadata ordination 
anoCC = anosim(m_cc, cc_metadata1$CST, distance = "bray", permutations = 9999)
anoCC
#Plot for NMDS Biological Processes vs Metadata ordination 
CST.anoCC <- with(cc_metadata, anosim(m_cc, cc_metadata1$CST, distance = "bray", permutations = 9999, strata = NULL, parallel = getOption("mc.cores")))
summary(CST.anoCC) #ANOSIM 0.4655 p <0.0001
plot(CST.anoCC) #Shows distance ranking per CST Group 

#Creating plot of envfit for large meta-data set against BP nmds plot 
cc_ordination = ggplot(data.scoresCC, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = data.scoresCC$CST)) + 
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  annotate("text", size = 3.0, x=-1.0, y=0.2, label= " Stress = 0.058") +
  annotate("text", size = 3.0, x=-1.0, y=0.1, label= " ANOSIM R = 0.4666
           p <0.0001 ") + ##Add NMDS Stress ANOSIM Stat R 0.9967, p <0.0001, stress no. 0.06722932
  labs(x = "NMDS1", colour = "CST", y = "NMDS2")  +
  ggtitle( "Cell Cycle Functions", size = 12) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  stat_ellipse(type = "t", level = 0.95, aes(color = data.scoresBP$CST))
cc_ordination


#Envfit of BP Functional Data against metadata covariates  
nmds_cc
envfit_cc <- envfit(nmds_cc, cc_metadata, perm=1000, na.rm=TRUE)
envfit_cc   
str(envfit_cc)
plot(envfit_cc)

#to annotate envfit_BP plot
envfit_CC_coord_cont = as.data.frame(scores(envfit_cc, "vectors")) * ordiArrowMul(envfit_cc)
envfit_CC_coord_cat = as.data.frame(scores(envfit_cc, "factors")) * ordiArrowMul(envfit_cc)

significant_cc_vectors <- names(envfit_cc$vectors$pvals)[envfit_cc$vectors$pvals < 0.05]

ggCCEnvfit <- ggplot(data = data.scoresCC, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = CST), size = 3, alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = as.data.frame(scores(envfit_cc, "vectors")),
               size = 2, alpha = 0.5, colour = "grey30") +
  geom_point(data = as.data.frame(scores(envfit_cc, "vectors")),
             aes(x = NMDS1, y = NMDS2), shape = "diamond", size = 2, alpha = 0.6,
             colour = "grey30") +
  geom_text(data = as.data.frame(scores(envfit_cc, "vectors"))[significant_cc_vectors, ],
            aes(x = NMDS1, y = NMDS2 + 0.04), size = 3.5,
            label = significant_cc_vectors,
            colour = "grey30", fontface = "bold") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 10),
        legend.text = element_text(size = 12, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  annotate("text", size = 3.0, x=-1.1, y=0.2, label= " CST (Adj-R 0.212)
           Smoking (Adj-R 0.053) 
           Prev VD (Adj-R 0.044)
           Vit K1 (Adj-R 0.031)
           Phenylalanine (Adj-R 0.056)") +
  labs(x = "NMDS1", colour = "CST", y = "NMDS2")
ggCCEnvfit

#to Create horizontal barchart

envfit_cc  #Generated above 
#Make tibble of R and p values for Vectors then Factors 
envfit_cc_Vectors <- tibble((envfit_cc$vectors)$r, (envfit_cc$vectors)$pvals)
envfit_cc_Factors <- tibble((envfit_cc$factors)$r, (envfit_cc$factors)$pvals)

#Rename columns as R-Squared and p value
names(envfit_cc_Vectors)[1] <- "R2"
names(envfit_cc_Vectors)[2] <- "pvalue"
View(envfit_cc_Vectors)

#Repeat for factors 
names(envfit_cc_Factors)[1] <- "R2"
names(envfit_cc_Factors)[2] <- "pvalue"
View(envfit_cc_Factors)

#Create this table as an object of all R2 and p Values for covariates in species envfit model 
EnvFitCCMetadata <- (rbind(envfit_cc_Vectors, envfit_cc_Factors)) #Covariate Row names lost before this 
View(EnvFitCCMetadata)

#compare both envfit results and your R2/pval dataframe to ensure correct
envfit_cc
View(EnvFitCCMetadata)

#get list from running enShort to display results 
covariate_cc_categories <- data.frame("Covariate" =c('Maternal Age',
                                                     'Education', 
                                                     'HP Index', 
                                                     'BMI', 
                                                     'MET per week', 
                                                     'WHO Wellbeing Score', 
                                                     'Fat',
                                                     'Starch',
                                                     'Sugars',
                                                     'NSP',
                                                     'Vitamin K1',
                                                     'Niacin',
                                                     'Selenium',
                                                     'Leucine',
                                                     'Lysine',
                                                     'Methionine',
                                                     'Phenylalanine',
                                                     'Valine', 
                                                     'CST',
                                                     'Ethnicity', 
                                                     'Prev BF', 
                                                     'Prev VD', 
                                                     'Smoking'),
                                      "Covariate type" = c('Demography',
                                                           'Demography',
                                                           'Demography',
                                                           'Demography',
                                                           'Lifestyle',
                                                           'Lifestyle',
                                                           'Macronutrient',
                                                           'Macronutrient',
                                                           'Macronutrient',
                                                           'Macronutrient',
                                                           'Micronutrient',
                                                           'Micronutrient',
                                                           'Micronutrient',
                                                           'Amino Acid',
                                                           'Amino Acid',
                                                           'Amino Acid',
                                                           'Amino Acid',
                                                           'Amino Acid',
                                                           'Species',
                                                           'Demography',
                                                           'Obstetric', 
                                                           'Obstetric', 
                                                           'Lifestyle'))

envfit_cc_results <- cbind(covariate_cc_categories,EnvFitCCMetadata) 
View(envfit_cc_results)
#Confirm these are correct corresponding covariate names and r/p by looking at enShort result 
envfit_cc
#Relabel R2 as Vaginal beta diversity for labelling 
colnames(envfit_cc_results)<- c('Covariate',
                                'Covariate_type',
                                'Cell Cycle', 
                                'p-value')
names(envfit_cc_results)
envfitR2_longCC <- gather(envfit_cc_results, Sample_type, R2, 'Cell Cycle', factor_key=TRUE)
envfitR2_longCC$Covariate_type = factor(envfitR2_longCC$Covariate_type, levels=c('Demography', 'Lifestyle', 'Obstetric', 'Macronutrient','Micronutrient', 'Amino Acid', 'Species'))
envfitR2_longCC$Sample_type = factor(envfitR2_longCC$Sample_type, levels=c('Cell Cycle'))
head(envfitR2_longCC)

#Horizontal Barchart ggplot for BP Function and metadata 
ggEnvfitCCChart <- ggplot(envfitR2_longCC, aes(x = fct_reorder(Covariate, +R2, max), y = R2, fill = Covariate_type, label = ifelse(envfitR2_longCC$`p-value` < 0.05, "*", "")))+
  facet_wrap(~Sample_type, nrow = 2) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0) +
  scale_color_brewer(palette = "Set2")+
  theme (axis.text.x = element_text(size = 7, colour = "black", angle = 45)) +
  theme (axis.text.y = element_text(size = 7, colour = "black"))+
  theme(axis.title = element_text(size = 15)) +
  theme(legend.position = "bottom",
        legend.text = element_text(color = "black", size = 8)) +
  theme(strip.text = element_text(face="bold", size=10),
        strip.background = element_rect(linewidth=2)) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        linewidth = 0.25, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white")) +
  xlab('')+ #Removed Axis label 
  ylab(bquote(''~R^2)) +
  labs(color=NULL) + 
  coord_flip()+ 
  theme(plot.title = element_textbox(size =15, hjust = 0.5, margin = margin(t = 5, b = 5)))
ggEnvfitCCChart

#RDA of BP Function

#Functional metadata set to include CST 
fun_metadataRDA <- MM_Selected_Metadata_RDA_030224
fun_metadataRDA <- fun_metadataRDA[, 2:24] #Remove Sample ID, keep CST

specCC <- CC_GO_Data
specCC <- specCC[, -1] #Remove non-species data columns (Sample ID)

#Explore species dataframe
names(specCC)
dim(specCC)  # dataset dimensions
head(specCC)  # look at first 5 rows
str(specCC)  # structure of objects in dataset
summary(specCC)  # summary statistics for all objects (min, mean, max, etc.)

# Count number of species frequencies in each abundance class
abCC <- table(unlist(specCC))
# Plot distribution of species frequencies
barplot(abBP, las = 1, # make axis labels perpendicular to axis
        xlab = "Abundance class", ylab = "Frequency", # label axes
        col = grey(5:0/5)) # 5-colour gradient for the bars
# Count the number of zeros in the dataset
sum(specCC == 0)
# Calculate proportion of zeros in the dataset = 68.0%
sum(specCC == 0)/(nrow(specCC) * ncol(specCC))
# Apply Hellinger transformation to correct for the double zero problem (The Hellinger transformation expresses abundances as the square-root of their relative abundance at each site (Borcard, Gillet, and Legendre 2011), solving the issue with double zeros)
specCC.hel <- decostand(specCC, method = "hellinger")

#Explore environmental dataset 
names(fun_metadataRDA)
dim(fun_metadataRDA)
head(fun_metadataRDA)
str(fun_metadataRDA)
summary(fun_metadataRDA)

# good idea to check for correlations between variables, as the constrained ordination methods we will be using are highly sensitive to collinearities in the explanatory matrix.

# We can visually look for colinearity between variables - HEAT MAP 
heatmap(abs(cor(fun_metadataRDA)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topleft", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
#From this you can observe colinearity  

#Standardise environmental variable 
# Scale and center variables
fun_metadataRDA.z <- decostand(fun_metadataRDA, method = "standardize")
# Variables are now centered around a mean of 0
round(apply(fun_metadataRDA.z, 2, mean), 1)
# and scaled to have a standard deviation of 1
apply(fun_metadataRDA.z, 2, sd)
#Now we have centred standardised our  variables. 

# Model the effect of all environmental variables on bp community composition
cc.rda <- rda(specCC.hel ~ ., data = fun_metadataRDA.z)
sqrt(vif.cca(cc.rda)) #Examine degree of colinearity - aim <2 for each 
summary(cc.rda)

# global R2 of the Model 
global_r2cc <- RsquareAdj(cc.rda)$adj.r.squared
global_r2cc #r2 -0.0277

###Examining just covariates identified in envfit - CST, BF, VD, Smoking, Vitamin K 
cc_rda <- rda(specCC.hel ~ fun_metadataRDA$CST + fun_metadataRDA$Prev_VD + fun_metadataRDA$Smoking + fun_metadataRDA$`Vitamin K1 (Phylloquinone)` + fun_metadataRDA$`Niacin (preformed)` + fun_metadataRDA$Phenylalanine)
summary(cc_rda)
plot(cc_rda)
RsquareAdj(cc_rda) ### -0.000962
anova.cca(cc_rda) #p 0.456 - do not perform forward selection
anova.cca(cc_rda, by = "axis") # Nil significant
anova.cca(cc_rda, by = "terms") # Nil significant 
sqrt(vif.cca(cc_rda)) #No colinearity all <2 

#### Creating bp RDA Plot 
# To plot data via ggplot2, information on the locations of sites (circles in the ordiplot) is obtained via function sites.long. Information on species is extracted by function species.long.
cc_rda_plot1 <- ordiplot(cc_rda, scaling = 1, type = "text")
sites.longCC <- sites.long(cc_rda_plot1, env.data=fun_metadataRDA)
head(sites.longCC)

species.longCC <- species.long(cc_rda_plot1)
species.longCC

#Information on the labeling of the axes is obtained with function axis.long. This information is extracted from the ordination model and not the ordiplot, hence it is important to select the same axes via argument ‘choices’. The extracted information includes information on the percentage of explained variation. I suggest that you cross-check with the summary of the redundancy analysis, where information on proportion explained is given.
Ordination.modelCC <- rda(specCC.hel ~ ., data = fun_metadataRDA, scaling = "species")
axis.longCC <- axis.long(Ordination.modelCC, choices=c(1, 2))
axis.longCC ### Gives x axis RDA 1 11% and y axis RDA 2 3.3% 

library(BiodiversityR)
library(ggThemeAssist) #Find code for BioR Theme at: Page 7 of https://cran.r-project.org/web/packages/BiodiversityR/BiodiversityR.pdf
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

## bp RDA Plot - many black spots - can this be fixed by using different colour dataframe/CST for geom_point? 
plotRDA_CC <- ggplot() + 
  xlab(axis.longCC[1, "label"]) +
  ylab(axis.longCC[2, "label"]) +  
  annotate("text", x=-0.3, y=0.2, label= "p = 0.456") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  geom_point(data=sites.longCC, 
             aes(x=axis1, y=axis2, colour=df2$CST), 
             size=5) +
  geom_point(data=species.longCC,color = "#3300CC", 
             aes(x=axis1, y=axis2), size=5)+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") + 
  geom_segment(data=species.longCC, aes(x=0, xend=axis1, y=0, yend=axis2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) 
plotRDA_CC

#Summary of cc plots 
cc_ordination
ggCCenvfit 
ggEnvfitCCChart
plotRDA_CC

#Assimilating all cc ordination, envfit and RDA plots 
Figure3b <- ggarrange(cc_ordination, ggCCenvfit, plotRDA_CC, ggEnvfitCCChart, 
                      labels = c("E", "F", "G", "H"),
                      common.legend = TRUE, legend = "left",
                      ncol = 4, nrow = 1)
Figure3b


# Metabolic Function Analysis ---------------------------------------------
#Create mf functional dataframe
MF_GO_Data <- read_excel("Desktop/Documents/PhD/Low Risk Vaginal Microbial Profile - MicroMOM Analysis/MM Datasets/MM_Early_Vaginal_Jan24_Analysis/MF_GO_Data.xlsx")
View(MF_GO_Data)
mf = MF_GO_Data 
mf = mf[,2:1959] #Remove Sample ID from BP data

#Create bp metadata dataframe
mf_metadata1 <- MM_Selected_Metadata_030224 #All 22 covariates with Sample ID and CST
mf_metadata = mf_metadata1[, 2:24] #Remove Sample ID, Keep CST

#Perform bp nmds ordination 
#Remove any columns where all samples have zero of that function
mf <- mf %>% 
  select_if(negate(function(col) is.numeric(col) && sum(col) < 0.000001))

m_mf = as.matrix(mf) #Create bp Matrix to run NMDS
set.seed(123)
nmds_mf = metaMDS(m_mf, distance = "bray", autotransform = FALSE)
nmds_mf
plot(nmds_mf) #Record stress no. 0.0979

data.scoresMF = as.data.frame(scores(nmds_mf)$sites) #Creates object of BP-nmds coordination points 
data.scoresMF$CST = mf_metadata1$CST #Adds CST labels 

#If you want the legend ordered
data.scoresMF$CST <- factor(data.scoresMF$CST, levels = c("CST I", "CST II", 
                                                          "CST III", "CST IV", "CST V", "CST VIII", 
                                                          "CST IX"))
View(data.scoresMF) #Ensure no NA Values 

#Stats for NMDS Biological Processes vs Metadata ordination 
anoMF = anosim(m_mf, mf_metadata1$CST, distance = "bray", permutations = 9999)
anoMF
#Plot for NMDS Biological Processes vs Metadata ordination 
CST.anoMF <- with(mf_metadata, anosim(m_mf, mf_metadata1$CST, distance = "bray", permutations = 9999, strata = NULL, parallel = getOption("mc.cores")))
summary(CST.anoMF) #ANOSIM 0.8601 p <0.0001
plot(CST.anoMF) #Shows distance ranking per CST Group 

#Creating plot of envfit for large meta-data set against BP nmds plot 
mf_ordination = ggplot(data.scoresMF, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(colour = data.scoresMF$CST)) + 
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  annotate("text", size = 3.0, x=1.1, y=0.35, label= " Stress = 0.0979") +
  annotate("text", size = 3.0, x=1.1, y=0.25, label= " ANOSIM R = 0.8601
           p <0.0001 ") +
  ggtitle( "Metabolic Functions", size = 12) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x = "NMDS1", colour = "CST", y = "NMDS2")  +
  stat_ellipse(type = "t", level = 0.95, aes(color = data.scoresMF$CST))
mf_ordination


#Envfit of BP Functional Data against metadata covariates  
nmds_mf
envfit_mf <- envfit(nmds_mf, mf_metadata, perm=1000, na.rm=TRUE)
envfit_mf
str(envfit_mf)
plot(envfit_mf)

#to annotate envfit_BP plot
envfit_MF_coord_cont = as.data.frame(scores(envfit_mf, "vectors")) * ordiArrowMul(envfit_mf)
envfit_MF_coord_cat = as.data.frame(scores(envfit_mf, "factors")) * ordiArrowMul(envfit_mf)

significant_mf_vectors <- names(envfit_mf$vectors$pvals)[envfit_mf$vectors$pvals < 0.05]

ggMFEnvfit <- ggplot(data = data.scoresMF, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = CST), size = 3, alpha = 0.5) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = as.data.frame(scores(envfit_mf, "vectors")),
               size = 2, alpha = 0.5, colour = "grey30") +
  geom_point(data = as.data.frame(scores(envfit_mf, "vectors")),
             aes(x = NMDS1, y = NMDS2), shape = "diamond", size = 2, alpha = 0.6,
             colour = "grey30") +
  geom_text(data = as.data.frame(scores(envfit_mf, "vectors"))[significant_mf_vectors, ],
            aes(x = NMDS1, y = NMDS2 + 0.04), size = 3.5,
            label = significant_mf_vectors,
            colour = "grey30", fontface = "bold") +
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
        axis.text.x = element_text(colour = "black", face = "bold", size = 10),
        legend.text = element_text(size = 10, face ="bold", colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10),
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) +
  annotate("text", size = 3.0, x=1.0, y=0.2, label= " CST (Adj-R 0.483)
           Smoking (Adj-R 0.011) 
           Prev VD (Adj-R 0.047)
           Prev BF (Adj-R 0.031)
           Vit K1 (Adj-R 0.050)") +
  labs(x = "NMDS1", colour = "CST", y = "NMDS2")
ggMFEnvfit


#to Create horizontal barchart

envfit_mf  #Generated above 


#Make tibble of R and p values for Vectors then Factors 
envfit_mf_Vectors <- tibble((envfit_mf$vectors)$r, (envfit_mf$vectors)$pvals)
envfit_mf_Factors <- tibble((envfit_mf$factors)$r, (envfit_mf$factors)$pvals)

##How to keep covariate names for these - answer: matched below 

#Rename columns as R-Squared and p value
names(envfit_mf_Vectors)[1] <- "R2"
names(envfit_mf_Vectors)[2] <- "pvalue"
View(envfit_mf_Vectors)

#Repeat for factors 
names(envfit_mf_Factors)[1] <- "R2"
names(envfit_mf_Factors)[2] <- "pvalue"
View(envfit_mf_Factors)

#Create this table as an object of all R2 and p Values for covariates in species envfit model 
EnvFitMFMetadata <- (rbind(envfit_mf_Vectors, envfit_mf_Factors)) #Covariate Row names lost before this 
View(EnvFitMFMetadata)

#compare both envfit results and your R2/pval dataframe to ensure correct
envfit_mf
View(EnvFitMFMetadata)

#get list from running enShort to display results 
covariate_mf_categories <- data.frame("Covariate" =c('Maternal Age',
                                                     'Education', 
                                                     'HP Index', 
                                                     'BMI', 
                                                     'MET per week', 
                                                     'WHO Wellbeing Score', 
                                                     'Fat',
                                                     'Starch',
                                                     'Sugars',
                                                     'NSP',
                                                     'Vitamin K1',
                                                     'Niacin',
                                                     'Selenium',
                                                     'Leucine',
                                                     'Lysine',
                                                     'Methionine',
                                                     'Phenylalanine',
                                                     'Valine', 
                                                     'CST',
                                                     'Ethnicity', 
                                                     'Prev BF', 
                                                     'Prev VD', 
                                                     'Smoking'),
                                      "Covariate type" = c('Demography',
                                                           'Demography',
                                                           'Demography',
                                                           'Demography',
                                                           'Lifestyle',
                                                           'Lifestyle',
                                                           'Macronutrient',
                                                           'Macronutrient',
                                                           'Macronutrient',
                                                           'Macronutrient',
                                                           'Micronutrient',
                                                           'Micronutrient',
                                                           'Micronutrient',
                                                           'Amino Acid',
                                                           'Amino Acid',
                                                           'Amino Acid',
                                                           'Amino Acid',
                                                           'Amino Acid',
                                                           'Species',
                                                           'Demography',
                                                           'Obstetric', 
                                                           'Obstetric', 
                                                           'Lifestyle'))

envfit_mf_results <- cbind(covariate_mf_categories,EnvFitMFMetadata) 
View(envfit_mf_results)
#Confirm these are correct corresponding covariate names and r/p by looking at enShort result 
envfit_mf
#Relabel R2 as Vaginal beta diversity for labelling 
colnames(envfit_mf_results)<- c('Covariate',
                                'Covariate_type',
                                'Metabolic Function', 
                                'p-value')
names(envfit_mf_results)
envfitR2_longMF <- gather(envfit_mf_results, Sample_type, R2, 'Metabolic Function', factor_key=TRUE)
envfitR2_longMF$Covariate_type = factor(envfitR2_longMF$Covariate_type, levels=c('Demography', 'Lifestyle', 'Obstetric', 'Macronutrient','Micronutrient', 'Amino Acid', 'Species'))
envfitR2_longMF$Sample_type = factor(envfitR2_longMF$Sample_type, levels=c('Metabolic Function'))
head(envfitR2_longMF)

#Horizontal Barchart ggplot for BP Function and metadata 
ggEnvfitMFChart <- ggplot(envfitR2_longMF, aes(x = fct_reorder(Covariate, +R2, max), y = R2, fill = Covariate_type, label = ifelse(envfitR2_longMF$`p-value` < 0.05, "*", "")))+
  facet_wrap(~Sample_type, nrow = 2) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0) +
  scale_color_brewer(palette = "Set2")+
  theme (axis.text.x = element_text(size = 7, colour = "black", angle = 45)) +
  theme (axis.text.y = element_text(size = 7, colour = "black"))+
  theme(axis.title = element_text(size = 15)) +
  theme(legend.position = "bottom",
        legend.text = element_text(color = "black", size = 8)) +
  theme(strip.text = element_text(face="bold", size=10),
        strip.background = element_rect(linewidth=2)) +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        linewidth = 0.25, linetype = "solid"),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"), 
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white")) +
  xlab('')+ #Removed Axis label 
  ylab(bquote(''~R^2)) +
  labs(color=NULL) + 
  coord_flip()+ 
  theme(plot.title = element_textbox(size =15, hjust = 0.5, margin = margin(t = 5, b = 5)))
ggEnvfitMFChart

#RDA of BP Function

#Functional metadata set to include CST 
fun_metadataRDA <- MM_Selected_Metadata_RDA_030224
fun_metadataRDA <- fun_metadataRDA[, 2:24] #Remove Sample ID, keep CST

specMF <- MF_GO_Data
specMF <- specMF[, -1] #Remove non-species data columns (Sample ID)

#Explore species dataframe
names(specMF)
dim(specMF)  # dataset dimensions
head(specMF)  # look at first 5 rows
str(specMF)  # structure of objects in dataset
summary(specMF)  # summary statistics for all objects (min, mean, max, etc.)

# Count number of species frequencies in each abundance class
abMF <- table(unlist(specMF))
# Plot distribution of species frequencies
barplot(abMF, las = 1, # make axis labels perpendicular to axis
        xlab = "Abundance class", ylab = "Frequency", # label axes
        col = grey(5:0/5)) # 5-colour gradient for the bars
# Count the number of zeros in the dataset
sum(specMF == 0)
# Calculate proportion of zeros in the dataset = 71.2%
sum(specMF == 0)/(nrow(specMF) * ncol(specMF))
# Apply Hellinger transformation to correct for the double zero problem (The Hellinger transformation expresses abundances as the square-root of their relative abundance at each site (Borcard, Gillet, and Legendre 2011), solving the issue with double zeros)
specMF.hel <- decostand(specMF, method = "hellinger")

#Explore environmental dataset 
names(fun_metadataRDA)
dim(fun_metadataRDA)
head(fun_metadataRDA)
str(fun_metadataRDA)
summary(fun_metadataRDA)

# good idea to check for correlations between variables, as the constrained ordination methods we will be using are highly sensitive to collinearities in the explanatory matrix.

# We can visually look for colinearity between variables - HEAT MAP 
heatmap(abs(cor(fun_metadataRDA)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topleft", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
#From this you can observe colinearity  

#Standardise environmental variable 
# Scale and center variables
fun_metadataRDA.z <- decostand(fun_metadataRDA, method = "standardize")
# Variables are now centered around a mean of 0
round(apply(fun_metadataRDA.z, 2, mean), 1)
# and scaled to have a standard deviation of 1
apply(fun_metadataRDA.z, 2, sd)
#Now we have centred standardised our  variables. 

# Model the effect of all environmental variables on bp community composition
mf.rda <- rda(specMF.hel ~ ., data = fun_metadataRDA.z)
sqrt(vif.cca(mf.rda)) #Examine degree of colinearity - aim <2 for each 
summary(mf.rda)

# global R2 of the Model 
global_r2mf <- RsquareAdj(mf.rda)$adj.r.squared
global_r2mf #r2 -0.00235

###Examining just covariates identified in envfit - CST, BF, VD, Smoking, Vitamin K 
mf_rda <- rda(specMF.hel ~ fun_metadataRDA$CST + fun_metadataRDA$Prev_BF + fun_metadataRDA$Prev_VD + fun_metadataRDA$Smoking)
plot(mf_rda)
RsquareAdj(mf_rda) ### -0.0149
anova.cca(mf_rda) #p 0.911 - do not perform forward selection
anova.cca(mf_rda, by = "axis") # Nil significant
anova.cca(mf_rda, by = "terms") # Nil significant 
sqrt(vif.cca(mf_rda)) #No colinearity all <2 

#### Creating bp RDA Plot 
# To plot data via ggplot2, information on the locations of sites (circles in the ordiplot) is obtained via function sites.long. Information on species is extracted by function species.long.
mf_rda_plot1 <- ordiplot(mf_rda, scaling = 1, type = "text")
sites.longMF <- sites.long(mf_rda_plot1, env.data=fun_metadataRDA)
head(sites.longMF)

species.longMF <- species.long(mf_rda_plot1)
species.longMF

#Information on the labeling of the axes is obtained with function axis.long. This information is extracted from the ordination model and not the ordiplot, hence it is important to select the same axes via argument ‘choices’. The extracted information includes information on the percentage of explained variation. I suggest that you cross-check with the summary of the redundancy analysis, where information on proportion explained is given.
Ordination.modelMF <- rda(specMF.hel ~ ., data = fun_metadataRDA, scaling = "species")
axis.longMF <- axis.long(Ordination.modelMF, choices=c(1, 2))
axis.longMF ### Gives x axis RDA 1 10.3% and y axis RDA 2 3.5% 

library(BiodiversityR)
library(ggThemeAssist) #Find code for BioR Theme at: Page 7 of https://cran.r-project.org/web/packages/BiodiversityR/BiodiversityR.pdf
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12, family="Arial"),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

## MF RDA Plot - many black spots - can this be fixed by using different colour dataframe/CST for geom_point? 
plotRDA_MF <- ggplot() + 
  xlab(axis.longMF[1, "label"]) +
  ylab(axis.longMF[2, "label"]) +  
  labs(colour = "CST")  +
  annotate("text", x=0.11, y=-1.0, label= " p = 0.911") +
  theme(axis.text.y = element_text(colour = "black", size = 10, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 10), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 10), 
        axis.title.x = element_text(face = "bold", size = 10, colour = "black"), 
        legend.title = element_text(size = 10, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key = element_blank()) + 
  geom_point(data=sites.longBP, 
             aes(x=axis1, y=axis2, colour=df2$CST), 
             size=5) +
  geom_point(data=species.longBP,color = "#3300CC", 
             aes(x=axis1, y=axis2), size=5)+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") + 
  geom_segment(data=species.longBP, aes(x=0, xend=axis1, y=0, yend=axis2), 
               color="black", arrow=arrow(length=unit(0.01,"npc"))) 
plotRDA_MF

#Summary of mf plots 
mf_ordination
ggMFenvfit 
ggEnvfitMFChart
plotRDA_MF

#Assimilating all species ordination, envfit and RDA plots 
Figure3c <- ggarrange(mf_ordination, ggMFenvfit, plotRDA_MF, ggEnvfitMFChart, 
                      labels = c("I", "J", "K", "L"),
                      common.legend = TRUE, legend = "left",
                      ncol = 4, nrow = 1)
Figure3c


# Figure 3 - Combined Functional Analysis Plots ----------------------------------------------------------

#Summary of bp plots 
bp_ordination
ggBPenvfit 
ggEnvfitBPChart
plotRDA_BP

#Summary of cc plots 
cc_ordination
ggCCenvfit 
ggEnvfitCCChart
plotRDA_CC

#Summary of mf plots 
mf_ordination
ggMFenvfit 
ggEnvfitMFChart
plotRDA_MF


#Trial as a 3 x 4 Table instead 
#Figure 3i nmds ordination 
Figure3i <- ggarrange(bp_ordination, cc_ordination,mf_ordination, 
                      labels = c("A", "B", "C"),
                      common.legend = TRUE, legend = "left", 
                      ncol = 3, nrow = 1)
Figure3i
#envfit plots 
Figure3ii <- ggarrange(ggBPEnvfit, ggCCEnvfit, ggMFEnvfit,
                       labels = c("D", "E", "F"),
                       legend = "none", 
                       ncol = 3, nrow = 1)
Figure3ii

#RDA Plots

Figure3iii <- ggarrange(plotRDA_BP, plotRDA_CC, plotRDA_MF, 
                        labels = c("G", "H", "I"),
                        legend = "none", 
                        ncol = 3, nrow = 1)
Figure3iii

#envfit barcharts 
Figure3iv <- ggarrange(ggEnvfitBPChart, ggEnvfitCCChart, ggEnvfitMFChart,
                       labels = c("J", "K", "L"),
                        common.legend = TRUE, legend = "bottom", 
                        ncol = 3, nrow = 1)
Figure3iv

#Vertical plot
Figure3V2 <- ggarrange(Figure3i, Figure3ii, Figure3iii, Figure3iv,
                       common.legend = FALSE, legend = "bottom", 
                       ncol = 1, nrow = 4)
Figure3V2


# Macronutrient Diversity Regression  -------------------------------------

MM_Metadata_Macronutrient_Analysis_040224 <- read_excel("Desktop/Documents/PhD/Low Risk Vaginal Microbial Profile - MicroMOM Analysis/MM Datasets/MM_Early_Vaginal_Jan24_Analysis/MM_Metadata_Macronutrient_Analysis 040224.xlsx")
View(MM_Metadata_Macronutrient_Analysis_040224)

#Clean regression Metadata Dataframe
regr_Metadata <-MM_Metadata_Macronutrient_Analysis_040224
#Ensure all confounding variables are in of correct class - factor or vector
class(regr_Metadata$Caucasian)
class(regr_Metadata$Education)
class(regr_Metadata$CST)
#convert each variable using following and recheck the class - should now be factors: 
regr_Metadata$Caucasian <- as.factor(regr_Metadata$Caucasian)
regr_Metadata$Education <- as.factor(regr_Metadata$Education)
regr_Metadata$CST <- as.factor(regr_Metadata$CST)
#Deal with values in diversity of 0 - either log transform or add 0.0001
regr_Metadata$Shannon_Diversity <- regr_Metadata$Shannon_Diversity - 0.0001     #Add 0.0001 - check changes values in vaginalmetadata dataframe 
View(regr_Metadata)

#Visualise each covariate individually, then plot against Diversity, create regression object, then get summary and regression line
#MACRONUTRIENTS - MOSTLY NORMAL DISTRIBUTION 
hist(regr_Metadata$Fat)
plot(regr_Metadata$Fat, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Fibre)
plot(regr_Metadata$Fibre, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$`NSP (Englyst)`)
plot(regr_Metadata$`NSP (Englyst)`, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Carbohydrate)
plot(regr_Metadata$Carbohydrate, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Starch)
plot(regr_Metadata$Starch, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Maltose)
plot(regr_Metadata$Maltose, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$GlyLoad)
plot(regr_Metadata$GlyLoad, regr_Metadata$Shannon_Diversity)

#AMINO ACIDS - NORMAL DISTRIBUTION 
hist(regr_Metadata$Protein)
plot(regr_Metadata$Protein, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Lysine)
plot(regr_Metadata$Lysine, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Leucine)
plot(regr_Metadata$Leucine, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Valine)
plot(regr_Metadata$Valine, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Alanine)
plot(regr_Metadata$Alanine, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$Phenylalanine)
plot(regr_Metadata$Phenylalanine, regr_Metadata$Shannon_Diversity)

#MICRONUTRIENTS - NON-NORMAL LEFT SKEWED DISTRIBUTION 
hist(regr_Metadata$`Niacin (preformed)`)
plot(regr_Metadata$`Niacin (preformed)`, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$`Vitamin K1 (Phylloquinone)`)
plot(regr_Metadata$`Vitamin K1 (Phylloquinone)`, regr_Metadata$Shannon_Diversity)

hist(regr_Metadata$`Selenium (Se)`)
plot(regr_Metadata$`Selenium (Se)`, regr_Metadata$Shannon_Diversity)


#Create objects of adjusted linear regression models - adjusted for BMI, education, ethnicity, age, parity and HP Index 
#MACRONUTRIENTS

#Fat 
ModelFat <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Fat + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelFat)
plot(regr_Metadata$Fat, regr_Metadata$Shannon_Diversity, main = "Fat (Adj-R-Squared 0.015, p=0.025)", xlab = "Fat g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelFat, col=3, lwd=2)

Regr_Fat_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Fat, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="chartreuse3", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "#339933") + 
  labs(x="Fat g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Fat Intake", subtitle = "Diversity -0.114 (+002/g Fat); p=0.025
          Adj R2 0.016")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Fat_Diversity

#Fibre
ModelFibre <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Fibre + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelFibre)
plot(regr_Metadata$Fibre, regr_Metadata$Shannon_Diversity, main = "Fibre (Adj-R-Squared -0.031, p=0.923)", xlab = "Fibre g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelFibre, col=3, lwd=2)

Regr_Fibre_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Fibre, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="chartreuse3", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Fibre g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Fibre Intake", subtitle = "Diversity -0.144 (+000/g Fat); p=0.923
          Adj R2 -0.031")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Fibre_Diversity

#NSP
ModelNSP <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$`NSP (Englyst)` + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelNSP)
plot(regr_Metadata$`NSP (Englyst)`, regr_Metadata$Shannon_Diversity, main = "NSP (Adj-R-Squared -0.026, p=0.464)", xlab = "NSP g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelNSP, col=3, lwd=2)

Regr_NSP_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$`NSP (Englyst)`, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="chartreuse3", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="NSP g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary NSP Intake", subtitle = "Diversity -0.137 (-0.004/g NSP); p=0.464
          Adj R2 -0.026")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_NSP_Diversity

#Carb - not plotting as signif overlap with starch 
ModelCarb <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Carbohydrate + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelCarb)
plot(regr_Metadata$Carbohydrate, regr_Metadata$Shannon_Diversity, main = "Carbohydrate (Adj-R-Squared 0.001, p=0.064)", xlab = "Carbohydrate g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelCarb, col=3, lwd=2)

#Starch
ModelStarch <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Starch + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelStarch)
plot(regr_Metadata$Starch, regr_Metadata$Shannon_Diversity, main = "Starch (Adj-R-Squared 0.007, p=0.0425)", xlab = "Starch g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelStarch, col=3, lwd=2)

Regr_Starch_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Starch, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="chartreuse3", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "#339933") + 
  labs(x="Starch g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Starch Intake", subtitle = "Diversity -0.123 (+002/g Starch); p=0.043
          Adj R2 0.007")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Starch_Diversity

#Maltose
ModelMaltose <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Maltose + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelMaltose)
plot(regr_Metadata$Maltose, regr_Metadata$Shannon_Diversity, main = "Maltose (Adj-R-Squared 0.026, p=0.0129)", xlab = "Maltose g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelMaltose, col=3, lwd=2)

Regr_Maltose_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Maltose, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="chartreuse3", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "#339933") + 
  labs(x="Maltose g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Maltose Intake", subtitle = "Diversity -0.110 (+044/g Maltose); p=0.013
          Adj R2 0.026")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Maltose_Diversity

#Glycaemic Load 
ModelGlyLoad <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$GlyLoad + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelGlyLoad)
plot(regr_Metadata$GlyLoad, regr_Metadata$Shannon_Diversity, main = "Glycaemic Load (Adj-R-Squared -0.005, p=0.010)", xlab = "Gly Load u/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelGlyLoad, col=3, lwd=2)

Regr_GlyLoad_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$GlyLoad, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="chartreuse3", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Glycaemic Load units/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Glycaemic Load", subtitle = "Diversity -0.133 (+001/u GlyLoad); p=0.0957
          Adj R2 0.016")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_GlyLoad_Diversity

#Macronutrient Figure 
Regr_Fat_Diversity #signif
Regr_Starch_Diversity #signif
Regr_Maltose_Diversity #Signif 

Regr_Fibre_Diversity #nonsignif
Regr_NSP_Diversity #nonsignif
Regr_GlyLoad_Diversity #nonsignif

Figure4a <- ggarrange(Regr_Fat_Diversity, Regr_Starch_Diversity, Regr_Maltose_Diversity, Regr_Fibre_Diversity, Regr_NSP_Diversity, Regr_GlyLoad_Diversity, 
                      labels = c("A", "B", "C", "D", "E", "F"),
                      ncol = 3, nrow = 2)
Figure4a

# Boxplot 3xMacronutrients and CST  ---------------------------------------

#Create boxplot on r studio not spss for 

#ANOVA Results from Table 1 Supplementary Data - Compared only CST I, II and IV where n>=10 
bp_Fat_CST <- ggplot(regr_Metadata, aes(x=regr_Metadata$CST, y=regr_Metadata$Fat, fill = regr_Metadata$CST, show.legend = FALSE) ) + 
  geom_boxplot(aes(fill=regr_Metadata$CST), show.legend = FALSE) + 
  scale_color_viridis(option = "D")+
  xlab('Vaginal Community State Type')+
  ylab('Dietary Fat Intake')+
  ggtitle(label = "Fat vs Vaginal Community State Type", subtitle = "No difference on ANOVA, p = 0.283")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
bp_Fat_CST

bp_Starch_CST <- ggplot(regr_Metadata, aes(x=regr_Metadata$CST, y=regr_Metadata$Starch, fill = regr_Metadata$CST, show.legend = FALSE) ) + 
  geom_boxplot(aes(fill=regr_Metadata$CST), show.legend = FALSE) + 
  scale_color_viridis(option = "D")+
  xlab('Vaginal Community State Type')+
  ylab('Dietary Starch Intake')+
  ggtitle(label = "Starch vs Vaginal Community State Type", subtitle = "Starch highest in CST IV Group, p = 0.002")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
bp_Starch_CST

bp_Maltose_CST <- ggplot(regr_Metadata, aes(x=regr_Metadata$CST, y=regr_Metadata$Maltose, fill = regr_Metadata$CST, show.legend = FALSE) ) + 
  geom_boxplot(aes(fill=regr_Metadata$CST), show.legend = FALSE) + 
  scale_color_viridis(option = "D")+
  xlab('Vaginal Community State Type')+
  ylab('Dietary Maltose Intake')+
  ggtitle(label = "Maltose vs Vaginal Community State Type", subtitle = "No difference on ANOVA, p = 0.306")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
bp_Maltose_CST

#Boxplots of CST vs Macronutrients - Add ANOVA Results 
Figure4b <- ggarrange(bp_Fat_CST, bp_Starch_CST, bp_Maltose_CST, 
                      labels = c("G", "H", "I"),
                      ncol = 3, nrow = 1)
Figure4b 

# Secretor Status Analysis  -----------------------------------------------
#Secretor stratification of signif macronutrients by Secretor status
#First Create Secretor Pos/Neg dataframes - then run regression models of both pos/neg against nutrient 
SecretorPositive <- subset(regr_Metadata, Secretor_Status == "Positive") #n=53
SecretorNegative <- subset(regr_Metadata, Secretor_Status == "Negative") #n=18 - all same education (third level complete) and ethnicity (caucasian) so these were removed from lm model 

#Fat - Only signif relationship w Secretor Positive 
ModelFatSecretorPos <- lm(SecretorPositive$Shannon_Diversity ~ SecretorPositive$Fat + SecretorPositive$BMI + SecretorPositive$Education + SecretorPositive$Caucasian + SecretorPositive$HP_Index + SecretorPositive$Maternal_age)
summary(ModelFatSecretorPos) #Diversity -0.131 (+0.003/ g Fat, p = 0.050, Adj R=0.029)
ModelFatSecretorNeg <- lm(SecretorNegative$Shannon_Diversity ~ SecretorNegative$Fat + SecretorNegative$BMI + SecretorNegative$HP_Index + SecretorNegative$Maternal_age) #Remove Education/Caucasian as all were same value - third level, caucasian
summary(ModelFatSecretorNeg) #Diversity -0.085 (-0.000/ g Fat, p = 0.766, Adj R= -0.170)

#Starch - signif between starch&diversity only for Secretor Negative women
ModelStarchSecretorPos <- lm(SecretorPositive$Shannon_Diversity ~ SecretorPositive$Starch + SecretorPositive$BMI + SecretorPositive$Education + SecretorPositive$Caucasian + SecretorPositive$HP_Index + SecretorPositive$Maternal_age)
summary(ModelStarchSecretorPos) #Diversity -0.145 (0.002/ g Starch, p = 0.086, Adj R= 0.010)
ModelStarchSecretorNeg <- lm(SecretorNegative$Shannon_Diversity ~ SecretorNegative$Starch + SecretorNegative$BMI + SecretorNegative$HP_Index + SecretorNegative$Maternal_age) #Remove Education/Caucasian as all were same value - third level, caucasian
summary(ModelStarchSecretorNeg) #Diversity -0.020 (-0.005/ g Starch, p = 0.038, Adj R = 0.164)

#Maltose - secretor status does not impact relationship between maltose & diversity 
ModelMaltoseSecretorPos <- lm(SecretorPositive$Shannon_Diversity ~ SecretorPositive$Maltose + SecretorPositive$BMI + SecretorPositive$Education + SecretorPositive$Caucasian + SecretorPositive$HP_Index + SecretorPositive$Maternal_age)
summary(ModelMaltoseSecretorPos) #Diversity -0.120 (0.033/ g Maltose, p = 0.137, Adj R= -0.006)
ModelMaltoseSecretorNeg <- lm(SecretorNegative$Shannon_Diversity ~ SecretorNegative$Maltose + SecretorNegative$BMI + SecretorNegative$HP_Index + SecretorNegative$Maternal_age) #Remove Education/Caucasian as all were same value - third level, caucasian
summary(ModelMaltoseSecretorNeg) #Diversity -0.079 (-0.003/ g Maltose, p = 0.639, Adj R = -0.158)

#Secretor status impacts relationship between Macronutrients and Diversity - 
#Secretor Positive status sees increase in diversity with increase in macronutrient - does not occur in secretor negative group 
#starch & diversity (secretor negative does not see increase in diversity)
#Fat and diversity (Secretor positive sees increase in diversity, not secretor neg)


# VD and BF and Smoking with Diversity ---------------------------

#Vaginal Delivery and Regression 
group_by(regr_Metadata,Prev_VD) %>%
  summarise(
    count = n(),
    median = median(Shannon_Diversity, na.rm = TRUE),
    IQR = IQR(Shannon_Diversity, na.rm = TRUE))

Wilcox_VD_Diversity <- wilcox.test(Shannon_Diversity ~ Prev_VD, 
                                   data = regr_Metadata,
                                   exact = FALSE)
Wilcox_VD_Diversity

bp_VD_SD <- ggplot(regr_Metadata, aes(x=regr_Metadata$Prev_VD, y=regr_Metadata$Shannon_Diversity, fill = regr_Metadata$Prev_VD, show.legend = FALSE) ) + 
  geom_boxplot(aes(fill=regr_Metadata$Prev_VD), show.legend = FALSE) + 
  scale_fill_manual(values = c("#ccccff", "#9944cc")) +
  xlab('Previous Vaginal Delivery')+
  ylab('Shannon Diversity')+
  ggtitle(label = "Vaginal Delivery", subtitle = "MWU U = 1543.5, p = 0.9246")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
bp_VD_SD
#Prev VD median(iqr) diversity index = 0.0125 (0.428)
#No Prev VD median(iqr) diversity index = 0.0178 (0.185)
#Wilcox rank based test between them showed no signif diff "MWU U = 1543.5, p = 0.9246"

#Breastfeeding history and Regression 
group_by(regr_Metadata,Prev_BF) %>%
  summarise(
    count = n(),
    median = median(Shannon_Diversity, na.rm = TRUE),
    IQR = IQR(Shannon_Diversity, na.rm = TRUE))

Wilcox_BF_Diversity <- wilcox.test(Shannon_Diversity ~ Prev_BF, 
                                   data = regr_Metadata,
                                   exact = FALSE)
Wilcox_BF_Diversity

bp_BF_SD <- ggplot(regr_Metadata, aes(x=regr_Metadata$Prev_BF, y=regr_Metadata$Shannon_Diversity, fill = regr_Metadata$Prev_BF, show.legend = FALSE) ) + 
  geom_boxplot(aes(fill=regr_Metadata$Prev_BF), show.legend = FALSE) + 
  scale_fill_manual(values = c("#ccccff", "#9944cc")) +
  xlab('Previous Breastfeeding')+
  ylab('Shannon Diversity')+
  ggtitle(label = "Prior Breastfeeding", subtitle = "MWU U = 1454, p = 0.466")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
bp_BF_SD
#Prev BF median(iqr) diversity index = 0.029 (0.413)
#No Prev BF median(iqr) diversity index = 0.0118 (0.184)
#Wilcox rank based test between them showed no signif diff "MWU U = 1454, p = 0.466"

#Smoking history and Regression 
group_by(regr_Metadata,Smoking) %>%
  summarise(
    count = n(),
    median = median(Shannon_Diversity, na.rm = TRUE),
    IQR = IQR(Shannon_Diversity, na.rm = TRUE))

Wilcox_Smoking_Diversity <- wilcox.test(Shannon_Diversity ~ Smoking, 
                                        data = regr_Metadata,
                                        exact = FALSE)
Wilcox_Smoking_Diversity

bp_Smoking_SD <- ggplot(regr_Metadata, aes(x=regr_Metadata$Smoking, y=regr_Metadata$Shannon_Diversity, fill = regr_Metadata$Smoking, show.legend = FALSE) ) + 
  geom_boxplot(aes(fill=regr_Metadata$Smoking), show.legend = FALSE) + 
  scale_fill_manual(values = c("#ccccff", "#9944cc")) +
  xlab('Smoking History')+
  ylab('Shannon Diversity')+
  ggtitle(label = "Smoking History", subtitle = "MWU U = 1474, p = 0.374")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
bp_Smoking_SD
#Prev Smoking median(iqr) diversity index = 0.0174 (0.463)
#No Prev Smoking median(iqr) diversity index = 0.0138 (0.185)
#Wilcox rank based test between them showed no signif diff "MWU U = 1474, p = 0.374"

#Obstetric/Lifestyle and Diversity 
bp_VD_SD
bp_BF_SD
bp_Smoking_SD

Figure5a <- ggarrange(bp_VD_SD, bp_BF_SD, bp_Smoking_SD, 
                      labels = c("A", "B", "C"),
                      ncol = 3, nrow = 1)
Figure5a


# Supp Figure - Amino Acids and Regression --------------------------------
#AMINO ACIDS AND DIVERSITY 
ModelProtein <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Protein + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelProtein)
plot(regr_Metadata$Protein, regr_Metadata$Shannon_Diversity, main = "Protein (Adj-R-Squared -0.020, p=0.271)", xlab = "Protein g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelProtein, col=3, lwd=2)

Regr_Protein_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Protein, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkgoldenrod1", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Protein g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Protein Intake", subtitle = "Diversity -0.137 (0.000/g NSP); p=0.271
          Adj R2 -0.020")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Protein_Diversity

#Lysine
ModelLysine <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Lysine + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelLysine)
plot(regr_Metadata$Lysine, regr_Metadata$Shannon_Diversity, main = "Lysine (Adj-R-Squared -0.029, p=0.622)", xlab = "Lysine g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelLysine, col=3, lwd=2)

Regr_Lysine_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Lysine, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkgoldenrod1", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Lysine g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Lysine Intake", subtitle = "Diversity -0.147 (0.000/g Lysine); p=0.622
          Adj R2 -0.029")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Lysine_Diversity

ModelLeucine <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Leucine + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelLeucine)
plot(regr_Metadata$Leucine, regr_Metadata$Shannon_Diversity, main = "Leucine (Adj-R-Squared -0.019, p=0.266)", xlab = "Leucine g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelLeucine, col=3, lwd=2)

Regr_Leucine_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Leucine, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkgoldenrod1", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Leucine g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Leucine Intake", subtitle = "Diversity -0.143 (0.000/g Leucine); p=0.266
          Adj R2 -0.019")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Leucine_Diversity

ModelValine <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Valine + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelValine)
plot(regr_Metadata$Valine, regr_Metadata$Shannon_Diversity, main = "Valine (Adj-R-Squared -0.023, p=0.360)", xlab = "Valine g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelValine, col=3, lwd=2)

Regr_Valine_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Valine, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkgoldenrod1", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Protein g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Valine Intake", subtitle = "Diversity -0.147 (0.000/g Valine); p=0.360
          Adj R2 -0.023")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Valine_Diversity

ModelAlanine <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Alanine + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelAlanine)
plot(regr_Metadata$Alanine, regr_Metadata$Shannon_Diversity, main = "Alanine (Adj-R-Squared -0.023, p=0.367)", xlab = "Alanine g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelAlanine, col=3, lwd=2)

Regr_Alanine_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Alanine, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkgoldenrod1", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Alanine g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Alanine Intake", subtitle = "Diversity -0.144 (0.000/g Alanine); p=0.367
          Adj R2 -0.023")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Alanine_Diversity

ModelPhenylalanine <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$Phenylalanine + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelPhenylalanine)
plot(regr_Metadata$Phenylalanine, regr_Metadata$Shannon_Diversity, main = "Phenylalanine (Adj-R-Squared -0.017, p=0.223)", xlab = "Phenylalanine g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelPhenylalanine, col=3, lwd=2)

Regr_Phenylalanine_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$Phenylalanine, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkgoldenrod1", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Phenylalanine g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Phenylalanine Intake", subtitle = "Diversity -0.143 (0.000/g Phenylalanine); p=0.223
          Adj R2 -0.020")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Phenylalanine_Diversity

#Amino Acids Figure 
Regr_Protein_Diversity
Regr_Lysine_Diversity
Regr_Leucine_Diversity
Regr_Valine_Diversity
Regr_Alanine_Diversity
Regr_Phenylalanine_Diversity

SuppFig_AminoAcids_Diversity <- ggarrange(Regr_Protein_Diversity, Regr_Lysine_Diversity, Regr_Leucine_Diversity, Regr_Valine_Diversity, Regr_Alanine_Diversity, Regr_Phenylalanine_Diversity, 
                                          labels = c("A", "B", "C", "D", "E", "F"),
                                          ncol = 3, nrow = 2)
SuppFig_AminoAcids_Diversity


# Micronutrients and Diversity --------------------------------------------

ModelNiacin <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$`Niacin (preformed)` + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelNiacin)
plot(regr_Metadata$`Niacin (preformed)`, regr_Metadata$Shannon_Diversity, main = "Niacin (Adj-R-Squared -0.029, p=0.631)", xlab = "Niacin g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelNiacin, col=3, lwd=2)

Regr_Niacin_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$`Niacin (preformed)`, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkmagenta", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Niacin g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Niacin Intake", subtitle = "Diversity -0.146 (-0.003/g Niacin); p=0.631
          Adj R2 -0.029")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Niacin_Diversity

ModelVitK1 <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$`Vitamin K1 (Phylloquinone)` + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelVitK1)
plot(regr_Metadata$`Vitamin K1 (Phylloquinone)`, regr_Metadata$Shannon_Diversity, main = "Vitamin K1 (Adj-R-Squared -0.030, p=0.749)", xlab = "Maltose g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelVitK1, col=3, lwd=2)

Regr_VitK1_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$`Vitamin K1 (Phylloquinone)`, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkmagenta", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Vit K1 g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary VitK1 Intake", subtitle = "Diversity -0.145 (0.000/g VitK1); p=0.749
          Adj R2 -0.030")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_VitK1_Diversity

ModelSelenium <- lm(regr_Metadata$Shannon_Diversity ~ regr_Metadata$`Selenium (Se)` + regr_Metadata$BMI + regr_Metadata$Education + regr_Metadata$Caucasian + regr_Metadata$HP_Index + regr_Metadata$Maternal_age + regr_Metadata$Prev_VD)
summary(ModelSelenium)
plot(regr_Metadata$`Selenium (Se)`, regr_Metadata$Shannon_Diversity, main = "Selenium (Adj-R-Squared -0.027, p=0.518)", xlab = "Selenium g/day", ylab = "Shannon Diversity", col=5, pch=19)
abline(ModelSelenium, col=3, lwd=2)

Regr_Selenium_Diversity <- ggplot(regr_Metadata, aes(x=regr_Metadata$`Selenium (Se)`, y=regr_Metadata$Shannon_Diversity)) +
  geom_point(color="darkmagenta", size=2, alpha=0.8)+
  geom_smooth(method = "lm",colour = "red3") + 
  labs(x="Selenium g/day", y="Shannon Diversity")+
  ggtitle(label = "Dietary Selenium Intake", subtitle = "Diversity -0.132 (0.000/g Selenium); p=0.518
          Adj R2 -0.027")+
  theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(t = 5, b = 5)))+
  theme(plot.subtitle = element_text(size=8, hjust=0.5, margin = margin(t = 5, b = 5)))
Regr_Selenium_Diversity

#Create Micronutrient vs Diversity Plot - no relationship between them. 
Regr_Niacin_Diversity
Regr_VitK1_Diversity
Regr_Selenium_Diversity

SuppFig_Micronutrients_Diversity <- ggarrange(Regr_Niacin_Diversity, Regr_VitK1_Diversity, Regr_Selenium_Diversity, 
                                              labels = c("A", "B", "C"),
                                              ncol = 3, nrow = 1)
SuppFig_Micronutrients_Diversity




# ANOSIM Plots - Species and 3 Function -----------------------------------
par(mfrow = c(2, 2))
plot(CST.anoSpec, main = " A. Species Ordination")
plot(CST.anoBP, main = "B. BP Ordination")
plot(CST.anoCC, main = "C. CC Ordination")
plot(CST.anoMF, main = "D. MF Ordination")
