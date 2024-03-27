ant_full<-read.csv("Riparian Community Final.csv")
Ant_community<-ant_full[,6:23]
ant_full$Habitat.Type<-as.factor(ant_full$Habitat.Type)
ant_full$Season<-as.factor(ant_full$Season)

library(vegan)
library(tidyverse)
library(ggtext)
library(FactoMineR)
library(corrplot)
library(factoextra)

###sampling a random number to set a seed but remain random in analysis starting point
sample(1:999, 1)
set.seed(257) ###257 was the number selected in our case

###NMDS Ordination Analysis
ant.nmds <- metaMDS(Ant_community, k=2, model = "global", distance = "bray")
ant.nmds$stress

###Permutational ANOVA by Microhabitat
ant.dist<- vegdist(Ant_community, method = "bray") ###calculates the Bray-Curtis dissimilarity between survey locations as was used in the NMDS ordination
set.seed(257)
PERMANOVA<-adonis2(ant.dist~ Habitat.Type, data=ant_full, strata = ant_full$Location,
                    permutations = 1000)
PERMANOVA

###envfit of vegetation influences on differentiation in ant communities
Total.enviro<-cbind(ant_full[,3:5], ant_full[,25:81]) ###subsetting environmental data from full dataset
sample(300, 3) ###sampling 3 random numbers as starting points to ensure results are not an aberration
set.seed(252) ###seed selected for data reporting
ef<-envfit(ant.nmds, Total.enviro, permu=999, na.rm= TRUE )
ef$factors
plot(ant.nmds) +
  ordiellipse(ant.nmds, groups = ant_full$Habitat.Type, draw = "polygon", 
              label = TRUE, lty = 1)
plot(ef, p.max = 0.05)

###Stacked barplot by species (Figure 3)
Ant_community_long<- Ant_community %>%
  mutate(Habitat=ant_full$Habitat.Type) %>%
  pivot_longer(!Habitat, names_to = "species", values_to = "count")

ant_c_l_sum <- Ant_community_long %>%
  group_by(Habitat, species) %>%
  summarise(across(everything(),sum))

ant_c_l_sum$Habitat <- factor(ant_c_l_sum$Habitat, levels=c("Arundo", "R. Forest", "Revegetated",
                                                            "F. Agriculture", "R. Channel", "R. Scrub"))

colvec<-rev(tmaptools::get_brewer_pal("Spectral", n = 18)) ###Getting the color palette for the plot

ggplot(data = ant_c_l_sum, aes(x= Habitat, y =count, fill=species)) +
  geom_bar(position = "fill", stat = "identity") +
  theme_classic()+
  labs(x='Habitat', y='Proportion of Total Encounters', fill = 'Species', face='bold') +
  scale_fill_manual(values=colvec,
                    labels = c("*Cardiocondyla mauritanica*",
                               "*Crematogaster hespera*",
                               "*Dorymyrmex insanus*",
                               "*Forelius pruinosus*",
                               "*Formica francoeuri*",
                               "*Hypoponera opaciceps*",
                               "*Linepithema humile*",
                               "*Liometopum occidentale*",
                               "*Monomorium ergatogyna*",
                               "*Pheidole hyatti*",
                               "*Pogonomyrmex californicus*",
                               "*Solenopsis molesta*",
                               "*Solenopsis validiuscula*",
                               "*Solenopsis xyloni*",
                               "*Stenamma punctatoventre*",
                               "*Tapinoma sessile*",
                               "*Temnothorax andrei*",
                               "*Temnothorax nitens*")) +
  theme(legend.text = element_markdown()) +
  theme(axis.text = element_text(size = 10))

###NMDS Ordination (Figure 4)
pchvec<-c(15,17)
library(RColorBrewer)
brewer.pal(6, "Set1")
colvec<-c("#E41A1C","#FFFF33", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
plot(ant.nmds)
orditorp(ant.nmds, display = 'species')
with(ant_full,
     points(ant.nmds, cex=1.0,
            col = colvec[Habitat.Type],
            pch = pchvec[Season]))
mtext(text = "Stress = 0.092", side = 3, line = -1.5, adj = 0.95, 
      cex = 1.0)
legend("bottomright", c("Arundo", "F. Agriculture", "R. Channel", "R. Forest", "R. Scrub", "Revegetated"), title = "Microhabitat", 
       col= c("#E41A1C","#FFFF33", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00"), pch=15, cex=0.8)
ordiellipse(ant.nmds, groups = ant_full$Habitat.Type, draw = "polygon", 
            label = F, lty = 1, col = colvec)
orditorp(ant.nmds, display = 'species')
ordiellipse(ant.nmds, groups = Ant_community$cluster, draw = "polygon")

###Vectoring Significant Environmental Drivers/Plants (Figure 5)
library(dplyr)
library(tibble)
  ###creating a dataframe of the analyses for ggplot
plot_df <- scores(ant.nmds, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("rowname") %>% 
  full_join(rownames_to_column(ant_full[,3:4])) %>% 
  dplyr::select(-rowname) %>% 
  mutate(Habitat.Type = factor(Habitat.Type))

fit_pvals <- ef$vectors$pvals %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  dplyr::rename("pvals" = ".")

fit_spp <- ef %>% 
  scores(., display = "vectors") %>% 
  as.data.frame() %>% 
  rownames_to_column("species") %>% 
  full_join(., fit_pvals, by = "species") %>% 
  filter(pvals <= 0.05) %>% 
  mutate(label = case_when(species == "Litter" ~ "Leaf Litter",
                           species == "Soil" ~ "Exposed Soil",
                           species == "Mulch" ~ "Mulch Cover",
                           species == "canopy.cover" ~ "Canopy Cover",
                           species == "soil.moisture" ~ "Soil Moisture",
                           species == "percent.clay" ~ "Clay Content",
                           species == "percent.sand" ~ "Sand Content",
                           species == "Arundo.donax" ~ "Arundo donax",
                           species == "Baccharis.pilularis" ~ "Baccharis pilularis",
                           species == "Bromus.madritensis.ssp..rubens" ~ "Bromus rubens",
                           species == "Carduus.pycnocephalus" ~ "Carduus pycnocephalus",
                           species == "Curcubita.foetidissima" ~ "Curcubita foetidissima",
                           species == "Erigeron.canadensis" ~ "Erigeron canadensis",
                           species == "Heliotropium.curassavicum" ~ "Heliotropium curassavicum",
                           species == "Hirschfeldia.incana" ~ "Hirschfeldia incana",
                           species == "Malva.parviflora" ~ "Malva parviflora",
                           species == "Populous.trichocarpa" ~ "Populous trichocarpa",
                           species == "Rubus.ursinus" ~ "Rubus ursinus",
                           species == "Salix.lasiolepis" ~ "Salix lasiolepis",
                           species == "Salsola.tragus" ~ "Salsola tragus",
                           species == "Sisyimbrium.irio" ~ "Sisyimbrium irio",
                           species == "Urtica.dioica" ~ "Urtica dioica",
                           species == "Total.native.Cover" ~ "Total Native Cover",
                           species == "Total.Non.native.Cover" ~ "Total Non-native Cover"))
env_spp<-fit_spp[1:7,]
pla_spp<-fit_spp[8:20,]
  ###the function the Vegan package uses for plotting ellipses
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}
  
df_ell <- data.frame()
for(g in levels(plot_df$Habitat.Type)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(plot_df[plot_df$Habitat.Type==g,],
                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),wt=rep(1/length(NMDS1),length(NMDS1)))$cov,center=c(mean(NMDS1),mean(NMDS2)))))
                                ,Habitat.Type=g))
}

library(ggrepel)
p_plant<- ggplot(plot_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Habitat.Type, shape = Season), cex = 3, stroke = 2) +
  scale_shape_manual(values = c(15,17), name = "Season") +
  scale_color_manual(values = c("Arundo" = "#E41A1C", "F. Agriculture" = "#FFFF33", "R. Channel" = "#377EB8", 
                                "R. Forest" = "#4DAF4A", "R. Scrub" = "#984EA3", "Revegetated" = "#FF7F00"), name = "Habitat Type") +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=Habitat.Type), lty = 1, size=1.5) +
  geom_segment(data = pla_spp, aes(x = 0, xend = NMDS1*2, y = 0, yend = NMDS2*2), size=0.8,
               arrow = arrow(length = unit(0.5, "cm")), lty = 1,
               col = "black") +
  theme_classic()+
  theme(legend.text = element_text(size=15), legend.title = element_text(size=15))+
  geom_label_repel(data = pla_spp, aes(label = label, x = NMDS1*2,y = NMDS2*2), size=5, fontface = "italic", color = "black")


p_enviro<- ggplot(plot_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Habitat.Type, shape = Season), cex = 3, stroke = 2) +
  scale_shape_manual(values = c(15,17), name = "Season") +
  scale_color_manual(values = c("Arundo" = "#E41A1C", "F. Agriculture" = "#FFFF33", "R. Channel" = "#377EB8", 
                                "R. Forest" = "#4DAF4A", "R. Scrub" = "#984EA3", "Revegetated" = "#FF7F00"), name = "Habitat Type") +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=Habitat.Type), lty = 1, size=1.5) +
  geom_segment(data = env_spp, aes(x = 0, xend = NMDS1*2, y = 0, yend = NMDS2*2), size=0.8,
               arrow = arrow(length = unit(0.5, "cm")), lty = 1,
               col = "black") +
  theme_classic()+
  theme(legend.text = element_text(size=15), legend.title = element_text(size=15))+
  geom_label_repel(data = env_spp, aes(label = label, x = NMDS1*2,y = NMDS2*2), size=5)

library(ggpubr)
ggarrange(p_plant, p_enviro, nrow = 1, common.legend = T, legend = "right", align = "h", labels = c("A","B"))

###Variation Partitioning (Figure 6)
varpar_table_no_NA<-na.omit(ant_full)
sigplants<-cbind(varpar_table_no_NA$Arundo.donax, varpar_table_no_NA$Baccharis.pilularis, varpar_table_no_NA$Bromus.madritensis.ssp..rubens,
                 varpar_table_no_NA$Curcubita.foetidissima, varpar_table_no_NA$Heliotropium.curassavicum, varpar_table_no_NA$Hirschfeldia.incana,
                 varpar_table_no_NA$Malva.parviflora, varpar_table_no_NA$Populous.trichocarpa, varpar_table_no_NA$Rubus.ursinus,
                 varpar_table_no_NA$Salix.lasiolepis, varpar_table_no_NA$Salsola.tragus, varpar_table_no_NA$Urtica.dioica)
var_dist<-vegdist(varpar_table_no_NA[,6:23], method = "bray")
v.part<-varpart(var_dist, sigplants, varpar_table_no_NA[,25:32])        
plot (v.part, digits = 2, Xnames = c('Sig. Plants', 'Other Env.'), bg = c('forestgreen', 'blue3'))

        
###summary stats calculations for supplemental table
Enviro_means<- Total.enviro %>%
group_by(`Habitat Type`) %>%
summarise_at(vars(Litter:`Total Non-native Cover`), mean, na.rm=T) 
Enviro_sd<- Total.enviro %>%
group_by(`Habitat Type`)%>%
summarise_at(vars(Litter:`Total Non-native Cover`), sd, na.rm=T)

###testing biogeography and ecological diversity null hypotheses (Supplemental Figure 1)
library(rareNMtests)

###Ecological null model
sample(1:999, 1)
set.seed(44)
cbecoq0<-EcoTest.sample(Ant_community, by=ant_full$Habitat.Type, MARGIN = 1, niter=200, method = "coverage", q=0)
windows();plot(cbecoq0)###sig diff, combined with non-significant bio null test means it is Cayuela et al 2015 
### scenario 2 of figure 1 (Differing community composition) for overall sample set.

###Biogeographic null model
set.seed(792)
cbbiogq0<-BiogTest.sample(Ant_community, by=ant_full$Habitat.Type, MARGIN = 1, niter=200, q=0, method = "coverage")
windows();plot(cbbiogq0)###not sig

