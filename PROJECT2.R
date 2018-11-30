setwd('~/OneDrive - SAP SE/SaanichInlet_Project')
library(tidyr)
library(dplyr)
library(pathview)
library(RColorBrewer)
library(knitr)
library(ggplot2)
library(tidyverse)
library(cowplot)

#ar and bac tsv
arc_class <- read.table("gtdbtk.ar122.classification_pplacer.tsv", sep="\t")
bac_class <- read.table("gtdbtk.bac120.classification_pplacer.tsv", sep="\t")

gtdb_dat <- rbind(arc_class, bac_class) %>% 
  dplyr::rename(mag = V1) %>% 
  separate(V2, sep=';', into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

#KO
ko <- read.table("SaanichInlet_MAGx_ORFs_ko.cleaned.txt") %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(ko = V2)

#checkm
checkm_dat <- read.table("MetaBAT2_SaanichInlet_10m_min1500_checkM_stdout.tsv",
                         header=TRUE,
                         sep="\t",
                         comment.char = '') %>% 
  dplyr::rename(mag = Bin.Id) %>% 
  dplyr::select(mag, Completeness, Contamination) %>% 
  mutate(mag = as.character(mag))

#metag rpkm
metag_rpkm <- read.table("SaanichInlet_10m_binned.rpkm.csv", header=T, sep=',') %>% 
  mutate(Sequence = gsub('m_', 'm.', Sequence)) %>% 
  mutate(Sequence = gsub('Inlet_', 'Inlet.', Sequence)) %>% 
  separate(col=Sequence, into=c("mag", "contig"), sep='_', extra="merge") %>% 
  group_by(Sample, mag) %>% 
  summarise(g_rpkm = sum(RPKM)) %>% 
  mutate(mag = gsub('Inlet.', 'Inlet_', mag))

#KO log 
KO_log_test <- read.table(file="KO_log.csv", header=F, sep=',') %>% 
  separate(V1, into = c("1", "2", "3", "4", "5", "6")) %>% 
  dplyr::select("1", "3", "4") %>%
  dplyr::rename(metabolism = "1", ko = "3", gene = "4")

#prokka 
prokka_mag_map <- read.table("Prokka_MAG_map_1.csv", header=F, sep=',') %>% 
  dplyr::rename(prokka_id = V1) %>% 
  dplyr::rename(mag = V2)

geo_dat <- read_csv("Saanich_Data.csv", col_names=TRUE)
geo_ts_dat <- read_csv("Saanich_TimeSeries_Chemical_DATA.csv", col_names=TRUE)


#Investgating MAG 227
#042 Cruise
rpkm042 <- read.table("SI042.rpkm.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm42 = V2)

#048 Cruise
rpkm048 <- read.table("SI048.rpkm.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm48 = V2)

#072 Cruise
rpkm072 <- read.table("SI072.rpkm.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm72 = V2)

#073 Cruise
rpkm073 <- read.table("SI073.rpkm.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm73 = V2)

#074 Cruise
rpkm074 <- read.table("SI074.rpkm.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm74 = V2)

#075 Cruise
rpkm075 <- read.table("SI075.rpkm.csv", sep=',') %>% 
  dplyr::rename(orf = V1) %>% 
  dplyr::rename(rpkm75 = V2)


SI042_rpkm_dat <- read.table(file="SI042.rpkm.csv", sep=',') %>%
  dplyr::rename(orf=V1, SI042_rpkm=V2)  #sequence ID = orf
SI048_rpkm_dat <- read.table(file = "SI048.rpkm.csv", sep=',') %>%
  dplyr::rename(orf=V1, SI048_rpkm=V2)
SI072_rpkm_dat <- read.table(file = "SI072.rpkm.csv", sep=',') %>%
  dplyr::rename(orf=V1, SI072_rpkm=V2)
SI073_rpkm_dat <- read.table(file = "SI073.rpkm.csv", sep=',') %>%
  dplyr::rename(orf=V1, SI073_rpkm=V2)
SI074_rpkm_dat <- read.table(file = "SI074.rpkm.csv", sep=',') %>%
  dplyr::rename(orf=V1, SI074_rpkm=V2)
SI075_rpkm_dat <- read.table(file = "SI075.rpkm.csv", sep=',') %>%
  dplyr::rename(orf=V1, SI075_rpkm=V2)

## Geo plots
carbondat <-
  geo_ts_dat %>%
  dplyr::select(Cruise, Date, Depth, Mean_co2) %>%
  filter(!is.na(Mean_co2)) %>%
  filter(Depth == 10) %>%
  dplyr::rename(CO2=Mean_co2)


carbondat %>%
  ggplot() +
  geom_point(aes(x=Date, y=CO2), colour = 'purple', size = 3) + 
  labs(y = expression('CO'[2]*' (uM)'))


tempdat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, Temperature) %>%
  filter(!is.na(Temperature)) %>%
  mutate(Depth_m=Depth*1000) %>%
  filter(Depth_m == 10)

tempdat %>%
  ggplot() +
  geom_point(aes(x=Date, y=Temperature), colour = 'red', size = 3) + 
  labs(y = expression('Temperature (' ^o*'C)'))


saldat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, Salinity) %>%
  filter(!is.na(Salinity)) %>%
  mutate(Depth_m=Depth*1000) %>%
  filter(Depth_m == 10)
saldat %>%
  ggplot() +
  geom_point(aes(x=Date, y=Salinity), colour = 'blue', size = 3) +
  labs(y = expression('Salinity (psu)'))

phosdat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, WS_PO4) %>%
  dplyr::rename(PO4=WS_PO4) %>%
  filter(!is.na(PO4)) %>%
  mutate(Depth_m=Depth*1000) %>%
  filter(Depth_m == 10) %>%
  dplyr::select(Date, PO4)

sidat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, SI) %>%
  filter(!is.na(SI)) %>%
  mutate(Depth_m=Depth*1000) %>%
  filter(Depth_m == 10) %>%
  dplyr::select(Date, SI)
phos_si_merge <- Reduce(function(x, y) merge(x,y, all =TRUE), list(phosdat, sidat))

phos_si_merge %>%
  dplyr::select(Date, PO4, SI) %>% 
  gather(key="Nutrients", value="Concentration", -Date) %>% 
  
  ggplot(aes(x=Date, y=Concentration, shape=Nutrients, color=Nutrients, size = 3)) +
  geom_point() +
  facet_wrap(~Nutrients, scales="free") +
  labs(y = expression('Concentration (uM)'))


nitratedat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, WS_NO3) %>%
  dplyr::rename(NO3_uM=WS_NO3) %>%
  filter(!is.na(NO3_uM)) %>%
  mutate(Depth_m=Depth*1000) %>%
  dplyr::select(Date, NO3_uM) 

nitritedat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, Mean_NO2) %>%
  dplyr::rename(NO2_uM=Mean_NO2) %>%
  filter(!is.na(NO2_uM)) %>%
  mutate(Depth_m=Depth*1000) %>%
  dplyr::select(Date, NO2_uM)

ammodat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, Mean_NH4) %>%
  dplyr::rename(NH4_uM=Mean_NH4) %>%
  filter(!is.na(NH4_uM)) %>%
  mutate(Depth_m=Depth*1000) %>%
  dplyr::select(Date, NH4_uM)

nitrousdat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, Mean_N2O) %>%
  dplyr::rename(N2O_uM=Mean_N2O) %>%
  filter(!is.na(N2O_uM)) %>%
  mutate(Depth_m=Depth*1000) %>%
  dplyr::select(Date, N2O_uM)

nitcomp_merge1 <- Reduce(function(x, y) merge(x,y, all =TRUE), list(ammodat, nitratedat))
nitcomp_merge2 <- Reduce(function(x, y) merge(x,y, all =TRUE), list(nitritedat, nitrousdat))

nitcomp_merge1 %>%
  dplyr::select(Date, NO3_uM, NH4_uM) %>% 
  dplyr::rename(NO3=NO3_uM) %>%
  dplyr::rename(NH4=NH4_uM) %>%
  gather(key="Nitrogen_Compounds", value="Concentration", -Date) %>% 
  ggplot(aes(x=Date, y=Concentration, shape=Nitrogen_Compounds, color=Nitrogen_Compounds)) +
  geom_point() +
  facet_wrap(~Nitrogen_Compounds, scales="free") +
  labs(y = expression('Concentration (uM)')) +
  scale_color_manual(values=c("#00CED1", "#32CD32"))

nitcomp_merge2 %>%
  dplyr::select(Date, NO2_uM, N2O_uM) %>% 
  dplyr::rename(NO2=NO2_uM) %>%
  dplyr::rename(N2O=N2O_uM) %>%
  gather(key="Nitrogen_Compounds", value="Concentration", -Date) %>% 
  ggplot(aes(x=Date, y=Concentration, shape=Nitrogen_Compounds, color=Nitrogen_Compounds)) +
  geom_point() +
  facet_wrap(~Nitrogen_Compounds, scales="free") +
  labs(y = expression('Concentration (uM)')) +
  scale_color_manual(values=c("#FF69B4", "#87CEEB"))
oxydat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, WS_O2) %>%
  dplyr::rename(O2_uM=WS_O2) %>%
  filter(!is.na(O2_uM)) %>%
  mutate(Depth_m=Depth*1000) %>%
  filter(Depth_m == 10)

oxydat %>%
  ggplot() +
  geom_point(aes(x=Date, y=O2_uM), colour = 'dark green', size = 3) +
  labs(y = expression('O'[2]* ' (uM)'))

sulfidedat <-
  geo_dat %>%
  dplyr::select(Cruise, Date, Depth, WS_H2S) %>%
  dplyr::rename(H2S_uM=WS_H2S) %>%
  filter(!is.na(H2S_uM)) %>%
  mutate(Depth_m=Depth*1000) %>%
  filter(Depth_m == 10)


sulfidedat %>%
  ggplot() +
  geom_point(aes(x=Date, y=H2S_uM), colour = 'black', size = 3) +
  labs(y = expression('H'[2]*'S (uM)'))

##END GEO
## Bubble plot

Nitrogen_KO_log <- subset(KO_log_test, metabolism == "Nitrogen",
                          select=c(ko, gene))

Sulfur_KO_log <- subset(KO_log_test, metabolism == "Sulfur",
                        select=c(ko, gene))

SI_ALL <- Reduce(function(x, y) left_join(x,y, all =TRUE), list(SI042_rpkm_dat, SI048_rpkm_dat, SI072_rpkm_dat, SI073_rpkm_dat, SI074_rpkm_dat, SI075_rpkm_dat))

#Calculate the total RPKM across all cruises
rpkm_total <- mutate(SI_ALL, rpkm_sum = (apply(SI_ALL[2:7], 1, sum))) %>%
  dplyr::select(orf, rpkm_sum) 

SI_ALL_ko <- left_join(SI_ALL, ko, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id"))  

SI_ALL_NS_ko <- merge(SI_ALL_ko, KO_log_test, by="ko") %>%
  dplyr::select(prokka_id, ko, gene, SI042_rpkm, SI048_rpkm, SI072_rpkm, SI073_rpkm, SI074_rpkm, SI075_rpkm) #mag = prokka_id

SI_rpkm_ko <- mutate(SI_ALL_NS_ko, total = (apply(SI_ALL_NS_ko[4:9], 1, sum)))

All_data_test <- left_join(SI_rpkm_ko, prokka_mag_map, by="prokka_id")

All_data_final <- left_join(All_data_test, gtdb_dat, by="mag")

All_data_final_f <- left_join(All_data_final, KO_log_test)

Bubble <- left_join(All_data_final_f, checkm_dat, by="mag") %>% 
  filter(Completeness > 50 & Contamination < 10)# %>% #good quality MAGs
#  dplyr::select(mag) %>%
#  str(data_frame(mag))
#  reorder(dplyr::Bubble$mag, Bubble$Phylum, sort(Bubble$mag, decreasing = FALSE))
#  dplyr::mutate(mag)
# mag  = reorder(mag,Phylum, sort)
#  reorder(mag,total, sort(mag, decreasing = FALSE)) %>%
#  dplyr::mutate(mag)
#  dplyr::mutate(mag = reorder(mag,Phylum, sort(mag, decreasing = FALSE))) # sort by their taxonomic Order so everything shows up together
# filter(Phylum == "p__Proteobacteria") %>% # Proteobacteria only; this can easily be changed


n_genes <- Bubble %>%
  mutate(total = if_else(total > 400, 400, total)) %>% 
  filter(metabolism == "Nitrogen")

s_genes <- Bubble %>%
  mutate(total = if_else(total > 400, 400, total)) %>%
  filter(metabolism == "Sulfur")

ggplot(n_genes, aes(x=gene, y=Class, col=Phylum)) +
  geom_point(aes(size=total)) +
  xlab("Nitrogen Cycle Metabolic Genes") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))

ggplot(s_genes, aes(x=gene, y=Class, col=Phylum)) +
  geom_point(aes(size=total)) +
  xlab("Sulfur Cycle Metabolic Genes") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(angle = 0),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))

#Contamination vs Completeness Bubble Plot

rpkm_dat <- left_join(metag_rpkm, checkm_dat, by="mag") %>% 
  left_join(gtdb_dat, by="mag") %>% 
  group_by(mag, Kingdom, Phylum, Class, Completeness, Contamination) %>% 
  summarise(g_rpkm = mean(g_rpkm))

ggplot(rpkm_dat, aes(x=Completeness, y=Contamination, col=Class)) +
  geom_point(aes(size=g_rpkm)) +
  scale_size(range=c(1,10)) +
  xlim(c(50,100)) +
  ylim(c(0,30)) +
  geom_hline(mapping=NULL, data=NULL, yintercept=5)+
  geom_hline(mapping=NULL, data=NULL, yintercept=10)+
  geom_hline(mapping=NULL, data=NULL, yintercept=15)+
  geom_vline(mapping=NULL, data=NULL, xintercept=5)+
  geom_hline(mapping=NULL, data=NULL, yintercept=10)+
  geom_hline(mapping=NULL, data=NULL, yintercept=15)+  
  theme(panel.background = element_blank(),
        panel.grid.major = element_line(colour = "#bdbdbd", linetype = "dotted"),
        panel.grid.minor = element_blank())

## ENDBUBBLE
##BEGIN PATHVIEW

merge
#Merging all cruises
rpkm_all <- Reduce(function(x, y) merge(x,y, all =TRUE), list(rpkm042, rpkm048, rpkm072, rpkm073, rpkm074, rpkm075))
#Summing RPKM across cruises
rpkm_total <- mutate(rpkm_all, rpkm_sum = (apply(rpkm_all[2:7], 1, sum))) %>%
  dplyr::select(orf, rpkm_sum)

rpkm_avg <- mutate(rpkm_total, rpkm_avg=rpkm_sum/6) %>%
  dplyr::select(orf,rpkm_avg)

ko_rpkm_total <- left_join(ko, rpkm_total, by="orf") %>% 
  separate(orf, into=c("prokka_id", "orf_id")) %>% # Split the Prokka ORF names into MAG identifier and ORF number for joining
  left_join(prokka_mag_map, by="prokka_id") %>% 
  left_join(gtdb_dat, by="mag")

t_rpkm_all <- ko_rpkm_total %>% 
  filter(mag == "SaanichInlet_10m.227") %>% 
  group_by(mag, ko) %>% 
  summarise(total = sum(rpkm_sum)) %>% 
  spread(key = mag, value = total)


pv_mat_all <- dplyr::select(t_rpkm_all, -ko)
rownames(pv_mat_all) <- t_rpkm_all$ko

pv.outall_N <- pathview(gene.data = pv_mat_all,
                        limit = list(gene = c(0,4)),
                        low = list(gene = "#91bfdb"),
                        mid = list(gene = "#ffffbf"),
                        high = list(gene = "#fc8d59"),
                        species = "ko",
                        pathway.id="00910",
                        kegg.dir = "~/OneDrive - SAP SE/SaanichInlet_Project")

pv.outall_S <- pathview(gene.data = pv_mat_all,
                        limit = list(gene = c(0,4)),
                        low = list(gene = "#91bfdb"),
                        mid = list(gene = "#ffffbf"),
                        high = list(gene = "#fc8d59"),
                        species = "ko",
                        pathway.id="00920",
                        kegg.dir = "~/OneDrive - SAP SE/SaanichInlet_Project")

pv.outall_CH4 <- pathview(gene.data = pv_mat_all,
                          limit = list(gene = c(0,4)),
                          low = list(gene = "#91bfdb"),
                          mid = list(gene = "#ffffbf"),
                          high = list(gene = "#fc8d59"),
                          species = "ko",
                          pathway.id="00680",
                          kegg.dir = "~/OneDrive - SAP SE/SaanichInlet_Project")

pv.outall_oxy <- pathview(gene.data = pv_mat_all,
                          limit = list(gene = c(0,4)),
                          low = list(gene = "#91bfdb"),
                          mid = list(gene = "#ffffbf"),
                          high = list(gene = "#fc8d59"),
                          species = "ko",
                          pathway.id="00190",
                          kegg.dir = "~/OneDrive - SAP SE/SaanichInlet_Project")

pv.outall_glyco <- pathview(gene.data = pv_mat_all,
                            limit = list(gene = c(0,4)),
                            low = list(gene = "#91bfdb"),
                            mid = list(gene = "#ffffbf"),
                            high = list(gene = "#fc8d59"),
                            species = "ko",
                            pathway.id="00010",
                            kegg.dir = "~/OneDrive - SAP SE/SaanichInlet_Project")

##END PATHVIEW

