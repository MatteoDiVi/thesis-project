rm(list=ls()[! ls() %in% c("maf_filtered","lung","lung_summarised", "maf_filtered_cellularity")])

library(maftools)
library(dplyr)
library(tidyverse)
library(mclust)
library(ggplot2)
library(reshape2)
library(ggpubr)

lung <- filter(maf_filtered_cellularity, Tumor_Type == "Lung", Allele_Frequency >= 0.05)
unique(lung$Tumor_Type)
unique(lung$Hugo_Symbol)

lung_summarised <- group_by(lung, Tumor_Sample_Barcode, Hugo_Symbol, Tumor_Type) %>%
  summarise(Allele_Frequency = mean(Allele_Frequency))%>%
  mutate(Tumor_Sample_Barcode=as.character(Tumor_Sample_Barcode))



# density plot ------------------------------------------------------------


# nras <- filter(lung, Hugo_Symbol == "NRAS", Allele_Frequency >= 0.05)
# 
# NRAS <- nras %>%
#   select_if(function(col) length(unique(col))>1) 
# unique(NRAS$Tumor_Sample_Barcode)
# unique(maf_filtered_cellularity$Tumor_Sample_Barcode)
# NRAS_Het <- ggplot(NRAS, aes(x=Allele_Frequency)) + 
#   geom_density()
# NRAS_Het

lung_het <- ggplot(lung, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
lung_het

# mutated in cbioportal for lung == TP53 = 50%   PIK3CA = 5%  FBXW7 = 9% (all missense mutations)



first3 <- filter(lung, Hugo_Symbol %in%  c("AKT1","ALK","BRAF"))
four_to_6 <- filter(lung, Hugo_Symbol %in%  c("CTNNB1","DDR2","EGFR"))
seven_to_9 <- filter(lung, Hugo_Symbol %in%  c("ERBB4","FBXW7","FGFR1"))
ten_to_12 <- filter(lung, Hugo_Symbol %in%  c("FGFR2","FGFR3","KRAS"))
thirteen_to_15 <- filter(lung, Hugo_Symbol %in%  c("NOTCH1","NRAS","PIK3CA"))
sixteen_to_18 <- filter(lung, Hugo_Symbol %in%  c("MET","ERBB2","TP53"))
last_3 <- filter(lung, Hugo_Symbol %in%  c("PTEN","SMAD4","STK11"))

first3_het <- ggplot(first3, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
first3_het

four_to_6_het <- ggplot(four_to_6, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
four_to_6_het

seven_to_9_het <- ggplot(seven_to_9, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
seven_to_9_het

ten_to_12_het <- ggplot(ten_to_12, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
ten_to_12_het

thirteen_to_15_het <- ggplot(thirteen_to_15, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
thirteen_to_15_het

sixteen_to_18_het <- ggplot(sixteen_to_18, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
sixteen_to_18_het

last_3_het <- ggplot(last_3, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
last_3_het




# cbioportal --------------------------------------------------------------
#c-bioportal calculations
setwd("C:/Users/matte/Desktop/Uni/tesi/cbioportal/lung/")
library(tidyverse)

lung_cbioportal_DDR2 <- read_tsv("DDR2_cbioportal_lung.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

lung_cbioportal_EGFR <- read_tsv("EGFR_cbioportal_lung.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

lung_cbioportal_FGFR1 <- read_tsv("FGFR1_cbioportal_lung.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

lung_cbioportal_KRAS <- read_tsv("KRAS_cbioportal_lung.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

lung_cbioportal_PIK3CA <- read_tsv("PIK3CA_cbioportal_lung.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

lung_cbioportal_STK11 <- read_tsv("STK11_cbioportal_lung.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

lung_cbioportal_TP53 <- read_tsv("TP53_cbioportal_lung.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

lung_cbioportal <- dplyr::bind_rows(list(DDR2=lung_cbioportal_DDR2,
                                         EGFR=lung_cbioportal_EGFR,
                                         FGFR1=lung_cbioportal_FGFR1,
                                         KRAS=lung_cbioportal_KRAS,
                                         PIK3CA=lung_cbioportal_PIK3CA,
                                         STK11=lung_cbioportal_STK11,
                                         TP53=lung_cbioportal_TP53),
                                         .id = 'Hugo_Symbol')

# unfiltered --------------------------------------------------------------
lung_unfiltered <- filter(maf_filtered, Tumor_Type == "Lung")
lung_summarised_unfiltered <- group_by(lung_unfiltered, Tumor_Sample_Barcode, Hugo_Symbol, Tumor_Type) %>%
  summarise(Allele_Frequency = mean(Old_Allele_Frequency))%>%
  mutate(Tumor_Sample_Barcode=as.character(Tumor_Sample_Barcode))

lung_summarised_unfiltered_graph <- filter(lung_summarised_unfiltered, Hugo_Symbol %in% c("DDR2","EGFR","FGFR1","KRAS","PIK3CA","STK11","TP53"))

lungtotunfiltered_vs_lungcbioportal <- dplyr::bind_rows(list(cbioportal=lung_cbioportal,
                                                             "Ion-Ampliseq Colon Lung"=lung_summarised_unfiltered_graph),
                                              .id = 'Dataset')

het_plot <- ggplot() +
  geom_density(alpha=0.4, data = lungtotunfiltered_vs_lungcbioportal, aes(x=Allele_Frequency, group= Dataset, fill=Dataset)) +
  ggtitle("Lung") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
het_plot+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

wcx <-wilcox.test(lung_cbioportal$Allele_Frequency,lung_summarised_unfiltered_graph$Allele_Frequency)




ggplot(data = lungtotunfiltered_vs_lungcbioportal, aes(y = Allele_Frequency,
                                             x = Dataset)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Dataset,
                                                  scales="free") #faceting


ggplot(data = lungtotunfiltered_vs_lungcbioportal, aes(y = Allele_Frequency,
                                             x = Hugo_Symbol)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Hugo_Symbol,
                                                  scales="free") +
  ggtitle("Lung") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  stat_compare_means(aes(group = Dataset), method = "t.test", label = "p.format", vjust = 1)


# total -------------------------------------------------------------------
lung_summarised_graph <- filter(lung_summarised, Hugo_Symbol %in% c("DDR2","EGFR","FGFR1","KRAS","PIK3CA","STK11","TP53"))

lungtot_vs_lungcbioportal <- dplyr::bind_rows(list(cbioportal=lung_cbioportal,
                                                   "Ion-Ampliseq Colon Lung"=lung_summarised_graph),
                                         .id = 'Dataset')

het_plot <- ggplot() +
  geom_density(alpha=0.4, data = lungtot_vs_lungcbioportal, aes(x=Allele_Frequency, group= Dataset, fill=Dataset)) +
  ggtitle("Lung") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
het_plot+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

ggplot(data = lungtot_vs_lungcbioportal, aes(y = Allele_Frequency,
                                        x = Dataset)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Dataset,
                                                  scales="free") #faceting

ggplot(data = lungtot_vs_lungcbioportal, aes(y = Allele_Frequency,
                                             x = Hugo_Symbol)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Hugo_Symbol,
                                                  scales="free") +
  ggtitle("Lung") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  stat_compare_means(aes(group = Dataset), method = "t.test", label = "p.format", vjust = 1)


# ERBB4 -------------------------------------------------------------------


# mean(lung_cbioportal_ERBB4$`Allele Freq (T)`)
#prova
prova <- filter(lung_summarised, Tumor_Sample_Barcode == "105", Hugo_Symbol == "TP53")%>%
  select(Tumor_Sample_Barcode, Hugo_Symbol, Allele_Frequency)

#normal
# lung_ERBB4 <- filter(lung, Hugo_Symbol == "ERBB4")
# filter(lung, Hugo_Symbol == "ERBB4") %>%

lung_TP53 <- filter(lung_summarised, Hugo_Symbol == "TP53")%>%
  mutate(Tumor_Sample_Barcode=as.character(Tumor_Sample_Barcode))
# mean(mean_sample_ERBB4$Allele_Frequency)


normal_vs_cbioportal <- dplyr::bind_rows(list(cbioportal=lung_cbioportal_TP53,
                                              normal=lung_TP53),
                                         .id = 'Dataset')


het_plot <- ggplot() +
  geom_density(alpha=0.4, data = normal_vs_cbioportal, aes(x=Allele_Frequency, group= Dataset, color=Dataset, fill=Dataset)) +
  labs(title="TP53")
het_plot     


#boxplot
# ggplot(data = normal_vs_cbioportal, aes(x = Allele_Frequency,
#                                    y = Dataset)) + 
#   geom_boxplot(aes(color = Dataset)) #plotting both tcga and geo

ggplot(data = normal_vs_cbioportal, aes(y = Allele_Frequency,
                                   x = Dataset)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Dataset,
                                                  scales="free") #faceting


# oncoplot ----------------------------------------------------------------

maf_lung <- read.maf(maf = lung)
oncoplot(maf=maf_lung, top = 10, pwLineWd = 0.0000001, borderCol = NA, fill = F, bgCol = "grey100")
