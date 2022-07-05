rm(list=ls()[! ls() %in% c("maf_filtered","colon", "maf_filtered_cellularity")])

library(maftools)
library(dplyr)
library(tidyverse)
library(mclust)
library(ggplot2)

colon <- filter(maf_filtered_cellularity, Tumor_Type == "Colon", Allele_Frequency >= 0.05)
colon_summarised <- group_by(colon, Tumor_Sample_Barcode, Hugo_Symbol, Tumor_Type) %>%
  summarise(Allele_Frequency = mean(Allele_Frequency))%>%
  mutate(Tumor_Sample_Barcode=as.character(Tumor_Sample_Barcode))
unique(colon$Tumor_Type)
unique(colon$Hugo_Symbol)



# cbioportal --------------------------------------------------------------
#c-bioportal calculations
setwd("C:/Users/matte/Desktop/Uni/tesi/cbioportal/colon/")
library(tidyverse)
colon_cbioportal_ERBB4 <- read_tsv("ERBB4_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_SMAD4 <- read_tsv("SMAD4_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_NOTCH1 <- read_tsv("NOTCH1_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_KRAS <- read_tsv("KRAS_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_PIK3CA <- read_tsv("PIK3CA_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_ALK <- read_tsv("ALK_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_TP53 <- read_tsv("TP53_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_PTEN <- read_tsv("PTEN_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_BRAF <- read_tsv("BRAF_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal_FBXW7 <- read_tsv("FBXW7_cbioportal_colon.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

colon_cbioportal <- dplyr::bind_rows(list(ERBB4=colon_cbioportal_ERBB4,
                                         ALK=colon_cbioportal_ALK,
                                         NOTCH1=colon_cbioportal_NOTCH1,
                                         KRAS=colon_cbioportal_KRAS,
                                         PIK3CA=colon_cbioportal_PIK3CA,
                                         SMAD4=colon_cbioportal_SMAD4,
                                         TP53=colon_cbioportal_TP53,
                                         FBXW7=colon_cbioportal_FBXW7,
                                         BRAF=colon_cbioportal_BRAF,
                                         PTEN=colon_cbioportal_PTEN),
                                    .id = 'Hugo_Symbol')

# unfiltered --------------------------------------------------------------
colon_unfiltered <- filter(maf_filtered, Tumor_Type == "Colon")
colon_summarised_unfiltered <- group_by(colon_unfiltered, Tumor_Sample_Barcode, Hugo_Symbol, Tumor_Type) %>%
  summarise(Allele_Frequency = mean(Old_Allele_Frequency))%>%
  mutate(Tumor_Sample_Barcode=as.character(Tumor_Sample_Barcode))

colon_summarised_unfiltered_graph <- filter(colon_summarised_unfiltered, Hugo_Symbol %in% c("ERBB4","BRAF","KRAS","PIK3CA","FBXW7","TP53", "SMAD4", "PTEN"))

colontotunfiltered_vs_coloncbioportal <- dplyr::bind_rows(list(cbioportal=colon_cbioportal,
                                                             "Ion-Ampliseq Colon Lung"=colon_summarised_unfiltered_graph),
                                                        .id = 'Dataset')

het_plot <- ggplot() +
  geom_density(alpha=0.4, data = colontotunfiltered_vs_coloncbioportal, aes(x=Allele_Frequency, group= Dataset, fill=Dataset)) +
  ggtitle("Colon") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
het_plot+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

ggplot(data = colontotunfiltered_vs_coloncbioportal, aes(y = Allele_Frequency,
                                                       x = Dataset)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Dataset,
                                                  scales="free") #faceting

ggplot(data = colontotunfiltered_vs_coloncbioportal, aes(y = Allele_Frequency,
                                                       x = Hugo_Symbol)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Hugo_Symbol,
                                                  scales="free") +
  ggtitle("Colon") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  stat_compare_means(aes(group = Dataset), method = "t.test", label = "p.format", vjust = 1)


# total -------------------------------------------------------------------
colon_summarised_graph <- filter(colon_summarised, Hugo_Symbol %in% c("ERBB4","BRAF","KRAS","PIK3CA","FBXW7","TP53", "SMAD4", "PTEN"))

colontot_vs_coloncbioportal <- dplyr::bind_rows(list(cbioportal=colon_cbioportal,
                                                   "Ion-Ampliseq Colon Lung"=colon_summarised_graph),
                                              .id = 'Dataset')

het_plot <- ggplot() +
  geom_density(alpha=0.4, data = colontot_vs_coloncbioportal, aes(x=Allele_Frequency, group= Dataset, fill=Dataset)) +
  ggtitle("Colon") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
het_plot+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

ggplot(data = colontot_vs_coloncbioportal, aes(y = Allele_Frequency,
                                             x = Dataset)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Dataset,
                                                  scales="free") #faceting

ggplot(data = colontot_vs_coloncbioportal, aes(y = Allele_Frequency,
                                             x = Hugo_Symbol)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Hugo_Symbol,
                                                  scales="free") +
  ggtitle("Colon") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  stat_compare_means(aes(group = Dataset), method = "t.test", label = "p.format", vjust = 1)

# density plot ------------------------------------------------------------

colon_het <- ggplot(colon, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
colon_het

unknown <- filter(colon, Hugo_Symbol == "Unknown")
unknown_het <- ggplot(unknown, aes(x=Allele_Frequency)) +
  geom_density()
unknown_het

first3 <- filter(colon, Hugo_Symbol %in%  c("ERBB4","FGFR3","EGFR"))
four_to_6 <- filter(colon, Hugo_Symbol %in%  c("PTEN","KRAS","TP53"))
seven_to_9 <- filter(colon, Hugo_Symbol %in%  c("STK11","PIK3CA","MET"))
ten_to_12 <- filter(colon, Hugo_Symbol %in%  c("DDR2","SMAD4","BRAF"))
thirteen_to_15 <- filter(colon, Hugo_Symbol %in%  c("FBXW7","NRAS","ALK"))
sixteen_to_18 <- filter(colon, Hugo_Symbol %in%  c("MAP2K1","FGFR1","ERBB2"))
eighteen_to_20 <- filter(colon, Hugo_Symbol %in%  c("FGFR2","AKT1"))
last2 <- filter(colon, Hugo_Symbol %in%  c("NOTCH1","CTNNB1"))

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

eighteen_to_20_het <- ggplot(eighteen_to_20, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
eighteen_to_20_het

last_2_het <- ggplot(last2, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
last_2_het
