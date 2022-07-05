rm(list=ls()[! ls() %in% c("maf_filtered","melanoma", "maf_filtered_cellularity")])

library(maftools)
library(dplyr)
library(tidyverse)
library(mclust)
library(ggplot2)

melanoma <- filter(maf_filtered_cellularity, Tumor_Type == "Melanoma", Allele_Frequency >= 0.05)
unique(melanoma$Tumor_Type)
unique(melanoma$Hugo_Symbol)
unique(maf_filtered$Tumor_Type)

melanoma_summarised <- group_by(melanoma, Tumor_Sample_Barcode, Hugo_Symbol, Tumor_Type) %>%
  summarise(Allele_Frequency = mean(Allele_Frequency))%>%
  mutate(Tumor_Sample_Barcode=as.character(Tumor_Sample_Barcode))

# cbioportal --------------------------------------------------------------
#c-bioportal calculations
setwd("C:/Users/matte/Desktop/Uni/tesi/cbioportal/melanoma/")
library(tidyverse)
melanoma_cbioportal_ERBB4 <- read_tsv("ERB44_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_MET <- read_tsv("MET_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_EGFR <- read_tsv("EGFR_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_FGFR2 <- read_tsv("FGFR2_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_NRAS <- read_tsv("NRAS_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_PTEN <- read_tsv("PTEN_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_BRAF <- read_tsv("BRAF_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_TP53 <- read_tsv("TP53_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_ALK <- read_tsv("ALK_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal_CTNNB1 <- read_tsv("CTNNB1_cbioportal_melanoma.tsv", col_names = T) %>%
  select( `Allele Freq (T)`,`Sample ID`) %>%
  filter(`Allele Freq (T)` >= '0.05')%>%
  rename(Allele_Frequency = `Allele Freq (T)`, Tumor_Sample_Barcode = `Sample ID`)

melanoma_cbioportal <- dplyr::bind_rows(list(ERBB4=melanoma_cbioportal_ERBB4,
                                         BRAF=melanoma_cbioportal_BRAF,
                                         EGFR=melanoma_cbioportal_EGFR,
                                         FGFR2=melanoma_cbioportal_FGFR2,
                                         NRAS=melanoma_cbioportal_NRAS,
                                         PTEN=melanoma_cbioportal_PTEN,
                                         ALK=melanoma_cbioportal_ALK,
                                         TP53=melanoma_cbioportal_TP53,
                                         CTNNB1=melanoma_cbioportal_CTNNB1,
                                         MET=melanoma_cbioportal_MET),
                                    .id = 'Hugo_Symbol')

# density plot ------------------------------------------------------------

melanoma_het <- ggplot(melanoma, aes(x=Allele_Frequency, color=Hugo_Symbol)) +
  geom_density()
melanoma_het

first3 <- filter(melanoma, Hugo_Symbol %in%  c("ERBB4","FGFR3","EGFR"))
four_to_6 <- filter(melanoma, Hugo_Symbol %in%  c( "BRAF","PTEN","NRAS"))
seven_to_9 <- filter(melanoma, Hugo_Symbol %in%  c("DDR2","ALK","CTNNB1"))
ten_to_12 <- filter(melanoma, Hugo_Symbol %in%  c("PIK3CA","FBXW7","MET"))
thirteen_to_15 <- filter(melanoma, Hugo_Symbol %in%  c("FGFR1","FGFR2","AKT1"))
sixteen_to_18 <- filter(melanoma, Hugo_Symbol %in%  c("TP53","ERBB2","SMAD4"))
eighteen_to_20 <- filter(melanoma, Hugo_Symbol %in%  c("STK11","MAP2K1"))
last2 <- filter(melanoma, Hugo_Symbol %in%  c("KRAS","NOTCH1"))

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
# unfiltered --------------------------------------------------------------
melanoma_unfiltered <- filter(maf_filtered, Tumor_Type == "Melanoma")
melanoma_summarised_unfiltered <- group_by(melanoma_unfiltered, Tumor_Sample_Barcode, Hugo_Symbol, Tumor_Type) %>%
  summarise(Allele_Frequency = mean(Old_Allele_Frequency))%>%
  mutate(Tumor_Sample_Barcode=as.character(Tumor_Sample_Barcode))

melanoma_summarised_unfiltered_graph <- filter(melanoma_summarised_unfiltered, Hugo_Symbol %in% c("ERBB4","BRAF","EGFR","FGFR2","NRAS","ALK","MET","TP53","PTEN","CTNNB1"))

melanomatotunfiltered_vs_melanomacbioportal <- dplyr::bind_rows(list(cbioportal=melanoma_cbioportal,
                                                                     "Ion-Ampliseq Colon Lung"=melanoma_summarised_unfiltered_graph),
                                                        .id = 'Dataset')

het_plot <- ggplot() +
  geom_density(alpha=0.4, data = melanomatotunfiltered_vs_melanomacbioportal, aes(x=Allele_Frequency, group= Dataset, fill=Dataset)) +
  ggtitle("Melanoma") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
het_plot+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

ggplot(data = melanomatotunfiltered_vs_melanomacbioportal, aes(y = Allele_Frequency,
                                                       x = Dataset)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Dataset,
                                                  scales="free") #faceting

ggplot(data = melanomatotunfiltered_vs_melanomacbioportal, aes(y = Allele_Frequency,
                                                       x = Hugo_Symbol)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Hugo_Symbol,
                                                  scales="free") +
  ggtitle("Melanoma") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  stat_compare_means(aes(group = Dataset), method = "t.test", label = "p.format", vjust = 1)


# total -------------------------------------------------------------------
melanoma_summarised_graph <- filter(melanoma_summarised, Hugo_Symbol %in% c("ERBB4","BRAF","EGFR","FGFR2","NRAS","ALK","MET","TP53","PTEN","CTNNB1"))

melanomatot_vs_melanomacbioportal <- dplyr::bind_rows(list(cbioportal=melanoma_cbioportal,
                                                           "Ion-Ampliseq Colon Lung"=melanoma_summarised_graph),
                                              .id = 'Dataset')

het_plot <- ggplot() +
  geom_density(alpha=0.4, data = melanomatot_vs_melanomacbioportal, aes(x=Allele_Frequency, group= Dataset, fill=Dataset)) +
  ggtitle("Melanoma") +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")
het_plot+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

ggplot(data = melanomatot_vs_melanomacbioportal, aes(y = Allele_Frequency,
                                             x = Dataset)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Dataset,
                                                  scales="free") #faceting

ggplot(data = melanomatot_vs_melanomacbioportal, aes(y = Allele_Frequency,
                                             x = Hugo_Symbol)) + 
  geom_boxplot(aes(color = Dataset)) + facet_wrap(~Hugo_Symbol,
                                                  scales="free") +
  ggtitle("Melanoma") + 
  theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9"))+
  stat_compare_means(aes(group = Dataset), method = "t.test", label = "p.format", vjust = 1)



