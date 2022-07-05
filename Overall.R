rm(list=ls()[! ls() %in% c("maf_filtered", "maf_filtered_nosnps", "dbsnp_info", "maf_filtered_cellularity")])

library(maftools)
library(dplyr)
library(tidyverse)
library(mclust)


# ERBB4 -------------------------------------------------------------------

ERBB4 <- filter(maf_filtered, Hugo_Symbol == "ERBB4", Allele_Frequency >= 0.05)

mean_per_sample <- select(ERBB4, Hugo_Symbol, Tumor_Sample_Barcode, Allele_Frequency) %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise_each(funs(mean)) 






# rsnps -------------------------------------------------------------------

library(rsnps)
library(dslabs)
library(purrr)
library(tidyverse)

maf_filtered$dbSNP_RS
snps_table <- filter(maf_filtered, !dbSNP_RS %in% c("","novel") & Old_Allele_Frequency != "NaN") %>%
  select(dbSNP_RS, Tumor_Sample_Barcode, Hugo_Symbol, Old_Allele_Frequency)

snps <- unique(snps_table$dbSNP_RS) 
(dbsnp_info <- ncbi_snp_query(snps))

snps_maf_mean <- map(.x=dbsnp_info$maf_population, .f = ~mean(.x$MAF, na.rm = T))


dbsnp_info <- dbsnp_info %>%
  mutate(snps_maf_mean = map(.x=dbsnp_info$maf_population, .f = ~mean(.x$MAF, na.rm = T))) %>%
  filter(snps_maf_mean != "NaN")

filtered_snps <- dbsnp_info %>%
  filter(snps_maf_mean >= 0.02)

maf_filtered_0.2_nosnps <- maf_filtered %>%
  filter(!dbSNP_RS %in% filtered_snps$query)

maf_filtered_nosnps <- maf_filtered %>%
  filter(!dbSNP_RS %in% dbsnp_info$query)%>%
  filter(!dbSNP_RS %in% c("", "novel") & !Old_Allele_Frequency %in% (0.45:0.55)| !Allele_Frequency > 0.95)

# Cellularity -------------------------------------------------------------
maf_filtered <- maf_filtered %>%
  filter(Old_Allele_Frequency >= "0.01", Old_Allele_Frequency != "NaN")
maf_no_cellularity <- maf_filtered_nosnps %>%
  group_by(Hugo_Symbol,Tumor_Sample_Barcode, Tumor_Type) %>%
  summarise(Old_Allele_Frequency = mean(Old_Allele_Frequency)) %>%
  filter(Old_Allele_Frequency != "NaN")


vediamo <- filter(maf_filtered_cellularity, Allele_Frequency < 0.05)


het_plot <- ggplot() +
  geom_density(alpha=0.4, data = maf_filtered_nosnps, aes(x=Old_Allele_Frequency, group= Tumor_Sample_Barcode, color=Tumor_Sample_Barcode, fill=Tumor_Sample_Barcode))
het_plot     


library(ggpmisc)
data <- select(maf_filtered_nosnps, Tumor_Sample_Barcode, Old_Allele_Frequency, Hugo_Symbol, Tumor_Type)

maxi <- data %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(max = max(Old_Allele_Frequency, na.rm=TRUE)) 

# data_grouped <- data %>%
#   group_by(Tumor_Sample_Barcode)
# 
# x = 0 : 1
# y = data_grouped$Old_Allele_Frequency
# x[ggpmisc:::find_peaks(y)]
# 
# d <- density(data_grouped$Old_Allele_Frequency)
# peakx <- d$x[d$y==max(d$y)]

maf_filtered_nosnps <- select(maf_filtered_nosnps, - Allele_Frequency)

maf_filtered_cellularity <-  merge(x=maxi, y=maf_filtered_nosnps, by="Tumor_Sample_Barcode")%>%
  mutate(Allele_Frequency = Old_Allele_Frequency / (max+0.01))
  
maf_filtered_cellularity <- maf_filtered_cellularity %>%
  filter(!(Allele_Frequency > 0.45 & Allele_Frequency < 0.55))%>%
  filter(Allele_Frequency < 0.95)

prova <- maf_filtered$Tumor_Type[maf_filtered$Tumor_Type == "Breast"] <-  "Melanoma"

library(GenVisR)

for_lolli <- maf_filtered_cellularity %>%
  filter(VARIANT_CLASS == "SNV")%>%
  rename(amino_acid_change = HGVSp_Short) %>%
  rename(gene = SYMBOL) %>%
  rename(transcript_name = Transcript_ID)

NRAS <- for_lolli %>%
  select(amino_acid_change, gene, transcript_name)%>%
  filter(gene == "NRAS")
library (randomcoloR)
palette_NRAS = randomColor(4)

lolliplot(NRAS, labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

TP53 <- for_lolli %>%
  filter(gene == "TP53", amino_acid_change != "", Tumor_Type == "Colon")%>%
  select(amino_acid_change, gene, transcript_name)

lolliplot(TP53,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)


PTEN <- for_lolli %>%
  select(amino_acid_change, gene, transcript_name)%>%
  filter(gene == "PTEN", amino_acid_change != "")
lolliplot(PTEN,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)


EGFR <- for_lolli %>%
  filter(gene == "EGFR", amino_acid_change != "", Tumor_Type == "Lung")%>%
  select(amino_acid_change, gene, transcript_name)
  
lolliplot(EGFR,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

KRAS <- for_lolli %>%
  filter(gene == "KRAS", amino_acid_change != "", Tumor_Type == "Lung")%>%
  select(amino_acid_change, gene, transcript_name)
lolliplot(KRAS,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

NRAS <- for_lolli %>%
  filter(gene == "NRAS", amino_acid_change != "", Tumor_Type == "Melanoma")%>%
  select(amino_acid_change, gene, transcript_name)
lolliplot(NRAS,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

KRAS_colon <- for_lolli %>%
  filter(gene == "KRAS", amino_acid_change != "", Tumor_Type == "Colon")%>%
  select(amino_acid_change, gene, transcript_name)
lolliplot(KRAS_colon,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

data <- brcaMAF[brcaMAF$Hugo_Symbol == 'TP53',c('Hugo_Symbol', 'amino_acid_change_WU')]
data <- as.data.frame(cbind(data, 'ENST00000269305'))
colnames(data) <- c('gene', 'amino_acid_change', 'transcript_name')

BRAF_colon <- for_lolli %>%
  filter(gene == "BRAF", amino_acid_change != "", Tumor_Type == "Colon")%>%
  select(amino_acid_change, gene, transcript_name)

lolliplot(BRAF_colon,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

BRAF_melanoma <- for_lolli %>%
  filter(gene == "BRAF", amino_acid_change != "", Tumor_Type == "Melanoma")%>%
  select(amino_acid_change, gene, transcript_name)

lolliplot(BRAF_melanoma,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

BRAF_lung <- for_lolli %>%
  filter(gene == "BRAF", amino_acid_change != "", Tumor_Type == "Lung")%>%
  select(amino_acid_change, gene, transcript_name)

lolliplot(BRAF_lung,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

CTNNB1 <- for_lolli %>%
  filter(gene == "CTNNB1", amino_acid_change != "", Tumor_Type == "Melanoma")%>%
  select(amino_acid_change, gene, transcript_name)
lolliplot(CTNNB1,labelCol = "amino_acid_change", txtSize = 2, pntSize = 2)

# Call lolliplot
lolliplot(data) 

table(for_lolli$VARIANT_CLASS)

table(maf_filtered_cellularity$VARIANT_CLASS)
table(maf_filtered$VARIANT_CLASS)

oncoplot <- maf_filtered_cellularity %>%
  filter(!Hugo_Symbol %in% c("TP53", "PTEN"))
oncomaf <- read.maf(maf = oncoplot)
oncoplot(maf=oncomaf, top = 10, pwLineWd = 0.0000001, borderCol = NA, draw_titv = TRUE, fill = F, bgCol = "grey90")

oncoplot2 <- maf_filtered_cellularity %>%
  filter(Variant_Type == "SNP", Tumor_Type != "Tyroid")
oncomaf2 <- read.maf(oncoplot2)
oncoplot(maf=oncomaf2, top = 15, pwLineWd = 0.0000001, borderCol = NA, draw_titv = TRUE, fill = F, bgCol = "grey90")
