#

# Import libraries####
library(tidyverse)
library(DRomics)
library(readxl)
library(writexl)
library(tidyverse)
library(janitor)
library(wesanderson)
library(cowplot)
library(mclust)

source("util.R")

#Set params####
path <- "./2_imports/"
tpm_filename <- "ring_tpm.xlsx"
rc_filename <- "rc.xlsx"

#Import data####
ring_rawtpm <- read_excel(paste(path, tpm_filename, sep = "")) #Change per stage as needed
tpm00 <- ring_rawtpm #Change this part during operation per stage as needed
rc_00 <- read_excel(rc_filename)

#Filter low counts####
tpm01 <- tpm00 %>%
  rowwise(geneID) %>%
  summarise(
    min.one = sum(across(-1, ~. >= 1)),
    zeros = sum(across(-1, ~. == 0))
  )

p1 <- tpm01 %>% 
ggplot(aes(x = min.one)) + 
  geom_histogram(bins = 10) +
  labs(
    title = "Detected genes",
    x = "Samples with at least 1 TPM",
    y = "Gene count") +
  theme_cowplot() 
p1

p2 <- tpm01 %>% 
ggplot(aes(x = zeros)) + 
  geom_histogram(bins = 10) +
  labs(
    title = "Detected genes",
    x = "Samples with 0 TPM",
    y = "Gene count") +
  theme_cowplot() 
p2

tpm02 <- tpm01 %>% filter(min.one > 40)

p3 <- tpm02 %>%
ggplot(aes(x = min.one)) +
geom_histogram(bins = 5) +
labs(
  title = "Detected genes post-filtering",
  x = "Samples with at least 1 TPM",
  y = "Gene count"
) +
theme_cowplot()
p3

tpm03 <- tpm00 %>%
  filter(geneID %in% tpm02$geneID) %>%
  column_to_rownames(var = "geneID") %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column(var = "sample") %>%
  mutate(sample_grouped = rep(c(1:19), each = 3)) %>%
  group_by(sample_grouped) %>%
  mutate_at(vars(-c("sample", "sample_grouped")), mean) %>%
  t() %>% 
  as.data.frame() %>%
  row_to_names(row_number = 1) %>%
  select(ends_with('1')) %>%
  mutate(across(everything(), as.numeric)) %>%
  mutate(across(everything(), ~ log2(. + 1)))

tpm04 <- tpm03 %>%
  rownames_to_column(var = "geneID") %>%
  pivot_longer(cols = -1, names_to = "sample", values_to = "logtpm") %>%
  mutate("tp" = substr(sample, 2,2)) %>%
  mutate(tp = as.factor(tp)) %>%
  mutate("dose" = substr(sample, 4, 4)) %>%
  mutate(dose = as.factor(dose))

p4 <- ggplot(tpm04, aes(x = sample, y = logtpm, fill = tp)) +
  geom_violin() +
  scale_fill_manual(values = wes_palette("Royal2")) +
  ylim(0,18) +
  stat_summary(
    fun = "mean",
    geom = "crossbar",
    width = 0.5,
    color = "black"
  ) +
  labs(
    title = "Distribution of Mean TPMs",
    x = "Sample",
    y = "log2(meanTPM + 1)",
    fill = "Timepoint"
  ) +
  scale_x_discrete(guide = guide_axis(angle = 60)) +
  theme_cowplot()
p4

#Dose-dependent genes ####
rc <- read_excel(paste(path, rc_filename, sep = ""))
dro.in01 <- rc %>%
  filter(geneID %in% tpm02$geneID) %>%
  column_to_rownames(var = "geneID") %>%
  mutate(across(everything(), as.numeric)) %>%
  rownames_to_column(var = "geneID")

names <- c("geneID", 150, 150, 150, 50, 50, 50, 16.7, 16.7, 16.7, 5.6, 5.6, 5.6, 1.9, 1.9, 1.9, 0, 0, 0)

dro.in02 <- rbind(names, dror.in01)
dro.in02 <- dro.in02 %>%
  mutate_at(-1, as.numeric)
colnames(dro.in02) <- names

dro.in03 <- RNAseqdata(as.data.frame(dro.in02))
dro.out01 <- itemselect(dro.in03, select.method = "quadratic", FDR = 0.05)
dro.out01
dro.out02 <- drcfit(dro.out01, progressbar = TRUE)
dro.out03 <- bmdcalc(dro.out02)
dro.out03

#Save results####
dror6.02 <- dro.out02     #CHANGE THE OBJECT NAMES
dror6.03 <- dro.out03$res #CHANGE THE OBJECT NAMES
plotfit2pdf(dro.out02, plot.type = "dose_fitted", 
            BMDoutput = dro.out03, path2figs = "./3_exports") 

#Explore results####
bmd.hist <- ggplot(data = dror4.03, aes(x = BMD.zSD)) +
  geom_histogram(bins = 20, color = "#446455", fill = "#d3dddc") +
  scale_x_continuous(breaks = c(seq(0, max(dros6.03$BMD.zSD) + 10, 10)))+
  labs(
    title = "Ring 4h BMDs", #Adjust accordingly
    y = "n Genes",
    x = "BMD (nM DHA)") +
  theme_cowplot()
bmd.hist

#Gaussian Mixture Modeling
gmm_model <- Mclust(drot1.03$BMD.zSD, G = 2) #Adjust (before $BMD.zSD portion) accordingly
clusters <- predict(gmm_model)$classification
plot(drot1.03$BMD.zSD, col = clusters, main = "Troph 1h GMM Clustering") #Adjust accordingly
points(gmm_model$parameters$mean, col = 1:3, pch = 8, cex = 2)
drot1.03 <- cbind(drot1.03, "clust" = clusters) #Adjust accordingly

#For unclustered sets
bmd_cut <- 25
drot4.03 <- drot4.03 %>%
  mutate(clust = case_when(
    BMD.zSD <= 25 ~ 1,
    BMD.zSD > 25 ~ 2
  ))
drot6.03 <- drot6.03 %>%
  mutate(clust = case_when(
    BMD.zSD <= 25 ~ 1,
    BMD.zSD > 25 ~ 2
  ))

#Compile results
sets <- c("dror1.03", "dror4.03", "dror6.03", "drot1.03", "drot4.03", "drot6.03", "dros1.03", "dros4.03", "dros6.03")

dror1.03 <- dror1.03 %>% mutate("set" = "ring 1h")
dror4.03 <- dror4.03 %>% mutate("set" = "ring 4h")
dror6.03 <- dror6.03 %>% mutate("set" = "ring 6h")
drot1.03 <- drot1.03 %>% mutate("set" = "troph 1h")
drot4.03 <- drot4.03 %>% mutate("set" = "troph 4h")
drot6.03 <- drot6.03 %>% mutate("set" = "troph 6h")
dros1.03 <- dros1.03 %>% mutate("set" = "schiz 1h")
dros4.03 <- dros4.03 %>% mutate("set" = "schiz 4h")
dros6.03 <- dros6.03 %>% mutate("set" = "schiz 6h")
deg.full0 <- rbind(
  dror1.03, dror4.03, dror6.03,
  drot1.03, drot4.03, drot6.03, 
  dros1.03, dros4.03, dros6.03)
deg.full1 <- deg.full0 %>%
  mutate_at("set", as.factor)
  mutate("shift" = case_when(
    yatdosemax > y0 ~ "up",
    yatdosemax < y0 ~ "down"
  )) %>%
  group_by(set)
# vv Optional cleanup vv
# rm(dror1.03, dror4.03, dror6.03, drot1.03, drot4.03, drot6.03, dros1.03, dros4.03, dros6.03)

#Bar chart summary of typology
q1 <- deg.full1 %>% 
  filter(set %in% c("schiz 1h", "schiz 4h", "schiz 6h")) %>% #Adjust accordingly
  ggplot(aes(x = shift, fill = model)) +
    geom_bar(position = "stack") +
    facet_grid(~set, switch = "x") +
    scale_fill_manual(values = wes_palette(name = "Royal2", n = 5)) +
    labs(y = "n DEGs", x = "", fill = "Model", title = "Schiz DEG summary") +
    theme_cowplot(font_size = 16) +
    theme(strip.placement = "outside",
          strip.background = element_rect(fill = NA, color = "white"),
          panel.spacing = unit(-0.1, "cm"))
q1

write_xlsx(deg.full1, "./3_exports/deg.full.xlsx") 
