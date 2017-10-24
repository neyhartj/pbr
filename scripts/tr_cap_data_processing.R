## Process example data from T3
##
## This is data on the two-row CAP
##
library(tidyverse)
library(stringr)
library(purrrlyr)
library(lme4)

# Read in line information
line_data <- read_csv("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/BOPA/2R_CAP/2R_CAP_line_details.csv")


# Read in geno information
genos <- read_tsv("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/BOPA/2R_CAP/2R_CAP_T3_download/genotype.hmp.txt")

# Parse the line names
line_name <- genos %>%
  select(-rs:-pos) %>%
  names()

## Phenotype data
phenos <- read_tsv(file = "other_data/2R_cap_trait_data/traits.txt")

# Parse the line names and trait names
phenos1 <- phenos %>%
  mutate(line = str_replace_all(line, "'", "")) %>%
  rename_at(vars(`grain yield`:`kernel weight`), str_replace_all, " ", "")

# Get only the lines of interest
phenos2 <- phenos1 %>%
  filter(line %in% line_name)

## Find the proportion of missingness per trial/trait
obs_per_trial <- phenos2 %>%
  gather(trait, value, -line, -trial) %>%
  complete(line, trial, trait) %>%
  group_by(trial, trait) %>%
  summarize(prop_na = mean(is.na(value))) %>%
  arrange(prop_na)

# Pick the top 6 dryland environment-traits
selected_trials <- obs_per_trial %>%
  filter(str_detect(trial, "Dry"),
         trait %in% c("grainyield", "plantheight")) %>%
  group_by(trait) %>%
  top_n(n = -6, prop_na)

# Filter for dryland CAP trials
tr_cap_phenos_met <- selected_trials %>%
  select(trial, trait) %>%
  left_join(., gather(phenos2, trait, value, -line, -trial)) %>%
  ungroup()


# Calculate means per line per trait
tr_cap_phenos <- phenos2 %>%
  group_by(line) %>%
  summarize_at(vars(grainyield:kernelweight), mean, na.rm = TRUE)

# Extract the genotypes for those lines with phenotype data
genos_anno <- genos %>%
  select(marker = rs, alleles:pos, which(names(.) %in% tr_cap_phenos$line))


# Merge the calls with the metadata
tr_cap_genos <- genos_anno %>%
  # Filter unmapped
  filter(chrom != "UNK") %>%
  mutate_at(vars(-marker:-pos), as.numeric)

# Convert to a matrix
genos1 <- t(subset(tr_cap_genos, select = -1:-4))

tr_cap_genos_mat <- genos1 %>%
  structure(dimnames = list(names(tr_cap_genos)[-c(1:4)], tr_cap_genos$marker))

# Filter
# All missing
to_keep <- colMeans(is.na(tr_cap_genos_mat)) < 1

mat1 <- tr_cap_genos_mat[,to_keep]

to_keep <- rowMeans(is.na(tr_cap_genos_mat)) < 1

mat2 <- mat1[to_keep,]

## MAF
af <- colMeans(mat2 + 1, na.rm = TRUE) / 2
maf <- pmin(af, 1-af)

to_keep <- maf > 0.01

mat3 <- mat2[,to_keep]

## Missingness
to_keep <- colMeans(is.na(mat3)) <= 0.10

mat4 <- mat3[,to_keep]

to_keep <- rowMeans(is.na(mat4)) <= 0.10

mat5 <- mat4[to_keep,]

## Extract marker metadata
tr_cap_genos_hmp <- tr_cap_genos %>%
  filter(marker %in% colnames(mat5)) %>%
  select(marker:pos, row.names(mat5)) %>%
  mutate(pos = pos / 1000)

tr_cap_genos_mat <- mat5[,colnames(mat5) %in% tr_cap_genos_hmp$marker]




# Edit the phenos
tr_cap_phenos <- tr_cap_phenos %>%
  filter(line %in% row.names(tr_cap_genos_mat))

tr_cap_phenos_met <- tr_cap_phenos_met %>%
  filter(line %in% row.names(tr_cap_genos_mat))

## Save
devtools::use_data(tr_cap_genos_hmp, tr_cap_genos_mat, tr_cap_phenos, tr_cap_phenos_met,
                   overwrite = T)




