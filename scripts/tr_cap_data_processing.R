## Process example data from T3
##
## This is data on the two-row CAP
##
library(tidyverse)
library(stringr)
library(purrrlyr)
library(lme4)

# Read in line information
line_data <- read_csv("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Genotypic Data/BOPA Markers/2R_CAP/2R_CAP_line_details.csv")


# Read in geno information
genos <- read_tsv("C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Genotypic Data/BOPA Markers/2R_CAP/2R_CAP_T3_download/genotype.hmp.txt")

# Parse the line names
line_name <- genos %>%
  select(-`rs#`:-QCcode) %>% names()

## Phenotype data
phenos <- read_tsv(file = "data/download_NDZK/download_NDZK/traits.txt")

# Parse the line names and trait names
phenos1 <- phenos %>%
  mutate(line = str_replace_all(line, "'", "")) %>%
  rename_at(vars(`grain yield`:`kernel weight`), str_replace_all, " ", "")

# Get only the lines of interest
phenos2 <- phenos1 %>%
  filter(line %in% line_name)

# Calculate means per line per trait
tr_cap_phenos <- phenos2 %>%
  group_by(line) %>%
  summarize_at(vars(grainyield:kernelweight), mean, na.rm = TRUE)

# Extract the genotypes for those lines with phenotype data
genos_anno <- genos %>%
  select(marker = `rs#`, alleles:pos, which(names(.) %in% tr_cap_phenos$line))

genos1 <- genos_anno%>%
  apply(X = ., MARGIN = 1, FUN = function(i) {

    alleles <- i[2] %>%
      str_split(pattern = "/") %>%
      unlist()

    convert = c("-1" = str_c(alleles[1], alleles[1]),
                "0" = str_c(alleles[1], alleles[2]),
                "1" = str_c(alleles[2], alleles[2])) %>%
      structure(as.character(names(.)), names = .)

    calls <- unlist(i[-c(1:4)])

    str_replace_all(calls, convert) %>% str_replace_all("NN", "NA") %>% parse_number()

  })


# Merge the calls with the metadata
tr_cap_genos <- bind_cols(select(genos_anno, marker:pos), as.data.frame(t(genos1))) %>%
  structure(names = names(genos_anno)) %>%
  # Filter unmapped
  filter(chrom != "UNK")

# Convert to a matrix
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
tr_cap_genos_map <- tr_cap_genos %>%
  select(1:4) %>%
  filter(marker %in% colnames(mat5)) %>%
  mutate(pos = pos / 1000)

tr_cap_genos <- mat5[,colnames(mat5) %in% tr_cap_genos_map$marker]


# Edit the phenos
tr_cap_phenos <- tr_cap_phenos %>%
  filter(line %in% row.names(tr_cap_genos))






