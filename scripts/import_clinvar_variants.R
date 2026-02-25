library(readr)
library(dplyr)
library(tidyr)
library(stringr)

clinvar_df <- read_tsv(
  "ref/clinvar_raw.tsv",
  col_names = c("chrom","pos","ref","alt","clnsig","geneinfo"),
  show_col_types = FALSE
) %>%
  separate(geneinfo, into = c("gene", "geneid"), sep = ":", fill = "right", extra = "merge") %>%
  mutate(clnsig = str_replace_all(clnsig, "_", " ")) %>%
  filter(clnsig %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")) %>%
  select(chrom, pos, ref, alt, clnsig, gene)

write_tsv(clinvar_df, "ref/clinvar_variants.tsv")