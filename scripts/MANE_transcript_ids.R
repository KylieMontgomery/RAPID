library(rtracklayer)
library(dplyr)

mane <- readGFF("MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz") %>%
  filter(type == "transcript", tag == "MANE_Select") %>%
  select(transcript_id, gene_name)

write.table(
  mane,
  "MANE_select_transcripts.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)