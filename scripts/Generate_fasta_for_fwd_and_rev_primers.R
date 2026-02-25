library(tidyverse)

# Example primer definitions
primer_df_fwd <- data.frame(
  name = c("HBB_forward_primer", "CFTR_forward_primer"),
  sequence = c(
    "GTGCACCTGACTCCTGAGGAGA",   # HBB example forward primer
    "TTGGCATGCTTTGATGACGCTT"    # CFTR example forward primer
  )
)

primer_df_rev <- data.frame(
  name = c("HBB_reverse_primer", "CFTR_reverse_primer"),
  sequence = c(
    "CCTTGATACCAACCTGCCCAG",    # HBB example reverse primer
    "AGTGGAGTTCTGAGCCCAGC"      # CFTR example reverse primer
  )
)

# Convert to FASTA format
fwd_fasta <- as.vector(rbind(
  paste0(">", primer_df_fwd$name),
  primer_df_fwd$sequence
))

rev_fasta <- as.vector(rbind(
  paste0(">", primer_df_rev$name),
  primer_df_rev$sequence
))

# Write files
writeLines(fwd_fasta, "ref/forward_primers.fasta")
writeLines(rev_fasta, "ref/reverse_primers.fasta")