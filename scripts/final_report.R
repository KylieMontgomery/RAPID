#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggtranscript)
  library(cowplot)
  library(rtracklayer)
  library(ggtext)
  library(base64enc)
})

options(stringsAsFactors = FALSE)

# -------------------- CLI --------------------
option_list <- list(
  make_option(c("--case_name"),      type = "character"),
  make_option(c("--case_sample"),    type = "character"),
  make_option(c("--case_control"),   type = "character"),
  make_option(c("--chrom"),          type = "character"),
  make_option(c("--start_pos"),      type = "numeric",   default = 0),
  make_option(c("--end_pos"),        type = "numeric",   default = 0),
  make_option(c("--gene"),           type = "character"),
  make_option(c("--mane_id"),        type = "character", default = NULL),
  make_option(c("--variant_name"),   type = "character", default = ""),
  make_option(c("--variant_pos"),    type = "character", default = ""),
  make_option(c("--merged_gtf"),     type = "character"),
  make_option(c("--output_html"),    type = "character", default = "report.html"),
  make_option(c("--MANE_transcripts"), type = "character"),
  make_option(c("--clinvar_variants"), type = "character"),
  make_option(c("--prop_threshold"), type = "double", default = 0.10,
              help = "Coverage proportion threshold, e.g. 0.10 or 0.05"),
  make_option(c("--gene_id"),        type = "character", default = "",
              help = "Optional gene_id (e.g. MSTRG.78) to force plotting")
) 
opt <- parse_args(OptionParser(option_list = option_list))

thresh <- if (!is.null(opt$prop_threshold)) as.numeric(opt$prop_threshold) else 0.05
target_gene_id <- if (!is.null(opt$gene_id) && nzchar(opt$gene_id)) opt$gene_id else NA_character_

# Parse comma-separated variant name/pos
parse_csv <- function(x) if (is.null(x) || x == "") character(0) else trimws(strsplit(x, ",")[[1]])
vnames <- parse_csv(opt$variant_name)
vpos   <- suppressWarnings(as.numeric(parse_csv(opt$variant_pos)))
if (length(vnames) == length(vpos) && length(vpos) > 0) {
  variant_labels <- tibble(pos = vpos, label = vnames, transcript_id = "Case Variants")
} else { 
  variant_labels <- tibble(pos = numeric(0), label = character(0), transcript_id = character(0))
} 

# -------------------- Paths --------------------
output_dir <- dirname(opt$output_html)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
transcripts_detected_plot_file <- file.path(output_dir, paste0(opt$case_name, "_transcripts_detected_plot.png"))
proportion_plot_file           <- file.path(output_dir, paste0(opt$case_name, "_proportion_plot.png"))
transcript_structure_plot_file <- file.path(output_dir, paste0(opt$case_name, "_transcript_structure_plot.png"))
additional_info_file           <- file.path(output_dir, paste0(opt$case_name, "_additional_info.csv"))

# -------------------- Load data --------------------
gtf_data <- rtracklayer::import(opt$merged_gtf)
df <- as.data.frame(gtf_data) %>%
  mutate(
    cov = as.numeric(cov),
    seqnames = as.character(seqnames),
    strand   = as.character(strand),
    source   = as.character(source),
    type     = as.character(type)
  ) 

# -------------------- Choose gene_id group --------------------
if (!is.na(target_gene_id)) {
  gene_group_ids <- target_gene_id
} else { 
  gene_group_ids <- df %>%
    filter(seqnames == opt$chrom,
           between(start, opt$start_pos, opt$end_pos),
           type %in% c("transcript","exon")) %>%
    filter(ref_gene_name == opt$gene | gene_name == opt$gene) %>%
    distinct(gene_id) %>% pull()
  if (length(gene_group_ids) == 0) {
    gene_group_ids <- df %>%
      filter(seqnames == opt$chrom,
             between(start, opt$start_pos, opt$end_pos),
             type %in% c("transcript","exon")) %>%
      group_by(gene_id) %>% summarise(tcov = sum(as.numeric(cov), na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(tcov)) %>% slice_head(n = 1) %>% pull(gene_id)
  }
} 
df <- df %>% filter(gene_id %in% gene_group_ids)

# -------------------- MANE + ClinVar --------------------
MANE_transcripts <- read.csv(opt$MANE_transcripts)
clinvar_df <- read.delim(opt$clinvar_variants, header = FALSE, sep = "\t",
                         col.names = c("chrom", "pos", "ref", "alt", "clnsig", "gene", "geneid")) %>%
  filter(clnsig %in% c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"))

clinvar_variants_case <- clinvar_df %>%
  filter(gene == opt$gene) %>%
  mutate(clnsig = factor(clnsig, levels = c("Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic")),
         pos = as.numeric(pos))

# -------------------- Build transcript table --------------------
transcripts <- df %>%
  filter(seqnames == opt$chrom,
         between(start, opt$start_pos, opt$end_pos),
         type == "transcript",
         sample_id %in% c(opt$case_sample, opt$case_control)) %>%
  mutate(sample_id = case_when(
    sample_id == opt$case_sample   ~ "Case",
    sample_id == opt$case_control  ~ "Control",
    TRUE ~ sample_id
  )) %>% 
  left_join(MANE_transcripts, by = "transcript_id")

# -------------------- 1) Detected-over-threshold curves --------------------
transcripts_detected <- transcripts %>%
  group_by(sample_id) %>%
  group_modify(~ {
    total_cov <- sum(.x$cov, na.rm = TRUE)
    if (total_cov <= 0 || nrow(.x) == 0) {
      return(tibble(threshold = double(0), transcripts = integer(0), sample_id2 = .y$sample_id))
    } 
    .x <- mutate(.x, NFLRT = cov / total_cov)
    thresholds <- seq(0, max(.x$NFLRT, na.rm = TRUE), length.out = 100)
    tibble(threshold = thresholds,
           transcripts = sapply(thresholds, function(t) sum(.x$NFLRT >= t)),
           sample_id2 = .y$sample_id)
  }) %>% ungroup() 

plot_curve <- function(transcripts_detected, main_title, vline = NULL) {
  p <- ggplot(transcripts_detected, aes(x = threshold, y = transcripts, color = sample_id2)) +
    geom_line(size = 1) +
    labs(title = main_title, x = "Proportion threshold", y = "Number of transcripts detected", color = NULL) +
    xlim(0, 1) + theme_minimal()
  if (!is.null(vline)) p <- p + geom_vline(xintercept = vline, linetype = 2)
  p
} 
get_legend <- function(plot) cowplot::get_legend(plot + theme(legend.position = "bottom"))

p_for_legend <- plot_curve(transcripts_detected, "Combined", vline = thresh) +
  scale_color_manual(values = c("Case" = "midnightblue", "Control" = "cornflowerblue"))
shared_legend <- get_legend(p_for_legend)
p_case    <- plot_curve(filter(transcripts_detected, sample_id2 == "Case"), "Case") +
  scale_color_manual(values = "midnightblue") + theme(legend.position = "none")
p_control <- plot_curve(filter(transcripts_detected, sample_id2 == "Control"), "Control") +
  scale_color_manual(values = "cornflowerblue") + theme(legend.position = "none")
p_merged  <- p_for_legend + theme(legend.position = "none")

transcripts_detected_plot <- plot_grid(p_case, p_control, p_merged, ncol = 3, align = "v")
transcripts_detected_plot_with_legend <- plot_grid(transcripts_detected_plot, shared_legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave(transcripts_detected_plot_file, plot = transcripts_detected_plot_with_legend, width = 7, height = 4, dpi = 300)

# -------------------- 2) Proportion barplot--------------------
filtered_data_threshold <- transcripts %>%
  mutate(transcript_id = sub("^MSTRG", "Novel", transcript_id)) %>%
  group_by(sample_id) %>%
  mutate(total_cov = sum(cov, na.rm = TRUE),
         prop_cov  = ifelse(total_cov > 0, cov / total_cov, 0)) %>%
  ungroup()

mane_id <- opt$mane_id
filtered_data_threshold <- filtered_data_threshold %>%
  mutate(transcript_id = factor(transcript_id,
                                levels = if (!is.null(mane_id) && mane_id %in% transcript_id)
                                  c(setdiff(unique(transcript_id), mane_id), mane_id)
                                else unique(transcript_id)))

custom_labels <- levels(filtered_data_threshold$transcript_id)
if (!is.null(mane_id) && mane_id %in% custom_labels) {
  custom_labels <- ifelse(custom_labels == mane_id,
                          paste0("<span style='color:#006400;'><b>", custom_labels, "</b></span>"),
                          custom_labels)
} 
names(custom_labels) <- levels(filtered_data_threshold$transcript_id)

# keep a transcript if it reaches the threshold in ANY sample
transcripts_to_keep <- filtered_data_threshold %>%
  group_by(transcript_id) %>%
  filter(any(prop_cov >= thresh, na.rm = TRUE)) %>%
  distinct(transcript_id)

plot_data <- filtered_data_threshold %>%
  semi_join(transcripts_to_keep, by = "transcript_id") %>%
  group_by(transcript_id) %>%
  tidyr::complete(sample_id = c("Case", "Control"), fill = list(prop_cov = 0, total_cov = 0)) %>%
  ungroup() %>%
  mutate(prop_cov = ifelse(prop_cov < 0.05, 0, prop_cov))

# Wide per-transcript proportions (used later for the combined table)
prop_summary <- plot_data %>%
  dplyr::select(transcript_id, sample_id, prop_cov) %>%
  dplyr::mutate(transcript_id = as.character(transcript_id)) %>%
  tidyr::pivot_wider(
    names_from  = sample_id,
    values_from = prop_cov,
    values_fill = 0,
    names_prefix = "prop_"
  ) %>% 
  dplyr::mutate(
    max_prop          = pmax(dplyr::coalesce(prop_Case, 0), dplyr::coalesce(prop_Control, 0)),
    passes_threshold  = max_prop >= thresh,
    is_MANE           = if (!is.null(mane_id)) transcript_id == mane_id else FALSE
  ) %>%
  dplyr::arrange(dplyr::desc(max_prop))

# Plot
proportion_plot <- ggplot(plot_data, aes(x = transcript_id, y = prop_cov, fill = sample_id)) +
  geom_bar(stat = "identity", width = 0.7, position = position_dodge2(width = 0.8, preserve = "single")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, .25, .5, .75, 1)) +
  scale_x_discrete(labels = custom_labels) +
  scale_fill_manual(values = c("Case" = "midnightblue", "Control" = "cornflowerblue")) +
  labs(title = "Transcript Proportions", x = NULL, y = "Proportion of Coverage") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = ggtext::element_markdown(size = 12),
        legend.position = "bottom", legend.title = element_blank())
ggsave(proportion_plot_file, plot = proportion_plot, width = 7, height = 4, dpi = 300)

# -------------------- 3) Structure plot + to_diff + combined CSV --------------------
# valid_transcripts computed from the SAME threshold as above
fdt_for_threshold <- filtered_data_threshold %>%
  mutate(transcript_id = as.character(transcript_id)) %>%
  filter(type == "transcript")

valid_transcripts <- fdt_for_threshold %>%
  group_by(transcript_id) %>%
  filter(any(prop_cov >= thresh, na.rm = TRUE)) %>%
  distinct(transcript_id) %>% pull() %>% as.character()

# ALWAYS add MANE if present in your data, even if it failed the threshold
if (!is.null(mane_id) && mane_id %in% fdt_for_threshold$transcript_id) {
  valid_transcripts <- union(valid_transcripts, mane_id)
}

# Fallback if nothing passes threshold
if (length(valid_transcripts) == 0) {
  if (!is.null(mane_id) && mane_id %in% fdt_for_threshold$transcript_id) {
    valid_transcripts <- c(valid_transcripts, mane_id)
  }
  topN <- fdt_for_threshold %>%
    group_by(transcript_id) %>%
    summarise(tcov = sum(cov, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(tcov)) %>% slice_head(n = min(5, n())) %>% pull(transcript_id)
  valid_transcripts <- unique(c(valid_transcripts, topN))
}

# Build exons
exons <- df %>%
  mutate(transcript_id = sub("^MSTRG", "Novel", transcript_id),
         start = as.numeric(start), end = as.numeric(end), strand = as.character(strand)) %>%
  filter(seqnames == opt$chrom,
         between(start, opt$start_pos, opt$end_pos),
         type == "exon",
         sample_id %in% c(opt$case_sample, opt$case_control),
         transcript_id %in% valid_transcripts) %>%
  distinct(seqnames, start, end, strand, transcript_id, .keep_all = FALSE)

make_blank_plot <- function(title = "No data to plot") {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = title) +
    theme_void() + xlim(0,1) + ylim(0,1)
}

if (nrow(exons) > 0) {
  x_min <- min(exons$start, na.rm = TRUE)
  x_max <- max(exons$end,   na.rm = TRUE)
  
  # MANE split (if present)
  mane <- if (!is.null(mane_id)) filter(exons, transcript_id == mane_id) else exons[0,]
  not_mane <- if (!is.null(mane_id)) filter(exons, transcript_id != mane_id) else exons
  
  # Compute differences
  df_diff <- tryCatch({
    if (!is.null(mane_id) && nrow(mane) > 0 && nrow(not_mane) > 0) {
      to_diff(exons = not_mane, ref_exons = mane, group_var = "transcript_id")
    } else tibble(start = numeric(), end = numeric(), transcript_id = character(), diff_type = character())
  }, error = function(e) tibble(start = numeric(), end = numeric(), transcript_id = character(), diff_type = character()))
  
  # -------------------- Normalize to_diff columns --------------------
  if (!"transcript_id" %in% names(df_diff) && "group" %in% names(df_diff)) {
    df_diff <- dplyr::rename(df_diff, transcript_id = group)
  }
  if ("xstart" %in% names(df_diff) && !"start" %in% names(df_diff)) {
    df_diff <- dplyr::rename(df_diff, start = xstart)
  }
  if ("xend" %in% names(df_diff) && !"end" %in% names(df_diff)) {
    df_diff <- dplyr::rename(df_diff, end = xend)
  }
  
  df_diff <- df_diff %>%
    dplyr::mutate(
      transcript_id = sub("^MSTRG", "Novel", as.character(transcript_id)),
      start = as.numeric(start),
      end   = as.numeric(end)
    ) %>%
    dplyr::filter(!is.na(transcript_id), !is.na(start), !is.na(end))
  
  diff_rows <- df_diff %>%
    dplyr::transmute(
      transcript_id,
      diff_start = start,
      diff_end   = end,
      diff_width = end - start + 1,
      diff_type
    )
  
  # -------------------- Single additional-info table (proportions + to_diff) --------------------
  additional_info <- prop_summary %>%
    dplyr::left_join(diff_rows, by = "transcript_id") %>%
    dplyr::arrange(dplyr::desc(max_prop), transcript_id, diff_start)
  write.csv(additional_info, additional_info_file, row.names = FALSE)
  
  # -------------------- Structure plot --------------------
  all_ids <- unique(exons$transcript_id)
  levels_top <- if (!is.null(mane_id) && mane_id %in% all_ids) c(mane_id, setdiff(all_ids, mane_id)) else all_ids
  include_variant_row <- nrow(variant_labels) > 0
  levels_top_with_variants <- if (include_variant_row) c("Case Variants", levels_top) else levels_top
  levels_for_scale <- rev(levels_top_with_variants)
  lbls <- setNames(levels_top_with_variants, levels_top_with_variants)
  if (!is.null(mane_id) && mane_id %in% lbls) {
    lbls[mane_id] <- sprintf("<span style='color:#006400'><b>%s</b></span>", mane_id)
  }
  if (include_variant_row) {
    lbls["Case Variants"] <- "<span style='color:red;'>Case Variants</span>"
  }
  
  p_transcripts <- ggplot(exons, aes(y = factor(transcript_id, levels = levels_top_with_variants))) +
    geom_range(aes(xstart = start, xend = end), color = "black", alpha = 0.7) +
    geom_intron(data = to_intron(exons, "transcript_id"),
                aes(xstart = start, xend = end, strand = strand)) +
    geom_range(data = df_diff, aes(xstart = start, xend = end, fill = diff_type)) +
    scale_x_continuous(limits = c(x_min, x_max)) + 
    scale_y_discrete(limits = levels_for_scale, labels = lbls, drop = FALSE, position = "right") +
    scale_fill_manual(values = c("in_ref" = "#FF8888", "not_in_ref" = "#55BCC2"),
                      labels = c("in_ref" = "In MANE", "not_in_ref" = "Not in MANE"), na.translate = FALSE) +
    coord_cartesian(clip = "off") +
    theme_bw() +
    theme(plot.margin = margin(20, 50, 20, 30),
          axis.title.y = element_blank(),
          axis.text.y  = ggtext::element_markdown(size = 9, margin = margin(l = 6)),
          legend.title = element_blank(),
          legend.position = "bottom") +
    labs(x = "Genomic coordinate", title = "Structural differences compared to MANE Select")
  
  p_transcripts <- p_transcripts +
    geom_text(
      data = df_diff,
      aes(x = (start + end)/2, y = factor(transcript_id, levels = levels_top_with_variants)),
      label = "*", nudge_y = 0.2, size = 5, color = "black", inherit.aes = FALSE
    )
  
  if (include_variant_row) {
    # Points
    p_transcripts <- p_transcripts +
      geom_point(
        data = variant_labels,
        aes(x = pos, y = factor("Case Variants", levels = levels_top_with_variants)),
        color = "red", size = 3, shape = 4, inherit.aes = FALSE
      ) +
      # Labels (basic overlap control; for many labels consider ggrepel)
      geom_text(
        data = variant_labels,
        aes(x = pos, y = factor("Case Variants", levels = levels_top_with_variants), label = label),
        color = "red", size = 1.5, vjust = 2, check_overlap = TRUE, inherit.aes = FALSE
      )
  }
  
  p_variants <- if (nrow(clinvar_variants_case) > 0) {
    ggplot(clinvar_variants_case, aes(x = pos, y = 1)) +
      geom_point(aes(shape = clnsig, fill = clnsig), size = 3, color = "black",
                 position = position_jitter(height = 0.3)) +
      scale_shape_manual(values = c("Pathogenic" = 21, "Likely pathogenic" = 21, "Pathogenic/Likely pathogenic" = 21),
                         name = "ClinVar Variants") +
      scale_fill_manual(values = c("Pathogenic" = "red", "Likely pathogenic" = "yellow", "Pathogenic/Likely pathogenic" = "orange"),
                        name = "ClinVar Variants") +
      scale_x_continuous(limits = c(x_min, x_max)) +
      theme_bw() +
      theme(axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.x = element_blank(),
            legend.position = "bottom",
            legend.title = element_text(size = 10, face = "bold"))
  } else {
    make_blank_plot("No ClinVar variants in gene window")
  }
  
  transcript_structure_plot <- plot_grid(p_transcripts, p_variants, ncol = 1, align = "v", rel_heights = c(4, 2))
  
} else {
  # No exons after filtering: still emit a valid additional_info table with NA diff columns
  additional_info <- prop_summary %>%
    dplyr::mutate(
      diff_start = NA_real_,
      diff_end   = NA_real_,
      diff_width = NA_real_,
      diff_type  = NA_character_
    )
  write.csv(additional_info, additional_info_file, row.names = FALSE)
  
  transcript_structure_plot <- make_blank_plot("No exons after filtering")
}

ggsave(transcript_structure_plot_file, plot = transcript_structure_plot, width = 7, height = 5, dpi = 300)

# -------------------- 4) HTML report --------------------
img1 <- dataURI(file = transcripts_detected_plot_file, mime = "image/png")
img2 <- dataURI(file = proportion_plot_file, mime = "image/png")
img3 <- dataURI(file = transcript_structure_plot_file, mime = "image/png")

html_content <- paste0(
  "<html><head><title>Final Report for ", opt$case_name, "</title></head><body>\n",
  "<h1>Final Report for ", opt$case_name, "</h1>\n",
  "<p>Gene: ", opt$gene, " (", opt$chrom, ":", opt$start_pos, "-", opt$end_pos, ")</p>\n",
  "<h2>1) Transcripts Detected Over Coverage Threshold</h2>\n",
  "<img src='", img1, "' width='1200' />\n",
  "<h2>2) Transcript Coverage Proportion</h2>\n",
  "<img src='", img2, "' width='1200' />\n",
  "<h2>3) Transcript Structure & Variants</h2>\n",
  "<img src='", img3, "' width='1400' />\n",
  "<h2>4) Downloadable Tables</h2>\n",
  "<ul>\n",
  "<li>Additional info CSV: ", basename(additional_info_file), "</li>\n",
  "</ul>\n",
  "</body></html>"
)
cat(html_content, file = opt$output_html)
