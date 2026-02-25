# RAPID: RNA Analysis Pipeline for Integrated Diagnostics

RAPID is a fully automated Snakemake workflow for the analysis of
targeted long-read RNA sequencing data generated using Oxford Nanopore
Technologies platforms. The pipeline performs quality control, primer
trimming, genome alignment, transcript assembly, quantification, and
generates case-control comparative transcript visualizations and
reports.

RAPID is designed for targeted transcript analysis in rare disease
diagnostics to identify aberrant splicing, novel transcripts, and
altered isoform usage.

Wet-lab protocol:\
<https://dx.doi.org/10.17504/protocols.io.8epv5ow45g1b/v1>

------------------------------------------------------------------------

# Installation

Clone the repository:

``` bash
git clone https://github.com/KylieMontgomery/RAPID.git
cd RAPID
```

Create the conda environment:

``` bash
mamba env create -f environment.yml
conda activate rapid_env
```

## Required input files Reference files

### Download GENCODE v47 reference genome and annotation:

Genome fasta:
<ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz>

Annotation GTF:
<ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz>

Place in:

```         
mkdir -p ref

wget -O ref/genome.fa.gz \
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

wget -O ref/annotation.v47.gtf.gz \
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

gunzip ref/genome.fa.gz
gunzip ref/annotation.v47.gtf.gz
```

This results in reference fasta and gtf files

```         
ref/genome.fa
ref/annotation.v47.gtf
```

### Primer fasta files

Generate using:

``` bash
Rscript scripts/Generate_fasta_for_fwd_and_rev_primers.R
```

Outputs:

```         
ref/forward_primers.fasta
ref/reverse_primers.fasta
```

Expected FASTA format Example:

```         
>HBB_forward_primer
GTGCACCTGACTCCTGAGGAGA

>HBB_reverse_primer
CCTTGATACCAACCTGCCCAG
```

### ClinVar variants

Generate using:

```         
## ClinVar (GRCh38) variants (Pathogenic / Likely pathogenic)

RAPID uses a preprocessed ClinVar TSV (chrom, pos, ref, alt, clnsig, gene) for variant annotation in plots.

### Quickstart

```bash
mkdir -p ref
wget -O ref/clinvar.vcf.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20250209.vcf.gz

Rscript scripts/import_clinvar_variants.R \
  --clinvar_vcf ref/clinvar.vcf.gz \
  --out_tsv ref/clinvar_variants.tsv
  
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%CLNSIG\t%INFO/GENEINFO\n' \
  ref/clinvar.vcf.gz > ref/clinvar_raw.tsv
```

Output:

```         
ref/clinvar_variants.tsv
```

### MANE select transcript ID's

Generate using:

```         
mkdir -p ref

wget -O ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
  https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz

Rscript scripts/MANE_transcript_ids.R \
  --mane_gtf ref/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
  --out_tsv ref/MANE_select_transcripts.tsv
```

Output:

```         
ref/MANE_select_transcripts.tsv
```

Each sample must have a single FASTQ file. If multiple FASTQ files exist
per sample, concatenate:

``` bash
mkdir -p data

cat sample_run*.fastq.gz > data/sample.fastq.gz
```

## Prepare the config.yaml

Edit config.yaml to define samples and analysis parameters.

Example:

``` bash
data_dir: "data"
output_dir: "results"

reference_fasta: "ref/genome.fa"
reference_gtf: "ref/annotation.v47.gtf"

fwd_primer: "ref/forward_primers.fasta"
rev_primer: "ref/reverse_primers.fasta"

clinvar_variants: "ref/clinvar_variants.tsv"
MANE_transcripts: "ref/MANE_transcripts.tsv"

threads: 8

cases:
  - name: Case1
    case_sample: sample_1
    case_control: control_1
    chrom: "chr11"
    start_pos: 5225464
    end_pos: 5227071
    gene: "HBB"
    variant_name: []
    variant_pos: []
    prop_threshold: 0.05
    mane_id: "ENST00000335295.4"
```

## Running the workflow

Dry run:

``` bash
snakemake -n
```

Run workflow:

``` bash
snakemake --cores 8
```

## Output files

Main outputs are written to:

```         
results/
```

Key files include:

Alignment files:

```         
results/aligned/\*.sorted.bam
```

Transcript assemblies:

```         
results/stringtie/\*.gtf results/stringtie/merged.gtf
```

Quantification:

```         
results/quant/\*.abundance.txt
```

Final transcript annotation:

```         
results/final_output/final_project_output.gtf
```

HTML reports:

```         
results/final_output/\*.html Workflow steps
```

The pipeline performs the following steps:

Quality control using NanoStat

Primer trimming using cutadapt

Alignment using minimap2

BAM sorting using samtools

Transcript assembly using StringTie in unguided mode

Transcript quantification using StringTie in merge mode

Case-control transcript comparison

Visualization and report generation using R

Report output

Each case generates an HTML report containing:

• Transcript detection plots • Case-control transcript proportion plots
• Transcript structure visualization • MANE transcript comparison •
ClinVar pathogenic variant annotation

These plots assist in identifying aberrant splicing and novel transcript
structures.

Citation

<https://www.medrxiv.org/content/10.64898/2025.12.30.25342835v1>

License

This project is licensed under the Apache License 2.0 --- see the LICENSE file
for details.

------------------------------------------------------------------------
