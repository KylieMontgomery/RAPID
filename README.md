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

Wet-lab protocol:
RAPID - RNA Analysis Pipeline for Integrated Diagnostics. protocols.io
https://dx.doi.org/10.17504/protocols.io.8epv5ow45g1b/v2 

------------------------------------------------------------------------

# Installation

Clone the repository:

``` bash
git clone https://github.com/KylieMontgomery/RAPID.git
cd RAPID
```

Create the conda environment:

``` bash
conda env create -f environment.yml
conda activate RAPID_env
```

## Docker image

A prebuilt Docker image with the software environment (R, Snakemake, minimap2,
samtools, StringTie, cutadapt, NanoStat, bcftools, and the required R packages
including ggtranscript) is available, so you do not need to build the conda
environment manually. The image contains the software environment only — the
pipeline code (Snakefile and scripts) is run from your cloned copy of this
repository, mounted into the container at runtime.

```bash
# 1. Clone this repository and enter it
git clone https://github.com/KylieMontgomery/RAPID.git
cd RAPID

# 2. Pull the image
docker pull kyliemontgomery/rapid:1.0.1

# 3. Run the pipeline, mounting the repo into the container
docker run --rm -v "$PWD":/data -w /data kyliemontgomery/rapid:1.0.1 \
  snakemake --cores 8
```

This is the recommended way to run RAPID reproducibly across different systems.

## Install additional R packages

One R package used for the transcript-structure plots is not available
through conda and must be installed separately after creating the environment:

With the `RAPID_env` environment activated:

```bash
conda activate RAPID_env

Rscript -e 'remotes::install_github("dzhang32/ggtranscript")'
```

## Test dataset

A small, fully synthetic test dataset is included in `test_data/` so you can verify
your installation and see what a successful run looks like before using your own data.
No patient data is used — reads are simulated with
[Badread](https://github.com/rrwick/Badread) (R10.4 / `nanopore2023` error model)
from the genomic HBB locus.

The dataset models a multi-isoform HBB profile so that the report plots resemble a
real run. The case sample contains a predominant novel exon-2-skip transcript
(~60%), the canonical MANE transcript (~30%), a low-abundance intron-retention
isoform (~10%), and a trace cryptic-junction "noise" transcript; the control sample
is dominated by the canonical transcript with the same low-level intron retention and
noise, and no exon-2 skip. This case-versus-control design exercises the full
pipeline — detection, quantification, case–control comparison, and the
transcript-structure and ClinVar plots — and illustrates the transcript detectability
plot: as the proportion threshold is raised, low-abundance noise and then minor
isoforms drop out, while the case-specific skip stands out against the control.

All transcripts are fabricated test structures, not real HBB biology; they exist only
to produce a realistic multi-transcript profile for testing.

The test config (`config.yaml`) and all required inputs ship in the repository, so the
test runs with no external downloads. From the repository root:

```bash
snakemake --cores 4 --configfile config.yaml
```

Outputs are written to `test_results/`; open the HTML report at
`test_results/final_output/example_HBB_case.html` to confirm the run succeeded.

## Required input files Reference files

### Download GENCODE v47 reference genome and annotation:

Genome fasta:
<ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz>

Annotation GTF:
<ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz>

Place in:

```         
mkdir -p ref

wget -O ref/genome.hg38.fa.gz \
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

wget -O ref/annotation.v47.gtf.gz \
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

gunzip ref/genome.hg38.fa.gz
gunzip ref/annotation.v47.gtf.gz
```

This results in reference fasta and gtf files

```         
ref/genome.hg38.fa
ref/annotation.v47.gtf
```

### Primer fasta files

> **Note:** Primer and amplicon *design* (target selection, MANE Select prioritisation,
> 1–3 kb amplicon sizing, multi-amplicon strategies for large genes, and primer
> verification) is part of library preparation and is documented in the
> [protocols.io protocol](https://dx.doi.org/10.17504/protocols.io.8epv5ow45g1b/v1).
> This section only covers formatting your designed primers for input to the pipeline.

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

RAPID uses a preprocessed ClinVar TSV (chrom, pos, ref, alt, clnsig, gene) for pathogenic / likely pathogenic variant annotation in the transcript-structure plots.

Generate it in three steps — download the VCF, extract the raw fields with bcftools, then filter and tidy with the R script:

```bash
mkdir -p ref

# 1. Download the latest ClinVar GRCh38 VCF.
#    NCBI updates this file weekly (Mondays) with a full monthly release on the
#    first Thursday, so the rolling 'clinvar.vcf.gz' always resolves.
#    For a reproducible analysis, pin a specific dated release instead, e.g.
#    clinvar_20251116.vcf.gz, and record that filename in your methods.
#    Available releases: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/
wget -O ref/clinvar.vcf.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz

# 2. Extract raw fields (chrom, pos, ref, alt, clnsig, geneinfo)
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%CLNSIG\t%INFO/GENEINFO\n' \
  ref/clinvar.vcf.gz > ref/clinvar_raw.tsv

# 3. Filter to Pathogenic / Likely pathogenic and write the final TSV
Rscript scripts/import_clinvar_variants.R
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
# Global config, e.g. references

threads: 8
data_dir: "data"
output_dir: "results"
reference_fasta: "ref/genome.hg38.fa"
reference_gtf: "ref/annotation.v47.gtf"
fwd_primer: "ref/forward_primers.fasta"
rev_primer: "ref/reverse_primers.fasta"
clinvar_variants: "ref/clinvar_variants.tsv"
MANE_transcripts: "ref/MANE_select_transcripts.tsv"

cases:
  - name: Case1
    case_sample: sample_1
    case_control: control_1
    chrom: "chr11"
    start_pos: 5225364 #Start position of the gene -100bp to increase visualisation window
    end_pos: 52271171 #End position of the gene +100bp to increase visualisation window
    gene: "HBB"
    variant_name: [] #Example format ["c.1234A>G", "c.111+222A>G"]
    variant_pos: [] #Genomic position of variant for visualising on the output figures, cannot be a range
    prop_threshold: 0.05 #This allows for manual resetting of the proportion threshold for complex transcript profiles e.g. 0.001
    mane_id: "ENST00000335295.4" #Or the most highly expressed transcript ID for the tissue utilised, for comparison
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

#### Transcript detectability plots

The first report output is a set of transcript detectability plots showing how the
abundance threshold affects transcript interpretation. Three line plots are produced:
one for the case sample, one for the pooled-control sample, and one overlay comparing
the two. The x-axis is the transcript-abundance threshold applied; the y-axis is the
number of transcript models retained above that threshold. As the threshold increases,
progressively fewer low-abundance transcript models are retained. These plots let you
assess whether the profile is dominated by a few abundant isoforms or whether many
low-abundance models are present, and help guide threshold choice.

The default (`prop_threshold: 0.05`) is set at 5% to simplify visualisation. **It is not
a biological or diagnostic cut-off.** Lower it in `config.yaml` where rare transcripts
are of interest, including transcripts potentially subject to nonsense-mediated decay
(NMD), which may fall below 5% of the steady-state RNA pool.

#### 2. Transcript proportion plot

A bar plot showing the expression of each transcript in the case and control samples,
expressed as a proportion of all transcripts generated from the target locus in that
sample type. Because values are within-sample proportions, this plot highlights shifts
in relative isoform usage between case and control — for example a novel transcript
dominating the case profile while being absent or minor in the pooled control.

#### 3. Transcript structure plot

Visualises transcript structures aligned against the relevant reference transcript
(usually MANE Select). Regions matching the reference are shown separately from regions
not present in it, such as retained intronic sequence or novel exonic sequence. Exonic
regions that differ structurally from the reference are marked with an asterisk (`*`).
Connecting lines between exons represent intronic sequence / splice junctions.

A separate lower **ClinVar track** displays pathogenic and likely pathogenic variants,
providing clinical context by highlighting regions of the gene known to be sensitive to
pathogenic variation.

Citation

<https://www.medrxiv.org/content/10.64898/2025.12.30.25342835v1>

License

This project is licensed under the Apache License 2.0 --- see the LICENSE file
for details.

------------------------------------------------------------------------
