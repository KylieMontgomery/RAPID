############################################
# Load the config file
############################################
configfile: "config.yaml"

CASES = config["cases"]
CASE_NAMES = [c["name"] for c in CASES]

ALL_SAMPLES = sorted(set(
    [c["case_sample"] for c in CASES] +
    [c["case_control"] for c in CASES]
))

DATA_DIR = config["data_dir"].rstrip("/")
OUTPUT_DIR = config["output_dir"].rstrip("/")
REFERENCE = config["reference_fasta"]
REFERENCE_GTF = config["reference_gtf"]
THREADS = config["threads"]

############################################
# Sample list
############################################
SAMPLES = ALL_SAMPLES

############################################
# Define the final rule "all"
############################################
rule all:
    input:
        # Sorted BAMs (global)
        expand(f"{OUTPUT_DIR}/all_samples/{{sample}}.sorted.bam", sample=SAMPLES),
        
        # NanoStat outputs (global)
        expand(f"{OUTPUT_DIR}/all_samples/{{sample}}_stats", sample=SAMPLES),

        # MergedQuant outputs (global)
        expand(f"{OUTPUT_DIR}/all_samples/{{sample}}_mergedQuant.gtf", sample=SAMPLES),
        expand(f"{OUTPUT_DIR}/all_samples/{{sample}}_mergedQuant.abundance.txt", sample=SAMPLES),

        # Final merged GTFs and HTML reports
        f"{OUTPUT_DIR}/merged.gtf",
        f"{OUTPUT_DIR}/final_output/final_project_output.gtf",
        expand(f"{OUTPUT_DIR}/final_output/{{case}}.html", case=CASE_NAMES)

############################################
# Per-sample Rules (deduplicated)
############################################
rule nanostat:
    input:
        fastq=lambda wc: f"{DATA_DIR}/{wc.sample}.fastq"
    output:
        directory(f"{OUTPUT_DIR}/all_samples/{{sample}}_stats")
    threads: THREADS
    params:
        outdir=lambda wc: f"{OUTPUT_DIR}/all_samples/{wc.sample}_stats",
        statname=lambda wc: f"{wc.sample}_stats"
    shell:
        """
        NanoStat --fastq {input.fastq} --outdir {params.outdir} --name {params.statname}
        """

rule cutadapt:
    input:
        fastq=lambda wc: f"{DATA_DIR}/{wc.sample}.fastq"
    output:
        trimmed=protected(f"{OUTPUT_DIR}/all_samples/{{sample}}.trimmed.fq")
    threads: THREADS
    params:
        fwd_primer=config.get("fwd_primer", ""),
        rev_primer=config.get("rev_primer", "")
    shell:
        """
        cutadapt -j {threads} -g file:{params.fwd_primer} -a file:{params.rev_primer} \
        --discard-untrimmed -o {output.trimmed} {input.fastq}
        """

rule minimap2:
    input:
        trimmed=f"{OUTPUT_DIR}/all_samples/{{sample}}.trimmed.fq"
    output:
        sam=protected(f"{OUTPUT_DIR}/all_samples/{{sample}}.sam")
    params:
        ref=REFERENCE
    shell:
        """
        minimap2 -ax splice -t {THREADS} -uf --secondary=no -C5 {params.ref} {input.trimmed} > {output.sam}
        """

rule sam_to_sorted_bam:
    input:
        sam=f"{OUTPUT_DIR}/all_samples/{{sample}}.sam"
    output:
        sorted_bam=f"{OUTPUT_DIR}/all_samples/{{sample}}.sorted.bam"
    threads: THREADS
    shell:
        """
        samtools view -@ {threads} -bS {input.sam} | samtools sort -@ {threads} -o {output.sorted_bam}
        """

rule stringtie_unguided:
    input:
        sorted_bam=f"{OUTPUT_DIR}/all_samples/{{sample}}.sorted.bam"
    output:
        gtf=protected(f"{OUTPUT_DIR}/all_samples/{{sample}}_unguided.gtf")
    threads: THREADS
    shell:
        """
        stringtie {input.sorted_bam} -o {output.gtf} -p {threads}
        """

rule gtf_list:
    input:
        expand(f"{OUTPUT_DIR}/all_samples/{{sample}}_unguided.gtf", sample=SAMPLES)
    output:
        protected(f"{OUTPUT_DIR}/gtf_filepath.txt")
    run:
        with open(output[0], "w") as out_handle:
            for gtf in input:
                out_handle.write(gtf + "\n")

rule stringtie_merge:
    input:
        gtf_filepath=f"{OUTPUT_DIR}/gtf_filepath.txt"
    output:
        merged_gtf=protected(f"{OUTPUT_DIR}/merged.gtf")
    params:
        ref=REFERENCE_GTF
    threads: THREADS
    shell:
        """
        stringtie --merge -o {output.merged_gtf} -G {params.ref} -p {threads} {input.gtf_filepath}
        """

rule stringtie_requant:
    input:
        sorted_bam=f"{OUTPUT_DIR}/all_samples/{{sample}}.sorted.bam",
        merged_gtf=f"{OUTPUT_DIR}/merged.gtf"
    output:
        gtf=f"{OUTPUT_DIR}/all_samples/{{sample}}_mergedQuant.gtf",
        abundance=f"{OUTPUT_DIR}/all_samples/{{sample}}_mergedQuant.abundance.txt"
    threads: THREADS
    shell:
        """
        stringtie {input.sorted_bam} -e -G {input.merged_gtf} -o {output.gtf} -A {output.abundance} -p {threads}
        """

rule naive_concat_gtf_cov:
    input:
        expand(f"{OUTPUT_DIR}/all_samples/{{sample}}_mergedQuant.gtf", sample=SAMPLES)
    output:
        f"{OUTPUT_DIR}/final_output/final_project_output.gtf"
    shell:
        """
        mkdir -p $(dirname {output})
        rm -f {output}
        for gtf in {input}; do
            sample_base=$(basename "$gtf" .gtf | sed 's/_mergedQuant$//')
            awk -v SAMPLE="$sample_base" -v FO="\\t" 'BEGIN {{FS=FO; OFS=FO}}
            /^#/ {{print; next}}
            {{
                if ($9 !~ /;$/) {{ $9 = $9 ";" }}
                $9 = $9 " sample_id \\"" SAMPLE "\\";"
                print
            }}' "$gtf" >> {output}
        done
        """

############################################
# Helper function for reports
############################################
def _samples_for_case(case_name):
    for c in CASES:
        if c["name"] == case_name:
            return [c["case_sample"], c["case_control"]]
    return []

rule final_html_report:
    input:
        merged_gtf=f"{OUTPUT_DIR}/final_output/final_project_output.gtf",
        sorted_bams=lambda wc: [
            f"{OUTPUT_DIR}/all_samples/{sample}.sorted.bam"
            for sample in _samples_for_case(wc.case)
        ]
    output:
        html=f"{OUTPUT_DIR}/final_output/{{case}}.html"
    params:
        mane_transcripts=config["MANE_transcripts"],
        clinvar_variants=config["clinvar_variants"]
    run:
        c = next(x for x in CASES if x["name"] == wildcards.case)
        variant_names_str = ",".join(c.get("variant_name", []))
        variant_pos_str = ",".join(str(pos) for pos in c.get("variant_pos", []))
        shell(f"""
            Rscript scripts/final_report.R \
              --case_name "{c['name']}" \
              --case_sample "{c['case_sample']}" \
              --case_control "{c['case_control']}" \
              --chrom "{c['chrom']}" \
              --start_pos {c['start_pos']} \
              --end_pos {c['end_pos']} \
              --gene "{c['gene']}" \
              --mane_id "{c['mane_id']}" \
              --variant_name "{variant_names_str}" \
              --variant_pos "{variant_pos_str}" \
              --merged_gtf "{input.merged_gtf}" \
              --MANE_transcripts "{params.mane_transcripts}" \
              --clinvar_variants "{params.clinvar_variants}" \
              --output_html "{output.html}"
        """)
