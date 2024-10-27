from pathlib import Path
import os


read_align_env = "../envs/variant_calling.yaml"
aln_file = config["mapping"]["outname"]
ref_stem = Path(config["assembly"]["ref"]).stem

rule bowtie2_build:
    input:
        ref=config["assembly"]["ref"]
    output:
        index=multiext(
            "output/bowtie2/ref",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        )
    params:
        lambda wildcards, output: os.path.commonprefix(output).rstrip(".")
    log:
        f"logs/bowtie2/bt.index.{aln_file}.log",
    conda:
        read_align_env
    threads: 8
    shell:
        "bowtie2-build {input.ref} {params}" 

rule bowtie_map_reads:
    input:
        reads=config["assembly"]["reads"],
        index=rules.bowtie2_build.output,
    output:    
        f"output/bowtie2/{aln_file}.sam"
    params:
        lambda wildcards, input: os.path.commonprefix(input.index).rstrip(".")
    conda:
        read_align_env
    threads: 8 
    shell:
        "bowtie2 -x {params} -1 {input.reads[0]} -2 {input.reads[1]} "
        "-S {output}"

rule samtools_faidx:
    input:
        ref=config["assembly"]["ref"]
    output:
        f"data/{ref_stem}.fasta.fai",
    log:
        f"logs/samtools/samtools.faidx.{ref_stem}.log",
    conda:
        read_align_env
    threads: 2
    shell:
        "samtools faidx {input.ref}"
    
rule gatk_sam_sort:
    input:
        rules.bowtie_map_reads.output
    output:
        f"output/gatk/{aln_file}.sorted.bam"
    conda:
        "bioinfo221"
    log:
        "logs/samtools/samtools.sort.log"
    threads: 4
    shell:
        "gatk SortSam INPUT={input} OUTPUT={output} "
        "SO=coordinate"