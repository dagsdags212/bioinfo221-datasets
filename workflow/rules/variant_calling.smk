variant_calling_env = "bioinfo221"

aln_file = config["mapping"]["outname"]
vcf_file = config["vc"]["outname"]

rule gatk_create_dict:
    input:
        ref=config["assembly"]["ref"]
    output:
        "output/gatk/ref.dict"
    conda:
        variant_calling_env
    log:
        "logs/gatk/create.dict.log"
    threads: 1
    shell:
        """
        gatk CreateSequenceDictionary REFERENCE={input.ref} OUTPUT={output}
        cp {output} data/
        """

rule gatk_preprocess:
    input:
        ref=config["assembly"]["ref"],
        sorted_bam=rules.gatk_sam_sort.output,
        ref_dict=rules.gatk_create_dict.output
    output:
        rg=f"output/gatk/{aln_file}.rg.bam",
        dedup=f"output/gatk/{aln_file}.dedup.bam",
        dedup_idx=f"output/gatk/{aln_file}.dedup.bam.bai"
    params:
        RGID="sam",
        RGLB="lib",
        RGPL="illumina",
        RGSM="Sample1",
        RGPU=1
    conda:
        variant_calling_env
    log:
        "logs/gatk/gatk.preprocess.log" 
    threads: 8
    script:
        "../scripts/gatk-preprocess.sh"

rule gatk_call_snps:
    input:
        ref=rules.gatk_preprocess.input.ref,
        dedup=rules.gatk_preprocess.output.dedup,
        ref_dict=rules.gatk_preprocess.input.ref_dict
    output:
        f"output/gatk/{vcf_file}.vcf"
    conda:
        variant_calling_env
    log:
        "logs/gatk/gatk.call.log"
    threads: 8
    shell:
        "gatk HaplotypeCaller -R {input.ref} -I {input.dedup} "
        "-O {output} --sequence-dictionary {input.ref_dict}"

rule snpeff_annotate_snps:
    input:
        vcf=rules.gatk_call_snps.output
    output:
        f"output/snpeff/{vcf_file}.ann.vcf"
    params:
        db=config["vc"]["db"]
    conda:
        variant_calling_env
    log:
        "logs/snpeff/snpeff.annotate.log"
    threads: 8
    shell:
        """
        snpEff {params.db} {input.vcf} > {output}
        mv snpEff* output/snpeff
        """