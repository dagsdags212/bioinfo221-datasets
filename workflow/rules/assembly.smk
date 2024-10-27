assembly_env = "../envs/assembly.yaml"

rule spades_assembly:
    input:
        reads=config["assembly"]["reads"]
    output:
        dirname=directory("output/spades"),
        contigs="output/spades/scaffolds.fasta"
    conda:
        assembly_env
    log:
        "logs/spades/spades.log"
    threads: 8
    message:
        "Spades assembly done. View contigs in {output}/contigs.fasta"
    shell:
        "spades.py -1 {input.reads[0]} -2 {input.reads[1]} -o {output.dirname}"

rule quast_qc:
    input:
        rules.spades_assembly.output.contigs
    output:
        directory("output/quast")
    params:
        ref=config["assembly"]["ref"]
    conda:
        assembly_env
    log:
        "logs/quast/quast.log"
    threads: 4
    shell:
        "quast -r {params.ref} {input} -o {output}"
