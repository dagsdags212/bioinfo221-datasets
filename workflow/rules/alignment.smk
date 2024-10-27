alignment_env = "../envs/alignment.yaml"

aln_file = config["alignment"]["outname"]

rule mafft_align:
    input:
        config["alignment"]["seqs"]
    output:
        f"output/mafft/{config["alignment"]["outname"]}.fasta"
    conda:
        alignment_env
    log:
        "logs/mafft/mafft.log"
    threads: 12
    shell:
        "mafft {input} > {output}"

rule raxml_tree:
    input:
        rules.mafft_align.output
    output:
        multiext(
            "output/raxml/RAxML",
            f"_bestTree.{aln_file}",
            f"_info.{aln_file}",
            f"_log.{aln_file}",
            f"_parsimonyTree.{aln_file}",
            f"_result.{aln_file}",
        )
    params:
        seed=config["raxml"]["seed"],
        model=config["raxml"]["model"],
        outdir=config["raxml"]["outdir"],
        outname=config["alignment"]["outname"]
    conda:
        alignment_env
    log:
        "logs/raxml/raxml.log"
    threads: 4 
    shell:
        "raxmlHPC -s {input} -n {params.outname} "
        "-p {params.seed} -m {params.model} -T {threads} -w {params.outdir}"