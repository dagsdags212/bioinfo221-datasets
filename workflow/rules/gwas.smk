env = "R"

rule run_gwas:
    output:
        em_model="output/gwas/em.out.RData",
        hk_model="output/gwas/hk.out.RData",
        fq_model="output/gwas/fq.model.RData",
        stepwise_model="output/gwas/stepwise.model.RData"
    params:
        dataset=config["gwas"]["dataset"]
    log:
        "logs/gwas/rqtl.log"
    conda:
        env
    threads: 4
    script:
        "../scripts/run-gwas.R"