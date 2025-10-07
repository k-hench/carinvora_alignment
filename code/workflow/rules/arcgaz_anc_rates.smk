"""
snakemake -n -R ancarcgaz

snakemake --jobs 50 \
    --latency-wait 30 \
    -p \
    --default-resources mem_mb=51200 threads=1 \
    --use-singularity \
    --singularity-args "--bind $CDATA" \
    --use-conda \
    --rerun-triggers mtime \
    --cluster '
      sbatch_delay \
        --export=ALL \
        -n {threads} \
        -e logs/{name}.{jobid}.err \
        -o logs/{name}.{jobid}.out \
        --mem={resources.mem_mb}' \
        --jn job_c.{name}.{jobid}.sh \
        -R ancarcgaz
"""

rule ancarcgaz:
    input:
        "../results/anc_allele/arcgaz_anc41_snp_pos.tsv.gz",  "../results/gerp/rs/gerp_rs.bed.gz.tbi"
#      tsv = "../results/anc_arcgaz_rs.tsv.gz"

rule anc_pos:
    input:
      tsv = "../results/anc_allele/arcgaz_anc41_snps.tsv.gz"
    output:
      tsv = "../results/anc_allele/arcgaz_anc41_snp_pos.tsv.gz"
    conda: "popgen_basics"
    shell:
      """
      zcat {input.tsv} | \
        awk '{{print $1"\t"$2"}}' | \
        bgzip > {output.tsv}

      tabix -s 1 -S 1 -b 2 -e 2 {output.tsv}
      """

# rule querry_anc_gerp:
#     input:
#       tsv = "../results/anc_allele/arcgaz_anc41_snp_pos.tsv.gz",
#       bed = "../results/gerp/rs/gerp_rs.bed.gz",
#       tbi = "../results/gerp/rs/gerp_rs.bed.gz.tbi"
#     output:
#       tsv = "../results/anc_allele/arcgaz_anc41_snp_rs.tsv.gz"
#     conda: "popgen_basics"
#     shell:
#       """
#       tabix -R {input.tsv} {input.bed} | \
#         bgzip > {output}
#       """

#      "../results/gerp/rs/gerp_rs.bed.gz"