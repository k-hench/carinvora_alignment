"""
snakemake -n --rerun-triggers mtime -R ancarcgaz

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
      tsv = "../results/anc_allele/arcgaz_anc41_snp_rs.tsv.gz"

rule anc_pos:
    input:
      tsv = "../results/anc_allele/arcgaz_anc41_snps.tsv.gz"
    output:
      tsv = "../results/anc_allele/arcgaz_anc41_snp_pos.tsv.gz"
    conda: "popgen_basics"
    shell:
      """
      zcat {input.tsv} | \
        cut -f 1,2 | \
        bgzip > {output.tsv}
      
      tabix -s 1 -S 1 -b 2 -e 2 {output.tsv}
      """

rule snp_pos_to_bed:
    input:
      tsv = "../results/anc_allele/arcgaz_anc41_snp_pos.tsv.gz",
      bed = "../results/gerp/rs/gerp_rs.bed.gz",
      tbi = "../results/gerp/rs/gerp_rs.bed.gz.tbi"
    output:
      bed = "../results/anc_allele/arcgaz_anc41_snp_pos.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      zcat {input.tsv} | \
        grep -v ^ref | \
         awk '{{print $1"\t"$2-1"\t"$2}}'  | \
        bgzip > {output.bed}
      """

rule attach_rs_to_bed:
    input:
      bed = "../results/gerp/rs/gerp_rs.bed.gz",
      bed_pos = "../results/anc_allele/arcgaz_anc41_snp_pos.bed.gz"
    output:
      tsv = "../results/anc_allele/arcgaz_anc41_snp_rs.tsv.gz"
    conda: "r_tidy"
      """
      # tabix -R {input.bed_pos} {input.bed} | bgzip > {output.pre_bed} 
      Rscript --vanilla code/R/instersect_pos_rs.R
      """
# 
# rule querry_anc_gerp:
#     input:
#       bed = "../results/anc_allele/arcgaz_anc41_snp_pos.bed.gz",
#       pre_bed = "../results/anc_allele/arcgaz_anc41_snp_rs.pre_bed.gz"
#     output:
#       tsv = "../results/anc_allele/arcgaz_anc41_snp_rs.tsv.gz"
#     conda: "popgen_basics"
#       """
#       bedtools intersect -a {input.pre_bed} -b {input.bed} -wa -wb | \
#         awk 'BEGIN{{print"chrom\tpos\trs"}}{{print $1"\t"$7"\t"$4}}' | \
#         bgzip > {output.tsv}
#       """