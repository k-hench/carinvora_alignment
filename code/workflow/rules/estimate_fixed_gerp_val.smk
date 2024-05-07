"""
snakemake -n --rerun-triggers mtime -R fixed_gerp

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
        -R fixed_gerp
"""
localrules: extract_ancestral_tree, compile_fixed_gerp_vcf, call_fixed_gerp

rule fixed_gerp:
    input:
      gerp = "../results/fixed_gerp/fixed_gerp.maf.rates"

rule extract_ancestral_tree:
    input:
      hal = '../results/cactus/{name}.hal'.format(name = P_NAME)
    output:
      tree = "../results/fixed_gerp/anc.tree"
    params:
      txt = '../results/cactus/{name}.txt'.format(name = P_NAME)
    shell:
      """
      head -n 1 {params.txt} > {output.tree}
      """

rule compile_fixed_gerp_vcf:
    input:
      tree = "../results/fixed_gerp/anc.tree"
    output:
      vcf = "../results/fixed_gerp/fixed_gerp.vcf",
      maf = "../results/fixed_gerp/fixed_gerp.maf",
      pdf = "../results/img/clades_tree.pdf"
    log: "logs/fixed_gerp_vcf.log"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/gerp_fixation_levels.R {input.tree} {output.pdf} {output.vcf} {output.maf} &>> {log}
      """

rule call_fixed_gerp:
    input:
      maf = "../results/fixed_gerp/fixed_gerp.maf",
      tree = "../results/neutral_tree/rerooted.tree"
    output:
      rates = "../results/fixed_gerp/fixed_gerp.maf.rates"
    params:
      refname = SPEC_REF 
    conda: "msa_phast"
    log: "logs/gerp_fixed.log"
    shell:
      """
      gerpcol -t {input.tree} -f {input.maf} -e {params.refname} -j -z -x ".rates" &> {log}
      """