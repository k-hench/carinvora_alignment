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
      gerp = "../results/fixed_gerp/fixed_gerp_annotated.vcf.gz"

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

rule fixed_gerp_beds:
    input:
      rates = "../results/fixed_gerp/fixed_gerp.maf.rates"
    output:
      bed = "../results/fixed_gerp/fixed_gerp.bed.gz"
    conda: "popgen_basics"
    params:
      col_idx = 2
    shell:
      """
      awk -v s="dummy" \
        '{{ print s"\t"NR-1"\t"NR"\t"${params.col_idx} }}' {input.rates} | \
        grep -v "\s0$" | \
        bgzip > {output.bed}
      """

rule bgzip_vcf:
    input:
      vcf = "../results/fixed_gerp/fixed_gerp.vcf"
    output:
      vcf = "../results/fixed_gerp/fixed_gerp.vcf.gz",
      tbi = "../results/fixed_gerp/fixed_gerp.vcf.gz.tbi"
    conda: "popgen_basics"
    shell:
      """
      bgzip {input.vcf}
      tabix -p vcf {output.vcf}
      """

rule add_gerp_values_to_vcf:
    input:
      vcf = "../results/fixed_gerp/fixed_gerp.vcf.gz",
      gerp = "../results/fixed_gerp/fixed_gerp.bed.gz"
    output:
      vcf_head = temp( "../results/fixed_gerp/new_info_field_rs.txt" ),
      vcf = "../results/fixed_gerp/fixed_gerp_annotated.vcf.gz"
    conda: "popgen_basics"
    shell:
      """
      echo '##INFO=<ID=RS,Number=1,Type=Float,Description="GERP RS score based on 67-way carnivora alignment (max. value = cummulative branch length of underlying phylogeny: 1.348829)">' > {output.vcf_head}

      bcftools annotate \
        -a {input.gerp} \
        -c "CHROM,FROM,TO,INFO/RS" \
        -h {output.vcf_head} \
        {input.vcf} | \
        bgzip > {output.vcf}
      
      tabix -p vcf {output.vcf}
      """