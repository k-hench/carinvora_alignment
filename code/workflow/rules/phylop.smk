"""
snakemake --rerun-triggers mtime -n all_phylop

snakemake --jobs 50 \
  --latency-wait 30 \
  -p \
  --default-resources mem_mb=51200 threads=1 \
  --use-singularity \
  --singularity-args "--bind $CDATA" \
  --use-conda \
  --rerun-triggers mtime \
  --cluster '
    sbatch \
      --export ALL \
      -n {threads} \
      -e logs/{name}.{jobid}.err \
      -o logs/{name}.{jobid}.out \
      --mem={resources.mem_mb}' \
      --jn job_c.{name}.{jobid}.sh \
      -R all_phylop
"""

localrules: merge_mafs, unzip_autosome_maf, extract_cds_by_scaff

PHYLOP_SCORES = [ "SPH", "LRT", "SCORE", "GERP" ]

rule all_phylop:
    input: 
      maf = "../results/neutral_tree/windows/autosomes.maf.gz",
      model = "../results/phylop/autosomes_neutral.mod",
      scores = expand( "../results/phylop/{score_type}/{score_type}_{mscaf_nr}.tsv.gz", mscaf_nr = MSCAFS, score_type = PHYLOP_SCORES )

rule merge_mafs:
    input:
      maf = "../results/neutral_tree/windows/{mscaf}.maf.gz".format( mscaf = SCFS[0] ),
      mafs = expand( "../results/neutral_tree/windows/{mscaf}.maf.gz", mscaf = SCFS[1:17] )
    output:
      maf = "../results/neutral_tree/windows/autosomes.maf.gz"
    params:
      prefix = "../results/neutral_tree/windows/autosomes.maf"
    shell:
      """
      zcat {input.maf}  > {params.prefix}

      for k in {input.mafs}; do zcat ${{k}} | tail -n +2 >> {params.prefix}; done

      gzip {params.prefix}
      """

rule unzip_autosome_maf:
    input:
       maf = "../results/neutral_tree/windows/autosomes.maf.gz"
    output:
      maf = temp( "../results/phylop/autosomes.maf" )
    shell:
      """
      zcat {input.maf} > {output.maf}
      """

# hal prefixes chrom names with species, this needs to be 
# added to the gff seq names too for msa_view to match
rule extract_cds_by_scaff:
    input:
      gff = GFF_FILE
    output:
      gff = "../results/phylop/cds/cds_{mscaf}.gff"
    params:
      scf = "mscaf_a1_{mscaf}",
      refspec = SPEC_REF
    shell:
      """
      grep "mscaf_a1_{wildcards.mscaf}" {input.gff} | \
        grep -w CDS | \
        sed 's/{params.scf}/{params.refspec}.{params.scf}/g' > {output.gff}
      """

rule phylofit_extract_codons:
    input:
      maf = "../results/maf/carnivora_set_{mscaf}.maf",
      gff = "../results/phylop/cds/cds_{mscaf}.gff"
    output:
      stats = "../results/phylop/sufficient_stats/mscaf_a1_{mscaf}_codons.ss"
    conda: "msa_phast"
    shell:
      """
      msa_view \
        {input.maf} \
         --in-format MAF \
         --4d \
         --features {input.gff} > {output.stats}
      """

rule phylofit_extract_3rd_codons:
    input:
      stats = "../results/phylop/sufficient_stats/{mscaf}_codons.ss"
    output:
      sites = "../results/phylop/sufficient_stats/{mscaf}_sites.ss"
    conda: "msa_phast"
    shell:
      """
      msa_view \
        {input.stats} \
        --in-format SS \
        --out-format SS \
        --tuple-size 1 > {output.sites}
      """

rule phylofit_aggregate_sites:
    input:
      sites = expand( "../results/phylop/sufficient_stats/{mscaf}_sites.ss", mscaf = SCFS[0:17] )
    output:
      sites = "../results/phylop/sufficient_stats/all_4d_sites.ss"
    params:
      species = ",".join(SPEC_ALL)
    conda: "msa_phast"
    shell:
      """
      msa_view \
        --unordered-ss \
        --out-format SS \
        --aggregate {params.species} \
        {input.sites} > {output.sites}
      """

rule phylop_model:
    input:
      sites = "../results/phylop/sufficient_stats/all_4d_sites.ss",
      tree = "../results/neutral_tree/rerooted.tree"
    output:
      model = "../results/phylop/autosomes_neutral.mod"
    log:
      "logs/phylofit.log"
    conda: "msa_phast"
    shell:
      """
      phyloFit \
        --tree {input.tree} \
        --EM \
        --subst-mod REV \
        --msa-format SS \
        --log {log} \
        {input.sites} > {output.model}
      """

rule call_phylop:
    input:
      maf = "../results/maf/carnivora_set_{mscaf_nr}.maf",
      model = "../results/phylop/autosomes_neutral.mod"
    output:
      txt = "../results/phylop/{score_type}/raw/{score_type}_{mscaf_nr}.txt.gz"
    log: "logs/phylop_{score_type}_{mscaf_nr}.log"
    params:
      mscaf = "mscaf_a1_{mscaf_nr}"
    conda: "msa_phast"
    shell:
      """
      phyloP \
         --base-by-base \
         --method {wildcards.score_type} \
         --mode CONACC \
         --msa-format MAF \
         --log {log} \
         --refidx 0 \
         --chrom {params.mscaf} \
         {input.model} {input.maf} | \
         gzip > {output.txt}
      """

rule pylop_to_tsv:
    input:
      txt = "../results/phylop/{score_type}/raw/{score_type}_{mscaf_nr}.txt.gz"
    output:
      tsv = "../results/phylop/{score_type}/{score_type}_{mscaf_nr}.tsv.gz"
    params:
      mscaf = "mscaf_a1_{mscaf_nr}"
    conda: "popgen_basics"
    shell:
      """
      py/phylotxt2tsv \
        --txt {input.txt} \
        --drop-zeros \
        --seqname {params.mscaf} | \
        bgzip > {output.tsv}
      
      tabix -p bed {output.tsv}
      """