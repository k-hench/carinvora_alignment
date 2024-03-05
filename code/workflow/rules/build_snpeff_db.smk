"""
snakemake --rerun-triggers mtime -n all_ml_snpeff
snakemake --dag  --rerun-triggers mtime -R all_ml_snpeff | dot -Tsvg > ../results/img/control/dag_snpeff_db.svg
snakemake --rulegraph all_ml_snpeff | dot -Tsvg > ../results/img/control/rules_snpeff_db.svg

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
  -R all_ml_snpeff
"""

GTF_FILE = "../data/genomes/annotation/arcGaz4_h1_annotation.gtf.gz"
GFF_FILE = "../data/genomes/annotation/arcGaz4_h1_annotation.gff3.gz"

rule all_ml_snpeff:
    input: 
      snpeff_gff = "../results/snp_eff/data/{ref}/genes.gtf.gz".format( ref = SPEC_REF ) 
      # vcf = expand( "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz", vcf_pre = "<vcf_base_name_here>" )

rule create_snpeff_config:
    output:
      conf = "../results/snp_eff/snpEff.config"
    shell:
      """
      echo "# Arctocephalus gazella genome, version arcGaz4_h1" > {output.conf}
      echo "arcgaz.4h1.genome : {SPEC_REF}" >> {output.conf}
      """

rule gunzip_fa:
    input:
      fa = REF_GENOME
    output:
      fa = temp( REF_GENOME.strip( ".gz" ) )
    script:
      """
      zcat {input.fa} > {output.fa}
      """

rule extract_cds:
    input:
      fa = REF_GENOME.strip( ".gz" ),
      gff = GFF_FILE
    output:
      cds = "../results/snp_eff/data/{ref}/cds.fa.gz".format( ref = SPEC_REF )
    params:
      cds_prefix = "../results/snp_eff/data/{ref}/".format( ref = SPEC_REF )
    conda: "gff3toolkit"
    shell:
      """
      gff3_to_fasta \
        -g {input.gff} \
        -f {input.fa} \
        -st cds \
        -d complete \
        -o {params.cds_prefix}/{SPEC_REF}

      mv {params.cds_prefix}/{SPEC_REF}_cds.fa {params.cds_prefix}/cds.fa
      gzip {params.cds_prefix}/cds.fa
      """

rule extract_prot:
    input:
      fa = REF_GENOME.strip( ".gz" ), 
      gff = GFF_FILE
    output:
      pep = "../results/snp_eff/data/{ref}/protein.fa.gz".format( ref = SPEC_REF )
    params:
      pep_prefix = "../results/snp_eff/data/{ref}".format( ref = SPEC_REF )
    conda: "gff3toolkit"
    shell:
      """
      gff3_to_fasta \
        -g {input.gff} \
        -f {input.fa} \
        -st pep \
        -d complete \
        -o {params.pep_prefix}/{SPEC_REF}
      
      mv {params.pep_prefix}/{SPEC_REF}_pep.fa {params.pep_prefix}/protein.fa
      gzip {params.pep_prefix}/protein.fa
      """

rule create_snpeff_db:
    input:
      fa = REF_GENOME.strip( ".gz" ),
      gtf = GTF_FILE,
      cds = "../results/snp_eff/data/{ref}/cds.fa.gz".format( ref = SPEC_REF ),
      prot = "../results/snp_eff/data/{ref}/protein.fa.gz".format( ref = SPEC_REF ),
      conf = "../results/snp_eff/snpEff.config"
    output:
      snp_fa = "../results/snp_eff/data/genomes/{ref}.fa".format( ref = SPEC_REF ),
      snp_gff = "../results/snp_eff/data/{ref}/genes.gtf.gz".format( ref = SPEC_REF )
    params:
      snpeff_path = "../results/snp_eff"
    resources:
      mem_mb=25600
    container: c_ml
    shell:
      """
      mkdir -p {params.snpeff_path}/data/{SPEC_REF} {params.snpeff_path}/data/genomes
      cd {code_dir}/{params.snpeff_path}/data/{SPEC_REF}
      ln -s {code_dir}/{input.gtf} ./genes.gtf.gz
      cd {code_dir}/{params.snpeff_path}/data/genomes
      ln -s {code_dir}/{input.fa} ./{SPEC_REF}.fa
      cd {code_dir}/{params.snpeff_path}
      snpEff build -Xmx24G -c {code_dir}/{input.conf} -dataDir $(pwd)/data -gtf22 -v {SPEC_REF}
      """

# # remaining workflow once there are genotypes to be annotated:
# rule snpeff_link_vcf:
#     input:
#       snpeff_gff = "../results/snp_eff/data/{ref}/genes.gtf.gz".format( ref = SPEC_REF ),
#       vcf = "../results/genotyping/filtered/{vcf_pre}.vcf.gz"
#     output:
#       vcf_ln = "../results/snp_eff/{vcf_pre}.vcf.gz"
#     params:
#       snpeff_path = "../results/snp_eff"
#     shell:
#       """
#       cd {code_dir}/{params.snpeff_path}
#       ln -s {code_dir}/{input.vcf} ./
#       """
# 
# rule run_snpeff:
#     input:
#       snpeff_gff = "../results/snp_eff/data/{ref}/genes.gtf.gz".format( ref = SPEC_REF ),
#       vcf = "../results/snp_eff/{vcf_pre}.vcf.gz"
#     output:
#       snpef_vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf",
#       report = "../results/snp_eff/{vcf_pre}_stats.html"
#     params:
#       snpeff_path = "../results/snp_eff"
#     resources:
#       mem_mb=25600
#     container: c_ml
#     shell:
#       """
#       cd {code_dir}/{params.snpeff_path}
#       snpEff ann -Xmx24G -stats {wildcards.vcf_pre}_stats.html \
#           -no-downstream \
#           -no-intergenic \
#           -no-intron \
#           -no-upstream \
#           -no-utr \
#           -v \
#           {SPEC_REF} {wildcards.vcf_pre}.vcf.gz > {code_dir}/{output.snpef_vcf}
#       """
# 
# rule bgzip_vcf:
#     input:
#       vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf"
#     output:
#       vcf = "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz",
#       vcf_idx = "../results/genotyping/annotated/{vcf_pre}_ann.vcf.gz.tbi"
#     conda: "popgen_basics"
#     shell:
#       """
#       bgzip {input.vcf}
#       tabix -p vcf {output.vcf}
#       """
