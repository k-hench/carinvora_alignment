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

localrules: merge_mafs

rule all_phylop:
    input: 
      maf = "../results/neutral_tree/windows/autosomes.maf.gz",
      model = "../results/phylop/autosomes_neutral.mod"

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

rule phylop_model:
    input:
      maf = "../results/neutral_tree/windows/autosomes.maf.gz",
      tree = "../results/neutral_tree/rerooted.tree"
    output:
      model = "../results/phylop/autosomes_neutral.mod"
    log:
      "logs/phylofit.log"
    params:
      overlay = config[ 'img_phylop' ],
      bind_paths = config[ 'singularity_bind_paths' ]
    shell:
      """
      PATH_UPDATE="PATH=/home/cactus/cactus_env/bin:/home/cactus/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/conda/envs/phast/bin"

      apptainer exec \
        --fakeroot --overlay {params.overlay}:ro \
        --bind {params.bind_paths} \
        --env "PATH=${{PATH_UPDATE}}" \
         {c_cactus} \
         phyloFit \
           --tree {input.tree} \
           --msa-format MAF \
           --log {log} \
           {input.maf} > {output.model}
      """
