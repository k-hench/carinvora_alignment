"""
snaekmake --rerun-triggers mtime -n all_phylop

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
      maf = "../results/neutral_tree/windows/autosomes.maf.gz"

rule merge_mafs:
    input:
      mafs = expand( "../results/neutral_tree/windows/{mscaf}.maf.gz", mscaf = SCFS[:17] )
    output:
      maf = "../results/neutral_tree/windows/autosomes.maf.gz"
    params:
      prefix = "../results/neutral_tree/windows/autosomes.maf"
    shell:
      """
      zcat {input.mafs[0]} | head -n 1 > {params.prefix}
      for k in {input.mafs}; do
        zcat ${{k}} | tail -n +2 >> {params.prefix}
      done

      gzip {params.prefix}
      """
