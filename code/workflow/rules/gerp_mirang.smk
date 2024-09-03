"""
snakemake --rerun-triggers mtime -n -R gerp_mirang

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
        -R gerp_mirang
"""

MIR_SCF_NAMESS = [ "NC_072369.1", "NC_072360.1", "NC_072372.1", "NC_072367.1",
                   "NC_072358.1", "NC_072359.1", "NC_072357.1", "NC_072371.1", 
                   "NC_072368.1", "NC_072362.1", "NC_072366.1", "NC_072370.1",
                   "NC_072363.1", "NC_072361.1", "NC_072364.1", "NC_072365.1" ]

MIR_SCF_NRS = [ str(x + 1).zfill(2) for x in range(16) ]

def get_mir_scaf_name(wildcards):
    return(MIR_SCF_NAMESS[int(wildcards.mscaf)-1])

rule gerp_mirang:
    input:
      ""

rule hal_to_maf_mirang:
    input:
      hal = '../results/cactus/{name}.hal'.format(name = P_NAME)
    output:
      maf = "../results/mirang/maf/mirang_scaf_{mscaf}.maf"
    log: "logs/mirang_hal_to_maf_{mscaf}.log"
    params:
      sif = c_cactus,
      js = "../results/cactus/scratch/{name}/".format(name = P_NAME),
      local_js = "js_mirang_{mscaf}",
      run = "run_mirang_{mscaf}",
      scaf_name = lambda wc: get_mir_scaf_name(wc)
    shell:
      """
      readonly CACTUS_IMAGE={params.sif} 
      readonly CACTUS_SCRATCH={params.js}\

      mkdir -p {params.js}{params.run}

      apptainer exec --cleanenv \
        --fakeroot --overlay ${{CACTUS_SCRATCH}} \
        --bind ${{CACTUS_SCRATCH}}/tmp:/tmp,{params.js}{params.run}:/run,$(pwd),{s_bind_paths} \
        --env PYTHONNOUSERSITE=1 \
        {params.sif} \
        cactus-hal2maf \
        /tmp/{params.local_js} \
        {input.hal} \
        {output.maf} \
        --refGenome mirang \
        --refSequence {params.scaf_name} \
        --dupeMode single \
        --filterGapCausingDupes \
        --chunkSize 1000000 \
        --noAncestors 2> {log}
      """

# we need to remove hits from scaffolds other than the target
# one that made it into the maf file
rule reformat_maf_mirang:
    input:
      maf = "../results/mirang/maf/mirang_scaf_{mscaf}.maf",
      genome = "../data/mirang/mirang_filt.fa.gz.fai"
    output:
      scf_bed = "../results/mirang/bed/mascaf_a1_{mscaf}.bed",
      maf = "../results/mirang/maf/mirang_no_own_hits_{mscaf}.maf" 
    params:
      ref_spec = "mirang",
      scaf_name = lambda wc: get_mir_scaf_name(wc)
    conda: "biopython"
    shell:
      """
      grep "{params.scaf_name}" {input.fai} | \
        awk '{{print $1"\t0\t"$2}}' > {output.scf_bed}
      
      py/intersect_maf_bed \
        --maf {input.maf} \
        --bed {output.scf_bed} \
        --ref {params.ref_spec} \
        --output /dev/stdout \
        --min_overlap_length 1 > {output.maf}
      """

rule call_gerp_mirang:
    input:
      maf = "../results/mirang/maf/mirang_no_own_hits_{mscaf}.maf",
      tree = "../results/neutral_tree/rerooted.tree"
    output:
      rates = "../results/mirang/gerp/mirang_{mscaf_nr}.maf.rates"
    params:
      refname = "mirang" 
    conda: "msa_phast"
    log: "logs/mirang_gerp_{mscaf_nr}.log"
    shell:
      """
      gerpcol -t {input.tree} -f {input.maf} -e {params.refname} -j -z -x ".rates" &> {log}

      mv {input.maf}.rates {output.rates}
      """

# bed coordinates need to be shifted by 1 bp,
# because bed is "left-open".
# gerp rates start at position 0 (-z parameter),
# but maf files also are also 0-indexed
rule parse_gerp_beds_mirang:
    input:
      rates = "../results/mirang/gerp/mirang_{mscaf}.maf.rates"
    output:
      bed = "../results/mirang/gerp/mirang_{mscaf}_gerp.bed.gz"
    conda: "popgen_basics"
    params:
      col_idx = 2,
      scaf_name = lambda wc: get_mir_scaf_name(wc)
    shell:
      """
      awk -v s="{params.scaf_name}" \
        '{{ print s"\t"NR-1"\t"NR"\t"${params.col_idx} }}' {input.rates} | \
        grep -v "\s0$" | \
        bgzip > {output.bed}
      """

# similarely to the coverage, the original bed is single bp elements,
# so we collapse continous chunks of equal gerp RS 
rule collapse_gerp_bed_mirang:
    input:
      bed = "../results/mirang/gerp/mirang_{mscaf}_gerp.bed.gz"
    output:
      bed = "../results/mirang/gerp/mirang_{mscaf}_gerp.collapsed.bed.gz"
    log: "logs/collapse_gerp_bed_{mscaf}.log"
    container: c_conda
    conda: "r_tidy"
    shell:
      """
      Rscript --vanilla R/collapse_bed_gerp.R {input.bed} {output.bed} &>> {log}
      """

rule merge_all_gerp_beds_mirang:
    input:
      beds = expand( "../results/mirang/gerp/mirang_{mscaf}_gerp.collapsed.bed.gz", mscaf = MIR_SCF_NRS )
    output:
      bed = "../results/mirang/gerp/mirang_RS.bed.gz"
    conda: "popgen_basics"
    shell:
      """
      zcat {input.beds} | bgzip > {output.bed}
      """