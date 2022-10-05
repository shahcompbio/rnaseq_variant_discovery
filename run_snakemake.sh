CLUSTER_CMD=("bsub -n {threads} -R {cluster.resources} -M {cluster.memory} -o {cluster.output} -J {cluster.name} -W {cluster.time}") 
snakemake --configfile config.yaml \
    --use-singularity \
    --singularity-args "--bind /juno --bind /rtsess01" \
    --jobs 30 \
    --skip-script-cleanup \
    --cluster-config ./cluster.yaml \
    --cluster "${CLUSTER_CMD}" \
    --cluster-cancel bkill \
    --allowed-rules vt_decompose_rna bam_readcount_rna vcf_readcount_annotator_snv_rna vcf_readcount_annotator_indel_rna \
    --rerun-incomplete #--dry-run
