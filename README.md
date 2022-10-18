# Bulk RNA-seq variant discovery pipeline
Follows Broad/GATK best practices as in [link](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-)

## Prerequisites
- You need `snakemake` to be installed
- Data will be stored in `./main_run`, but you can change the `Snakefile` and `rules/*.smk` accordingly to change your results path.

## Running on APOLLO cohort
```bash
bash run_snakemake.sh
```
