import os
import pandas as pd

cohort_file = config['cohort_file']
cohort = pd.read_table(cohort_file) # patient, sample, rna_sample w/o NA
#cohort = cohort[cohort['sample']=='ADT167R2']
   
vcfs=[f'main_run/{patient}/{sample}/outputs/haplotypecaller/{rna_sample}.vcf' 
            for rix, (patient, sample, rna_sample) in cohort.iterrows()]
print(vcfs)
rule all:
    input:
        #[f'main_run/{patient}/{sample}/outputs/pvacseq/MHC_Class_I/{sample}.filtered.tsv' 
        #    for rix, (patient, sample, rna_sample) in cohort.iterrows()] + \
        #[f'main_run/{patient}/{sample}/outputs/pvacseq/MHC_Class_I/{sample}.expression_filtered.tsv' 
        #    for rix, (patient, sample, rna_sample) in cohort.iterrows()] + \
        [f'main_run/{patient}/{sample}/outputs/haplotypecaller/{rna_sample}.vcf' 
            for rix, (patient, sample, rna_sample) in cohort.iterrows()]

include: "rules/pvacseq.smk"
include: "rules/rnaseq_variantcalling.smk"
