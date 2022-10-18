import os
import json
import pandas as pd

def _get_runinfo(rna_sample):
    paths_path = "/juno/work/shah/users/chois7/apollo/APOLLO.KALLISTO.tsv"
    assert os.path.exists(paths_path)
    paths = pd.read_table(paths_path)
    paths = paths[paths['isabl_sample_id']==rna_sample]
    paths = paths[ (paths["result_type"] == "run_info")]['result_filepath'].values
    assert len(paths) == 1, f'run_info paths length not 1 for {rna_sample}: {paths}'
    path = paths[0]
    assert os.path.exists(path), f'run_info path: {path} does not exist'
    runinfo = json.loads(open(path, 'r').read())
    return runinfo

def _get_R1(wildcards):
    if wildcards.sample == 'testRNA': return config['test_R1']
    runinfo = _get_runinfo(wildcards.sample)
    R1_list = list(filter(lambda x: ('_R1_' in x), runinfo['call'].split(' ')))
    assert len(R1_list) == 1, f'R1_list = {R1_list}'
    fastq = R1_list[0]
    return fastq
    
def _get_R2(wildcards):
    if wildcards.sample == 'testRNA': return config['test_R2']
    runinfo = _get_runinfo(wildcards.sample)
    R2_list = list(filter(lambda x: ('_R2_' in x), runinfo['call'].split(' ')))
    assert len(R2_list) == 1, f'R2_list = {R2_list}'
    fastq = R2_list[0]
    return fastq


rule star:
    input:
        R1=_get_R1, # forward
        R2=_get_R2, # reverse
    output:
        bam=temp('main_run/{patient}/{sample}/outputs/star/{rna_sample}.Aligned.out.bam'),
    params:
        star_ref_dir=config['star']['ref_dir'],
        threads=config['star']['threads'],
        out_prefix='main_run/{patient}/{sample}/outputs/star/{rna_sample}.',
    singularity:
        'docker://dceoy/star'
    threads: config['star']['threads'],
    shell:
        'touch {output.bam} && '
        'STAR '
        '--genomeDir {params.star_ref_dir} '
        '--readFilesIn {input.R1} {input.R2} ' # forward reverse
        '--outFileNamePrefix {params.out_prefix} '
        '--runThreadN {params.threads} ' 
        '--readFilesCommand "gunzip -c" ' # how to read input
        '--outSAMtype BAM Unsorted ' # output file format (unsorted BAM file)
        '--outSAMmapqUnique 60 ' # mapQ of unique mapped reads
        '--sjdbOverhang 100 ' # sequence length near junctions when making junction-DB
        '--twopassMode Basic ' # STAR 2-pass mode
        '--limitOutSJcollapsed 1000000 ' # max number of collapsed junctions

rule add_or_replace_rg:
    input:
        bam='main_run/{patient}/{sample}/outputs/star/{rna_sample}.Aligned.out.bam',
    output:
        bam=temp('main_run/{patient}/{sample}/outputs/add_or_replace_rg/{rna_sample}.sorted.bam'),
    params:
        tmp_dir=config['tmp_dir'],
    singularity:
        'docker://broadinstitute/gatk',
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
        "AddOrReplaceReadGroups "
        '-I {input.bam} -O {output.bam} '
        '--CREATE_INDEX true '
        '-SO coordinate ' # sort order
        '-PL illumina '
        '-LB {wildcards.rna_sample} '
        '-PU {wildcards.rna_sample} '
        '-SM {wildcards.rna_sample} ' # fixed post-hoc

rule mark_duplicates:
    input:
        bam='main_run/{patient}/{sample}/outputs/add_or_replace_rg/{rna_sample}.sorted.bam',
    output:
        bam=temp('main_run/{patient}/{sample}/outputs/mark_dupplicates/{rna_sample}.dedup.bam'),
        metrics='main_run/{patient}/{sample}/outputs/mark_dupplicates/{rna_sample}.dedup.bam.metrics',
    params:
        tmp_dir=config['tmp_dir'],
    singularity:
        'docker://broadinstitute/gatk',
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
        "MarkDuplicates "
        '-I {input.bam} -O {output.bam} -M {output.metrics} '
        '--CREATE_INDEX true '
        '--ASSUME_SORT_ORDER coordinate '

rule split_n_cigar_reads:
    input:
        bam='main_run/{patient}/{sample}/outputs/mark_dupplicates/{rna_sample}.dedup.bam',
    output:
        bam=temp('main_run/{patient}/{sample}/outputs/split_n_cigar_reads/{rna_sample}.splitncigar.bam'),
    params:
        tmp_dir=config['tmp_dir'],
        reference_fasta=config['reference_fasta'],
    singularity:
        'docker://broadinstitute/gatk',
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
        'SplitNCigarReads '
        '-R {params.reference_fasta} '
        '-I {input.bam} -O {output.bam} '

rule base_recalibrator:
    input:
        bam='main_run/{patient}/{sample}/outputs/split_n_cigar_reads/{rna_sample}.splitncigar.bam',
    output:
        table='main_run/{patient}/{sample}/outputs/recalibration/{rna_sample}.recal.table',
    params:
        tmp_dir=config['tmp_dir'],
        reference_fasta=config['reference_fasta'],
        known_snp_sites=config['known_snp_sites'],
    singularity:
        'docker://broadinstitute/gatk',
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
        'BaseRecalibrator '
        '-R {params.reference_fasta} '
        '-I {input.bam} -O {output.table} '
        '--known-sites {params.known_snp_sites} '

rule apply_bqsr:
    input:
        bam='main_run/{patient}/{sample}/outputs/split_n_cigar_reads/{rna_sample}.splitncigar.bam',
        table='main_run/{patient}/{sample}/outputs/recalibration/{rna_sample}.recal.table',
    output:
        bam='main_run/{patient}/{sample}/outputs/recalibration/{rna_sample}.bam',
    params:
        tmp_dir=config['tmp_dir'],
        reference_fasta=config['reference_fasta'],
    singularity:
        'docker://broadinstitute/gatk',
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
        'ApplyBQSR '
        '-R {params.reference_fasta} '
        '-I {input.bam} '
        '--bqsr-recal-file {input.table} '
        '-O {output.bam} '
    
rule haplotypecaller:
    input:
        bam='main_run/{patient}/{sample}/outputs/recalibration/{rna_sample}.bam',
    output:
        vcf='main_run/{patient}/{sample}/outputs/haplotypecaller/{rna_sample}.vcf',
    params:
        tmp_dir=config['tmp_dir'],
        reference_fasta=config['reference_fasta'],
    singularity:
        'docker://broadinstitute/gatk',
    shell:
        "gatk --java-options '-Xmx4g -Djava.io.tmpdir={params.tmp_dir}' "
        'HaplotypeCaller '
        '-R {params.reference_fasta} '
        '-I {input.bam} '
        '-O {output.vcf} '

rule vt_decompose_rna:
    input:
        vcf='main_run/{patient}/{sample}/outputs/haplotypecaller/{rna_sample}.vcf',
    output:
        vcf='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}.decomposed.vcf',
    params:
        vcf='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}.sample_fixed.vcf',
    singularity:
        "docker://zlskidmore/vt",
    shell:
        'sed "s/{wildcards.sample}/{wildcards.rna_sample}/g" {input.vcf} > {params.vcf} && '
        'vt decompose '
        '-s {params.vcf} ' # fix to input.vcf
        '-o {output.vcf}' 

rule bam_readcount_rna:
    input:
        vcf='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}.decomposed.vcf',
        bam='main_run/{patient}/{sample}/outputs/recalibration/{rna_sample}.bam',
        ref=config['reference_fasta'],
    output:
        snv_tsv='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}_bam_readcount_snv.tsv',
        indel_tsv='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}_bam_readcount_indel.tsv',
    params:
        outdir='main_run/{patient}/{sample}/outputs/rna_readcount'
    singularity:
        '/juno/work/shah/users/chois7/apollo/neoantigen/bam-readcount_helper_v0.0.3.sif',
    shell:
        'python /usr/bin/bam_readcount_helper.py ' 
        '{input.vcf} ' 
        '{wildcards.rna_sample} ' # wrongly mapped to DNA sample... searches sample with pyvcf
        '{input.ref} ' 
        '{input.bam} '
        '{params.outdir}'
    
rule vcf_readcount_annotator_snv_rna:
    input:
        vcf='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}.decomposed.vcf',
        snv_tsv='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}_bam_readcount_snv.tsv',
    output:
        vcf='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}.decomposed.snv_annotated.vcf',
    singularity:
        'docker://griffithlab/vatools',
    shell:
        'vcf-readcount-annotator '
        '{input.vcf} {input.snv_tsv} '
        'RNA -s {wildcards.rna_sample} '
        '-t snv -o {output.vcf}'

rule vcf_readcount_annotator_indel_rna:
    input:
        vcf='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}.decomposed.snv_annotated.vcf',
        indel_tsv='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}_bam_readcount_indel.tsv',
    output:
        vcf='main_run/{patient}/{sample}/outputs/rna_readcount/{rna_sample}.decomposed.both_annotated.vcf',
    singularity:
        'docker://griffithlab/vatools',
    shell:
        'vcf-readcount-annotator '
        '{input.vcf} {input.indel_tsv} '
        'RNA -s {wildcards.rna_sample} '
        '-t indel -o {output.vcf}'

