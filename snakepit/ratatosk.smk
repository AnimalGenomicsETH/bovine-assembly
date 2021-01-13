import pysam
import glob
from pathlib import Path,PurePath
import shutil

localrules: ratatosk_split_bins, ratatosk_merge_bin1, ratatosk_finish

OUT_PREFIX = 'rata_test'
Path(OUT_PREFIX).mkdir(exist_ok=True)

for _dir in ['segments','ratatosk']:
    Path(f'{OUT_PREFIX}/{_dir}').mkdir(exist_ok=True)

correct1_threads = 16
rule all:
    input:
        #OUT_PREFIX+'/bin_split'
        OUT_PREFIX+'/ratatosk/sample_corrected.fastq'

rule map_long_reads:
    input:
        reads = expand('data/offspring.read_R{N}.SR.fq.gz', N = (1,2)),
        asm = 'asm.scaffolds.fasta'
    output:
        'LR.bam'
    threads: 24
    resources:
        mem_mb = 6000,
        walltime = '4:00',
        disk_scratch = 250
    shell:
        'minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'


rule map_SR_reads:
    input:
        reads = expand('data/offspring.read_R{N}.SR.fq.gz', N = (1,2)),
        asm = 'asm.scaffolds.fasta'
    output:
         'SR.bam'
    threads: 24
    resources:
        mem_mb = 6000,
        walltime = '4:00',
        disk_scratch = 250
    shell:
        'minimap2 -ax sr -t {threads} {input.asm} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule index_bam:
    input:
        '{bam}.bam'
    output:
        '{bam}.bam.bai'
    threads: 8
    resources:
        mem_mb = 4000
    shell:
        'samtools index -@ {threads} {input}'

rule ratatosk_segment_bam:
    input:
        SR = 'SR.bam',
        SR_ind = 'SR.bam.bai',
        LR = 'LR.bam',
        LR_ind = 'LR.bam.bai'
    output:
        out = OUT_PREFIX+'/segments/sample.bin',
        long_unknown = OUT_PREFIX+'/segments/sample_lr_unknown.fq'
    params:
        lambda wildcards, output: PurePath(output['out']).with_suffix('')
    threads: 12
    resources:
        mem_mb = 6000
    shell:
        'python3 /cluster/work/pausch/alex/software/Ratatosk/script/reference_guiding/segmentBAM.py -t {threads} -s {input.SR} -l {input.LR} -o {params} > {output}'

checkpoint ratatosk_make_bins:
#rule ratatosk_make_bins:
    input:
        SR = 'SR.bam'
    output:
        script = directory(OUT_PREFIX+'/bin_split')
    params:
        threads = correct1_threads,
        seg_prefix = f'{OUT_PREFIX}/segments',
        unmapped = f'{OUT_PREFIX}/segments/sample_sr_unmapped.fa',
        BIN_SIZE = 5000000
    run:
        bamf = pysam.AlignmentFile(input.SR, "rb")
        Path(output.script).mkdir(exist_ok=True)
        for chr_name, chr_length in zip(bamf.references,bamf.lengths):
            for position in range(0,chr_length+1,params.BIN_SIZE):
                NAME_LR_IN_FILE = params.seg_prefix + f'/sample_lr_{chr_name}_{position}'
                NAME_SR_IN_FILE = params.seg_prefix + f'/sample_sr_{chr_name}_{position}.fa'
                NAME_LR_OUT_FILE = f'{NAME_LR_IN_FILE}_corrected'
                NAME_LR_IN_FILE += '.fq'

                if Path(NAME_LR_IN_FILE).exists() and Path(NAME_LR_IN_FILE).stat().st_size > 0:
                    if Path(NAME_SR_IN_FILE).exists() and Path(NAME_SR_IN_FILE).stat().st_size > 0:
                        with open(output.script+f'/bin_{chr_name}_{position}.sh','w') as fout:
                            fout.write(f'Ratatosk -v -c {params.threads} -s {NAME_SR_IN_FILE} -l {NAME_LR_IN_FILE} -u {params.unmapped} -o {NAME_LR_OUT_FILE}')
                    else:
                        shutil.copy(NAME_LR_IN_FILE,f'{NAME_LR_OUT_FILE}.fastq')

rule ratatosk_correct_bin1:
    input:
        OUT_PREFIX+'/bin_split/bin_{N}.sh'
    output:
        OUT_PREFIX+'/segments/sample_lr_{N}_corrected.fastq'
    threads: 4 #correct1_threads
    resources:
        mem_mb = 2000,
        walltime = '3:00'
    #group: 'correct1'
    shell:
        'bash {input}'


def aggregate_corrected_bin1(wildcards):
    checkpoint_output = checkpoints.ratatosk_make_bins.get(**wildcards).output[0]
    return expand(OUT_PREFIX+'/segments/sample_lr_{chunk}_corrected.fastq',fpath=checkpoint_output,chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('bin_{chunk}.sh')).chunk)

rule ratatosk_merge_bin1:
    input:
        aggregate_corrected_bin1
    output:
        OUT_PREFIX+'/segments/sample_lr_map.fastq'
    params:
        glob.glob(OUT_PREFIX+'/segments/sample_lr_*_corrected.fastq')
    shell:
        'cat {input} > {output}'

rule ratatosk_get_SR_fastq:
    input:
        SR = 'SR.bam'
    output:
        OUT_PREFIX+'/segments/sample_sr.fastq.gz'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'samtools bam2fq -@ {threads} -n {input.SR} > {output}'

rule ratatosk_correct_bin2:
    input:
        short_reads = OUT_PREFIX+'/segments/sample_sr.fastq.gz',
        long_unknown = OUT_PREFIX+'/segments/sample_lr_unknown.fq',
        long_mapped = OUT_PREFIX+'/segments/sample_lr_map.fastq'
    output:
        OUT_PREFIX+'/segments/sample_lr_unknown_corrected.fastq'
    params:
        lambda wildcards, output: PurePath(output[0]).with_suffix(''),
    threads: 36
    resources:
        mem_mb = 8000,
        walltime = '24:00'
    shell:
        'Ratatosk -v -c {threads} -s {input.short_reads} -l {input.long_unknown} -a {input.long_mapped} -o {params}'

rule ratatosk_finish:
    input:
        mapped = OUT_PREFIX+'/segments/sample_lr_map.fastq',
        unknown = OUT_PREFIX+'/segments/sample_lr_unknown_corrected.fastq'
    output:
        OUT_PREFIX+'/ratatosk/sample_corrected.fastq'
    shell:
        'cat {input} > {output}'
