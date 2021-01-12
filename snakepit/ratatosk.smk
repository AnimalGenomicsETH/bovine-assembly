import pysam
import glob

OUT_PREFIX = 'rata_test'
Path(OUT_PREFIX).mkdir(exist_ok=True)

for _dir in ['segments','ratatosk']:
    Path(f'{OUT_PREFIX}/{_dir}').mkdir(exist_ok=True)

rule all:
    input:
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
        OUT_PREFIX+'/segments'
    threads: 12
    resources:
        mem_mb = 6000
    shell:
        'python3 /cluster/work/pausch/alex/software/Ratatosk/script/reference_guiding/segmentBAM.py -t {threads} -s {input.SR} -l {input.LR} -o {output}/sample > {output}/sample.bin'


rule ratatosk_make_bins:
    input:
        SR = 'SR.bam'
    output:
        script = OUT_PREFIX+'/bin.sh'
    params:
        threads = 16,
        seg_prefix = f'{OUT_PREFIX}/segments',
        unmapped = f'{OUT_PREFIX}/segments/sample_sr_unmapped',
        BIN_SIZE = 5000000
    run:
        NAME_SR_UNMAPPED_IN_FILE = params.seg_prefix + '/sample_sr_unmapped'

        bamf = pysam.AlignmentFile(input.SR, "rb")
        with open(output.script,'w') as fout:
            for chr_name, chr_length in zip(bamf.references,bam.lengths):
                for position in range(0,chr_length+1,params.BIN_SIZE):
                    NAME_LR_IN_FILE = params.seg_prefix + f'/sample_lr_{chr_name}_{position}'
                    NAME_SR_IN_FILE = params.seg_prefix + f'/sample_sr_{chr_name}_{position}'
                    NAME_LR_OUT_FILE = f'{NAME_LR_IN_FILE}_corrected'

                    fout.write(f'if [ -f {NAME_LR_IN_FILE}.fq ] && [ -s {NAME_LR_IN_FILE}.fq ]; then if [ -f {NAME_SR_IN_FILE}.fa ] && [ -s {NAME_SR_IN_FILE}.fa ]; then')
                    fout.write(f'Ratatosk -v -c {params.threads} -s {NAME_SR_IN_FILE}.fa -l {NAME_LR_IN_FILE}.fq -u {params.unmapped}.fa -o {NAME_LR_OUT_FILE}')
                    fout.write(f'else cp {NAME_LR_IN_FILE}.fq {NAME_LR_OUT_FILE}.fastq; fi; fi;')

checkpoint ratatosk_split_bins:
    input:
        OUT_PREFIX+'/bin.sh'
    output:
        directory(OUT_PREFIX+'/bin_split')
    params:
        sharding = 100
    '''
    mkdir -p {output} && cd {output}
    split --number=l/{params.sharding} --additional-suffix=.sh {input} {output}
    '''


rule ratatosk_correct_bin1:
    input:
        'bin_split/chunk{N}.sh'
    output:
        'bin_split/chunk{N}.out'
    threads: 24
    resources:
        mem_mb = 400,
        walltime = '4:00'
    shell:
        './{input} > {output}'


def aggregate_corrected_bin1(wildcards):
    checkpoint_output = checkpoints.ratatosk_split_bins.get(**wildcards).output[0]
    return expand('{fpath}/{chunk}.out',fpath=checkpoint_output,chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('{chunk}.sh')).chunk)

rule ratatosk_merge_bin1:
    input:
        aggregate_corrected_bin1
    output:
        OUT_PREFIX+'/segments/sample_lr_map.fastq'
    params:
        OUT_PREFIX+'/segments'
    shell:
        'cat \$(ls -t {params}/sample_lr_*_corrected.fastq) > {output}'

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
        OUT_PREFIX+'/segments/sample_lr_unknown_corrected'
    threads: 36
    resources:
        mem_mb = 8000,
        walltime = '24:00'
    shell:
        'Ratatosk -v -c {threads} -s {input.short_reads} -l {input.long_unknown} -a {input.long_mapped} -o {output}'

rule ratatosk_finish:
    input:
        mapped = OUT_PREFIX+'/segments/sample_lr_map.fastq',
        unknown = OUT_PREFIX+'/segments/sample_lr_unknown_corrected.fastq'
    output:
        OUT_PREFIX+'/ratatosk/sample_corrected.fastq'
    shell:
        'cat {input} > {output}'
