import pysam
from pathlib import Path,PurePath
import shutil

localrules: ratatosk_make_bins, ratatosk_merge_bin1, ratatosk_finish

segmentBAM_path = '/cluster/work/pausch/alex/software/Ratatosk/script/reference_guiding/segmentBAM.py'
OUT_PREFIX = 'rata_test'
out_path = Path(OUT_PREFIX)

seg_path = out_path / 'segments'
res_path = out_path / 'ratatosk'

for path in (out_path,seg_path,res_path):
    path.mkdir(exist_ok=True)

correct1_threads = 16
rule all:
    input:
        #OUT_PREFIX+'/bin_split'
        res_path / 'sample_corrected.fastq'

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
        long_unknown = seg_path / 'sample_lr_unknown.fq'
    params:
        out = seg_path / 'sample',
        seg_script = segmentBAM_path
    threads: 12
    resources:
        mem_mb = 6000
    shell:
        'python3 {parmas.seg_script} -t {threads} -s {input.SR} -l {input.LR} -o {params.out} > {params.out}.bin'

checkpoint ratatosk_make_bins:
#rule ratatosk_make_bins:
    input:
        bam = 'SR.bam'
    output:
        script = directory(OUT_PREFIX+'/bin_split')
    params:
        threads = correct1_threads,
        unmapped = seg_path / 'sample_sr_unmapped.fa',
        BIN_SIZE = 5000000
    run:
        bamf = pysam.AlignmentFile(input.bam "rb")
        directory = Path(output.script)
        directory.mkdir(exist_ok=True)
        for chr_name, chr_length in zip(bamf.references,bamf.lengths):
            for position in range(0,chr_length+1,params.BIN_SIZE):
                short_read_in = seg_path / f'sample_sr_{chr_name}_{position}.fa'
                long_read_in = seg_path / f'sample_lr_{chr_name}_{position}'
                long_read_out = str(long_read_in) + '_corrected'
                long_read_in = long_read_in.with_suffix('.fq')
                short_read_in, long_read_in, long_read_out = (seg_path / f'sample_{read}_{chr_name}_{position}{ext}' for read,ext in zip(('sr','fa'),('lr','fq'),'lr','_corrected.fastq'))

                if long_read_in.exists() and long_read_in.stat().st_size > 0:
                    if short_read_in.exists() and short_read_in.stat().st_size > 0:
                        with open(directory / f'bin_{chr_name}_{position}.sh','w') as fout:
                            fout.write(f'Ratatosk -v -c {params.threads} -s {short_read_in} -l {long_read_in} -u {params.unmapped} -o {long_read_out.with_suffix("")}')
                    else:
                        shutil.copy(long_read_in,long_read_out)

rule ratatosk_correct_bin1:
    input:
        out_path / 'bin_split/bin_{N}.sh'
    output:
        temp(seg_path / 'sample_lr_{N}_corrected.fastq')
    threads: 4 #correct1_threads
    resources:
        mem_mb = 2000,
        walltime = '3:00'
    #group: 'correct1'
    shell:
        'bash {input}'


def aggregate_corrected_bin1(wildcards):
    checkpoint_output = checkpoints.ratatosk_make_bins.get(**wildcards).output[0]
    return expand(seg_path / 'sample_lr_{chunk}_corrected.fastq',fpath=checkpoint_output,chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('bin_{chunk}.sh')).chunk)

rule ratatosk_merge_bin1:
    input:
        aggregate_corrected_bin1
    output:
        seg_path / 'sample_lr_map.fastq'
    shell:
        'cat {input} > {output}'

rule ratatosk_get_SR_fastq:
    input:
        SR = 'SR.bam'
    output:
        temp(seg_path / 'sample_sr.fastq.gz')
    threads: 4
    resources:
        mem_mb = 10000,
        walltime = '1:00'
    shell:
        'samtools bam2fq -@ {threads} -n {input.SR} > {output}'

rule ratatosk_correct_bin2_p1:
    input:
        short_reads = seg_path / 'sample_sr.fastq.gz',
        long_unknown = seg_path / 'sample_lr_unknown.fq',
        long_mapped = seg_path / 'sample_lr_map.fastq'
    output:
        seg_path / 'sample_lr_unknown_corrected2.fastq'
    params:
        seg_path / 'sample_lr_unknown_corrected'
    threads: 36
    resources:
        mem_mb = 8000,
        walltime = '24:00'
    shell:
        'Ratatosk -1 -v -c {threads} -s {input.short_reads} -l {input.long_unknown} -a {input.long_mapped} -o {params}'

checkpoint ratatosk_shard_bin2:
    input:
        long_unknown = seg_path / 'sample_lr_unknown.fq',
    output:
        temp(directory(seg_path / 'unknown_shards'))
    params:
    shell:
        '''
        mkdir -p {output}
        split -a 2 -d -l {params} --additional_suffix=.fq {input} {output}/shard_
        '''

rule ratatosk_correct_bin2_p2:
    input:
        short_reads = seg_path / 'sample_sr.fastq.gz',
        long_unknown = seg_path / 'sample_lr_unknown.fq',
        long_mapped = seg_path / 'sample_lr_map.fastq',
        p1_correction = seg_path / 'sample_lr_unknown_corrected2.fastq',
        lr_shard = seg_path / 'unknown_shards/shard_{N}.fq'
    output:
        seg_path / 'shard_{N}_corrected.fastq'
    params:
        seg_path / 'sample_lr_unknown_corrected'
    threads: 36
    resources:
        mem_mb = 5750,
        walltime = '24:00'
    shell:
        'Ratatosk -2 -v -c {threads} -s {input.short_reads} -l {input.long_unknown} -a {input.long_mapped} -o {params} -L {input.lr_shard} -O {output}'


def aggregate_corrected_bin2(wildcards):
    checkpoint_output = checkpoints.ratatosk_shard_bin2.get(**wildcards).output[0]
    return expand(seg_path / 'shard_{chunk}_corrected.fastq',fpath=checkpoint_output,chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('shard_{chunk}.fq')).chunk)

rule ratatosk_merge_bin2:
    input:
        aggregate_corrected_bin2
    output:
        seg_path / 'sample_lr_unknown_corrected.fastq'
    shell:
        'cat {input} > {output}'

rule ratatosk_finish:
    input:
        mapped = seg_path / 'sample_lr_map.fastq',
        unknown = seg_path / 'sample_lr_unknown_corrected.fastq'
    output:
        res_path / 'sample_corrected.fastq'
    shell:
        'cat {input} > {output}'
