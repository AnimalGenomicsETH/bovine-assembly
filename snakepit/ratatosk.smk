from pathlib import Path,PurePath
import shutil

localrules: ratatosk_make_bins, ratatosk_merge_bin1, ratatosk_shard_bin2, ratatosk_merge_bin2, ratatosk_finish

if Path('config/ratatosk.yaml').exists():
    configfile: 'config/ratatosk.yaml'

#wildcard_constraints:
#    N = r'\d+'

class Default(dict):
    def __missing__(self, key):
        return '{'+key+'}'

def get_dir(base='segments',ext='', **kwargs):
    if base == 'result':
        base_dir = get_dir('main','corrected',**kwargs)
    elif base == 'main':
        base_dir = 'ratatosk_{animal}'
    elif base == 'segments':
        base_dir = get_dir('main','segments',**kwargs)
    else:
        raise Exception('Base not found')
    return str(Path(base_dir.format_map(Default(kwargs))) / ext)


for path in ('main','result','segments'):
    Path(get_dir(path,animal=config['animal'])).mkdir(exist_ok=True)

rule all:
    input:
        expand(get_dir('result','sample_corrected{ambiguous}.fastq',animal=config['animal']),ambiguous=('','_unambiguous'))

rule map_long_reads:
    input:
        reads = config['long_reads'],
        asm = config['reference']
    output:
        temp(get_dir('main','long_reads.bam'))
    threads: 12
    resources:
        mem_mb = 12000,
        walltime = '24:00',
        disk_scratch = 250
    shell:
        'minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} | samtools sort - -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule map_SR_reads:
    input:
        reads = config['short_reads'],
        asm = config['reference']
    output:
         temp(get_dir('main','short_reads.bam'))
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
        SR = multiext(get_dir('main','short_reads.bam'),'','.bai'),
        LR = multiext(get_dir('main','long_reads.bam'),'','.bai')
    output:
        get_dir('segments','sample_lr_unknown.fq'),
        get_dir('segments','sample_sr_unmapped.fa')
    params:
        out = lambda wildcards, output: PurePath(output[0]).with_name('sample')
    threads: 12
    resources:
        mem_mb = 6000
    shell:
        'python3 {config[segment_py]} -t {threads} -s {input.SR[0]} -l {input.LR[0]} -o {params.out} > {params.out}.bin'

checkpoint ratatosk_make_bins:
    input:
        SR = multiext(get_dir('main','short_reads.bam'),'','.bai'),
        LR = multiext(get_dir('main','long_reads.bam'),'','.bai'),
        unknown = get_dir('segments','sample_lr_unknown.fq'),
        unmapped = get_dir('segments','sample_sr_unmapped.fa')
    output:
        scripts = directory(get_dir('main','bin_split'))
    params:
        seg_path = lambda wildcards, input: Path(input.unknown).parent,
        BIN_SIZE = 5000000
    run:
        import pysam
        bamf = pysam.AlignmentFile(input.SR[0], "rb")
        directory = Path(output.scripts)
        directory.mkdir(exist_ok=True)
        for chr_name, chr_length in zip(bamf.references,bamf.lengths):
            for position in range(0,chr_length+1,params.BIN_SIZE):
                short_read_in, long_read_in, long_read_out = (params.seg_path / f'sample_{read}_{chr_name}_{position}{ext}' for read,ext in (('sr','.fa'),('lr','.fq'),('lr','_corrected.fastq')))
                if long_read_in.exists() and long_read_in.stat().st_size > 0:
                    if short_read_in.exists() and short_read_in.stat().st_size > 0:
                        with open(directory / f'bin_{chr_name}_{position}.sh','w') as fout:
                            fout.write(f'Ratatosk -v -c {config["phase1_threads"]} -s {short_read_in} -l {long_read_in} -u {input.unmapped} -o {long_read_out.with_suffix("")}')
                    else:
                        shutil.copy(long_read_in,long_read_out)

rule ratatosk_correct_bin1:
    input:
        get_dir('main','bin_split/bin_{N}.sh')
    output:
        temp(get_dir('segments','sample_lr_{N}_corrected.fastq'))
    threads: config['phase1_threads']
    resources:
        mem_mb = 2000,
        walltime = '4:00'
    shell:
        'bash {input}'

def aggregate_corrected_bin1(wildcards):
    checkpoint_output = PurePath(checkpoints.ratatosk_make_bins.get(**wildcards).output[0])
    return expand(get_dir('segments','sample_lr_{chunk}_corrected.fastq',**wildcards),chunk=glob_wildcards(checkpoint_output.joinpath('bin_{chunk}.sh')).chunk)

rule ratatosk_merge_bin1:
    input:
        aggregate_corrected_bin1
    output:
        temp(get_dir('segments','sample_lr_map.fastq'))
    shell:
        'cat {input} > {output}'

rule ratatosk_get_SR_fastq:
    input:
        SR = multiext(get_dir('main','short_reads.bam'),'','.bai')
    output:
        temp(get_dir('segments','sample_sr.fastq.gz'))
    threads: 4
    resources:
        mem_mb = 10000,
        walltime = '1:00'
    shell:
        'samtools bam2fq -@ {threads} -n {input.SR[0]} > {output}'

rule ratatosk_correct_bin2_p1:
    input:
        short_reads = get_dir('segments','sample_sr.fastq.gz'),
        long_unknown = get_dir('segments','sample_lr_unknown.fq'),
        long_mapped = get_dir('segments','sample_lr_map.fastq')
    output:
        temp(get_dir('segments','sample_lr_unknown_corrected2.fastq'))
    params:
        lambda wildcards, output: PurePath(output[0]).with_name('sample_lr_unknown_corrected')
    threads: 28
    resources:
        mem_mb = 9000,
        walltime = '24:00'
    shell:
        'Ratatosk -1 -v -c {threads} -s {input.short_reads} -l {input.long_unknown} -a {input.long_mapped} -o {params}'

rule ratatosk_build_graph:
    input:
        short_reads = get_dir('segments','sample_sr.fastq.gz'),
        long_unknown = get_dir('segments','sample_lr_unknown.fq'),
        long_mapped = get_dir('segments','sample_lr_map.fastq'),
        p1_correction = get_dir('segments','sample_lr_unknown_corrected2.fastq'),
    output:
        temp(get_dir('main','direct_write.gfa'))
    params:
        out = lambda wildcards, input: PurePath(input.p1_correction).with_name('sample_lr_unknown_corrected'),
        gfa = lambda wildcards, output: PurePath(output[0]).with_suffix('')
    threads: 24
    resources:
        mem_mb = 10000,
        walltime = '24:00'
    shell:
        'Ratatosk -2 -X {params.gfa} -v -c {threads} -s {input.short_reads} -l {input.p1_correction} -a {input.long_mapped} -o {params.out}'

checkpoint ratatosk_shard_bin2:
    input:
        long_unknown = get_dir('segments','sample_lr_unknown.fq'),
    output:
        temp(directory(get_dir('segments','unknown_shards')))
    shell:
        '''
        mkdir -p {output}
        split -a 2 -d -l {config[max_lines_per_phase2]} --additional-suffix=.fq {input} {output}/shard_
        '''

rule ratatosk_correct_bin2_p2:
    input:
        graph = get_dir('main','direct_write.gfa'),
        short_reads = get_dir('segments','sample_sr.fastq.gz'),
        long_unknown = get_dir('segments','sample_lr_unknown.fq'),
        long_mapped = get_dir('segments','sample_lr_map.fastq'),
        p1_correction = get_dir('segments','sample_lr_unknown_corrected2.fastq'),
        lr_shard = get_dir('segments','unknown_shards/shard_{N}.fq')
    output:
        temp(get_dir('segments','shard_{N}_corrected.fastq'))
    params:
        lambda wildcards, output: PurePath(output[0]).with_name('sample_lr_unknown_corrected')
    threads: 12
    resources:
        mem_mb = 5000,
        walltime = '24:00'
    shell:
        'Ratatosk -2 -Y {input.graph} -v -c {threads} -s {input.short_reads} -l {input.long_unknown} -a {input.long_mapped} -o {params} -L {input.lr_shard} -O {output}'

def aggregate_corrected_bin2(wildcards):
    checkpoint_output = checkpoints.ratatosk_shard_bin2.get(**wildcards).output[0]
    return expand(get_dir('segments','shard_{chunk}_corrected.fastq',animal=config['animal']),chunk=glob_wildcards(PurePath(checkpoint_output).joinpath('shard_{chunk}.fq')).chunk)

rule ratatosk_merge_bin2:
    input:
        aggregate_corrected_bin2
    output:
        temp(get_dir('segments','sample_lr_unknown_corrected.fastq'))
    shell:
        'cat {input} > {output}'

rule ratatosk_finish:
    input:
        mapped = get_dir('segments','sample_lr_map.fastq'),
        unknown = get_dir('segments','sample_lr_unknown_corrected.fastq')
    output:
        get_dir('result','sample_corrected.fastq')
    shell:
        'cat {input} > {output}'

rule replace_ambiguous_bases:
    input:
        get_dir('result','sample_corrected.fastq')
    output:
        get_dir('result','sample_corrected_unambiguous.fastq')
    run:
        import re, random
        IUPAC = {'R':['A','G'],'Y':['C','T'],'S':['G','C'],'W':['A','T'],'K':['G','T'],'M':['A','C'],
        'B':['C','G','T'],'D':['A','G','T'],'H':['A','C','T'],'V':['A','C','G'],'N':['A','C','G','T']}
        regex = re.compile(r'({})'.format('|'.join(IUPAC)))
        with open(input[0],'r') as fin, open(output[0],'w') as fout:
            for i, line in enumerate(fin):
                if i%4 == 1: # fasta sequence line
                    fout.write(regex.sub(lambda match: random.choice(IUPAC[match.group(1)]), line) + '\n')
                else:
                    fout.write(line + '\n')
