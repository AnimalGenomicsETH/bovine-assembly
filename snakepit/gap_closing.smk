
rule TGS_gapcloser_chromsome:
    input:
        scaffolds = get_dir('work','{haplotype}.scaffolds.fasta'),
        reads = config['TGS_reads']
    output:
        get_dir('work','{haplotype}.{chr}.filled.fasta'),
        temp = get_dir('work','{haplotype}_{chr}.fa')
    params:
        dir_ = lambda wildcards: get_dir('work',f'TGS_{wildcards.haplotype}_{wildcards.chr}'),
        scaffolds = lambda wildcards, input: '../' + PurePath(input['scaffolds']).name,
        reads = lambda wildcards, input: '../../' + input['reads']
    threads: 12
    resources:
        mem_mb = 10000,
        walltime = '4:00'
    shell:
        '''
        mkdir -p {params.dir_}
        seqtk subseq {input.scaffolds} <(echo {wildcards.chr})> {output.temp}
        (cd {params.dir_} && {config[tgs_root]}/TGS-GapCloser.sh --scaff {output.temp} --reads {params.reads} --output {wildcards.haplotype} --minmap_arg '-x ava-ont' --tgstype ont --racon --racon /cluster/work/pausch/alex/software/racon/build/bin/racon --thread {threads})
        cp {params.dir_}/{params.out}.scaff_seqs {output}
        '''

rule TGS_gapcloser_merge:
    input:
        expand(get_dir('work','{haplotype}.{chr}.filled.fasta'),chr=range(1,31)),
        get_dir('work','{haplotype}.X.filled.fasta')
    output:
        get_dir('work','{haplotype}.filled.fasta')
    shell:
        'cat {input} > {output}'
