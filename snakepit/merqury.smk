localrules: merqury_formatting, merqury_formatting_simple

### LOCAL EXECUTION TEMPORARY CODE ###

if 'get_dir' not in dir():
    from pathlib import PurePath

    wildcard_constraints:
        assembler = r'[^\W_]+',
        parent = r'dam|sire',
        individual = r'dam|sire|offspring',
        haplotype = r'asm|hap1|hap2|sire|dam|ref',
        hap = r'\w+',
        sample = r'\d+',
        data = r'[^\W_]+',
        modifier = r'\w+',
        read_t  = r'hifi|SR',
        mapper = r'mm2|wm2'

    class Default(dict):
        def __missing__(self, key):
            return '{'+key+'}'

    def get_dir(base,ext='',**kwargs):
        if base == 'work':
            base_dir = 'merqury_{haplotype}'
        elif base == 'result':
            base_dir = 'merqury_{haplotype}'
        elif base == 'data':
            base_dir = 'data'
        else:
            raise Exception('Base not found')
        if ext and ext[0] == '.':
            return f'{base_dir}{ext}'.format_map(Default(kwargs))
        return str(PurePath(base_dir.format_map(Default(kwargs))) / ext)

    config['k-mers']=21
    config['mem_adj']=1100
    config['merqury_root']='/cluster/work/pausch/alex/software/merqury'
    config['genome_est']=2.7

    rule all:
        input:
            get_dir('result','.merqury.full.stats',haplotype='asm')

for configuration in ('assembly','short_reads'):
    if configuration not in config:
        config[configuration] = None

rule count_asm_kmers:
    input:
        config['assembly'] or get_dir('work','{haplotype}.contigs.fasta')
    output:
        directory(get_dir('work','{haplotype}.contigs.meryl'))
    threads: 8
    resources:
        mem_mb = 2500
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule count_SR_reads:
    input:
        'data/{individual}.read_R{N}.SR.fq.gz' #config['short_reads'] or
    output:
        temp(directory('data/{individual}.read_R{N}.meryl'))
    threads: 24
    resources:
        mem_mb = 4000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl count k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule merge_SR_reads:
    input:
        expand('data/{{individual}}.read_R{N}.meryl', N = (1,2))
    output:
        directory('data/{individual}.meryl')
    threads: 12
    resources:
        mem_mb = 4000
    params:
        mem = lambda wildcards,resources,threads: resources['mem_mb']*threads/config['mem_adj']
    shell:
        'meryl union-sum k={config[k-mers]} memory={params.mem} threads={threads} output {output} {input}'

rule generate_hapmers:
    input:
        dam = 'data/dam.meryl',
        sire = 'data/sire.meryl',
        child = 'data/offspring.meryl'
    output:
        directory(expand('data/{parent}.hapmer.meryl', parent = ('dam', 'sire')))
    params:
        dam = lambda wildcards,input: PurePath(input['dam']).name,
        sire = lambda wildcards,input: PurePath(input['sire']).name,
        child = lambda wildcards,input: PurePath(input['child']).name
    threads: 12
    resources:
        mem_mb = 4000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd data
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/hapmers.sh {params.sire} {params.dam} {params.child}
        '''

rule meryl_lookup_exclude:
    input:
        reads = expand('data/offspring.read_R{N}.SR.fq.gz', N=(1,2)),
        #get opposite parental hapmers
        hapmer = lambda wildcards: f'data/{"dam" if wildcards.parent == "sire" else "sire"}.hapmer.meryl'
    output:
        reads = expand('data/offspring.{{parent}}.read_R{N}.SR.fq.gz', N=(1,2))
    shell:
        '''
        meryl-lookup -exclude -sequence {input.reads} -mers {input.hapmer} -threads {threads} -output {output.reads}
        '''

rule hap_blob:
    input:
        hapmers = expand('data/{parent}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: (get_dir('work',f'{X}.contigs.fasta') for X in (('asm',) if wildcards.hap == 'asm' else ('hap1','hap2')))
    output:
        get_dir('work','{hap}.hapmers.blob.png')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{hap}',
        hapmers = lambda wildcards, input: tuple('..' / PurePath(fname) for fname in input['hapmers']),
        asm = lambda wildcards, input: tuple(PurePath(fname).name for fname in input['asm'])
    threads: 8
    resources:
        mem_mb = 4000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/hap_blob.sh {params.hapmers} {params.asm} {params.out}
        '''

rule merqury_spectra_cn:
    input:
        read_db = lambda wildcards: 'data/offspring.meryl' if wildcards.hap in ('asm','trio') else 'data/{hap}.meryl',
        filt =  lambda wildcards: 'data/offspring.filt' if wildcards.hap in ('asm','trio') else 'data/{hap}.filt',
        asm = lambda wildcards: (get_dir('work',f'{X}.contigs.fasta') for X in (('asm',) if wildcards.hap == 'asm' else ('hap1','hap2'))),
        asm_dbs = lambda wildcards: (get_dir('work',f'{X}.contigs.meryl') for X in (('asm',) if wildcards.hap == 'asm' else ('hap1','hap2')))
    output:
        multiext(get_dir('work','{hap}'), '.qv', '.completeness.stats')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        base = lambda wildcards, input: PurePath(input['read_db']).stem,
        out = '{hap}',
        read_db = lambda wildcards, input: PurePath(input['read_db']).name,
        asm = lambda wildcards, input: tuple(PurePath(fname).name for fname in input['asm'])
    threads: 8
    resources:
        mem_mb = 6000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        ln -sfn  ../{input.filt}
        ln -sfn ../{input.read_db}
        find ../data/ -name "{params.base}.gt*.meryl" -exec ln -sfn ../data/{{}} . \;
        $MERQURY/eval/spectra-cn.sh {params.read_db} {params.asm} {params.out}
        '''

rule merqury_spectra_hap:
    input:
        reads = 'data/offspring.meryl',
        hapmers = expand('data/{parent}.hapmer.meryl',parent=('sire','dam')),
        asm = lambda wildcards: (get_dir('work',f'{X}.contigs.fasta') for X in (('asm',) if wildcards.hap == 'asm' else ('hap1','hap2'))),
        stats = get_dir('work','{hap}.completeness.stats')
    output:
        get_dir('work','{hap}.hap.completeness.stats')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{hap}.hap',
        hapmers = lambda wildcards, input: tuple('..' / PurePath(fname) for fname in input['hapmers']),
        asm = lambda wildcards, input: tuple(PurePath(fname).name for fname in input['asm'])
    threads: 12
    resources:
        mem_mb = 4000
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/spectra-hap.sh ../{input.reads} {params.hapmers} {params.asm} {params.out}
        '''

rule merqury_phase_block:
    input:
        asm = get_dir('work','{haplotype}.contigs.fasta'),
        hapmers = expand('data/{parent}.hapmer.meryl',parent=('sire','dam'))
    output:
        multiext(get_dir('work','{haplotype}.100_20000.'),'phased_block.bed','switch.bed', 'switches.txt','phased_block.stats')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{haplotype}',
        hapmers = lambda wildcards, input: tuple('..' / PurePath(fname) for fname in input['hapmers']),
        asm = lambda wildcards, input: PurePath(input['asm']).name
    threads: 12
    resources:
        mem_mb = 3500
    envmodules:
       'gcc/8.2.0',
       'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/phase_block.sh {params.asm} {params.hapmers} {params.out}
        '''

rule merqury_block_n_stats:
    input:
        asm_block = multiext(get_dir('work','{haplotype}'),'.contigs.fasta','.100_20000.phased_block.bed')
    output:
        multiext(get_dir('work','{haplotype}.{haplotype}.contigs'),'.continuity.NG.png','.block.NG.png')
    params:
        dir_ = lambda wildcards, output: PurePath(output[0]).parent,
        out = '{haplotype}',
        asm_block = lambda wildcards,input: tuple(PurePath(i).name for i in input['asm_block'])
    threads: 6
    resources:
        mem_mb = 2500
    envmodules:
        'gcc/8.2.0',
        'r/4.0.2'
    shell:
        '''
        cd {params.dir_}
        export MERQURY={config[merqury_root]}
        $MERQURY/trio/block_n_stats.sh {params.asm_block} {params.out} $( echo "({config[genome_est]}*1000000000)/1" | bc)
        '''

rule merqury_prep:
    input:
        'data/{individual}.meryl'
    output:
        hist = 'data/{individual}.hist',
        filt = 'data/{individual}.filt',
    threads: 12
    resources:
        mem_mb = 2000
    shell:
        '''
        meryl histogram {input} > {output.hist}
        java -jar -Xmx1g {config[merqury_root]}/eval/kmerHistToPloidyDepth.jar {output.hist} > {output.hist}.ploidy
        filt=`sed -n 2p {output.hist}.ploidy | awk '{{print $NF}}'`
        echo $filt > {output.filt}
        meryl greater-than $filt output data/{wildcards.individual}.gt$filt.meryl {input}
        '''

rule merqury_formatting:
    input:
        get_dir('work','{haplotype}.{haplotype}.contigs.continuity.NG.png'),
        stats = lambda wildcards: multiext(get_dir('work','asm' if wildcards.haplotype == 'asm' else 'trio'),'.hapmers.blob.png','.completeness.stats','.hap.completeness.stats','.qv'),
        switches = get_dir('work','{haplotype}.100_20000.switches.txt'),
        phase = get_dir('work','{haplotype}.100_20000.phased_block.stats'),
    output:
        get_dir('result','.merqury.full.stats')
    shell:
        '''
        awk '/{wildcards.haplotype}/ {{print "QV "$4}}' {input.stats[3]} > {output}
        awk '/{wildcards.haplotype}/ {{print "completeness "$5}}' {input.stats[1]} >> {output}
        awk '/{wildcards.haplotype}/ && /sire/ {{print "sire "$5}}' {input.stats[2]} >> {output}
        awk '/{wildcards.haplotype}/ && /dam/ {{print "dam "$5}}' {input.stats[2]} >> {output}
        awk '{{print "switches "$14+0}}' {input.switches} >> {output}
        awk '{{print "phased "$6}}' {input.phase} >> {output}
        '''

rule merqury_formatting_simple:
    input:
        stats = get_dir('work','{haplotype}.completeness.stats'),
        qv = get_dir('work','{haplotype}.qv'),
    output:
        get_dir('result','.merqury.simple.stats')
    shell:
        '''
        awk '/{wildcards.haplotype}/ {{print "QV "$4}}' {input.qv} > {output}
        awk '/{wildcards.haplotype}/ {{print "completeness "$5}}' {input.stats} >> {output}
        '''
