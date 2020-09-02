rule map_reads:
    input:
        asm = 'canu/{animal}.contigs_raw.fasta',
        reads = 'data/{animal}.hifi.fq.gz'
    output:
        temp('canu/{animal}_read_aln.paf')
    shell:
        'minimap2 -I 6G -xmap-pb -t {threads} {input.asm} {input.reads} > {output}'
 
rule cut_and_split:
    input:
        paf = 'canu/{animal}_read_aln.paf',
        contigs = 'canu/{animal}.contigs_raw.fasta'
    output:
        splits = 'canu/{animal}.split.fasta'
    params:
        config['pd_root']
    shell:
        '''
        {params}/bin/pbcstat {input.paf}
        {params}/bin/calcuts PB.stat > cutoffs
        {params}/bin/split_fa {input.contigs} > {output.splits}
        '''

rule map_splits:
    input:
        'canu/{animal}.split.fasta'
    output:
        'canu/{animal}_ctg-aln.paf'
    shell:
        'minimap2 -I 6G -xasm5 -DP {input} {input} -t {threads} > {output}'

rule purge_dups:
    input:
        paf = 'canu/{animal}_ctg-aln.paf',
        contigs = 'canu/{animal}.contigs.fasta'
    output:
        'test.txt'
    params:
        config['pd_root']
    shell:
        '''
        {params}/bin/purge_dups -2 -T cutoffs -c PB.base.cov {input.paf} > dups.bed
        {params}/bin/get_seqs dups.bed {input.contigs}
        touch {output}
        '''
