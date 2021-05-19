#intersect data
#bcftools isec --threads 4 -p filtered_all ../dedup_alignment/cohort.autosomes.vcf.gz      ../filter_impute/filtered_autosomes_gatk4.vcf.gz
#grep -v '^#' 0001.vcf | cut -f -5 > GATK.raw

#bcftools view -m2 -M2 -v snps 0001.vcf | awk '$1!~/#/ {print $1"\t"$2"\t"$4"\t"$5}' > GATK.ra
#awk '$1!~/#/ {print $1"\t"$2"\t"$4"\t"$5}' /cluster/work/pausch/alex/SNP_exclusive.vcf > TRUTH.ra
#LC_ALL=C fgrep -F -f
#parallel --pipepart -a shared.ra --block 10M LC_ALL=C fgrep -F -f TRUTH.ra > test.ra

wildcard_constraints:
    norm = r'normed|raw',
    collapse = r'none|all',
    cols = r'pos_only|bases',
    var = r'snp|indel'
    
localrules: format_record_tsv, summarise_matches, variant_density

rule all:
    input:
        expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm=('normed','raw'),collapse=('none','all'),cols=('pos_only','bases'),var='snp'),
        expand('intersection_{norm}_{collapse}.{cols}.{var}.summary.stats',norm=('normed','raw'),collapse=('none','all'),cols='pos_only',var='indel')
rule bcftools_norm:
    input:
        vcf = lambda wildcards: config['vcf'][wildcards.caller],
        ref = config['ref']
    output:
        '{caller}.normed.vcf.gz'
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '30'
    shell:
        'bcftools norm --threads {threads} -m- -f {input.ref} -Oz -o {output} {input.vcf}'

rule bcftools_isec:
    input:
        DV = lambda wildcards: 'DV.normed.vcf.gz' if wildcards.norm == 'normed' else config['vcf']['DV'],
        GATK = lambda wildcards: 'GATK.normed.vcf.gz' if wildcards.norm == 'normed' else config['vcf']['GATK']
    output:
        expand('intersection_{{norm}}_{{collapse}}/{N}.vcf',N=('DV','GATK','shared'))
    params:
        lambda wildcards, output: PurePath(output[0]).parent
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '30'
    shell:
        '''
        bcftools isec --threads {threads} -c {wildcards.collapse} -p {params} {input.DV} {input.GATK}
        mv {params}/0000.vcf {params/DV.vcf
        mv {params}/0001.vcf {params}/GATK.vcf
        mv {params}/0002.vcf {params}/shared.vcf
        '''

rule format_record_tsv:
    input:
        'intersection_{norm}_{collapse}/{callset}.vcf'
    output:
        'intersection_{norm}_{collapse}/{callset}.{cols}.{var}.tsv'
    params:
        target_cols = lambda wildcards: '$1"\t"$2"\t"$4"\t"$5' if wildcards.cols == 'bases' else '$1"\t"$2',
        v_type = lambda wildcards: '-m2 -M2 -v snps' if wildcards.var == 'snp' else '-m2 -M2 -v indels'
    shell:
        '''
        bcftools view {params.v_type} {input} | awk '$1!~/#/ {{print {params.target_cols}}}' > {output}
        '''

rule count_matches:
    input:
        variants = 'intersection_{norm}_{collapse}/{callset}.{cols}.{var}.tsv',
        truth = lambda wildcards: config['truth'][wildcards.cols]
    output:
        'intersection_{norm}_{collapse}/{callset}.{cols}.{var}.matches.count'
    threads: 8 #lambda wildcards, input: int(max(input.size_mb,1))
    resources:
        mem_mb = 3000,
        walltime = '15'
    shell:
        '''
        parallel --pipepart -a {input.variants} --jobs {threads} --block 1M LC_ALL=C fgrep -F -c -f {input.truth} | awk '{{c+=$1}}END{{print c}}' > {output}
        '''

#rule truth_only_positions:
#    'parallel --pipepart -a ../truth.pos_only.tsv --block 1M LC_ALL=C fgrep -F -v -f all.unique > truth_only.tsv'

rule summarise_matches:
    input:
        tsvs = expand('intersection_{{norm}}_{{collapse}}/{callset}.{{cols}}.{{var}}.tsv',callset=('DV','GATK','shared')),
        counts = expand('intersection_{{norm}}_{{collapse}}/{callset}.{{cols}}.{{var}}.matches.count',callset=('DV','GATK','shared'))
    output:
        'intersection_{norm}_{collapse}.{cols}.{var}.summary.stats'
    shell:
        '''
        echo $(wc -l {input.tsvs[0]}) $(cat {input.counts[0]}) | awk '{{print "DV "$3/$1}}' >> {output}
        echo $(wc -l {input.tsvs[1]}) $(cat {input.counts[1]}) | awk '{{print "GATK "$3/$1}}' >> {output}
        echo $(wc -l {input.tsvs[2]}) $(cat {input.counts[2]}) | awk '{{print "shared "$3/$1}}' >> {output}
        '''

rule variant_density:
    input:
        'intersection_{norm}_{collapse}/{callset}.pos_only.{var}.tsv'
    output:
        'intersection_{norm}_{collapse}/{callset}.pos_only.{var}.bed'
    run:
        from collections import defaultdict
        with open(output[0],'w') as fout, open(input[0],'r') as fin:
            fout.write('\t'.join(('Chr','Start','End','Value')) + '\n')
            window = 100000

            v_counts = defaultdict(int)
            for line in fin:
                chr, pos = line.rstrip().split()
                v_counts[f'{chr}_{int(pos)//window}'] += 1
            for key,value in v_counts.items():
                chr, pos = key.split('_')
                pos = int(pos) * window
                fout.write('\t'.join(map(str,(chr,pos,pos+window,value))) + '\n')
