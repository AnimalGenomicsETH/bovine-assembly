#intersect data
#bcftools isec --threads 4 -p filtered_all ../dedup_alignment/cohort.autosomes.vcf.gz      ../filter_impute/filtered_autosomes_gatk4.vcf.gz
#grep -v '^#' 0001.vcf | cut -f -5 > GATK.raw

#bcftools view -m2 -M2 -v snps 0001.vcf | awk '$1!~/#/ {print $1"\t"$2"\t"$4"\t"$5}' > GATK.ra
#awk '$1!~/#/ {print $1"\t"$2"\t"$4"\t"$5}' /cluster/work/pausch/alex/SNP_exclusive.vcf > TRUTH.ra
#LC_ALL=C fgrep -F -f
#parallel --pipepart -a shared.ra --block 10M LC_ALL=C fgrep -F -f TRUTH.ra > test.ra

rule bcftools_norm:
    input:
        vcf = lambda wildcards: config['vcf'][wildcards.caller],
        ref = config['ref']
    output:
        '{caller}.normed.vcf'
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '30'
    shell:
        'bcftools norm --threads {threads} -m- -f {input.ref} -o {output} {input.vcf}

rule bcftools_isec:
    input:
        DV = lambda wildcards: 'DV.normed.vcf' if wildcards.norm == 'normed' else config['vcf']['DV'],
        GATK = lambda wildcards: 'GATK.normed.vcf' if wildcards.norm == 'normed' else config['vcf']['GATK']
    output:
        directory('intersection_{norm}_{collapse}'),
        expand('intersection_{{norm}}_{{collapse}}/{N}.vcf',N=('DV','GATK','shared'))
    threads: 4
    resources:
        mem_mb = 2000,
        walltime = '30'
    shell:
        '''
        bcftools isec --threads {threads} -c {wildcards.collapse} -p {output} {input.DV} {input.GATK}
        mv {output}/0000.vcf {output}/DV.vcf
        mv {output}/0001.vcf {output}/GATK.vcf
        mv {output}/0002.vcf {output}/shared.vcf
        '''

rule format_record_tsv:
    input:
        'intersection_{norm}_{collapse}/{callset}.vcf'
    output:
        'intersection_{norm}_{collapse}/{callset}.{cols}.tsv'
    params:
        target_cols = lambda wildcards: '$1"\t"$2"\t"$4"\t"$5' if wildcards.cols == 'bases' else '$1"\t"$2'
    shell:
        '''
        awk '$1!~/#/ {{print {params.target_cols}}} {input} > {output}
        '''

rule count_matches:
    input:
        variants = 'intersection_{norm}_{collapse}/{callset}.{cols}.tsv',
        truth = config['truth']
    output:
        'intersection_{norm}_{collapse}/{callset}.{cols}.matches.count'
    threads: lambda wildcards, input: input.size_mb//10
    shell:
        '''
        parallel --pipepart -a {input.variants} --block 10M LC_ALL=C fgrep -F -c -f {input.truth} | awk '{c+=$1}END{print c}' > {output}
        '''

#rule truth_only_positions:
#    'parallel --pipepart -a ../truth.pos_only.tsv --block 1M LC_ALL=C fgrep -F -v -f all.unique > truth_only.tsv'

rule summarise_matches:
    input:
        tsvs = 'intersection_{norm}_{collapse}/{callset}.{cols}.tsv'
        counts = 'intersection_{norm}_{collapse}/{callset}.{cols}.matches.count'
    output:
        'intersection_{norm}_{collapse}_summary.stats'
    shell:
        '''
        echo $(wc -l {input.tsvs[0]}) $(cat {input.counts[0]}) | awk '{print "DV "$3/$1}' >> {output}
        echo $(wc -l {input.tsvs[1]}) $(cat {input.counts[1]}) | awk '{print "GATK "$3/$1}' >> {output}
        echo $(wc -l {input.tsvs[2]}) $(cat {input.counts[2]}) | awk '{print "shared "$3/$1}' >> {output}
