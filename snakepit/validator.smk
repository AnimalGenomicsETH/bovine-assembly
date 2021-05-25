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

rule average_quality:
    input:
        'intersection_{{norm}}_{{collapse}}/{N}.vcf'
    output:
        'intersection_{{norm}}_{{collapse}}/{N}.qc'
    shell:
        '''
        awk '!/#/ {{c+=$6;j+=1}} END {{print c/j}}' {input} > {output}
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

rule truth_only_positions:
    input:
        variants = expand('intersection_{{norm}}_{{collapse}}/{callset}.{{cols}}.snp.tsv',callset=('shared','DV','GATK')),
        truth = lambda wildcards: config['truth'][wildcards.cols]
    output:
        'intersection_{norm}_{collapse}/missing_truth.{cols}.snp.tsv'
    threads: 8
    resources:
        mem_mb = 3000,
        walltime = '60'
    shell:
        '''
        for i in {{1..29}}; do
          awk -v i=$i '$1 == i' {input.variants[0]} >> temp.var
          awk -v i=$i '$1 == i' {input.variants[1]} >> temp.var
          awk -v i=$i '$1 == i' {input.variants[2]} >> temp.var
          awk -v i=$i '$1 == i' {input.truth} > temp.truth
          parallel --pipepart -a temp.truth --jobs {threads} --block 100k LC_ALL=C fgrep -F -v -f temp.var >> {output}
        done
        '''

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

from collections import defaultdict

def get_density_counts(f_path,window):
    v_counts = defaultdict(int)
    with open(f_path,'r') as fin:
        for line in fin:
            chr, pos = line.rstrip().split()
            v_counts[(int(chr),int(pos)//window)] += 1
    return v_counts

rule variant_density:
    input:
        lambda wildcards: 'intersection_{norm}_{collapse}/{mode}.pos_only.{var}.tsv' if wildcards.mode != 'exlusive' else ('intersection_{norm}_{collapse}/DV.pos_only.{var}.tsv','intersection_{norm}_{collapse}/GATK.pos_only.{var}.tsv')
    output:
        'intersection_{norm}_{collapse}/{mode}.{var}.ideogram.tsv'
    run:
        window = 100000
        with open(output[0],'w') as fout:
            if wildcards.mode == 'shared':
                fout.write('\t'.join(('Chr','Start','End','Value','Color')) + '\n')
                counts = get_density_counts(input[0],window)
                for (chr,pos),value in counts.items():
                    pos *= window
                    fout.write('\t'.join(map(str,(chr,pos,pos+window,value,'17A589'))) + '\n')

            else:
                fout.write('\t'.join(('Chr','Start','End','Value_1','Color_1','Value_2','Color_2')) + '\n')
                counts_1 = get_density_counts(input[0],window)
                counts_2 = get_density_counts(input[1],window)
                for key in sorted(counts_1.keys() | counts_2.keys()):
                    chr,pos = key 
                    pos *= window
                    fout.write('\t'.join(map(str,(chr,pos,pos+window,counts_1[key],'2874A6',counts_2[key],'D4AC0D'))) + '\n')


rule picard_concordance:
    input:
        DV = 'DV.normed.vcf.gz',
        GATK = 'GATK.normed.vcf.gz'
    output:
        'HD_concordance/{sample}.{caller}.conc'
    shell:
        '''
        #picard GenotypeConcordance CALL_VCF={input.GATK} TRUTH_VCF={input.DV} CALL_SAMPLE={wildcards.sample} TRUTH_SAMPLE={wildcards.sample} O={wildcards.sample}
        tail -n 4 HD_concordance/{wildcards.sample}_{wildcards.caller}.genotype_concordance_summary_metrics | awk '$1~/SNP/ {{print $2"\t{wildcards.caller}\t"$10"\t"$12"\t"2*$4*$5/($4+$5)"\t"2*$7*$8/($7+$8)"\t"2*$10*$11/($10+$11)"\t"$13"\t"$14}}' > {output}
        '''

rule aggreate_concordance:
    input:
        expand('HD_concordance/{sample}.{caller}.conc',sample=('RM1894','RM1896','RM1899','BSWCHEF120071057962'),caller=('DV','GATK'))
    output:
        'HD_concordance/summary.txt'
    shell:
        '''
        echo -e "sample\tcaller\tSensitivity\tSpecificity\tHET F1\tHOMV F1\tVAR F1\tGC\tnR-GC" > {output}
        sort -k 13,13nr -k 11,11nr {input} >> {output}
        '''
