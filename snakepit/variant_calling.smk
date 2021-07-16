localrules: mumandco, mumandco_summarise, dnadiff, paf_variants, sam2delta_conversion

rule sam2delta_conversion:
    input:
        WORK_PATH + '{haplotype}_scaffolds_{reference}.wm2.sam'
    output:
        temp(WORK_PATH + '{haplotype}_scaffolds_{reference}.wm2.sam.delta')
    shell:
        'python {workflow.basedir}/scripts/sam2delta.py {input}'

rule mumandco:
    input:
        asm = WORK_PATH + '{haplotype}.scaffolds.fasta',
        ref = WORK_PATH + '{haplotype}_scaffolds_{reference}.wm2.sam.delta',
        qur = WORK_PATH + '{reference}_scaffolds_{haplotype}.wm2.sam.delta'
    output:
        dir_ = directory(WORK_PATH + '{haplotype}_{reference}_SV_output'),
        results = multiext(WORK_PATH + '{haplotype}_{reference}_SV_output/{haplotype}_{reference}_SV','.summary.txt','.SVs_all.tsv'),
    params:
        dir_ = lambda wildcards, input: PurePath(input['asm']).parent,
        out = '{haplotype}_{reference}_SV',
        asm = lambda wildcards, input: PurePath(input['asm']).name,
        ref = lambda wildcards: config['ref_genome'] if wildcards.reference == 'ref' else '{reference}.scaffolds.fasta',
        ref_delta = WORK_PATH + '{haplotype}_{reference}_SV_ref.delta',
        qur_delta = WORK_PATH + '{haplotype}_{reference}_SV_query.delta'
    shell:
        '''
        cp {input.ref} {params.ref_delta} && cp {input.qur} {params.qur_delta}
        cd {params.dir_}
        mumandco -r {params.ref} -q {params.asm} -g $( echo "({config[genome_est]}*1000000000)/1" | bc) -o {params.out}
        '''
rule mumandco_summarise:
    input:
        WORK_PATH + '{haplotype}_{reference}_SV_output/{haplotype}_{reference}_SV.SVs_all.tsv'
    output:
        RESULT_PATH + '.{reference}.mumSV.txt'
    shell:
        '''
        > {output}
        for SV_type in insertion deletion inversion transloc; do
            awk -v var=$SV_type '$6 ~ var {{sum += $5; count++}} END {{print "all "var" "count" "sum}}' {input} >> {output}
            awk -v var=$SV_type '$1 ~ /{config[ref_chrm]}/ && $2 ~ /{config[ref_chrm]}/ && $6 ~ var {{sum += $5; count++}} END {{print "chrm "var" "count" "sum}}' {input} >> {output}
        done
        '''

rule paf_variants:
    input:
        WORK_PATH + '{haplotype}_scaffolds_ref.{mapper}.paf'
    output:
        RESULT_PATH + '.{mapper}.vcf'
    shell:
        'sort {input} -k6,6 -k8,8n | paftools.js call -f {config[ref_genome]} - > {output}'

rule dipcall_variants:
    input:
        hifiasm = get_dir('work','hap2.contigs.fasta',assembler='hifiasm'),
        shasta = get_dir('work','asm.contigs.fasta',assembler='shasta')
    output:
        mak = temp(get_dir('comparison','{haplotype}.{sample}.hifiasm.shasta.dipcall.mak')),
        vcf = get_dir('comparison', '{haplotype}.{sample}.hifiasm.shasta.dipcall.dip.vcf.gz'),
        bed = get_dir('comparison', '{haplotype}.{sample}.hifiasm.shasta.dipcall.dip.bed')
    params:
        lambda wildcards, output: PurePath(output['mak']).with_suffix('')
    threads: 16
    resources:
        mem_mb = 5000
    shell:
        '''
        /cluster/work/pausch/alex/software/dipcall/run-dipcall -a -t {threads} {params} {config[ref_genome]} {input.hifiasm} {input.shasta} > {output.mak}
        make -j2 -f {output.mak}
        awk -F'\t' 'BEGIN{{SUM=0}}{{ SUM+=$3-$2 }}END{{print SUM}}' {output.bed}
        '''

rule nucmer:
    input:
        WORK_PATH + '{haplotype}.contigs.fasta'
    output:
        WORK_PATH + '{haplotype}.delta'
    threads: 24
    resources:
        mem_mb = 5500
    shell:
        'nucmer --maxmatch -l 100 -c 500 -t {threads} -p {wildcards.assembler}_{wildcards.sample}/{wildcards.haplotype} {config[ref_genome]} {input}'

rule dnadiff:
    input:
        WORK_PATH + '{haplotype}.delta'
    output:
        WORK_PATH + '{haplotype}.dnadiff.report'
    shell:
        'dnadiff -d {input} -p {wildcards.assembler}_{wildcards.sample}/.{wildcards.haplotype}.dnadiff'
