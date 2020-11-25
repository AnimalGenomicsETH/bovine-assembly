localrules: mumandco, dnadiff, paf_variants, sam2delta_conversion

rule sam2delta_conversion:
    input:
        '{assembler}_{sample}/{haplotype}_scaffolds_{reference}.wm2.sam'
    output:
        temp('{assembler}_{sample}/{haplotype}_scaffolds_{reference}.wm2.sam.delta')
    shell:
        'python {workflow.basedir}/scripts/sam2delta.py {input}'

rule mumandco:
    input:
        asm = '{assembler}_{sample}/{haplotype}.scaffolds.fasta',
        ref = '{assembler}_{sample}/{haplotype}_scaffolds_{reference}.wm2.sam.delta',
        qur = '{assembler}_{sample}/{reference}_scaffolds_{haplotype}.wm2.sam.delta'
    output:
        dir_ = directory('{assembler}_{sample}/{haplotype}_{reference}_SV_output'),
        results = multiext('{assembler}_{sample}/{haplotype}_{reference}_SV_output/{haplotype}_{reference}_SV','.summary.txt','.SVs_all.tsv'),
        summary = 'results/{haplotype}_{sample}_{assembler}.{reference}.mumSV.txt'
    params:
        dir_ = lambda wildcards, input: PurePath(input['asm']).parent,
        out = '{haplotype}_{reference}_SV',
        asm = lambda wildcards, input: PurePath(input['asm']).name,
        ref = lambda wildcards: config['ref_genome'] if wildcards.reference == 'ref' else '{reference}.scaffolds.fasta',
        ref_delta = '{assembler}_{sample}/{haplotype}_{reference}_SV_ref.delta',
        qur_delta = '{assembler}_{sample}/{haplotype}_{reference}_SV_query.delta'
    shell:
        '''
        cp {input.ref} {params.ref_delta} && cp {input.qur} {params.qur_delta}
        (cd {params.dir_} && mumandco -r {params.ref} -q {params.asm} -g $( echo "({config[genome_est]}*1000000000)/1" | bc) -o {params.out})
        for SV_type in insertion deletion inversion transloc; do
            awk -v var=$SV_type '$6 ~ var {{sum += $5; count++}} END {{print "all "var" "count" "sum}}' {output.results[1]} > {output.summary}
            awk -v var=$SV_type '$1 ~ /{config[ref_chrm]}/ && $2 ~ /{config[ref_chrm]}/ && $6 ~ var {{sum += $5; count++}} END {{print "chrm "var" "count" "sum}}' {output.results[1]} > {output.summary}
        done
        '''

rule paf_variants:
    input:
        '{assembler}_{sample}/{haplotype}_scaffolds_ref.{mapper}.paf'
    output:
        'results/{haplotype}_{sample}_{assembler}.{mapper}.vcf'
    shell:
        'sort {input} -k6,6 -k8,8n | paftools.js call -f {config[ref_genome]} - > {output}'

rule dipcall_variants:
    input:
        hapA = '{assembler}_{sample}/{hapA}.contigs.fasta',
        hapB = '{assembler}_{sample}/{hapB}.contigs.fasta'
    output:
        mak = temp('{assembler}_{sample}/{hapA}_{hapB}.mak'),
        vcf = '{assembler}_{sample}/{hapA}_{hapB}.dip.vcf.gz'
    params:
        '{hapA}_{hapB}'
    threads: 16
    resources:
        mem_mb = 3000
    shell:
        '''
        run-dipcall {params} {config[ref_genome]} {input.hapA} {input.hapB} > {output.mak}
        make -j2 -f {output.mak}
        '''

rule nucmer:
    input:
        '{assembler}_{sample}/{haplotype}.contigs.fasta'
    output:
        '{assembler}_{sample}/{haplotype}.delta'
    threads: 24
    resources:
        mem_mb = 5500
    shell:
        'nucmer --maxmatch -l 100 -c 500 -t {threads} -p {wildcards.assembler}_{wildcards.sample}/{wildcards.haplotype} {config[ref_genome]} {input}'

rule dnadiff:
    input:
        '{assembler}_{sample}/{haplotype}.delta'
    output:
        '{assembler}_{sample}/{haplotype}.dnadiff.report'
    shell:
        'dnadiff -d {input} -p {wildcards.assembler}_{wildcards.sample}/.{wildcards.haplotype}.dnadiff'
