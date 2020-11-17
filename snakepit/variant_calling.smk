localrules: mumandco, dnadiff, paf_variants

rule mumandco:
    input:
        '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        dir_ = directory('{assembler}_{sample}/{haplotype}_{reference}_SV_output'),
        results = multiext('{assembler}_{sample}/{haplotype}_{reference}_SV_output/{haplotype}_{reference}_SV','.summary.txt','.SVs_all.tsv'),
        summary = 'results/{haplotype}_{sample}_{assembler}.{reference}.mumSV.txt'
    params:
        dir_ = '{assembler}_{sample}/{haplotype}_{reference}_SV',
        ref = lambda wildcards: config['ref_genome'] if wildcards.reference == 'ref' else '{assembler}_{sample}/{reference}.scaffolds.fasta',
    shell:
        '''
        mumandco -r {params.ref} -q {input} -g $( echo "({config[genome_est]}*1000000000)/1" | bc) -o {params.dir_}
        for SV_type in insertion deletion inversion transloc; do
            awk -v var=$SV_type '$1 ~ /{config[ref_chrm]}/ && $2 ~ /{config[reg_chrm]}/ && $6 ~ var {sum += $5} END {print var" "sum}' {output.results[1]} > {outputs.summary}
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
        mak = temp('{hapA}_{hapB}.mak'),
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
