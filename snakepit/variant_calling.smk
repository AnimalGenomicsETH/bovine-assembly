localrules: mumandco, dnadiff

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
        cp {output.results[0]} {output.summary}
        '''

#rule summarise_variants:
#    input:
#        '{assembler}_{sample}/{haplotype}_{reference}_SV'

rule paf_variants:
    input:
        '{assembler}_{sample}/{haplotype}_scaffolds_ref.{mapper}.paf'
    output:
        'results/{haplotype}_{sample}_{assembler}.{mapper}.vcf'
    shell:
        'sort {input} -k6,6 -k8,8n | paftools.js call -f {config[ref_genome]} - > {output}'

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
