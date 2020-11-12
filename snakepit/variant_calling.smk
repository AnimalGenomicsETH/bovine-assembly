localrules: mumandco

rule mumandco:
    input:
        '{assembler}_{sample}/{haplotype}.scaffolds.fasta'
    output:
        dir_ = directory('{assembler}_{sample}/{haplotype}_{reference}_SV_output'),
        results = multiext('{assembler}_{sample}/{haplotype}_{reference}_SV_output/{haplotype}_{reference}_SV','.summary.txt','.SVs_all.tsv')
    params:
        dir_ = '{assembler}_{sample}/{haplotype}_{reference}_SV',
        ref = lambda wildcards: config['ref_genome'] if wildcards.reference == 'ref' else '{assembler}_{sample}/{reference}.scaffolds.fasta'
    shell:
        'mumandco -r {params.ref} -q {input} -g $( echo "({config[genome_est]}*1000000000)/1" | bc) -o {params.dir_}'

#rule summarise_variants:
#    input:
#        '{assembler}_{sample}/{haplotype}_{reference}_SV'
