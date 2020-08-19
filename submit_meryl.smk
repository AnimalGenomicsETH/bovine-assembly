#mkdir split
rule split:
    input:
        '/cluster/work/pausch/alex/assembly/BSWCHEF120152514636/BSWCHEF120152514636.fq.gz'
    output:
        'split/BSWCHEF120152514636.0000.fq'
    params:
        'split/BSWCHEF120152514636.'
    shell: 'zcat {input} | split -a 4 -d -C 10GB --additional-suffix=.fq - {params}'

rule rezip:
    input:
        'split/{animal}.{sample}.fq'
    output:
        'split/{animal}.'
    shell:
        '''
        pigz -p {threads} --recursive  split/{animal}.*.fq
        ls split/{animal}.[0-9][0-9][0-9][0-9].fq.gz > $FOFN.$tid
        '''
