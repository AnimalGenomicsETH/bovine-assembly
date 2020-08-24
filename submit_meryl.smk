
from pathlib import Path

rule diverge:
    input:
        hifi reads
    output:
        meryl
    run:
        if Path('{input}').stat().st_size() > 10*(1024**3):
            do_split
        else:
            count

rule split:
    input:
        '/cluster/work/pausch/alex/assembly/BSWCHEF120152514636/BSWCHEF120152514636.fq.gz'
    output:
        'split/BSWCHEF120152514636.0000.fq'
    params:
        'split/BSWCHEF120152514636.'
    shell: 'zcat {input} | split -a 4 -d -C 10GiB --additional-suffix=.fq - {params}'

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

rule merge:
    input:
        list to merge
    output:
        readdb
    shell:
        meryl sum

rule prep_for_count:
    input:
        reads
    output:
        something
    run:
        if file larger than 10GB -> split
        else dont

