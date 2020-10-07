localrules: sample_data, raw_QC

rule filter_hifi_data:
    input:
        'data/{animal}.raw.hifi.fq.gz'
    output:
        'data/{animal}.cleaned.hifi.fq.gz'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastp -i {input} -o {output} --average_qual {config[filtering][avg_qual]} --length_required {config[filtering][min_length]} --thread {threads} --html data/{wildcards.animal}.html --json /dev/null'

rule sample_data:
    input:
        'data/{animal}.cleaned.hifi.fq.gz'
    output:
        'data/{animal}.{sample}.hifi.fq.gz'
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    shell:
        '''
        if [ {wildcards.sample} -eq 100 ]
        then
            ln -s $(pwd)/{input} {output}
        else
            seqtk sample {input} $(bc <<<"scale=2;{wildcards.sample}/100") | pigz -p 4 > {output}
        fi
        '''

rule raw_QC:
    input:
        'data/{animal}.{read_t}.hifi.fq.gz'
    output:
        'data/{animal}.{read_t}.QC.txt'
    shell:
        '{workflow.basedir}/src/fasterqc {input} {output}'
