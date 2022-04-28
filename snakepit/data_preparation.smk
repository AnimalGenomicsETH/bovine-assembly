localrules: raw_merge_files, raw_QC, fastq_to_fasta

rule raw_read_conversion:
    input:
        lambda wildcards: f'{config["data"]["long_reads"][wildcards.individual]}{{read_name}}.ccs.bam'
    output:
        temp('data/{individual}_{read_name}.temp.fastq.gz')
    threads: 8
    resources:
        mem_mb = 3000
    shell:
        'samtools fastq -@ {threads} -c 6 -0 {output} {input}'

rule raw_merge_files:
    input:
        lambda wildcards: expand('data/{{individual}}_{read_name}.temp.fastq.gz',read_name=glob_wildcards(f'{config["data"]["long_reads"][wildcards.individual]}{{read_name}}.ccs.bam').read_name)
    output:
        protected('data/{individual}.raw.hifi.fq.gz')
    shell: 'cat {input} > {output}'

rule filter_hifi_data:
    input:
        'data/{individual}.raw.hifi.fq.gz'
        #'data/{individual}.raw.{read_tech}.fq.gz'
    output:
        'data/{individual}.cleaned.hifi.fq.gz'
        #'data/{individual}.cleaned.{read_tech}.fq.gz'
    params:
        html = 'data/{individual}.hifi.html'
        ##'data/{individual}.{read_tech}.html'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastp -i {input} -o {output} --average_qual {config[filtering][avg_qual]} --length_required {config[filtering][min_length]} --length_limit {config[filtering][max_length]} --thread {threads} --html {params.html} --json /dev/null'
        #can filter sepearately for hifi and pion reads

rule fastq_to_fasta:
    input:
        'data/{parent}.cleaned.hifi.fq.gz'
    output:
        temp('data/{parent}.hifi.fasta')
    shell:
        'seqtk seq -a {input} > {output}'

rule filter_SR_data:
    input:
        reads = lambda wildcards: expand(f'{config["data"]["short_reads"][wildcards.individual]}_R{{N}}.fastq.gz', N = (1,2))
    output:
        reads = expand('data/{{individual}}.read_R{N}.SR.fq.gz', N = (1,2))
    params:
        html = 'data/{individual}.SR.html'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastp -i {input.reads[0]} -I {input.reads[1]} -o {output.reads[0]} -O {output.reads[1]} -g --thread {threads} --html {params.html} --json /dev/null'

import random
rule sample_data:
    input:
        'data/offspring.cleaned.hifi.fq.gz'
    output:
        temp('data/offspring.{sample}.hifi.fa.gz')
    envmodules:
        'gcc/8.2.0',
        'pigz/2.4'
    threads: 4
    resources:
        mem_mb = 3000
    shell:
        '''
        if [ {wildcards.sample} -eq 100 ]
        then
            ln -s $(pwd)/{input} {output}
        else
            echo $(python -c 'import random; print(random.getrandbits(64))')
            seqtk seq -A -s $(python -c 'import random; print(random.getrandbits(64))') -f $(bc <<<"scale=2;{wildcards.sample}/100") {input} | pigz -p {threads} > {output}
        fi
        '''

rule raw_QC:
    input:
        'data/{individual}.{sample}.hifi.fq.gz'
    output:
        'data/{individual}.{sample}.QC.txt'
    shell:
        '{workflow.basedir}/src/fasterqc {input} {output}'
