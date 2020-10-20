localrules: sample_data, raw_QC

rule raw_read_conversion:
    input:
        f'{config["raw_data"]}/{config["animal"]}/ccs/{{read_name}}.ccs.bam'
    output:
        temp('data/{read_name}.temp.fastq.gz')
    threads: 8
    resources:
        mem_mb = 3000
    shell:
        'samtools fastq -@ {threads} -c 6 -0 {output} {input}'

rule raw_merge_files:
    input:
        expand('data/{read_name}.temp.fastq.gz',read_name=glob_wildcards(f'{config["raw_data"]}/{config["animal"]}/ccs/{{read_name}}.ccs.bam').read_name)
    output:
        protected('data/{animal}.raw.hifi.fq.gz')
    shell: 'cat {input} > {output}'

rule filter_hifi_data:
    input:
        'data/{animal}.raw.hifi.fq.gz'
    output:
        'data/{animal}.cleaned.hifi.fq.gz'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastp -i {input} -o {output} --average_qual {config[filtering][avg_qual]} --length_required {config[filtering][min_length]} --length_limit {config[filtering][max_length]} --thread {threads} --html data/{wildcards.animal}.html --json /dev/null'

rule filter_SR_data:
    input:
        R1 = lambda wildcards: f'{config["short_reads"]}/20200929.B-BOV_{config["trio"][wildcards.parent]}_R1.fastq.gz',
        R2 = lambda wildcards: f'{config["short_reads"]}/20200929.B-BOV_{config["trio"][wildcards.parent]}_R2.fastq.gz'
    output:
        R1 = 'data/{parent}_R1.fastq.gz',
        R2 = 'data/{parent}_R2.fastq.gz'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -g --thread {threads} --html data/{wildcards.parent}.html --json /dev/null'

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
