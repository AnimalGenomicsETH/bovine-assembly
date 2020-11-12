localrules: raw_merge_files, sample_data, raw_QC, fastq_to_fasta

rule raw_read_conversion:
    input:
        lambda wildcards: f'{config["data"][config["animal"]]["long_reads"][wildcards.individual]}{{read_name}}.ccs.bam'
    output:
        temp('data/{individual}_{read_name}.temp.fastq.gz')
    threads: 8
    resources:
        mem_mb = 3000
    shell:
        'samtools fastq -@ {threads} -c 6 -0 {output} {input}'

rule raw_merge_files:
    input:
        lambda wildcards: expand('data/{{individual}}_{read_name}.temp.fastq.gz',read_name=glob_wildcards(f'{config["data"][config["animal"]]["long_reads"][wildcards.individual]}{{read_name}}.ccs.bam').read_name)
    output:
        protected('data/{individual}.raw.hifi.fq.gz')
    shell: 'cat {input} > {output}'

rule filter_hifi_data:
    input:
        'data/{individual}.raw.hifi.fq.gz'
    output:
        'data/{individual}.cleaned.hifi.fq.gz'
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastp -i {input} -o {output} --average_qual {config[filtering][avg_qual]} --length_required {config[filtering][min_length]} --length_limit {config[filtering][max_length]} --thread {threads} --html data/{wildcards.individual}.html --json /dev/null'

rule fastq_to_fasta:
    input:
        'data/{parent}.cleaned.hifi.fq.gz'
    output:
        'data/{parent}.hifi.fasta'
    shell:
        'seqtk seq -a {input} > {output}'

rule filter_SR_data:
    input:
        reads = lambda wildcards: expand(f'{config["data"][config["animal"]]["short_reads"][wildcards.parent]}_R{{N}}.fastq.gz', N = (1,2))
    output:
        reads = expand('data/{{parent}}.read_R{N}.SR.fq.gz', N = (1,2))
    threads: 12
    resources:
        mem_mb = 4000
    shell:
        'fastp -i {input.reads[0]} -I {input.reads[1]} -o {output.reads[0]} -O {output.reads[1]} -g --thread {threads} --html data/{wildcards.parent}.html --json /dev/null'

rule sample_data:
    input:
        'data/offspring.cleaned.hifi.fq.gz'
    output:
        'data/offspring.{sample}.hifi.fq.gz'
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
        'data/offspring.{read_t}.hifi.fq.gz'
    output:
        'data/offspring.{read_t}.QC.txt'
    shell:
        '{workflow.basedir}/src/fasterqc {input} {output}'
