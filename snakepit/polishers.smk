
MODEL = '/cluster/work/pausch/alex/software/pepper/models/PEPPER_polish_haploid_guppy360.pkl'

rule all:
    input:
        'Assembly.contigs.fasta'

rule map_ONT_reads:
    input:
        reads = config['reads'],
        asm = '{haplotype}.fasta'
    output:
        temp('{haplotype}_ONT_reads.unsorted.bam')
    threads: 16
    resources:
        mem_mb = 4500,
        walltime = '24:00'
    shell:
        'minimap2 -ax map-ont -t {threads} {input.asm} {input.reads} | samtools view -hb -F 0x904 -o {output}'

rule sort_bam:
    input:
        '{haplotype}_ONT_reads.unsorted.bam'
    output:
        temp('{haplotype}_ONT_reads.sorted.bam')
    threads: 8
    resources:
        mem_mb = 6000,
        disk_scratch = 200
    shell:
        'samtools sort {input} -m 3000M -@ {threads} -T $TMPDIR -o {output}'

rule index_bam:
    input:
        '{haplotype}_ONT_reads.sorted.bam'
    output:
        temp('{haplotype}_ONT_reads.sorted.bam.bai')
    threads: 4
    resources:
        mem_mb = 4000
    shell:
        'samtools index -@ {threads} {input}'

rule pepper_make_images:
    input:
        asm = '{haplotype}.fasta',
        bam = '{haplotype}_ONT_reads.sorted.bam',
        bai = '{haplotype}_ONT_reads.sorted.bam.bai'
    output:
        temp(directory('pepper_images_{haplotype}'))
    threads: 12
    resources:
        mem_mb = 6000,
        walltime = '24:00'
    shell:
        'pepper make_images -b {input.bam} -f {input.asm} -o {output} -t {threads}'

rule pepper_call_consensus:
    input:
        images = 'pepper_images_{haplotype}'
    output:
        temp(directory('pepper_predictions_{haplotype}'))
    params:
        batch_size = 256,
        workers = 4,
        model = MODEL
    threads: 12
    resources:
        mem_mb = 7500,
        walltime = '24:00'
    shell:
        'pepper call_consensus -i {input.images} -bs {params.batch_size} -w {params.workers} -m {params.model} -t {threads} -o {output}'

rule pepper_stitch:
    input:
        'pepper_predictions_{haplotype}'
    output:
        '{haplotype}.contigs.fasta'
    threads: 6
    resources:
        mem_mb = 15000,
        walltime = '24:00'
    shell:
        'pepper stitch -i {input} -o {output} -t {threads}'



mkdir -p PEPPERv4;
mkdir -p PEPPERv4/logs;
mkdir -p PEPPERv4/intermediate_files;

rule pepper_snp:
    input:
        ''
    output:
        ''
    shell:
        '''
        pepper_snp call_variant -b bam -f fasta -t 24 -m /opt/pepper_models/PEPPER_SNP_R941_ONT_V4.pkl -o PEPPERv4/pepper_snp/ -s gaur -w 4 -bs 64 --ont 2>&1 | tee PEPPERv4/logs/1_pepper_snp.log
        '''

rule pepper_snp_post:
    input:
        ''
    output:
        ''
    shell:
        '''
        mv PEPPERv4/pepper_snp/*.vcf PEPPERv4/PEPPER_SNP_OUPUT.vcf;
        bgzip PEPPERv4/PEPPER_SNP_OUPUT.vcf;
        tabix -p vcf PEPPERv4/PEPPER_SNP_OUPUT.vcf.gz;
        rm -rf PEPPERv4/pepper_snp/;
        '''

rule pepper_margin:
    input:
        ''
    output:
        ''
    shell:
        '''
        margin phase bam fasta PEPPERv4/PEPPER_SNP_OUPUT.vcf.gz /opt/margin_dir/params/misc/allParams.ont_haplotag.json -t 24 -V -o PEPPERv4/MARGIN_PHASED.PEPPER_SNP_MARGIN 2>&1 | tee PEPPERv4/logs/2_margin_haplotag.log;
        '''

rule pepper_margin_post:
    input:
        ''
    output:
        ''
    shell:
        '''
        mv PEPPERv4/*.bam PEPPERv4/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam;
        samtools index -@24 PEPPERv4/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam
        '''

rule pepper_hp:
    input:
        ''
    output:
        ''
    shell:
        '''
        pepper_hp call_variant -b PEPPERv4/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam -f fasta -t 24 -m /opt/pepper_models/PEPPER_HP_R941_ONT_V4.pkl -o PEPPERv4/pepper_hp/ -s gaur -w 4 -bs 64 --ont_asm 2>&1 | tee PEPPERv4/logs/3_pepper_hp.log
        '''

rule pepper_hp_post:
    input:
        ''
    output:
        ''
    shell:
        '''

        mv PEPPERv4/pepper_hp/PEPPER_HP_OUTPUT_1.vcf PEPPERv4/PEPPER_HP_OUTPUT_1.vcf;
        mv PEPPERv4/pepper_hp/PEPPER_HP_OUTPUT_2.vcf PEPPERv4/PEPPER_HP_OUTPUT_2.vcf;
        bgzip PEPPERv4/PEPPER_HP_OUTPUT_1.vcf;
        bgzip PEPPERv4/PEPPER_HP_OUTPUT_2.vcf;
        tabix -p vcf PEPPERv4/PEPPER_HP_OUTPUT_1.vcf.gz;
        tabix -p vcf PEPPERv4/PEPPER_HP_OUTPUT_2.vcf.gz;
        rm -rf PEPPERv4/pepper_hp/;
        '''

rule pepper_DV:

mkdir -p PEPPERv4/dv_intermediate_outputs/;
echo "STARTING DEEPVARIANT";
time /opt/deepvariant/bin/run_deepvariant --model_type=WGS --customized_model=/opt/dv_models/202012_polish_nohp_rows/model.ckpt-28400 --ref=fasta --reads=PEPPERv4/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam --output_vcf=PEPPERv4/PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP1.vcf.gz --sample_name="gaur" --intermediate_results_dir=PEPPERv4/dv_intermediate_outputs/ --num_shards=24 --make_examples_extra_args="alt_aligned_pileup=rows,realign_reads=false,sort_by_haplotypes=true,parse_sam_aux_fields=true,add_hp_channel=false,hp_tag_for_assembly_polishing=1,variant_caller=vcf_candidate_importer,proposed_variants=PEPPERv4/PEPPER_HP_OUTPUT_1.vcf.gz"  2>&1 | tee PEPPERv4/logs/4_DeepVariant.log
mkdir -p PEPPERv4/dv_intermediate_outputs/;
echo "STARTING DEEPVARIANT";
    time /opt/deepvariant/bin/run_deepvariant --model_type=WGS --customized_model=/opt/dv_models/202012_polish_nohp_rows/model.ckpt-28400 --ref=fasta --reads=PEPPERv4/MARGIN_PHASED.PEPPER_SNP_MARGIN.haplotagged.bam --output_vcf=PEPPERv4/PEPPER_MARGIN_DEEPVARIANT_ASM_POLISHED_HAP2.vcf.gz --sample_name="gaur" --intermediate_results_dir=PEPPERv4/dv_intermediate_outputs/ --num_shards=24 --make_examples_extra_args="alt_aligned_pileup=rows,realign_reads=false,sort_by_haplotypes=true,parse_sam_aux_fields=true,add_hp_channel=false,hp_tag_for_assembly_polishing=2,variant_caller=vcf_candidate_importer,proposed_variants=PEPPERv4/PEPPER_HP_OUTPUT_2.vcf.gz"  2>&1 | tee PEPPERv4/logs/4_DeepVariant.log



# rule racon_polish:
#     input:
#         scaffolds = WORK_PATH + '{haplotype}.scaffolds.fasta',
#         aln = WORK_PATH + '{haplotype}_scaffolds_reads.sam',
#         reads = 'data/offspring.{sample}.hifi.fq.gz'
#     output:
#         WORK_PATH + '{haplotype}.polished.fasta'
#     threads: 16
#     resources:
#         mem_mb = 28000,
#         walltime = '2:30'
#     shell:
#         'racon -t {threads} {input.reads} {input.aln} {input.scaffolds} > {output}'
