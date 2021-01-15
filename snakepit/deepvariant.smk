BIN_VERSION='1.1.0'
INPUT='input'
OUTPUT='output'

rule deepvariant_make_examples:
    input:
        ''
    output:
        ''
    shell:
        '''
        time seq 0 23 | parallel -q --halt 2 --line-buffer /opt/deepvariant/bin/make_examples \
        --mode calling --ref "/input/asm.fasta" --reads "/input/hifi.bam" \
        --examples "/output/intermediate_results_dir/make_examples.tfrecord@24.gz" --add_hp_channel \
        --alt_aligned_pileup "diff_channels" --gvcf "/output/intermediate_results_dir/gvcf.tfrecord@24.gz" \
        --noparse_sam_aux_fields --norealign_reads --nosort_by_haplotypes --vsc_min_fraction_indels "0.12" --task {}
        '''

rule deepvariant_call_variants:
    input:
        ''
    output:
        ''
    shell:
        '''
        singularity run -B ${INPUT}:/input,${OUTPUT}:/output,${OUTPUT}/intermediate_results_dir:/output/intermediate_results_dir,$TMPDIR:$TMPDIR \
          deepvariant_1.1.0.sif \
          /opt/deepvariant/bin/call_variants \
          --outfile "/output/intermediate_results_dir/call_variants_output.tfrecord.gz" \
          --examples "/output/intermediate_results_dir/make_examples.tfrecord@24.gz" \
          --checkpoint "/opt/models/pacbio/model.ckpt"
        '''

rule deepvaraint_postprocess:



mkdir -p ${OUTPUT}
mkdir -p ${OUTPUT}/intermediate_results_dir
THREADS=16
