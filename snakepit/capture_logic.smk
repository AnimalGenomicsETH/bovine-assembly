def have_parental_data(read_t):
    for parent in ('dam', 'sire'):
        if parent not in config['data'][config['animal']][read_t]:
            return False
    return True

def capture_logic(wildcards):
    required_files = [f'data/offspring.{wildcards.sample}.QC.txt']
    for assembler in config['assemblers']:
        a_path = PurePath(f'{assembler}_{wildcards.sample}')
        r_path = PurePath('results')
        if 'asm' in config['haplotypes']:

            required_files.extend([a_path / 'asm.scaffolds.fasta'])
            for extension in config['target_metrics']['general'] + config['target_metrics']['primary']:
                required_files.append(r_path / f'asm_{wildcards.sample}_{assembler}.{extension}')

        if 'trio' in config['haplotypes'] and have_parental_data('short_reads'):
            for N in (1,2):
                required_files.extend([a_path / f'hap{N}.scaffolds.fasta'])

                for extension in config['target_metrics']['general'] + config['target_metrics']['trio']:
                    required_files.append(r_path / f'hap{N}_{wildcards.sample}_{assembler}.{extension}')

                if have_parental_data('long_reads'):
                    required_files.append(a_path / f'hap{N}.scaff_seq')

            required_files.append(r_path / f'hap2_100_{assembler}.hap1.mumSV.txt')

        if 'parents' in config['haplotypes'] and have_parental_data('long_reads') and wildcards.sample == 100:

            for parent in ('dam', 'sire'):
                required_files.extend([a_path / f'{parent}.scaffolds.fasta'])

                for extension in config['target_metrics']['general'] + config['target_metrics']['parents']:
                    required_files.append(r_path / f'{parent}_{wildcards.sample}_{assembler}.{extension}')


            if 'trio' in config['haplotypes']:
                required_files.extend([r_path / f'hap2_100_{assembler}.dam.mumSV.txt', r_path / f'hap1_100_{assembler}.sire.mumSV.txt'])
                required_files.extend([r_path / f'hap2_100_{assembler}.dam.mm2.dot.png', r_path / f'hap1_100_{assembler}.sire.mm2.dot.png'])
    return map(str,required_files)
