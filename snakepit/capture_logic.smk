def have_parental_data(read_t):
    for parent in ('dam', 'sire'):
        if parent not in config['data'][config['animal']][read_t]:
            return False
    return True

def capture_logic(wildcards):
    sample = 100
    required_files = [f'data/offspring.{sample}.QC.txt']

    for assembler in config['assemblers']:
        a_path = PurePath(WORK_PATH.format(assembler=assembler,sample=sample))
        r_path = PurePath(RESULT_PATH).parent

        if 'asm' in config['haplotypes']:
            for extension in config['target_metrics']['general'] + config['target_metrics']['primary']:
                required_files.append(r_path / f'asm_{sample}_{assembler}.{extension}')

        if 'trio' in config['haplotypes'] and have_parental_data('short_reads'):
            for N in (1,2):
                for extension in config['target_metrics']['general'] + config['target_metrics']['trio']:
                    required_files.append(r_path / f'hap{N}_{sample}_{assembler}.{extension}')

                if have_parental_data('long_reads'):
                    required_files.append(a_path / f'hap{N}.filled.fasta')

            required_files.append(r_path / f'hap2_100_{assembler}.hap1.mumSV.txt')

        if 'parents' in config['haplotypes'] and have_parental_data('long_reads'):
            for parent in ('dam', 'sire'):
                for extension in config['target_metrics']['general'] + config['target_metrics']['parents']:
                    required_files.append(r_path / f'{parent}_{sample}_{assembler}.{extension}')

            if 'trio' in config['haplotypes']:
                required_files.extend([r_path / f'hap2_100_{assembler}.dam.mumSV.txt', r_path / f'hap1_100_{assembler}.sire.mumSV.txt'])
                required_files.extend([r_path / f'hap2_100_{assembler}.dam.mm2.dot.png', r_path / f'hap1_100_{assembler}.sire.mm2.dot.png'])

        for sample_rate in config['sampling']:
            if sample_rate == 100:
                continue
            for extension in config['target_metrics']['sampling']:
                for hap in ('asm',):#,'hap1','hap2'):
                    required_files.append(r_path / f'{hap}_{sample_rate}_{assembler}.{extension}')

    return map(str,required_files)
