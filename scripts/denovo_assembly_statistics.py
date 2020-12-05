import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import datetime
import pandas
import seaborn
from itertools import cycle
from numpy import linspace, cumsum
from pathlib import PurePath

animal, haplotype, sample, assembler = '', '', '100', ''

def get_basename():
    return f'{haplotype}_{sample}_{assembler}'

def open_results(extension):
    return open(f'results/{get_basename()}.{extension}','r')

def save_figure(figure,path):
    figure.savefig(path)
    plt.close(figure)

def plot_chromosome_scaffolds():
    chromosome = 0
    fig, ax = plt.subplots()

    with open_results('gaps.txt') as file_in:
        for line in file_in:
            if 'CM' not in line:
                continue
            chromosome += 1
            parts = line.rstrip().split()
            scaffolds = [int(i) for i in parts[2].split(',')]
            h_sums = [0] + list(cumsum(scaffolds))

            colours = cycle(('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d','#666666'))
            for height, bottom in zip(scaffolds,h_sums):
                ax.bar(chromosome,height,bottom=bottom,width=0.5,color=next(colours))

    f_name = f'figures/{get_basename()}.chromosomes.png'
    save_figure(fig,f_name)
    return f_name

def load_auNCurves(d_type):
    auN_values, metrics = [[],[]], dict()

    with open_results(f'{d_type}.auN.txt') as file_in:
        for line in file_in:
            if line[:2] == 'NL':
                x, Nx, Lx  = (int(i) for i in line.rstrip().split()[1:])
                auN_values[0].append(Nx)
                auN_values[1].append(Lx)
                if x == 50:
                    metrics['N50'] = Nx
                    metrics['L50'] = Lx
            elif line[:2] != 'CC':
                metrics[line[:2]] = int(line.split()[-1])

    return auN_values, metrics

def plot_auNCurves():
    data, metrics = load_auNCurves('contigs')
    auN_data, aln_metrics = load_NGA()

    fig, (ax_N,ax_L) = plt.subplots(1,2,sharex=True,figsize=(6, 4))

    x_vals = linspace(0,100,len(data[0]))
    ax_N.plot(x_vals,data[0],'forestgreen',label='Nx')
    ax_N.plot(x_vals,auN_data[0],'darkorange',label='NGx')
    ax_N.plot(x_vals,auN_data[1],'darkmagenta',label='NGAx')
    ax_L.plot(x_vals,data[1],'forestgreen')

    for ax in (ax_N,ax_L):
        ax.set_yscale('log')
        ax.set_xlabel('x',fontsize=16)

    ax_N.set_title('Nx',fontsize=18)
    ax_N.set_ylabel('contig length',fontsize=14)
    ax_N.legend()

    ax_L.set_title('Lx',fontsize=18)
    ax_L.set_ylabel('number of contigs',fontsize=14)

    plt.tight_layout()
    save_path = f'figures/{get_basename()}_auN_curves.png'
    save_figure(fig,save_path)
    return metrics, aln_metrics, save_path

def load_NGA():
    auN_data, data = [], dict()
    with open_results('NGA50.txt') as file_in:
        for line in file_in:
            if line[:2] == 'NG':
                try:
                    auN_data.append(int(line.split()[1]))
                except ValueError:
                    auN_data.append(0)
                    print(f'Undefined result for {line.split()[0][2:]}')
            else:
                (key, value) = line.rstrip().split()
                data[key] = value
    return (auN_data[:len(auN_data)//2],auN_data[len(auN_data)//2:]), data

def plot_sampling_curves(df):
    pass


def load_key_pair_file(fname):
    return {key:value for (key,*value) in (line.rstrip().split() for line in open_results(fname))}

def kmer_QV(run_type='txt'):
    return {k:v[0] for (k,v) in load_key_pair_file(f'merqury.{run_type}.stats').items()}

def load_mumSV(reference='ref',key='chrm'):
    return {v[0]:v[1:] for (k,*v) in [line.rstrip().split() for line in open_results(f'{reference}.mumSV.txt')] if k == key}

def load_asmgene(threshold=97):
    table = [[],[],[]]
    with open_results(f'asmgene.{threshold}.txt') as file_in:
        for line in file_in:
            if line[0] != 'X':
                continue
            parts = line.rstrip().split()
            for i, p in enumerate(parts[1:]):
                table[i].append(p)
    return table

def busco_report():
    with open_results('BUSCO.txt') as file_in:
        for line in file_in:
            if 'lineage dataset' in line:
                LD_set = line.split()[5]
            if 'C:' in line:
                return LD_set, line.strip().rstrip()

def load_resource_benchmark(animal,assembler,sample):
    info_calls = {'CPU time':'cputime', 'Run time':'walltime', 'Max Memory':'max_mem', 'Average Memory':'mean_mem', 'Delta Memory':'delta_mem'}
    data = dict()
    reached_resources = False
    with open(f'logs/assembler_{assembler}/sample-{sample}.animal-{animal}.out','r') as benchmark:
        for line in benchmark:
            if not reached_resources:
                reached_resources = 'Resource usage summary:' in line
            elif len(data) == len(info_calls):
                return data
            elif len(line) > 1:
                code, raw_val = line.strip().split(' :')
                if code not in info_calls:
                    continue
                val = int(float(raw_val.strip().split()[0]))
                data[info_calls[code]] = val
    #default values
    data = {code:0 for code in info_calls.values()}
    data['walltime'] = 1
    return data

def img_sizer(width,dpi=96):
    full_width = (8.25-1) * dpi #A4 size - margin
    if isinstance(width,float):
        width = int(width*full_width)
    return f'style="object-fit:cover;width:{width}px;height:100%;"'

def IMAGE(path,scale):
    return f'<img src="{path}" {img_sizer(scale)} />'

def generate_markdown_reads():
    build_str = f'Animal ID: **{animal}**\n\n' \
                f'Sampled at: {sample}\n\n' \
                f'Time of report generation: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n'

    df = pandas.read_csv(f'data/offspring.{sample}.QC.txt')
    stats = df.describe().transpose()

    build_str += '## read metrics \n' \
                 f'Total reads: {len(df):,}\n\n'

    build_str += f'| | {" | ".join(stats.columns[1:])} |\n' + \
                 ' -- '.join('|'*(len(stats.columns)+1)) +'\n'

    for row in ('length',' quality'):
        build_str += f'| {row} | {" | ".join(map("{:.2f}".format,stats.loc[row][1:]))} |\n'
    build_str += '\n'

    seaborn.jointplot(data=df,x='length',y=' quality',kind='hex',joint_kws={'bins':'log'}).savefig(f'figures/offspring.{sample}.QC.png')

    build_str += IMAGE(f'figures/offspring.{sample}.QC.png',.6) + '\n\n'
    return build_str

def emph_haplotype(haplotype):
    if haplotype == 'asm':
        return f'<span style="color:blue"> **asm** </span>'
    else:
        return haplotype

def generate_markdown_string(summary_str,build_str=None):
    asm_metrics = load_auNCurves('contigs')[1]
    scaff_metrics = load_auNCurves('scaffolds')[1]
    lineage, busco_string = busco_report()
    kmer_stats = kmer_QV('full' if haplotype in ('asm','hap1','hap2') else 'simple')

    summary_str += f'| {assembler} | {emph_haplotype(haplotype) if build_str is not None else sample} | {asm_metrics["SZ"]/1e9:.2f} | {asm_metrics["NN"]:,} | ' \
                   f'{asm_metrics["N50"]/1e6:.2f} | {asm_metrics["L50"]} | {scaff_metrics["N50"]/1e6:.2f} | {busco_string[2:7]} | {float(kmer_stats["QV"]):.1f} |\n'

    if build_str is None:
        return summary_str, {'S50':scaff_metrics["N50"],'N50':asm_metrics["N50"],'QV':kmer_stats["QV"],'completness':kmer_stats['completness'],'BUSCO_C':busco_string[2:6],'BUSCO_S':busco_string[10:14]}

    build_str += '\n\n---\n\n' \
                f'# assembler: *{assembler}*, haplotype: {haplotype} \n'

    asm_metrics, aln_metrics, auN_plot = plot_auNCurves()
    build_str += '## assembly metrics\n'

    build_str += f'Genome length: {asm_metrics["SZ"]/1e9:.4f} gb\n\n' \
                 f'Total contigs: {asm_metrics["NN"]:,}\n\n'

    build_str += '### scaffolded chromosomes\n' + \
                 IMAGE(plot_chromosome_scaffolds(),.5) + '\n\n'

    build_str += f'Contig length and quantity\n\n' \
                 f'**N50**: {asm_metrics["N50"]/1e6:.2f} mb\n\n' \
                 f'auN value: {asm_metrics["AU"]}\n\n' + \
                 IMAGE(auN_plot,.8) + '\n\n'

    if assembler == 'canu':
        build_str += 'Purged coverage:\n\n' + \
                     IMAGE(f'{assembler}_{sample}/{haplotype}.contigs_raw.spectra.png',.4) + \
                     IMAGE(f'{assembler}_{sample}/{haplotype}.purged.spectra.png',.4) + '\n\n'

    build_str += '### reference metrics\n' \
                 '* Coverage\n' \
                 f'  * Rcov: {aln_metrics["Rcov"]}\n' \
                 f'  * Rdup: {aln_metrics["Rdup"]}\n' \
                 f'  * Qcov: {aln_metrics["Qcov"]}\n' \
                 '* Contigs\n' \
                 f'  * breaks: {aln_metrics["#breaks"]}\n' \
                 f'  * auNGA: {aln_metrics["AUNGA"]}\n\n' + \
                 IMAGE(f'results/{get_basename()}.dot.png',.7) + '\n\n'

    build_str += '## validation results\n\n'

    build_str += '### merqury k-mers\n' \
                 f'Coverage: {float(kmer_stats["completeness"])/100:.1%}\n\n' \
                 f'QV: {kmer_stats["QV"]}\n\n' + \
                 IMAGE(f'{assembler}_{sample}/{haplotype}.spectra-cn.ln.png',.4) + \
                 IMAGE(f'{assembler}_{sample}/{haplotype}.spectra-asm.ln.png',.4) + '\n\n'

    if haplotype not in ('dam','sire'):
        build_str += f'Switch error: {kmer_stats["switches"]}\n\n' \
                     f'dam: {kmer_stats["dam"]}\n\n' \
                     f'sire: {kmer_stats["sire"]}\n\n' + \
                     IMAGE(f'{assembler}_{sample}/{haplotype}.{haplotype}.contigs.block.NG.png',.45) + \
                     IMAGE(f'{assembler}_{sample}/{haplotype}.{haplotype}.contigs.continuity.NG.png',.45) + '\n\n'

    build_str += '### BUSCO \n' \
                 f'Lineage: **{lineage}**\n\nAnalysis: {busco_string}\n\n'

    mumsv = load_mumSV('ref')
    build_str += '### structural variants\n' + \
                 '\n\n'.join(f'variant: {variant}, count: {count}, total size: {int(size)/1e6:.1f} mB' for (variant, (count, size)) in mumsv.items()) + '\n\n'

    build_str += '### repeat content\n' \
                 'work coming soon\n\n'

    build_str += '### asmgene\n'

    gene_map = load_asmgene()
    build_str += f'| | {" | ".join(gene_map[0])} |\n' \
                 '| -- | -- | -- |\n'

    for row, values in zip(('ref','asm'),gene_map[1:]):
        build_str += f'| {row} | {" | ".join(values)} |\n'

    return summary_str, build_str

import argparse
from pathlib import Path
from markdown2 import markdown
from weasyprint import HTML, CSS
from shutil import rmtree
from itertools import product

def custom_PDF_writer(output,prepend_str,md_content,css):
    header = markdown(prepend_str,extras=['tables'])
    raw_html = markdown(md_content, extras={'tables':{},'header_ids':{},'toc':{'depth':2},'code-friendly':{},'cuddled-lists':{}})
    full_html = header + raw_html.toc_html + raw_html
    html = HTML(string=full_html,base_url=str(Path().cwd()))
    if css:
        html.write_pdf(output,stylesheets=[CSS(filename=css)])
    else:
        html.write_pdf(output)

def main(direct_input=None):
    Path('figures').mkdir(exist_ok=True)

    parser = argparse.ArgumentParser(description='Produce assembly report.')
    parser.add_argument('--animal', type=str, required=True)
    parser.add_argument('--samples', nargs='+', required=True)
    parser.add_argument('--input', nargs='+', required=True)
    parser.add_argument('--outfile', default='assembly_report.pdf', type=str)
    parser.add_argument('--keepfig', action='store_true')
    parser.add_argument('--css', default='report.css', type=str)

    args = parser.parse_args(direct_input)
    global animal
    animal = args.animal
    css_path = args.css if Path(args.css).is_file() else None

    prepend_str = generate_markdown_reads() + \
                  '# table of contents'

    md_string = ''
    summary_string = '# summary \n' \
                     '| assembler | haplotype | size | contigs | N50 | L50 | S50 | BUSCO | QV |\n' \
                     '| --------- | --------- | ---- | ------- | --- | --- | --- | ----- | -- |\n'

    sample_string = summary_string.replace('# summary ','# downsampling ').replace('haplotype','sample rate')

    contig_files = [PurePath(f).name.split('.')[0] for f in args.input if '.contigs.auN.txt' in f]
    haplotypes = sorted({ctg.split('_')[0] for ctg in contig_files})
    assemblers = sorted({ctg.split('_')[-1] for ctg in contig_files})

    global haplotype
    global assembler
    for haplotype_t, assembler_t in product(haplotypes,assemblers):
        haplotype = haplotype_t
        assembler = assembler_t

        summary_string, md_string = generate_markdown_string(summary_string,md_string)

    md_string = summary_string + md_string

    global sample
    row_data = []
    if len(args.samples) > 1:
        haplotype = 'asm'
        for sample_t, assembler_t in product(args.samples,assemblers):
            sample = sample_rate
            assembler = assembler_t
            summary_string, new_row = generate_markdown_string(summary_string)
            new_row.update({'sample':sample_t,'assembler':assembler_t})
            row_data.append(new_row)

        md_string += summary_string
        md_string += 'figure'

    custom_PDF_writer(args.outfile,prepend_str,md_string,css_path)
    if not args.keepfig:
        rmtree('figures')

if __name__ == "__main__":
    main()
