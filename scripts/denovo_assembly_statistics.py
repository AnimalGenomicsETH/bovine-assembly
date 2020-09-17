import matplotlib.pyplot as plt
import datetime
import pandas
import seaborn

def plot_chromosomes_scaffolds(animal,assembler):
    counts = []
    with open (f'results/{animal}_{assembler}.gaps.txt','r') as file_in:
        for line in file_in:
            if 'CH' not in line:
                continue
            parts = line.rstrip().split()
            counts.append([int(i) for i in parts[2].split(',')])

    length = max(map(len, counts))
    grid = numpy.array([xi+[None]*(length-len(xi)) for xi in counts]).T
    bottoms = grid.cumsum(axis=0)
    bottoms = np.insert(bottoms,0,np.zeros(len(bottoms[0])),axis=0)
    
    fig = plt.figure()

def load_auNCurves(animal,assembler):
    auN_values, metrics = [[],[]], dict()

    with open(f'results/{animal}_{assembler}.auN.txt','r') as file_in:
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

from numpy import linspace
def plot_auNCurves(animal,assembler):
    data, metrics = load_auNCurves(animal,assembler)
    auN_data, aln_metrics = load_NGA(animal,assembler)

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
    fig.savefig(f'{animal}_{assembler}_auN_curves.png')

    return metrics, aln_metrics

def load_NGA(animal,assembler):
    auN_data, data = [], dict()
    with open(f'results/{animal}_{assembler}.NGA50.txt','r') as file_in:
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

def kmer_QV(animal,assembler):
    data = dict()
    with open(f'results/{animal}_{assembler}.merqury.stats.txt','r') as file_in:
        for line in file_in:
            data[line[:2]] = line.rstrip().split()[1]
    return data

def load_asmgene(animal,assembler):
    table = [[],[],[]]
    with open(f'results/{animal}_{assembler}.asmgene.txt','r') as file_in:
        for line in file_in:
            if line[0] != 'X':
                continue
            parts = line.rstrip().split()
            for i, p in enumerate(parts[1:]):
                table[i].append(p)
    return table

def busco_report(animal,assembler):
    with open(f'results/{animal}_{assembler}.BUSCO.txt') as file_in:
        for line in file_in:
            if 'lineage dataset' in line:
                LD_set = line.split()[5]
            if 'C:' in line:
                return LD_set, line.strip().rstrip()

def load_resource_benchmark(animal,assembler):
    info_calls = {'CPU time':'cputime', 'Run time':'walltime', 'Max Memory':'max_mem', 'Average Memory':'mean_mem', 'Delta Memory':'delta_mem'}
    data = dict()
    reached_resources = False
    with open(f'logs/assembler_{assembler}/animal-{animal}.out','r') as benchmark:
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

def generate_markdown_reads(animal):
    build_str = f'Animal ID: **{animal}**\n\n' \
                 f'Time of report generation: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n'
    
    df = pandas.read_csv(f'data/{animal}.QC.txt')
    stats = df.describe().transpose()

    build_str += '## read metrics \n' \
                 f'Total reads: {len(df):,}\n\n'
  
    build_str += f'| | {" | ".join(stats.columns[1:])} |\n' + \
                 ' -- '.join('|'*(len(stats.columns)+1)) +'\n'                 

    for row in ('length','quality'):
        build_str += f'| {row} | {" | ".join(map("{:.2f}".format,stats.loc[row][1:]))} |\n'
    build_str += '\n'

    seaborn.jointplot(data=df,x='length',y='quality',kind='hex',joint_kws={'bins':'log'}).savefig('test.png')

    build_str += IMAGE('test.png',.6) + '\n\n'
    return build_str

def generate_markdown_string(animal,assembler,build_str,summary_str):
    build_str += '\n\n---\n\n' \
                f'# assembler: *{assembler}*\n'

    resource_stats = load_resource_benchmark(animal,assembler)

    build_str += '## resource details\n' \
                 f'Runtime (wall/cpu): {resource_stats["walltime"]} s / {resource_stats["cputime"]} s\n\n' \
                 f'>*average threading* = {resource_stats["cputime"]/resource_stats["walltime"]:.1f}\n\n' \
                 f'Memory (mean/max): {resource_stats["mean_mem"]} mb / {resource_stats["max_mem"]} mb\n\n'

    build_str += '## assembly metrics\n'

    asm_metrics, aln_metrics = plot_auNCurves(animal,assembler)
    build_str += f'Genome length: {asm_metrics["SZ"]/1e9:.4f} gb\n\n' \
                 f'Total contigs: {asm_metrics["NN"]:,}\n\n'

    build_str += f'Contig length and quantity\n\n' \
                 f'**N50**: {asm_metrics["N50"]/1e6:.2f} mb\n\n' \
                 f'auN value: {asm_metrics["AU"]}\n\n' + \
                 IMAGE(f'{animal}_{assembler}_auN_curves.png',.8) + '\n\n'

    if assembler == 'canu':
        build_str += 'Purged coverage:\n\n' + \
                     IMAGE(f'{assembler}/{animal}.contigs_raw.spectra.png',.4) + \
                     IMAGE(f'{assembler}/{animal}.purged.spectra.png',.4) + '\n\n'

        #              ' '.join(IMAGE(f'canu/{animal}.{sample}.spectra.png,',.4) for sample in ('contigs_raw','purged')) + '\n\n'
        
        print(' '.join(IMAGE(f'canu/{animal}.{sample}.spectra.png,',.4) for sample in ('contigs_raw','purged')))
        #'canu/{animal}.contigs_raw.spectra.png',.4) + ' ' + IMAGE('canu/{animal}.purged.spectra.png',.4) + '\n\n'

    build_str += '### reference metrics\n' \
                 '* Coverage\n' \
                 f'  * Rcov: {aln_metrics["Rcov"]}\n' \
                 f'  * Rdup: {aln_metrics["Rdup"]}\n' \
                 f'  * Qcov: {aln_metrics["Qcov"]}\n' \
                 '* Contigs\n' \
                 f'  * breaks: {aln_metrics["#breaks"]}\n' \
                 f'  * auNGA: {aln_metrics["AUNGA"]}\n\n' + \
                 IMAGE(f'results/{animal}_{assembler}.dot.png',.7) + '\n\n'

    build_str += '## validation results\n\n'

    kmer_stats = kmer_QV(animal,assembler)
    print(kmer_stats)
    build_str += '### merqury k-mers\n' \
                 f'Coverage: {float(kmer_stats["CV"])/100:.1%}\n\n' \
                 f'QV: {kmer_stats["QV"]}\n\n' + \
                 IMAGE(f'{assembler}/{animal}.spectra-asm.ln.png',.45) + '\n\n'

    lineage, busco_string = busco_report(animal,assembler)
    build_str += '### BUSCO \n' \
                 f'Lineage: **{lineage}**\n\nAnalysis: {busco_string}\n\n'

    build_str += '### asmgene\n'

    gene_map = load_asmgene(animal,assembler)
    build_str += f'| | {" | ".join(gene_map[0])} |\n' \
                 '| -- | -- | -- |\n'

    for row, values in zip(('ref','asm'),gene_map[1:]):
        build_str += f'| {row} | {" | ".join(values)} |\n'

    summary_str += f'| *{assembler}* | {asm_metrics["SZ"]/1e9:.2f} | {asm_metrics["NN"]:,} | ' \
                   f'{asm_metrics["N50"]/1e6:.2f} | {asm_metrics["L50"]} | {busco_string[2:7]} | {float(kmer_stats["QV"]):.1f} |\n'


    return build_str, summary_str

import argparse
from pathlib import Path
from markdown2 import markdown
from weasyprint import HTML, CSS

def custom_PDF_writer(output,prepend_str,md_content,css):
    header = markdown(prepend_str,extras=['tables'])
    raw_html = markdown(md_content, extras={'tables':{},'header_ids':{},'toc':{'depth':2},'code-friendly':{},'cuddled-lists':{}})
    full_html = header + raw_html.toc_html + raw_html
    print(full_html)
    html = HTML(string=full_html,base_url=str(Path().cwd()))
    html.write_pdf(output,stylesheets=[CSS(filename=css)])

def main(direct_input=None):
    parser = argparse.ArgumentParser(description='Produce assembly report.')
    parser.add_argument('--animal', type=str, required=True)
    parser.add_argument('--assemblers', nargs='+', required=True)
    parser.add_argument('--outfile', default='assembly_report.pdf', type=str)
    parser.add_argument('--keepfig', action='store_true')
    parser.add_argument('--css', default='', type=str)

    args = parser.parse_args(direct_input)
    css_path = args.css if Path(args.css).is_file() else ''

    prepend_str = generate_markdown_reads(args.animal) + \
                  '# table of contents'

    md_string = ''
    summary_string = '# summary \n' \
                     '| assembler | size | contigs | N50 | L50 | completeness | QV kmer |\n' \
                     '| :-------- | :--: | :-----: | :-: | :-: | :----------: | :-----: |\n'
    
    for assembler in args.assemblers:
        md_string, summary_string = generate_markdown_string(args.animal,assembler,md_string,summary_string)

    md_string = summary_string + md_string
    custom_PDF_writer(args.outfile,prepend_str,md_string,css_path)
    if not args.keepfig:
        for assembler in args.assemblers:
            Path(f'{args.animal}_{assembler}_auN_curves.png').unlink(missing_ok=True)

if __name__ == "__main__":
    main()
