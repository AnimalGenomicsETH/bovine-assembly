import matplotlib.pyplot as plt
import datetime

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
    with open(f'results/{animal}_{assembler}.asm-ccs.qv.txt','r') as file_in:
        for line in file_in:
            if line[:2] == 'CV':
                coverage = float(line.rstrip().split()[1])
            elif line[:2] == 'QV':
                raw_QV, adjusted_QV = (float(i) for i in line.rstrip().split()[1:])
    return (coverage, raw_QV, adjusted_QV)

def load_asmgene(animal,assembler):
    table = dict()
    with open(f'results/{animal}_{assembler}.asmgene.txt','r') as file_in:
        for line in file_in:
            if line[0] != 'X':
                continue
            parts = line.rstrip().split()
            table[parts[1]] = tuple(int(i) for i in parts[2:])
    return table

def busco_report(animal,assembler):
    with open(f'results/{animal}_{assembler}_busco_short_summary.txt') as file_in:
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

def generate_markdown_string(animal,assembler,build_str,summary_str):
    build_str += '\n\n---\n\n' \
                f'# assembler: *{assembler}*\n'

    resource_stats = load_resource_benchmark(animal,assembler)

    build_str += '## resource details\n' \
                 f'Runtime (wall/cpu): {resource_stats["walltime"]}s / {resource_stats["cputime"]}s\n\n' \
                 f'>*average threading* = {resource_stats["cputime"]/resource_stats["walltime"]:.1f}\n\n' \
                 f'Memory (mean/max): {resource_stats["mean_mem"]}mb / {resource_stats["max_mem"]}mb\n\n'

    build_str += '## assembly metrics\n'

    asm_metrics, aln_metrics = plot_auNCurves(animal,assembler)
    build_str += f'Genome length: {asm_metrics["SZ"]/1e9:.4f}gb\n\n' \
                 f'Total contigs: {asm_metrics["NN"]}\n\n'

    build_str += f'Contig length and quantity\n\n' \
                 f'**N50** {asm_metrics["N50"]/1e6:.2f}mb\n\n' \
                 f'![auN curve]({animal}_{assembler}_auN_curves.png)' \
                 f'auN value: {asm_metrics["AU"]}\n\n'

    if assembler == 'canu':
        build_str += 'Purged coverage:\n' \
                     f'![raw](canu/{animal}.contigs_raw.spectra.png) ![purged](canu/{animal}.purged.spectra.png)\n\n'

    build_str += '### reference metrics\n' \
                 '* Coverage\n' \
                 f'  * Rcov: {aln_metrics["Rcov"]}\n' \
                 f'  * Rdup: {aln_metrics["Rdup"]}\n' \
                 f'  * Qcov: {aln_metrics["Qcov"]}\n' \
                 '* Contigs\n' \
                 f'  * breaks: {aln_metrics["#breaks"]}\n' \
                 f'  * auNGA: {aln_metrics["AUNGA"]}\n\n'

    build_str += '## validation results\n\n'

    kmer_stats = kmer_QV(animal,assembler)
    build_str += '### k-mer spectra\n' \
                 f'Coverage: {kmer_stats[0]:.1%}\n\n' \
                 f'QV (raw/adjusted): {kmer_stats[1]} / {kmer_stats[2]}\n\n'

    lineage, busco_string = busco_report(animal,assembler)
    build_str += '### BUSCO \n' \
                 f'Lineage: **{lineage}**\n\nAnalysis: {busco_string}\n\n'

    build_str += '### asmgene\n'

    gene_map = load_asmgene(animal,assembler)
    build_str += '| | Ref | asm |\n' \
                 '| -- | -- | -- |\n'

    for row, values in gene_map.items():
        build_str += f'| {row} | {" | ".join(map(str, values))} |\n'

    summary_str += f'| **{assembler}** | {asm_metrics["SZ"]/1e9:.2f} | {asm_metrics["NN"]} | ' \
                   f'{asm_metrics["N50"]} | {asm_metrics["L50"]} | {busco_string[10:15]} | {kmer_stats[2]:.1f} |\n'

    return build_str, summary_str

import argparse
from pathlib import Path
from markdown2 import markdown
from weasyprint import HTML, CSS

def custom_PDF_writer(output,prepend_str,md_content,css):
    header = markdown(prepend_str)
    raw_html = markdown(md_content, extras=['tables','header_ids','toc','code-friendly','cuddled-lists'])
    full_html = header + raw_html.toc_html + raw_html
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

    prepend_str = f'Animal ID: **{args.animal}**\n\n' \
                  f'Time of report generation: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n\n' \
                  '# table of contents'

    md_string = ''
    summary_string = '# summary \n' \
                     '| assembler | size | contigs | N50 | L50 | completeness | QV |\n' \
                     '| :-------- | ---- | ------- | --- | --- | ------------ | -- |\n'

    for assembler in args.assemblers:
        md_string, summary_string = generate_markdown_string(args.animal,assembler,md_string,summary_string)

    md_string = summary_string + md_string
    custom_PDF_writer(args.outfile,prepend_str,md_string,css_path)
    if not args.keepfig:
        for assembler in args.assemblers:
            Path(f'{args.animal}_{assembler}_auN_curves.png').unlink(missing_ok=True)

if __name__ == "__main__":
    main()
