import matplotlib.pyplot as plt
import datetime

def load_auNCurve(animal):
    auN_values, metrics = [[],[]], dict()
    
    with open(f'{animal}.asm.auN.txt','r') as file_in:
        for line in file_in:
            if line[:2] == 'NL':
                x, Nx, Lx  = (int(i) for i in line.rstrip().split()[1:])
                auN_values[0].append(Nx)
                auN_values[1].append(Lx)
                #[x//10] = (Nx,Lx)
            elif line[:2] != 'CC':
                metrics[line[:2]] = int(line.split()[-1])
                 
    return auN_values, metrics
    
            
def plot_auNCurve(animal):
    data, metrics = load_auNCurve(animal)
    
    fig, (ax_N,ax_L) = plt.subplots(2,sharex=True,figsize=(7, 4), dpi=300)
    
    x_vals = list(range(0,101,10))
    ax_N.plot(x_vals,data[0],'r')
    ax_L.plot(x_vals,data[1],'b')

    for ax in (ax_N,ax_L):
        ax.set_yscale('log')
        ax.set_xlabel('x',fontsize=16)

    ax_N.set_title('Nx',fontsize=18)
    ax_N.set_ylabel('contig length',fontsize=14)
    
    ax_L.set_title('Lx',fontsize=18)
    ax_L.set_ylabel('number of contigs',fontsize=14)

    plt.tight_layout()
    fig.savefig('auN_curve.png')

    return metrics
    #plt.show(block=False)

def load_NGA(animal):
    data = dict()
    with open(f'{animal}.NGA50.txt','r') as file_in:
        for line in file_in:
            (key, value) = line.rstrip().split()
            data[key] = value
    return data

def kmer_QV(animal):
    with open(f'{animal}.asm-ccs.qv.txt','r') as file_in:
        for line in file_in:
            if line[:2] == 'CV':
                coverage = float(line.rstrip().split()[1])
            elif line[:2] == 'QV':
                raw_QV, adjusted_QV = (float(i) for i in line.rstrip().split()[1:])
    return (coverage, raw_QV, adjusted_QV)

import glob
def busco_report():
    for filename in glob.glob('short_summary.*.busco_results.txt'):
        with open(filename) as file_in:
            for line in file_in:
                if 'lineage dataset' in line:
                    LD_set = line.split()[5]
                if 'C:' in line:
                    return LD_set, line.strip().rstrip()

            
#pip install md2pdf
import os
def generate_markdown_string(args):
    build_str = '# Assembly Report\n'

    build_str += f'Animal ID: **{args.animal}**\n\n' \
                 f'Time of report generation: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n'

    build_str += '## Assembler details\n' \
                 f'Assembler: **{args.assembler}**\n\n' \
                 f'Runtime (wall/cpu): {args.walltime} / {args.cputime}\n\n' \
                 f'>_average threading_ = {args.cputime/args.walltime}\n\n' \
                 f'Memory (mean/max): {args.mean_mem} / {args.max_mem}\n\n'

    build_str += '## Assembly metrics\n'

    asm_metrics = plot_auNCurve(args.animal)
    build_str += f'Genome length: {asm_metrics["SZ"]/1e9:.2f}gb\n\n' \
                 f'Total contigs: {asm_metrics["NN"]}\n\n'
    
    build_str += f'Contig length and quantity\n\n' \
                 '![alt text](auN_curve.png)\n\n' \
                 f'auN value: {asm_metrics["AU"]}\n\n'

    NGA_data = load_NGA(args.animal)

    build_str += '---\n\n'
    build_str += 'Reference metrics\n' \
                 '+ Coverage\n' \
                 f'  + Rcov: {NGA_data["Rcov"]}\n' \
                 f'  + Qcov: {NGA_data["Qcov"]}\n' \
                 '+ Contig lengths\n' \
                 f'  + NG50: {NGA_data["NG50"]}\n' \
                 f'  + NGA50: {NGA_data["NGA50"]}\n' \
                 f'  + auNGA: {NGA_data["AUNGA"]}\n' \

    kmer_stats = kmer_QV(args.animal)
    build_str += '### K-mer validation\n' \
                 f'Coverage: {kmer_stats[0]:.1%}\n\n' \
                 f'QV (raw/adjusted): {kmer_stats[1]} / {kmer_stats[2]}\n\n'

    lineage, busco_string = busco_report()
    build_str += '### BUSCO analysis\n' \
                 f'Lineage: {lineage}\n\nAnalysis: {busco_string}\n'

    

    with open('assembly_report.md','w') as file_out:
        file_out.write(build_str)
    

import argparse
from md2pdf.core import md2pdf

def main():
    parser = argparse.ArgumentParser(description='Produce assembly report.')
    parser.add_argument('--animal', default='', type=str)
    parser.add_argument('--assembler', default='', type=str)
    parser.add_argument('--walltime', default=1, type=int)
    parser.add_argument('--cputime', default=0, type=int)
    parser.add_argument('--max_mem', default=0, type=int)
    parser.add_argument('--mean_mem', default=0, type=int)

    args = parser.parse_args()
    generate_markdown_string(args)

    css_path = 'github.css' if os.path.isfile('github.css') else ''
    md2pdf('out2.pdf',md_file_path='assembly_report.md',css_file_path=css_path,base_url=os.getcwd())
    
if __name__ == "__main__":
    main()

    
