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
    with open('busco_short_summary.txt') as file_in:
        for line in file_in:
            if 'lineage dataset' in line:
                LD_set = line.split()[5]
            if 'C:' in line:
                return LD_set, line.strip().rstrip()

def load_resource_benchmark(jobname):
    info_calls = {'CPU time':'cputime', 'Run time':'walltime', 'Max Memory':'max_mem', 'Average Memory':'mean_mem', 'Delta Memory':'delta_mem'}
    data = dict()
    reached_resources = False 
    with open(jobname,'r') as benchmark:
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
            

import os
def generate_markdown_string(args):
    build_str = '# Assembly Report\n'

    build_str += f'Animal ID: **{args.animal}**\n\n' \
                 f'Time of report generation: {datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")}\n'

    resource_stats = load_resource_benchmark(args.jobname)

    build_str += '## Assembler details\n' \
                 f'Assembler: **{args.assembler}**\n\n' \
                 f'Runtime (wall/cpu): {resource_stats["walltime"]}s / {resource_stats["cputime"]}s\n\n' \
                 f'>_average threading_ = {resource_stats["cputime"]/resource_stats["walltime"]:.3f}\n\n' \
                 f'Memory (mean/max): {resource_stats["mean_mem"]}mb / {resource_stats["max_mem"]}mb\n\n'

    build_str += '## Assembly metrics\n'

    asm_metrics = plot_auNCurve(args.animal)
    build_str += f'Genome length: {asm_metrics["SZ"]/1e9:.2f}gb\n\n' \
                 f'Total contigs: {asm_metrics["NN"]}\n\n'

    build_str += f'Contig length and quantity\n\n' \
                 '![alt text](auN_curve.png)' \
                 'apples on a stick\n\n' \
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

    #with open('assembly_report.md','w') as file_out:
    return build_str
    

import argparse
from md2pdf.core import md2pdf

def main():
    parser = argparse.ArgumentParser(description='Produce assembly report.')
    parser.add_argument('--animal', default='', type=str)
    parser.add_argument('--assembler', default='', type=str)
    parser.add_argument('--outfile', default='assembly_report.pdf', type=str)
    parser.add_argument('--jobname', default='', type=str)
    parser.add_argument('--keepfig', default=False, type=bool)
    

    args = parser.parse_args()
    css_path = 'github.css' if os.path.isfile('github.css') else ''
    

    md_string = generate_markdown_string(args)
    md2pdf(args.outfile,md_content=md_string,css_file_path=css_path,base_url=os.getcwd())
    if not args.keepfig:
        os.remove('auN_curve.png')
    
if __name__ == "__main__":
    main()



                
                
            
    
    
