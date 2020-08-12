import matplotlib.pyplot as plt
import numpy as np

def load_auNCurve():
    auN_values = np.empty((11,2))
    
    with open('BSWCHEF120152514636.asm.auN.txt','r') as file_in:
        for line in file_in:
            if line[:2] == 'NL':
                x, Nx, Lx  = (int(i) for i in line.rstrip().split()[1:])
                auN_values[x//10] = (Nx,Lx)
    return auN_values.T
    
            
def plot_auNCurve(data=None):
    data = data or load_auNCurve()
    
    fig, (ax_N,ax_L) = plt.subplots(1,2,sharex=True)
    
    x_vals = np.linspace(0,100,data.shape[1])
    ax_N.plot(x_vals,data[0],'r')
    ax_L.plot(x_vals,data[1],'b')

    for ax in (ax_N,ax_L):
        ax.set_yscale('log')
        ax.set_xlabel('x')

    ax_N.set_title('Nx')
    ax_L.set_title('Lx')

    plt.show(block=False)


def kmer_QV():
    with open('BSWCHEF120152514636.asm-ccs.qv.txt','r') as file_in:
        for line in file_in:
            if line[:2] == 'CV':
                coverage = float(line.rstrip().split()[1])
            elif line[:2] == 'QV':
                raw_QV, adjusted_QV = (float(i) for i in line.rstrip().split()[1:])

    
def busco_report():
    with open('short_summary.specific.cetartiodactyla_odb10.busco_results.txt') as file_in:
        for line in file_in:
            if 'C:' in line:
                return line.strip().rstrip()

            
#pip install md2pdf
def generate_markdown_string():
    build_str = '#ASSEMBLY\n'

    build_str += 'apple sauce\n'

    build_str += '##test 2'

    with open('assembly_report.md','w') as file_out:
        file_out.write(build_str)
    
