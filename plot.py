import matplotlib.pylab as plot
import os, os.path

# Count the number of output files produced by the hydro code.
DIR = './output'
N = len([name for name in os.listdir(DIR) if os.path.isfile(os.path.join(DIR, name))]) - 1

# Generate and save a 'frame' for each output file.
for i in range(N):
    f = open(f"./output/data{i}.txt",'r')

    x   = []
    rho = []
    mo   = []
    en = []

    for row in f:
        row = row.split(' ')
        x.append(float(row[0]))
        rho.append(float(row[1]))
        mo.append(float(row[2]))
        en.append(float(row[3]))

    fig, (ax1, ax2, ax3) = plot.subplots(1, 3)
    ax1.scatter(x, rho, s=5, c='red')
    ax2.scatter(x, mo, s=5, c='green')
    ax3.scatter(x, en, s=5, c='blue')
    ax1.title.set_text('Mass Density')
    ax2.title.set_text('Momentum Density')
    ax3.title.set_text('Energy Density')
    #ax1.set_ylim(0.06,1.04)
    #ax2.set_ylim(-0.018,0.42)
    #ax3.set_ylim(0.13,1.56)
    fig.set_size_inches(15, 5)
    fig.tight_layout(pad=2.5)
    fig.savefig(f'./frames/{i}.png')

    f.close()
