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
    ax1.scatter(x, rho, s=4, c='red')
    ax2.scatter(x, mo, s=4, c='green')
    ax3.scatter(x, en, s=4, c='blue')
    ax1.title.set_text('Mass Density')
    ax2.title.set_text('Velocity')
    ax3.title.set_text('Pressure')
    ax1.set_xlim(0,5)
    ax1.set_ylim(0,5)
    ax2.set_xlim(0,5)
    ax3.set_xlim(0,5)
    fig.set_size_inches(15, 5)
    fig.tight_layout(pad=2.5)
    fig.savefig(f'./frames/{i}.png')

    f.close()
