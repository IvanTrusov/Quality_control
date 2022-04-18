import seaborn as sns
import numpy as np
import pyfastx
import pandas as pd
import matplotlib.pyplot as plt
import os
from collections import Counter
from scipy.stats import shapiro, norm
from statsmodels.graphics.gofplots import qqplot
import random

def GC_distribution_over_all_sequences(path):

    fastq_file = pyfastx.Fastq(path, build_index=False)

    seq = []
    n = 0

    # loading sequences

    for name, sequence, qual in fastq_file:
        seq.append(sequence)
        n += 1
        if n == 100000:
            break

    divide_n = n / 100

    # graph styles

    sns.set_theme(style='ticks', rc={'axes.facecolor': (0.9994925028835063, 0.9192618223760093, 0.6061361014994233),
                                     'figure.facecolor': 'white', 'axes.grid': False, 'grid.color': 'white',
                                     'grid.linestyle': '-'})

    # calculating the average GC content in reads

    gc_list = [round(((i.count('G') + i.count('C')) / len(i)) * 100, 2) for i in seq]

    gc = pd.DataFrame(gc_list)

    gc_gr = Counter(gc_list)

    # GC distribution plot

    fig = plt.figure(figsize=(14, 7))

    grid = plt.GridSpec(4, 3, wspace=0.25)
    ax_main = fig.add_subplot(grid[:, :-1])
    ax_right = fig.add_subplot(grid[:, -1])

    x = sorted(gc_gr.keys())
    y = [round(gc_gr[i] / divide_n, 2) for i in x]

    ax_main.bar(x, y, color="blue", width=1)
    ax_main.grid(axis='both')

    ax_main.set(ylabel=('Percentage'))
    ax_main.set(xlabel=('Mean of GC content(%)'))

    # Normal Distribution Plot

    x = np.arange(1, 100, 0.1)
    f = [i * 100 for i in norm.pdf(x, gc.mean(), gc.std())]

    ax_main.plot(x, f, label='Normal Distribution', color='red')
    ax_main.legend(facecolor='white')

    ax_main.set(xlim=(0, 100))
    ax_main.set(xticks=(np.arange(0, 101, 4)))
    ax_main.set(title=('GC distribution over all sequences (' + os.path.basename(path) + ')'))

    # Probability Plot
    gc_short_list = gc_list[:5000]

    shapiro_test = shapiro(random.sample(gc_short_list, len(gc_short_list)))

    res = qqplot(np.array(random.sample(gc_short_list[:1000], len(gc_short_list[:1000]))), line='s', ax=ax_right)

    ax_right.get_lines()[0].set_markerfacecolor('blue')
    ax_right.get_lines()[0].set_markersize(5.0)
    ax_right.get_lines()[1].set_linewidth(2.0)
    ax_right.get_lines()[1].set_color('black')

    ax_right.annotate('Shapiro test: p-value = ' + '%.3e' % shapiro_test[1], (0.28, 0.02), size='small',
                      weight='roman',
                      xycoords='axes fraction', va='center')

    ax_right.grid(axis='both')
    ax_right.set(title=('GC Q-Q Plot ' + '(' + os.path.basename(path) + ')'))

    plt.savefig(r"core/PDF/GC distribution over all sequences.pdf")

    return