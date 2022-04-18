from collections import Counter
import seaborn as sns
import numpy as np
import os
import pyfastx
import matplotlib.pyplot as plt

def Per_sequence_quality_scores(path):

    fastq_file = pyfastx.Fastq(path, build_index=False)

    quality = []
    n = 0

    # loading quality

    for name, sequence, qual in fastq_file:
        quality.append(qual)
        n += 1
        if n == 100000:
            break

    divide_n = n / 100

    # decoding quality from ASCII and calculating the average content

    translate_data_mean = [int(sum([ord(j) - 33 for j in i]) / len(i)) for i in
                           quality]

    translate_data_mean = Counter(translate_data_mean)

    # graph styles

    sns.set_theme(style='ticks', rc={'axes.facecolor': (0.9994925028835063, 0.9192618223760093, 0.6061361014994233),
                                     'figure.facecolor': 'white', 'axes.grid': False, 'grid.color': 'white',
                                     'grid.linestyle': '-'})

    # creating the axes of the graph and its construction

    plt.figure(figsize=(14, 7))
    x = sorted(translate_data_mean.keys())
    plt.bar(x, [round(translate_data_mean[i] / divide_n, 2) for i in x], color="blue")

    plt.xticks(np.arange(0, max(sorted(translate_data_mean.keys())) + 2, 2))
    plt.grid(axis='both')
    plt.ylabel('Percentage')
    plt.xlabel('Quality')
    plt.title('Per sequence quality scores (' + os.path.basename(path) + ')')

    plt.savefig(r"core/PDF/Per sequence quality scores.pdf")

    return