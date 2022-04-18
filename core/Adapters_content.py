import matplotlib.pyplot as plt
from collections import Counter
import re
import os
import pyfastx
from core.FASTA import FASTA
import seaborn as sns

def Adapters_content(path):

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

    # loading adapters for search

    with open(r'core/adapters_list/adapters.txt', 'r') as motifs:
        adapters = FASTA(motifs)

    adapters_content, adapters_position = {}, []

    # search for adapters

    for position, motif in enumerate(adapters.values()):
        for s in seq:
            if s.find(motif) >= 0:
                [adapters_position.append(i.start() + 1) for i in re.finditer("(?=" + motif + ")", s)]

        adapters_content[list(adapters.keys())[position]] = Counter(adapters_position)
        adapters_position.clear()

    # graph styles

    sns.set_theme(style='ticks', rc={'axes.facecolor': (0.9994925028835063, 0.9192618223760093, 0.6061361014994233),
                                     'figure.facecolor': 'white', 'axes.grid': False, 'grid.color': 'white',
                                     'grid.linestyle': '-'})

    colors = [(0.8, 0.4, 0), (0.9, 0.6, 0), (0, 0.49, 0.7), (0.35, 0.7, 0.9), (0, 0.6, 0.5)]
    plt.figure(figsize=(14, 7))
    xt = []
    count = 0

    # creating the axes of the graph and its construction

    for color, adapter_content in enumerate(adapters_content.values()):
        if len(adapter_content) > 0:

            x = sorted(adapter_content.keys())
            y = [round(adapter_content[i] / divide_n, 2) for i in x]

            xt.append(max(x))

            for position in range(len(y)):

                if position == 0:
                    pass
                else:
                    y[position] += y[position - 1]

            if max(y) > 1:
                count += 1

                if len(y) == 1:
                    plt.bar(x, y, color=colors[color], label=list(adapters_content.keys())[color], width=0.1)

                else:
                    plt.plot(x, y, color=colors[color], label=list(adapters_content.keys())[color])

                plt.legend(facecolor='white')

    # If adapters are present then two graphs are plotted

    if count > 0:
        if max(xt) // 20 <= 0:
            k = 1
        else:
            k = max(xt) // 20

        plt.xticks(range(1, max(xt) + k, k))
        plt.grid(axis='both')
        plt.title('Adapter content with cumulative filling(' + os.path.basename(path) + ')')
        plt.xlabel('Position in read')
        plt.ylabel('Percentage')
        plt.yticks(range(0, 110, 10))

        plt.savefig(r"core/PDF/Adapters content with cumulative filling.pdf")

        plt.figure(figsize=(14, 7))

        xt, yt = [], []

        for color, adapter_content in enumerate(adapters_content.values()):
            if len(adapter_content) > 0:

                x = sorted(adapter_content.keys())
                y = [round(adapter_content[i] / divide_n, 2) for i in x]

                xt.append(max(x))
                yt.append(max(y))
                yt.append(min(y))

                if max(y) > 0.09:

                    # This is the part for writing information about the content of adapters to a separate file
                    # for the preprocessing

                    # with open(os.path.join(os.path.dirname(path), 'adapters_finded.txt'), 'w') as adapt:
                    #     adapt.write(list(adapters_content.keys())[color] + '\n' +
                    #                 adapters[list(adapters_content.keys())[color]])

                    if len(y) == 1:
                        plt.bar(x, y, color=colors[color], label=list(adapters_content.keys())[color], width=0.1)
                        plt.yticks(range(1, int(y[0] + 2)))

                    else:
                        plt.plot(x, y, color=colors[color], label=list(adapters_content.keys())[color])
                        plt.ylim(min(yt), max(yt) + max(yt) / 10)

                    plt.legend(facecolor='white')

        plt.xticks(range(1, max(xt) + k, k))

        plt.title('Adapter content(' + os.path.basename(path) + ')')
        plt.xlabel('Position in read')
        plt.ylabel('Percentage')
        plt.grid(axis='both')

        plt.savefig(r"core/PDF/Adapters content.pdf")

        # if there are no adapters

    else:
        plt.figure(figsize=(14, 3))
        plt.title('There are no Adapters in ' + os.path.basename(path), fontsize=14, pad=-10)
        plt.axis('off')

        plt.savefig(r"core/PDF/Adapters content.pdf")

    return
