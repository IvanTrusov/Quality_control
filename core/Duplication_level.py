from collections import Counter
import seaborn as sns
import pyfastx
import pandas as pd
import matplotlib.pyplot as plt
import os

def Duplication_level(path):

    fastq_file = pyfastx.Fastq(path, build_index=False)

    sequences = []
    n = 0

    # loading sequences

    for name, seq, qul in fastq_file:
        sequences.append(seq)
        n += 1
        if n == 100000:
            break

    divide_n = n / 100

    # graph styles

    sns.set_theme(style='ticks', rc={'axes.facecolor': (0.9994925028835063, 0.9192618223760093, 0.6061361014994233),
                                     'figure.facecolor': 'white', 'axes.grid': False, 'grid.color': 'white',
                                     'grid.linestyle': '-'})

    S = Counter(sequences)
    percentage = [0] * 17

    # creating the x-axis and counting the number of duplicates

    for i in Counter(S.values()).items():

        if i[0] < 10:
            percentage[i[0]] += (i[0] * i[1]) / divide_n

        elif i[0] > 10 and i[0] <= 50:
            percentage[10] += (i[0] * i[1]) / divide_n

        elif i[0] > 50 and i[0] <= 100:
            percentage[11] += (i[0] * i[1]) / divide_n

        elif i[0] > 100 and i[0] <= 500:
            percentage[12] += (i[0] * i[1]) / divide_n

        elif i[0] > 500 and i[0] <= 1000:
            percentage[13] += (i[0] * i[1]) / divide_n

        elif i[0] > 1000 and i[0] <= 5000:
            percentage[14] += (i[0] * i[1]) / divide_n

        elif i[0] > 5000 and i[0] <= 10000:
            percentage[15] += (i[0] * i[1]) / divide_n

        else:
            percentage[16] += (i[0] * i[1]) / divide_n

    percentage.pop(0)

    col_val = [str(i) for i in range(1, 10)] + ['>10', '>50', '>100', '>500', '>1k', '>5k', '>10k']

    # preparing two axes in a dataframe

    graphic = pd.DataFrame({'Percentage': percentage, 'Duplicate level': col_val})

    # constructing the graph

    graphic.plot(kind='bar', x='Duplicate level', y='Percentage', figsize=(14, 7), color='blue', legend=False)

    plt.grid(axis='both')

    if int(max(percentage)) // 10 > 0:
        k = int(max(percentage)) // 10
    else:
        k = 1

    plt.yticks(range(0, int(max(percentage)) + k, k))

    plt.xlabel('Sequence duplication level')
    plt.xticks(rotation=0)
    plt.ylabel('Percentage')
    plt.title('Sequence Duplication Levels (' + os.path.basename(path) + ')')

    plt.savefig(r"core/PDF/Duplication level.pdf")

    return