import seaborn as sns
import numpy as np
import pyfastx
import pandas as pd
import matplotlib.pyplot as plt
import os

def Base_content(path):

    fastq_file = pyfastx.Fastq(path, build_index=False)

    seq = []
    n = 0

    # loading sequences

    for name, sequence, qual in fastq_file:
        seq.append(list(sequence))
        n += 1
        if n == 100000:
            break

    df = pd.DataFrame(seq)
    seq.clear()

    df = df.T
    df = df.apply("".join, axis=1)

    A, C, G, T, N = [], [], [], [], []

    # calculation of the average content of nucleotides in a sequence

    for s in df:
        A.append((s.count('A') / len(s)) * 100)
        C.append((s.count('C') / len(s)) * 100)
        G.append((s.count('G') / len(s)) * 100)
        T.append((s.count('T') / len(s)) * 100)
        N.append((s.count('N') / len(s)) * 100)

    content = pd.DataFrame({'A': A, 'C': C, 'G': G, 'T': T, 'N': N})

    content.index += 1

    # data compression for the x-axis

    x_ax, a, c, g, t = [], int(content.iloc[0, 0]), int(content.iloc[0, 1]), int(content.iloc[0, 2]), int(
        content.iloc[0, 3])
    s = 1

    for i in range(1, len(content['A'])):

        if int(content.iloc[i, 0]) in {a, a + 1, a - 1} \
                and int(content.iloc[i, 1]) in {c, c + 1, c - 1} \
                and int(content.iloc[i, 2]) in {g, g + 1, g - 1} \
                and int(content.iloc[i, 3]) in {t, t + 1, t - 1}:

            s += 1

        else:
            x_ax.append(s)
            s = 1

            a, c, g, t = int(content.iloc[i, 0]), int(content.iloc[i, 1]), int(content.iloc[i, 2]), int(
                content.iloc[i, 3])

    x_ax.append(s)

    s, k = 0, 5
    X = []

    if len(x_ax) <= 12:
        k = 8

    if len(x_ax) < 50:

        for i in x_ax:
            s += i

            if i < k:
                for j in range(1, i + 1):
                    X.append(str(s - i + j))

            else:
                for j in range(1, i):
                    content.drop(index=[s - j], inplace=True)
                line = str(s - i + 1) + '-' + str(s)
                X.append(line)

    else:

        for i in x_ax:
            s += i

            if i == 1:
                X.append(str(s))

            elif i < k:
                for j in range(1, i):
                    content.drop(index=[s - j], inplace=True)

                X.append(str(s - i + 1))

            else:
                for j in range(1, i):
                    content.drop(index=[s - j], inplace=True)
                line = str(s - i + 1) + '-' + str(s)
                X.append(line)

    if X[-1].count(str(s)) < 1:
        X[-1] = X[-1] + '-' + str(s)

    # graph styles

    sns.set_theme(style='ticks', rc={'axes.facecolor': 'black', 'figure.facecolor': 'white', 'grid.linestyle': '-'})

    content['X'] = X

    # constructing the graph

    content[['A', 'C', 'G', 'T', 'X', 'N']].plot(kind='bar', x='X', stacked=True, figsize=(14, 7), rot=0, lw=0,
                                                 color=[(0.8, 0.4, 0), (0.9, 0.6, 0), (0.95, 0.9, 0.25),
                                                        (0.35, 0.7, 0.9), (0, 0.6, 0.5)]) \
        .legend(bbox_to_anchor=(1, 1), facecolor='white')

    plt.grid(axis='both')
    plt.xticks(range(0, len(X), 1 if sum(x_ax) // 40 <= 1 else sum(x_ax) // 40))
    plt.yticks(np.arange(0, 110, 10))
    plt.xlabel('position in read(bp)')
    plt.ylabel('Percentage')
    plt.title('Base content (' + os.path.basename(path) + ')')

    plt.savefig(r"core/PDF/Base content.pdf")

    return