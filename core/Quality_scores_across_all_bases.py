import seaborn as sns
import numpy as np
import pandas as pd
import os
import pyfastx
import matplotlib.pyplot as plt

def Quality_scores_across_all_bases(path):

    fastq_file = pyfastx.Fastq(path, build_index=False)

    quality = []
    n = 0

    # loading quality

    for name, sequence, qual in fastq_file:
        quality.append(qual)
        n += 1
        if n == 100000:
            break

    # decoding quality from ASCII

    translate_data = pd.DataFrame([(ord(j) - 33) for j in i] for i in quality)

    l = len(translate_data.T) + 1

    translate_data.columns = [i for i in range(1, l)]

    # data compression for the x-axis

    q3, q1 = np.percentile(translate_data[1].tolist(), 75), np.percentile(translate_data[1].tolist(), 25)
    x_ax = []
    s = 1

    for i in range(2, l):
        q3_, q1_ = np.percentile(translate_data[i].tolist(), 75), np.percentile(translate_data[i].tolist(), 25)

        if q3_ == q3 and q1_ == q1 and s < 15:
            s += 1
        else:
            q3 = q3_
            q1 = q1_
            x_ax.append(s)
            s = 1

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
                    translate_data.drop(columns=[s - j], inplace=True)
                line = str(s - i + 1) + '-' + str(s)
                X.append(line)

    else:

        for i in x_ax:
            s += i

            if i == 1:
                X.append(str(s))

            elif i < k:
                for j in range(1, i):
                    translate_data.drop(index=[s - j], inplace=True)

                X.append(str(s - i + 1))

            else:
                for j in range(1, i):
                    translate_data.drop(index=[s - j], inplace=True)
                line = str(s - i + 1) + '-' + str(s)
                X.append(line)

    if X[-1].count(str(s)) < 1:
        X[-1] = X[-1] + '-' + str(s)

    translate_data.columns = X

    mean = translate_data.mean()
    mean.index = [i for i in range(0, len(mean))]

    # graph styles

    sns.set_theme(style='ticks', rc={'axes.facecolor': (0.9994925028835063, 0.9192618223760093, 0.6061361014994233),
                                     'figure.facecolor': 'white', 'axes.grid': False, 'grid.color': 'white',
                                     'grid.linestyle': '-'})

    # creating the axes of the graph and its construction

    ax = sns.catplot(data=translate_data, kind='box', showfliers=False, color="blue")
    ax.fig.set_size_inches(14, 7)
    plt.grid(axis='both')

    sns.lineplot(data=mean, color='black')

    plt.yticks(np.arange(0, max(translate_data.max()) + 2, 2))
    plt.ylim((0, max(translate_data.max())))
    plt.xticks(range(0, len(X), 1 if sum(x_ax) // 40 <= 1 else sum(x_ax) // 40))
    plt.xlim(-1, len(X))

    plt.xlabel('Position in read(bp)')
    plt.ylabel('Quality')
    plt.title('Quality scores across all bases (' + os.path.basename(path) + ')')

    axis = [i for i in range(-1, len(X) + 1, 1)]
    plt.plot(axis, [20 for i in axis], color="r", linewidth=2)
    plt.plot(axis, [28 for i in axis], color=(0, 0.6, 0.5), linewidth=2)

    plt.savefig(r"core/PDF/Quality scores across all bases.pdf")

    return