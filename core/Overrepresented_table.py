from collections import Counter
import pandas as pd
import os
import pyfastx
import matplotlib.pyplot as plt

def Overrepresented_table(path):

    fastq_file = pyfastx.Fastq(path, build_index=False)

    sequences = []
    n = 0

    # loading sequences

    for name, seq, qul in fastq_file:
        sequences.append(seq[:51])
        n += 1
        if n == 100000:
            break

    # counting overrepresented sequences

    overrep = Counter(sequences)
    all = len(sequences)
    sequences.clear()

    #creating an information dataframe

    df = pd.DataFrame({'Sequence': overrep.keys(), 'Count': overrep.values(),
                       'Percentage': [round((i / all) * 100, 12) for i in overrep.values()]})

    df = df.sort_values('Count', ascending=False)

    #threshold of displayed values

    df_over_1 = df[df['Percentage'] >= 0.1]

    if len(df_over_1) > 0:

        df_over = pd.DataFrame({'Sequence': df_over_1['Sequence'],
                                'Percentage(%)': df_over_1['Percentage']})

        # ability to save sequences for preprocessing

        # df_over.to_csv(os.path.join(os.path.dirname(path), os.path.basename(path) + '_over' + '.csv'))

        # constructing the graph

        plt.figure(figsize=(14, (len(df_over.values) / 4) + 1))
        ax = plt.subplot()
        ax.axis('off')
        table = ax.table(colLabels=df_over.columns, colColours=['white'] * 2, cellText=df_over.values,
                         cellColours=[['lightgray'] * 2] * len(df_over.values), cellLoc='center', loc='center',
                         bbox=[-0.12, 0, 1.2, 1], colWidths=[8 for i in df_over.columns])

        table.set_fontsize(11)
        plt.title('Overrepresented sequences (' + os.path.basename(path) + ')', fontsize=14)

        plt.savefig(r"core/PDF/Overrepresented table.pdf")

    else:

        plt.figure(figsize=(14, 3))
        plt.title('There are no overrepresented sequences in ' + os.path.basename(path), fontsize=14, pad=-10)
        plt.axis('off')

        plt.savefig(r"core/PDF/Overrepresented table.pdf")

    return