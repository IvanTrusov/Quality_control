import multiprocessing, os, PyPDF2, time
from core.Adapters_content import Adapters_content
from core.Base_content import Base_content
from core.Duplication_level import Duplication_level
from core.GC_distribution_over_all_sequences import GC_distribution_over_all_sequences
from core.Overrepresented_table import Overrepresented_table
from core.Per_sequence_quality_scores import Per_sequence_quality_scores
from core.Quality_scores_across_all_bases import Quality_scores_across_all_bases
from core.welcome_file import start_work

if __name__ == '__main__':

    # specifying the path to the file

    path = start_work()

    start_time = time.time()

    files = os.listdir(r'core/PDF')

    if len(files) > 1:
        for file in files:
            os.remove(r'core/PDF/' + file)

    graphics = [Adapters_content,
                Base_content,
                Duplication_level,
                GC_distribution_over_all_sequences,
                Overrepresented_table,
                Per_sequence_quality_scores,
                Quality_scores_across_all_bases,
                ]

    with multiprocessing.Pool(processes=len(graphics)) as pool:

        multiple_results = [pool.apply_async(graphic, (path,)) for graphic in graphics]
        [res.get() for res in multiple_results]

    # combining separate pdf files into one

    merger = PyPDF2.PdfFileMerger()

    for filename in files:
        merger.append(fileobj=open(os.path.join(r'core/PDF', filename), 'rb'))

    merger.write(open(os.path.join('output', os.path.basename(path) + '.pdf'), 'wb'))

    print("It takes: %s seconds" % round((time.time() - start_time), 3))

    print('\n' + 'The processing was successful. The finished file is waiting for you in the "output" folder!')

