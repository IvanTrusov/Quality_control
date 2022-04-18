import os


def start_work():

    if 'FASTQ_files' or 'output' not in os.listdir():

        if 'FASTQ_files' not in os.listdir():
            os.mkdir('FASTQ_files')

        if 'output' not in os.listdir():
            os.mkdir('output')

        print('Missing folders have been added to your working directory, please check it.')

    print('\n'+ 'Hi, before you start working, make sure that you uploaded your files to the FASTQ_folder.' + '\n'
          + 'To continue, press Enter.')
    input()

    fastq_folder = os.listdir(r'FASTQ_files')

    if len(fastq_folder) < 1:

        while len(fastq_folder) < 1:

            if len(fastq_folder) < 1:
                print('Sorry, but the file folder is empty, please copy your files to it')

                print('Print e/t to exit or try again')

                inp = str(input())

                if inp == 'e':
                    raise Exception('The file folder is empty')

                elif inp == 't':
                    pass

                else:
                    raise Exception('Incorrect input')

            fastq_folder = os.listdir(r'FASTQ_files')

    available_files = []

    for file in fastq_folder:
        if os.path.basename(file).endswith(".gz") or os.path.basename(file).endswith(".fastq"):
            available_files.append(file)

    if len(available_files) < 1:

        while len(available_files) < 1:
            print('Sorry, but the folder must, contain fastq or fastq.gz files')

            print('Print e/t to exit or try again')

            inp = str(input())

            if inp == 'e':
                raise Exception('The file folder does not contain any available files')

            elif inp == 't':
                pass

            else:
                raise Exception('Incorrect input')

            available_files = []
            fastq_folder = os.listdir(r'FASTQ_files')

            for file in fastq_folder:
                if os.path.basename(file).endswith(".gz") or os.path.basename(file).endswith(".fastq"):
                    available_files.append(file)


    print('We found the following available files:')

    for file in available_files:
        print(file)

    print('\n' + 'Print the name of the file you want to process and press Enter')

    file = str(input())

    if file not in available_files:

        while file not in available_files:

            if file not in available_files:
                print('The file must be one of the available!')

                print('Print e/t to exit or try again')

                inp = str(input())

                if inp == 'e':
                    raise Exception('The file is not available')

                elif inp == 't':
                    pass

                else:
                    raise Exception('Incorrect input')

            print('\n' + 'Print the name of the file you want to process and press Enter')
            file = str(input())

    path = r'FASTQ_files/' + file

    print('\n' + 'Wait a few seconds, processing has already started...')

    return path
