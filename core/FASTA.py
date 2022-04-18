def FASTA(file):

    dy={}

    for line in file:
        line = line.strip()

        if line[0]=='>':
            key = line[1:]
            value=[]
        else:
            value.append(line)
            dy[key] = ''.join(value)

    return dy

