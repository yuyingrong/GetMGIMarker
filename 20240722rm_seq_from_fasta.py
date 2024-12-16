# 20231212
# remove seqs from .fasta by header name

# input: original .fasta
# output: trimmed .fasta

import sys

infile = sys.argv[1]
outfile = sys.argv[2]
seqname = sys.argv[3]
# seqname: list of seq names for removal

def rm_seq(infile, outfile, seqname):
    x = 0
    counter = 0
    # read seqname into python list
    with open(seqname, 'r') as fr:
        headers=[''.join(['>', name]) for name in fr.readlines()]
    # open outfile for writing
    fw = open(outfile, 'w+')
    # read infile
    with open(infile, 'r') as fr:
        print('Skipping:')
        for row in fr:
            if x == 1:
                x = 0
                # skip seq line after skipping header line
                continue
            if row in headers:
                print(row.strip('>').strip('\n'))
                counter += 1
                x = 1
                # skip header line
                continue
            #print(row)
            fw.write(row)
    print(f'Removed {counter} sequence(s)!')
    return 0


rm_seq(infile, outfile, seqname)

print('done')

