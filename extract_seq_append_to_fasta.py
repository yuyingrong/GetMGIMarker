# 20231212
# extract seqs from infile.fasta by header name, then append to outfile.fasta

# input: original .fasta
# output: extended .fasta

import sys

infile = sys.argv[1]
outfile = sys.argv[2]
seqname = sys.argv[3]
# seqname: list of seq names for extracting and appending

##### WILL MODIFY OUTFILE #####
##### APPENDING DIRECTLY TO OUTFILE #####
##### MUST FIRST MAKE A COPY OF OUTFILE #####

def extract_and_add_seq(infile, outfile, seqname):
    # infile is the file from which seqs are extracted
    # outfile is the file to which extracted seqs are appended
    x = 0
    counter = 0
    # read seqname into python list
    with open(seqname, 'r') as fr:
        lines = fr.readlines()
        if not lines[0].startswith('>'):
            headers=[''.join(['>', name]) for name in lines]
        elif lines[0].startswith('>'):
            headers=[name for name in lines]
    # open outfile for appending
    fw = open(outfile, 'a')
    
    # read infile 
    with open(infile, 'r') as fr:
        print('Appending:')
        for row in fr:
            if x == 1:
                # write seq to outfile
                fw.write(row)
                # restore x
                x = 0
                # skip lower loop
                continue
            if row in headers:
                # append seq header to outfile
                fw.write(row)
                print(row.strip('>').strip('\n'))
                counter += 1
                # enable upper loop
                x = 1
    print(f'Appended {counter} sequence(s)!')
    return 0


extract_and_add_seq(infile, outfile, seqname)

print('done')

