
#!/dss/dsshome1/08/ra65mav/.conda/envs/py3.11/bin/python

'''
Input 1: cellranger output fragments.tsv, after
- subsetting for wanted barcodes, which also removes headers starting with '#'
- removing unwanted cols
- sorting rows by barcodes

Output: txt tile containing max_depth, mean_depth, median_depth, and breadth
'''

import sys

import csv
import numpy as np

# sys.argv[0] is the name of the python script
#barcodes_file_path = sys.argv[1]
fragments_file_path = sys.argv[1]
stats_file_path = sys.argv[2]

# Example file paths
# barcodes_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/tb1000_barcodes.txt'
# fragments_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/data/A/frags1000.tsv'
# stats_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/stats_tb1000.csv'

#print(f"barcodes_file_path = '{barcodes_file_path}';")
print(f"fragments_file_path = '{fragments_file_path}';")
print(f"stats_file_path = '{stats_file_path}'.")



def calc_stats_per_barcode(fragments_file_path, stats_file_path):
    # open output file for writing
    fw = open(stats_file_path, 'w+')
    fw.write('barcode,max,mean,median,breadth\n')
    counter = 0
    # initiate barcode
    barcode = 'barcode'
    # input file contains only valid barcodes and is already ordered by barcode
    with open(fragments_file_path) as f:
        # read one row at a time; instead of reading in the whole file
        rows = csv.reader(f, delimiter="\t", quotechar="'")

        for row in rows:
            # format of a row in frags.tsv: 7\t43\tGGGTGTCAGCGCATTT-1\n
            
            # the start/end indices are not nucleotide positions, but (nucleotide position - 1)
            start = int(row[0])-4-1
            end = int(row[1])+5-1
            # convert start/end indices into readpair indicies
            read1 = np.arange(start, min(start+50, end), dtype='int32')
            read2 = np.arange(max(start, end-49), end, dtype='int32')
            
            current_barcode = row[2]
            
            if barcode == current_barcode:
                # update the arr
                arr[np.union1d(read1, read2)]+=1
            elif barcode != current_barcode:
                # this means that either all frags for this barcode have been iterated, since rows are ordered by barcode
                # or the first row is being iterated
                try: 
                    # write stats for this barcode to output file
                    fw.write(f'{barcode},{np.max(arr)},{round(np.average(arr[arr!=0]),2)},{np.median(arr[arr!=0])},{arr[arr!=0].shape[0]}\n')
                    counter += 1
                except NameError:
                    # this should only happen when the first row is iterated, so arr has not been declared
                    print(barcode)### should print 'barcode' ###
                    pass
                # update barcode
                barcode = current_barcode
                # (re-)initiate arr of all genomic positions
                arr = np.zeros(shape=(85145900,), dtype='int32')
                # update the arr
                arr[np.union1d(read1, read2)]+=1
        # write the last line to file
        fw.write(f'{barcode},{np.max(arr)},{round(np.average(arr[arr!=0]),2)},{np.median(arr[arr!=0])},{arr[arr!=0].shape[0]}\n')
        counter += 1
        fw.close()
        print(f'Wrote {counter} lines to file.')
    
    return 0



calc_stats_per_barcode(fragments_file_path, stats_file_path)
