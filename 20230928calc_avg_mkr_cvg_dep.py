#!/dss/dsshome1/08/ra65mav/.conda/envs/py3.11/bin/python
'''
This script generates avg depth with and without removal of barcode multiplets
# input: unmodified cellsnp output
# part 1: create defaultdict with d[barcode_index]((position, count)); output the dict
# part 2: reduce dict with duplicate vals (barcode multiplets)
# part 3: write avg depth per barcode_index, write to file
# to contain only the keys with unique values; output the cleaned dict
'''

import sys
from collections import defaultdict


# sys.argv[0] is the name of the python script
cellsnp_file_path = sys.argv[1]
output1_file_path = sys.argv[2]
#output2_file_path = sys.argv[3]

# Example file paths
# cellsnp_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/cellsnp_lite_results/B_nigra_paternal_5733_B/cellSNP.tag.DP.mtx'
# output1_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/B_avgs.csv'
# output2_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/B_avgs_dup_rm.csv'

print(f"cellsnp_file_path = '{cellsnp_file_path}';")
print(f"output1_file_path (all cells) = '{output1_file_path}';")
#print(f"output2_file_path (after barcode multiplet removal) = '{output2_file_path}'.")



def create_dict(in_file):
    # create dict of sets for each barcode_id
    d = defaultdict(set)
    
    with open(in_file, 'r') as f:
        # skip first three rows
        next(f)
        next(f)
        next(f)
        # format of in_file: position_id \t barcode_id \t count \n
        
        for line in f:
            line = line.strip('\n').split('\t')
            # add (position_id,count) to each barcode_id
            # using a set means no (position_id,count) duplicates will be added to the same barcode
            d[line[1]].add((line[0], int(line[2])))
    
    f.close()
    print(f'Number of cells before dup and zero removal: {len(d)}')
    
    return d


# remove duplicates based on 1) marker positions and 2) per marker position mapping counts
# if barcode_indices have identical 1) and 2) then they are considered barcode multiplets and removed

def remove_duplicates_and_zeros(d):
    # empty duplicates and retain unique
    count_dup = 0#
    for i,value in enumerate(d.values()):
        # search for (position_id,count) duplicates between barcodes
        # if counted more than one pair between barcodes
        if list(d.values()).count(value) > 1:
            # update counter
            count_dup += 1
            # remove the content of the current barcode and preserve the content of the other barcode(s)
            # this will remove the content of all barcodes with duplicate contents, until one remain
            d[list(d.keys())[i]] = ()
    
    # how many cells had multiple barcodes?
    print(f'Number of duplicate keys: {count_dup}')#
    
    # create tuple to contain barcodes with no content
    # either due to no content or removal in prev step
    zero_keys = ()
    for i,value in enumerate(d.values()):
        # if a key's value is empty, add key to the tuple
        # the only way to grab key from evaluating value is by enumerating the keys/values
        if not value:
            zero_keys += (list(d.keys())[i],)
            
    print(f'Number of zero keys after emptying duplicates: {len(zero_keys)}')#
    
    # rm zero_keys from output dict
    for key in zero_keys:
        d.pop(key)
    
    print(f'Number of cells after dup and zero removal: {len(d)}')#
    
    return d


# write output to file
d = create_dict(cellsnp_file_path)
avgs = [(int(key), sum([tup[1] for tup in d[key]])/len(d[key])) for key in d]
with open(output1_file_path, 'w+') as fw:
    fw.write('barcode_index,avg_depth\n')
    for pair in avgs:
        fw.write(f'{pair[0]},{round(pair[1],2)}\n')
print('done 1/2')
'''
d = remove_duplicates_and_zeros(d)
avgs = [(int(key), sum([tup[1] for tup in d[key]])/len(d[key])) for key in d]
with open(output2_file_path, 'w+') as fw:
    fw.write('barcode_index,avg_depth\n')
    for pair in avgs:
        fw.write(f'{pair[0]},{round(pair[1],2)}\n')
print('done 2/2')
'''

