
#!/dss/dsshome1/08/ra65mav/.conda/envs/py3.11/bin/python

'''
input file: modified cellsnp output table "m" from R
input arg: file path
input arg: sample

Output: csv file with barcode and distances between adjacent markers for each barcode
csv format: barcode, distance
'''

import sys

import pandas as pd
from collections import defaultdict


# sys.argv[0] is the name of the python script
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

#input_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/data/A/m.csv'
#output_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/res.txt'


def table_to_interval_dict(input_file_path, output_file_path):
    # read input file to table
    t = pd.read_table(input_file_path, sep=',')
    # first order table by barcode, then by position
    t = t.sort_values(by=['barcode_index', 'position'])
    # reset indices so enumerate will start from 0, not original indices
    t = t.reset_index(drop=True)
    # open output file
    with open(output_file_path, 'w+') as fw:
        # write header
        fw.write('barcode_index,marker_distance\n')
        for i,barcode_index in enumerate(t['barcode_index']):
            try:
                if t['barcode_index'][i+1]==barcode_index:
                    # calc distance for each pair of positions within the same barcode
                    distance = t['position'][i+1]-t['position'][i]
                    fw.write(f'{barcode_index},{distance}\n')
            except KeyError:
                # when t['barcode_index'][i+1] iterates out of range
                print('done')

    return 0


table_to_interval_dict(input_file_path, output_file_path)
