
#!/dss/dsshome1/08/ra65mav/.conda/envs/py3.11/bin/python

'''
input file: modified cellsnp output table "m" from R
input arg: file path
input arg: sample

Output: csv file with barcode and distances between adjacent markers for each barcode
csv format: barcode, distance
'''

import sys

import csv


# sys.argv[0] is the name of the python script
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

#input_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/data/A/m.csv'
#output_file_path = '/dss/dsslegfs01/pn29fi/pn29fi-dss-0008/Yuying/out/res.txt'


def calc_distance(input_file_path, output_file_path):
    # open output file
    with open(output_file_path, 'w+') as fw:
        # write header
        fw.write('barcode_index,marker_distance\n')
        
        # initiate barcode_index
        prev_barcode_index = ''
        prev_position = ''
        counter = 0
        # input file contains only valid barcodes, and is already ordered by barcode then by position
        with open(input_file_path) as f:
            # skip header
            next(f)
            # read one row at a time; instead of reading in the whole file
            rows = csv.reader(f, delimiter=',', quotechar='"')
            
            for row in rows:
                barcode_index = row[0]
                position = int(row[1])
                if barcode_index == prev_barcode_index:
                    distance = position - prev_position
                    fw.write(f'{barcode_index},{distance}\n')
                    prev_position = position
                elif barcode_index != prev_barcode_index:
                    # if the first row is being iterated
                    # or if all rows for prev_barcode_index has been iterated
                    prev_barcode_index = barcode_index
                    prev_position = position
                    counter += 1
                    continue
            print(f'Iterated {counter} barcodes.')
        print('done')

    return 0


calc_distance(input_file_path, output_file_path)
