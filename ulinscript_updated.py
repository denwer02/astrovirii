import argparse
import copy
import os
import pandas as pd
import re
from Bio import SeqIO
from parser_gb import read_csv

def orf_coord(input_file, orf_map):
    '''
    Retrieves the coordinates of ORFs
    
    Input:
        input_file - file with nucleotide sequences in genbank-format
        orf_map - csv file annotation of orfs and their codes
    Output:
        coord_file - file with coordinates
    '''

    # dictionary with annotations of ORFs
    orf_dict = read_csv(orf_map)
    # possible ORFs
    orf_types = ['1a', '1b', '1ab', '1ab_orf', '2']
    # list of ORFs in final table
    orf_types_final = ['1a', '1b', '2']


    out_file_name = os.path.splitext(input_file)[0] + '_orf.txt'
    print(out_file_name)
    out_file = open(out_file_name, 'w')

    out_file.write('id' + ',' + ','.join(orf_types_final) + '\n')
    

    dict_coord = {}

    with open(input_file) as handle:
        records = list(SeqIO.parse(handle, 'gb'))
        for rec in records:
            dict_coord[rec.id] = {}
            # counter for polymerase A and B genes
            # if counter==0, haven't met A gene
            pol_count = 0
            print(rec.id)
            for feature in rec.features:
                if 'codon_start' in feature.qualifiers.keys():
                    cod_start = int(feature.qualifiers['codon_start'][0]) - 1
                    if cod_start<0:
                        print('here')
                else:
                    cod_start = 0
                if feature.type == 'CDS':
                    if 'product' in feature.qualifiers.keys():
                        product = map_feature(feature.qualifiers['product'][0], orf_dict)
                        #print(feature.qualifiers['product'][0], product)
                        
                        if product not in dict_coord[rec.id].keys():
                            
                            if product in orf_types_final:
                                print(product, orf_types_final)
                                dict_coord[rec.id][product] = [int(feature.location._start) + cod_start, int(feature.location._end)]
                            elif product == '1ab':
                                if pol_count == 1:
                                    continue
                                else:
                                    print(feature.location.parts)
                                    print([int(feature.location.parts[0]._start) + cod_start, int(feature.location.parts[0]._end)])
                                    dict_coord[rec.id]['1a'] = [int(feature.location.parts[0]._start) + cod_start, int(feature.location.parts[0]._end)]
                                    dict_coord[rec.id]['1b'] = [int(feature.location.parts[1]._start) + cod_start, int(feature.location.parts[1]._end)]
                                    pol_count += 1
                            elif product == '1ab_orf':
                                if pol_count == 0:
                                    dict_coord[rec.id]['1a'] = [int(feature.location._start) + cod_start, int(feature.location._end)]
                                    pol_count += 1
                                elif pol_count == 1:
                                    dict_coord[rec.id]['1b'] = [int(feature.location._start) + cod_start, int(feature.location._end)]
                                else:
                                    print('Couldn\'t find annotation \'product\' qualifier for {}'.format(rec.id))
                        
                    elif 'note' in feature.qualifiers.keys():
                        note = map_feature(feature.qualifiers['note'][0], orf_dict)
                        if note not in dict_coord[rec.id].keys():
                            if note in orf_types_final:
                                dict_coord[rec.id][note] = [int(feature.location._start) + cod_start, int(feature.location._end)]
                            elif note == '1ab':
                                if pol_count == 1:
                                    continue
                                else:
                                    dict_coord[rec.id]['1a'] = [int(feature.location.parts[0]._start) + cod_start, int(feature.location.parts[0]._end)]
                                    dict_coord[rec.id]['1b'] = [int(feature.location.parts[1]._start) + cod_start, int(feature.location.parts[1]._end)]
                                    pol_count += 1
                            elif note == '1ab_orf':
                                if pol_count == 0:
                                    dict_coord[rec.id]['1a'] = [int(feature.location._start) + cod_start, int(feature.location._end)]
                                    pol_count += 1
                                elif pol_count == 1:
                                    dict_coord[rec.id]['1b'] = [int(feature.location._start) + cod_start, int(feature.location._end)]
                                else:
                                    print('Couldn\'t find annotation in \'gene\' qualifier for {}'.format(rec.id))
    for id in dict_coord.keys():
        # string to write to out_file
        s = id 
        for orf in orf_types_final:
            if orf in dict_coord[id].keys():
                st = dict_coord[id][orf][0]
                e = dict_coord[id][orf][1]
                s = s + ',' + str(st) + '-' + str(e)
            else:
                s = s + ',NA-NA'
        s += '\n'
        out_file.write(s)
    #print(dict_coord)
    out_file.close()


def map_feature(feature, feature_map):
    '''
    feature_map - dictionary, e.g. feature_map['Italia']='ITA'
    feature - possible key from feature_map
    return value if key is in feature_map
    '''
    for k, v in feature_map.items():
        if feature == k:
            return v
    return feature

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-orf_map", "--orf_map_file", type=str,
                        help="Csv-file with short codes for ORFs", required=True)
    args = parser.parse_args()

    orf_coord(args.input_file, args.orf_map_file)
