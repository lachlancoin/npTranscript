#!/usr/bin/env python3
'''
npTranscript utility for automated retrievel of Coronavirus ORF co-ordinates from NCBI

author: Daniel Rawlinson, The Peter Doherty Institute for Infection and Immunity
email: daniel.rawlinson@unimelb.edu.au


'''
'''
TODO:
    -enable saving of genome sequence
    '''


#imports
import requests
import argparse
import xml.etree.ElementTree as ET
import csv
import sys

#constants
target_feats = ['5\'UTR', 'gene', '3\'UTR']
write_dicts = [] #list to contain dictionaries of features for writing
feature_keys = ['Name','Type','Minimum','Maximum','Length', 'Direction','gene']

#args
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--genome', required=True, metavar='STRING', type=str, help='NCBI accession number for virus genome')
#parser.add_argument('--keep-seq', action= 'store_true', help='Flag to store sequence retrieved in fasta format')

opts = parser.parse_args()
#run
URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'

gff_params = {'db':'nuccore','id':[opts.genome],'rettype':'gb','idtype':'acc','retmode':'xml'}
#fa_params = {'db':'nucleotide','id':[opts.genome],'rettype':'fasta','idtype':'acc'}

print('Downloading Data for genome accesion#{}'.format(opts.genome))
gff = requests.get(URL, params = gff_params)
#fa = requests.get(URL, params= fa_params)

#evaluate download
if gff.status_code != 200:
    sys.exit('Error: failed to download. Received HTTP error {}'.format(gff.status_code))
else:
    print('Data successfully downloaded')

print('Parsing Data')
record = ET.ElementTree(ET.fromstring(gff.text))
root = record.getroot()
for r in root.iter('GBSeq_feature-table'):
    for f in r.findall('GBFeature'):
        if f[0].text in target_feats:
            typ = f[0].text
            start = int(f[2][0][0].text)
            stop = int(f[2][0][1].text)
            if typ == target_feats[1]:
                gene = f[3][0][1].text
                name = gene+' '+typ
            else: 
                gene = typ.replace('\'', '') 
                name = gene
            
            feature_values = [name, typ, start, stop, stop - start , 'forward',gene]
            
            write_dicts.append(dict(zip(feature_keys, feature_values)))

with open('Coordinates.csv', 'w',  newline='') as f:
    dwrite = csv.DictWriter(f, fieldnames = feature_keys)
    dwrite.writeheader()
    for feature in write_dicts:
        dwrite.writerow(feature)
    