<!-- ------------------ HEADER ------------------ -->
<!-- Developed and maintained by Genomics Division
<!-- of the Institute of Technology an Renewable Energy (ITER)
<!-- Tenerife, Canary Islands, SPAIN
<!-- See the "Contact us" section to collaborate with us to growth
<!-- this repository. ;=)

<!-- ------------------ SECTION ------------------ -->
<p align="left">
  <a href="https://github.com/genomicsITER/monkeypox" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/logos_GH.png" width="auto" /> 
      </a>
</p>

# Monkeypox
A public repository of resources related to MonkeyPox maintained by ITER.

This is the result of an ongoing joint-effort of the following institutions and laboratories:
<ul>
 <li><a href="https://www3.gobiernodecanarias.org/sanidad/scs/organica.jsp?idCarpeta=10b3ea46-541b-11de-9665-998e1388f7ed">Servicio de Microbiología, Hospital Universitario Ntra. Sra. de Candelaria</a>, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li><a href="https://fciisc.org/">Fundación Canaria Instituto de Investigación Sanitaria de Canarias</a> at the Research Unit, Hospital Universitario Ntra. Sra. de Candelaria</a>, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li><a href="https://portalciencia.ull.es/grupos/6361/detalle">Laboratorio de Inmunología Celular y Viral</a>, Unidad de Farmacología, Facultad de Medicina, Universidad de La Laguna, 38200 San Cristóbal de La Laguna, Spain.</li>
 <li><a href="https://www.iter.es/areas/area-genomica/">Genomics Division, Instituto Tecnológico y de Energías Renovables</a>, 38600 Santa Cruz de Tenerife, Spain

<hr>
<!-- ------------------ SECTION ------------------ -->

## Python Code to programmatically retrieve MPXV sequences from GenBank ##
<!-- ## Table of contents ## -->
<!-- <ul> -->
<!-- <li><a href="#1">1. ...</a></li> -->
<!-- </ul> -->

<!-- <hr> -->

<a name="1"></a>
####  Retrieve MPXV sequences of length > 190.000 bases from GenBank
  
```Python

#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Version: June 28, 2022
	
from pathlib import Path
import csv
import json
import re
import subprocess

# Usage: conda run --no-capture-output -n entrez-direct python3 download-monkeypox-sequences.py

MIN_SEQ_LEN = 190000
QUERY = 'esearch -db nuccore -query "monkeypox" | efilter -organism "Monkeypox virus" -molecule "genomic" | efetch -format docsum -mode json'
METADATA_FILENAME = 'metadata.tsv'
LOCAL_ACCESSIONS = 'accessions.tsv'
QUERY_LOCAL = 'efetch -input "{}" -db nuccore -format docsum -mode json'.format(LOCAL_ACCESSIONS)

to_csv = []
accessions = []

local_file = Path(LOCAL_ACCESSIONS)
if local_file.is_file():
  print("[1] Performing query from local accessions: {}".format(LOCAL_ACCESSIONS))
  ps = subprocess.Popen(QUERY_LOCAL,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  output = ps.communicate()[0]
  if ps.returncode == 0:
    print("Query was OK!")
	
  print("[1] Performing query: {}".format(QUERY_LOCAL))
  ps = subprocess.Popen(QUERY_LOCAL,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  output = ps.communicate()[0]
  if ps.returncode == 0:
    print("Query was OK!")
	
  print("[2] Processing query result...")
  data=json.loads(output)
  record_uids=data['result']['uids']
  print("Query has returned {} records".format(len(record_uids)))

  print("[3] Parsing result records...")
  for uid in record_uids:
    if data['result'][uid]['slen'] >= MIN_SEQ_LEN:
      data_dict = {
        'title': data['result'][uid]['title'],
        'filename': '{}.fasta'.format(data['result'][uid]['accessionversion']),
        'caption': data['result'][uid]['caption'],
        'accession': data['result'][uid]['accessionversion'],
        'creation_date': data['result'][uid]['createdate'],
        'update_date': data['result'][uid]['updatedate'],
        'ncbi_link': 'https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(data['result'][uid]['caption']),
        'completeness': data['result'][uid]['completeness'],
        'sequence_length_in_bp': data['result'][uid]['slen'],
        'extra_data': data['result'][uid]['extra'],
        'subtype_format': data['result'][uid]['subtype'],
        'subtype_data': data['result'][uid]['subname']
      }
      to_csv.append(data_dict)
      accessions.append(data['result'][uid]['caption'])
  print("Query has found {} valid records".format(len(to_csv)))
  
print("[1] Performing query: {}".format(QUERY))
ps = subprocess.Popen(QUERY,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
output = ps.communicate()[0]
if ps.returncode == 0:
  print("Query was OK!")

print("[2] Processing query result...")
data=json.loads(output)
record_uids=data['result']['uids']
print("Query has returned {} records".format(len(record_uids)))

print("[3] Parsing result records...")
to_csv = []
for uid in record_uids:
  if data['result'][uid]['slen'] >= MIN_SEQ_LEN and data['result'][uid]['caption'] not in accessions:
    data_dict = {
      'title': data['result'][uid]['title'],
      'filename': '{}.fasta'.format(data['result'][uid]['accessionversion']),
      'caption': data['result'][uid]['caption'],
      'accession': data['result'][uid]['accessionversion'],
      'creation_date': data['result'][uid]['createdate'],
      'update_date': data['result'][uid]['updatedate'],
      'ncbi_link': 'https://www.ncbi.nlm.nih.gov/nuccore/{}'.format(data['result'][uid]['caption']),
      'completeness': data['result'][uid]['completeness'],
      'sequence_length_in_bp': data['result'][uid]['slen'],
      'extra_data': data['result'][uid]['extra'],
      'subtype_format': data['result'][uid]['subtype'],
      'subtype_data': data['result'][uid]['subname']
    }
    to_csv.append(data_dict)
    accessions.append(data['result'][uid]['caption'])
print("Query has found {} valid records".format(len(to_csv)))

keys = to_csv[0].keys()
print("[4] Generating metadata file...")
with open(METADATA_FILENAME, 'w', encoding='utf8', newline='') as output_file:
  dict_writer = csv.DictWriter(output_file, keys, delimiter='\t')
  dict_writer.writeheader()
  dict_writer.writerows(to_csv)
print("Metadata written to {}".format(METADATA_FILENAME))

print("[5] Downloading sequences...")
for item in to_csv:
  print("Downloading {}...".format(item['filename']))
  cmd = 'efetch -db nuccore -id "{}" -format fasta > {}'.format(item['accession'], item['filename'])
  ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  output = ps.communicate()[0]
```
  
