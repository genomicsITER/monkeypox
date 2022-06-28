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

import csv
import json
import re
import subprocess

MIN_SEQ_LEN = 190000
CONDA_ENV_NAME = 'entrez-direct'
QUERY = 'esearch -db nuccore -query "monkeypox" | efilter -organism "Monkeypox virus" -molecule "genomic" | efetch -format docsum -mode json'

print("[1] Loading conda env...")

print("[2] Performing query: {}".format(QUERY))
ps = subprocess.Popen(QUERY,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
output = ps.communicate()[0]

print("Query result: {}".format(output))

print("[3] Processing query result...")
data=json.loads(output)
record_uids=data['result']['uids']
print("Query has returned: {} records".format(len(record_uids)))

print("[4] Parsing result records...")
to_csv = []
for uid in record_uids:
  if data['result'][uid]['slen'] >= MIN_SEQ_LEN:
    data_dict = {
      'title': data['result'][uid]['title'],
      'filename': re.sub(r'(?u)[\\\\/:*\?"<>|]', '', data['result'][uid]['title']).replace(', complete genome', '').replace(', partial genome', '').replace(" ", "_") + '.fasta',
      'caption': data['result'][uid]['caption'],
      'accession': data['result'][uid]['caption'],
      'creation_date': data['result'][uid]['createdate'],
      'update_date': data['result'][uid]['updatedate'],
      'completeness': data['result'][uid]['completeness'],
      'sequence_length_in_bp': data['result'][uid]['slen'],
      'extra_data': data['result'][uid]['extra'],
      'subtype_format': data['result'][uid]['subtype'],
      'subtype_data': data['result'][uid]['subname']
    }
    to_csv.append(data_dict)

keys = to_csv[0].keys()
print("[5] Generating metadata file...")
with open('metadata.tsv', 'w', encoding='utf8', newline='') as output_file:
  dict_writer = csv.DictWriter(output_file, keys, delimiter='\t')
  dict_writer.writeheader()
  dict_writer.writerows(to_csv)

print("[6] Downloading sequences...")
for item in to_csv:
  print("Downloading {} FASTA sequence...".format(item['filename']))
  cmd = 'efetch -db nuccore -id "{}" -format fasta > {}'.format(item['accession'], item['filename'])
  ps = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
  output = ps.communicate()[0]
```
  
