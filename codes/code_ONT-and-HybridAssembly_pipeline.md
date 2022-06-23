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

# Code for Nanopore long-reads processing and hybrid *de novo* assembly #
## Table of contents ##
<ul>
<li><a href="#1">1. Preprocessing of Nanopore reads</a></li>
<li><a href="#2">2. Assessment of human reads with Kraken 2</a></li>
<li><a href="#3">3. Hybrid de novo assembly with Unicycler</a></li>
<li><a href="#4">4. Benchmarking and QC of assembly with QUAST</a></li>
<li><a href="#5">5. Filter non-monkeypox contigs using Kraken 2 + PlusPF database</a></li>
<li><a href="#6">6. Create final consensus FASTA using reference genome</a></li>
</ul>

<hr>
  
<a name="1"></a>
#### 1. Preprocessing of Nanopore reads:
```Bash
# Filter short and low quality reads using FASTP:
infile=/path/to/ont/all_reads.fastq
outfile=all_reads.filtered.fastq
fastp -w 48 -i ${infile} -l 500 -q 10 -o ${outfile}
```

<a name="2"></a>
#### 2. Assessment of human reads with Kraken 2:
```Bash
# Example of a database (i.e. used in Viral-Recon)
# Also check here: https://benlangmead.github.io/aws-indexes/k2
database="kraken2_human"

# Define input and output files:
infile="all_reads.filtered.fastq"
unclassified="non_human_reads.fastq"
classified="only_human_reads.fastq"
report="kraken2_report.txt"

# Assess the origin of reads using Kraken 2:
kraken2 --db ${database} --report ${report} --unclassified-out ${unclassified} --classified-out ${classified} ${infile}
```

<a name="3"></a>
#### 3. Hybrid de novo assembly with Unicycler:
```Bash
# In case you install Unicycler through conda:
conda activate unicycler

# ILMN non-human reads:
r1="unclassified_1.fastq"
r2="unclassified_2.fastq"

# ONT non-human reads:
ont_reads="non_human_reads.fastq"

# Define output directory:
outdir="/path/to/unicycler/results"
unicycler -1 ${r1} -2 ${r2} -l ${ont_reads} -o ${outdir}
```

<a name="4"></a>
#### 4. Benchmarking and QC of assembly with QUAST:
```Bash
# In case you install QUAST through conda:
conda activate quast

assembly="/path/to/unicycler/results/assembly.fasta"
outdir="/path/to/quast/results"
quast.py --threads 32 --circos --output-dir ${outdir} --reference ${reference} --labels raw_assembly ${assembly}
```

<a name="5"></a>
#### 5. Filter non-monkeypox contigs using Kraken 2 + PlusPF database:
```Bash
# First, download PlusPF database:
wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20210517.tar.gz

# Then untar:
tar -xf k2_pluspf_20210517.tar.gz

# In case you install Kraken2 through conda:
conda activate kraken2

database="/path/to/k2_pluspf"
infile="assembly.fasta"
report="assembly.report.txt"
classified="assembly.classified.fasta"

# Run Kraken 2 with PlusPF database on the raw assembly:
kraken2 --db ${database} --report ${report} --classified-out ${classified} ${infile}

# Check assembly.report.txt and manually remove non-monkeypox contigs from the raw assembly.
```

<a name="6"></a>
#### 6. Create final consensus FASTA using reference genome:
```Bash
# In case you install Minimap2 through conda:
conda activate minimap2
assembly="assembly.onlyMPXV.fasta"
sam="aligned_contigs.sam"

# Map
minimap2 -ax asm5 -o ${output} ${reference} ${assembly}

# Sort BAM
bam="aligned_contigs.bam"
samtools view -bS ${sam} | samtools sort - -o ${bam}

# In case you install BCFtools through conda:
conda activate bcftools
fastq_cns="onlyMPXV.cns.fastq"

# Create a pileup with SAMtools, pipe to BCFtools, and convert to FASTA
samtools mpileup -uf ${reference} ${bam} | bcftools call -c | vcfutils.pl vcf2fq > ${cns_fastq}

# In case you install Seqtk through conda:
conda activate seqtk
cns_fasta="assembly.onlyMPXV.cns.fastq"
seqtk seq -aQ64 -U -n N ${cns_fastq} > ${cns_fasta}
```

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
