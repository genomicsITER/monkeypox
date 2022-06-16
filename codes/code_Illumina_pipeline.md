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

# Code for Illumina short-reads processing #
## Table of contents ##
<ul>
<li><a href="#1">1. FASTQ Aggregation and Interleaving with BBMap</a></li>
<li><a href="#2">2. Human reads Scrubbing</a></li>
<ul>
<li><a href="#2.1">2.1. With NCBI SRA Human Scrubber</a></li>
<li><a href="#2.2">2.2. With Kraken 2</a></li>
</ul>
<li><a href="#3">3. Mapping non-human reads to the reference sequence</a></li>
<ul>
<li><a href="#3.1">3.1. With BWA</a></li>
<li><a href="#3.2">3.2. With Minimap2</a></li>
<li><a href="#3.3">3.3. With Bowtie2</a></li>
</ul>
<li><a href="#4">4. Sort and index aligned reads</a></li>
<li><a href="#5">5. Variant Calling</a></li>
<ul>
<li><a href="#5.1">5.1. With iVar</a></li>
<li><a href="#5.2">5.2. With LoFreq</a></li>
</ul>
<li><a href="#6">6. Consensus FASTA Generation</a></li>
<li><a href="#7">7. Multisample Alignment and Phylogenetic Analysis</a></li>
  <ul>
<li><a href="#7.1">7.1. with MAFFT and IQ-TREE</a></li>
<li><a href="#7.2">7.2. With Nextstraing-monkeypox</a></li>
</ul>
</ul>

<hr>

<a name="1"></a>
#### 1. FASTQ Aggregation and Interleaving with BBMap:
```Bash
# Aggregate R1 (forward) FASTQ files:
zcat *_R1_001.fastq.gz > sample_r1.fastq.gz

# Aggregate R2 (reverse) FASTQ files:
zcat *_R2_001.fastq.gz > sample_r2.fastq.gz

# Path to BBMap installation path:
bbmap="/path/to/bbmap"
  
# Aggregated FASTQ:
r1="sample_r1.fastq.gz"
r2="sample_r2.fastq.gz"
  
# Interleaved FASTQ file:
outfile="interleaved.fastq"
  
# Run BBMap reformat command:
${bbmap}/reformat.sh n1="${r1}" in2="${r2}" out="${outfile}"
```

<a name="2"></a>
#### 2. Human reads Scrubbing:
<a name="2.1"></a>
##### 2.1 With NCBI SRA Human Scrubber:
```Bash
# Path to sra-human-scrubber installation path:
sraScrubber="/path/to/sra-human-scrubber/scripts/scrub.sh"

# Path to human database:
db="/path/to/sra-human-scrubber/data/human_filter.db"
  
# Interleaved FASTQ:
infile="interleaved.fastq"

# Scrubbed FASTQ:
outfile="interleaved.scrubbed.fastq"

# Run sra-human-scrubber command:
${sraScrubber} -d ${db} -i ${infile} -o ${outfile}
```

<a name="2.2"></a>
##### 2.2 With Kraken 2:
```Bash
# In case you install Kraken 2 through conda:
conda activate kraken2
  
# Path to Kraken 2 Human database (check this Kraken 2 manual section to learn how to download and build your own database.
# See: https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#custom-databases
database="/path/to/kraken2-human-db"

# Aggregated FASTQ files (interleaved reads are not necessary in this case):
r1="sample_r1.fastq.gz"
r2="sample_r2.fastq.gz"

# Output FASTQ with unclassified reads (# symbol is for paired-end reads):
unclassified="unclassified#.fastq"
  
# Output FASTQ with classified reads (# symbol is for paired-end reads):
classified="classified#.fastq"
  
# Kraken 2 classification report:
report="report.txt"
  
# Run Kraken 2 command:
kraken2 --db ${database} --unclassified-out ${unclassified} --classified-out ${classified} --report ${report} --paired ${infile}
```

<a name="3"></a>
#### 3. Mapping non-human reads to the reference sequence:
<a name="3.1"></a>  
##### 3.1 With BWA:
```Bash
# In case you install BWA through conda:
conda activate bwa
  
# Define Monkeypox reference genome:
reference="MT903344.1.fasta"
  
# Define FASTQ files (only ONE interlaved FASTQ in case you use SRA-human-scrubber or TWO if you use Kraken 2)
# Uncomment the corresponding code line/s
#interleaved="interleaved.scrubbed.fastq"
r1="unclassified_1.fastq"
r2="unclassified_2.fastq"
  
# Define aligned SAM file:
outfile="aligned.sam"
  
# Run BWA-MEM command:
bwa mem -Y ${reference}  ${r1} ${r2} > ${outfile}
```
<a name="3.2"></a> 
##### 3.2 With Minimap2:
```Bash
# In case you install Minimap2 through conda:
conda activate minimap2
  
# Define Monkeypox reference genome:
reference="MT903344.1.fasta"
  
# Define FASTQ files (only ONE interlaved FASTQ in case you use SRA-human-scrubber or TWO if you use Kraken 2)
# Uncomment the corresponding code line/s
#interleaved="interleaved.scrubbed.fastq"
r1="unclassified_1.fastq"
r2="unclassified_2.fastq"
  
# Define aligned SAM file:
outfile="aligned.sam"

# Run Minimap2 command:
minimap2 -ax sr ${reference} ${r1} ${r2} > ${outfile}
```
<a name="3.3"></a> 
##### 3.3 With Bowtie2:
```Bash
# In case you install Bowtie2 through conda:
conda activate bowtie2
  
# Define Monkeypox reference genome:
reference="MT903344.1.fasta"

# Define FASTQ files (only ONE interlaved FASTQ in case you use SRA-human-scrubber or TWO if you use Kraken 2)
# Uncomment the corresponding code line/s
#interleaved="interleaved.scrubbed.fastq"
r1="unclassified_1.fastq"
r2="unclassified_2.fastq"
  
# Define aligned SAM file:
outfile="aligned.sam"

# Before aligning reads, indexing the reference genome is needed.
# You need to run this command in the reference directory:
bowtie2-build ${reference} MT903344.1

# Run Bowtie2 command:
index_prefix=/path/to/index/files
bowtie2 -x ${index_prefix} -1 ${r1} -2 ${r2} -S ${outfile}
```

<a name="4"></a> 
#### 4. Sort and index aligned reads:
```Bash
# Sort SAM:
infile="aligned.sam"
outfile="aligned.sorted.bam"
samtools sort ${infile} > ${outfile}

# Discard unmapped reads from sorted BAM file:
infile="aligned.sorted.bam"
outfile="aligned.sorted.mapped.bam"
samtools view -F 0x04 -b ${infile} > ${outfile}

# Index BAM file:
samtools index ${outfile}
```

<a name="5"></a> 
#### 5. Variant Calling:
<a name="5.1"></a> 
##### 5.1 With iVar:
```Bash
# In case you install iVar through conda:
conda activate ivar

# Make a pileup and pipe to iVar to call variants:
infile="aligned.sorted.mapped.bam"
prefix="out_variants"
samtools mpileup --reference ${reference} ${infile} | ivar variants -r ${reference} -p ${prefix}
```
<a name="5.2"></a> 
##### 5.2 With LoFreq:
```Bash
# In case you install LoFreq through conda:
conda activate lofreq

# Call variants
infile="aligned.sorted.mapped.bam"
outfile="variants.vcf"
lofreq call -f ${reference} -o ${outfile} ${infile}
```

<a name="6"></a> 
#### 6. Consensus FASTA Generation:
```Bash
#In preparation...
```

<a name="7"></a> 
#### 7. Multisample Alignment and Phylogenetic Analysis:
```Bash
#In preparation...
```

<a name="7.1"></a>  
##### 7.1. With MAFFT and IQ-TREE:
```Bash
#In preparation...
  ```

<a name="7.2"></a>  
##### 7.2. With Nextstraing-monkeypox:
```Bash
#In preparation...
  ```
  
<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
