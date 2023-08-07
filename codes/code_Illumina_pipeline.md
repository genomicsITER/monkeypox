<!-- ------------------ HEADER ------------------ -->
<!-- Developed and maintained by Genomics Division -->
<!-- of the Institute of Technology an Renewable Energy (ITER) -->
<!-- Tenerife, Canary Islands, SPAIN -->
<!-- See the "Contact us" section to collaborate with us to growth -->
<!-- this repository. ;=) -->

<!-- ------------------ SECTION ------------------ -->
<p align="left">
  <a href="https://github.com/genomicsITER/influenza" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/logos_GH.png" width="auto" /> 
      </a>
</p>

# Influenza
A public repository of resources related to Influenza analysis maintained by ITER.

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
<li><a href="#1">1. Quality Control of raw reads with FastQC</a></li>
<li><a href="#2">2. Adapter trimming with fastp</a></li>
<li><a href="#3">3. Quality Control of trimmed reads with FastQC</a></li>
<li><a href="#4">4. Remove host reads using Kraken2 with humanDB database</a></li>
<li><a href="#5">5. Preliminary assembly using Unicycler</a></li>
<li><a href="#6">6. Detect hits using BLASTn with NCIB Influenza Virus Database</a></li>
<ul>
<li><a href="#6.1">6.1. Download and build NCBI Influenza Virus database</a></li>
<li><a href="#6.2">6.2. Run BLASTn to get the top hit from the NCBI Influenza Virus Database</a></li>
<li><a href="#6.3">6.3. Parse BLASTn results to detect the strain to use as reference</a></li>
</ul>
<li><a href="#7">7. Align non-host reads to a reference genome</a></li>
<ul>
<li><a href="#7.1">7.1. Select the reference based on BLASTn results</a></li>
<li><a href="#7.2">7.2. Align reads by segment using BWA</a></li>
<li><a href="#7.3">7.3. Discard unmapped reads and sort BAMs file using SAMtools</a></li>
<li><a href="#7.4">7.4. Quality Control of sorted BAMs using SAMtools and MosDepth</a></li>
</ul>
<li><a href="#8">8. Variant Calling with iVar</a></li>
<li><a href="#9">9. Consensus FASTA generation with iVar</a></li>
<li><a href="#10">10. Lineage classification of segment 4 (HA) with NextClade</a></li>
<li><a href="#11">11. Multisample Alignment and Phylogenetic Analysis</a></li>
</ul>

<hr>

<a name="1"></a>
#### 1. Quality Control of raw reads with FastQC:
```Bash
# In case you install FastQC through conda:
conda activate fastqc

r1="sample_R1_001.fastq.gz"
r2="sample_R2_001.fastq.gz"
outdir="./fastqc_results"

# For each sample run following commands:
fastqc ${r1} --outdir ${outdir}
fastqc ${r2} --outdir ${outdir}

conda deactivate
```

<a name="2"></a>
#### 2. Adapter trimming with fastp:
```Bash
conda activate fastp

threads=16

adapters="Influenza.adapters.fasta"

r1="sample_R1_001.fastq.gz"
r2="sample_R2_001.fastq.gz"

out1="sample_R1_001.fastp.fastq.gz"
out2="sample_R2_001.fastp.fastq.gz"

json="sample.fastp.json"
html="sample.fastp.html"

fastp --thread ${threads}\
  --in1 ${r1} \
  --in2 ${r2} \
  --out1 ${out1} \
  --out2 ${out2} \
  --json ${json} \
  --html ${html} \
  --adapter_fasta ${adapters}

conda deactivate
```

<a name="3"></a>
#### 3. Quality Control of trimmed reads with FastQC:
```Bash
# In case you install FastQC through conda:
conda activate fastqc

trimmed_r1="sample_R1_001.fastp.fastq.gz"
trimmed_r2="sample_R2_001.fastp.fastq.gz"
outdir="./fastqc_trimmed_results"

fastqc ${trimmed_r1} --outdir ${outdir}
fastqc ${trimmed_r2} --outdir ${outdir}

conda deactivate
```

<a name="4"></a>
#### 4. Remove host reads using Kraken2 with humanDB database:
```Bash
# In case you install Kraken2 through conda:
conda activate kraken2

database="/path/to/kraken2/databases/kraken2-human-db"

trimmed_r1="sample_R1_001.fastp.fastq.gz"
trimmed_r2="sample_R2_001.fastp.fastq.gz"

classified="sample.classified#.fastq"
unclassified="sample.unclassified#.fastq"

report="sample.kraken2-report.txt"

kraken2 --db ${database} \
  --unclassified-out ${unclassified} \
  --classified-out ${classified} \
  --report ${report} \
  --paired ${r1} ${r2}

conda deactivate
```

<a name="5"></a>
#### 5. Preliminary assembly using Unicycler:
```Bash
# In case you install Unicycler through conda:
conda activate unicycler

threads=16

r1="sample.unclassified_1.fastq"
r2="sample.unclassified_2.fastq"

outdir="./unicycler_results"

unicycler --threads ${threads} \
  -1 ${r1} \
  -2 ${r2} \
  --out ${outdir}

mv ${outdir}/assembly.fasta sample.unicycler.fasta

conda deactivate
```

<a name="6"></a>
#### 6. Detect hits using BLASTn with NCIB Influenza Virus Database:
<a name="6.1"></a>  
##### 6.1 Download and build NCBI Influenza Virus database:
```Bash
# 1. Download latest NCBI Influenza DB sequences and metadata
wget https://ftp.ncbi.nih.gov/genomes/INFLUENZA/influenza.fna.gz

# 2. Gunzip NCBI FLU FASTA
# Replace FASTA headers >gi|{gi}|gb|{accession}|{description} with >{accession}|{description} for easier parsing and processing
zcat influenza.fna.gz | sed -E 's/^>gi\|[0-9]+\|gb\|(\w+)\|/>/' > influenza.fna

# 3. Make BLASTDB:
# In case you install BLAST through conda:
conda activate blast

makeblastdb -in influenza.fna

mkdir blast_db

mv influenza.fna.* blast_db/

conda deactivate
```
<a name="6.2"></a> 
##### 6.2 Run BLASTn to get the top hit from the NCBI Influenza Virus Database:
```Bash
# In case you install BLAST through conda:
conda activate blast

database="/path/to/NCBI_Influenza_Virus_Database/blast_db/influenza.fna"

threads=16

asm="sample.unicycler.fasta"
out="sample.unicycler.blastn.txt"

blastn -num_threads ${threads} \
  -db ${database} \
  -num_alignments 1 \
  -outfmt "6 stitle" \
  -query ${asm} \
  -out ${out}

conda deactivate
```
<a name="6.3"></a> 
##### 6.3 Parse BLASTn results to detect the strain to use as reference:
```Bash
# This code creates a file with the strain to use as reference, as "A/H1N1", "A/H3N2" or "B":

infile="sample.unicycler.blastn.txt"
outfile="sample.strain.txt"

touch ${outfile}

strain=$( cat ${infile} | cut -d"(" -f2 | cut -d"/" -f1 | uniq )

if [[ -n ${strain} ]]; then
  if [[ ${strain} == "A" ]]; then
    substrain=$( cat "{infile}" | cut -d"(" -f3 | cut -d")" -f1 | uniq )
    echo -e ${strain}"/"${substrain} >> ${outfile}
  else
    # If ${strain}=="B", only save the strain:
    echo -e ${strain} >> ${outfile}
  fi
fi
```

<a name="7"></a> 
#### 7. Align non-host reads to a reference genome:
<a name="7.1"></a>  
##### 7.1 Select the reference based on BLASTn results:
```Bash

infile="sample.strain.txt"

strain=$( cat ${infile} )

# Path to multiple Influenza strains reference genomes separated by segments:
infA_h1n1_dir=/path/to/reference/Influenza_A_virus_H1N1_California/separated_segments
infA_h3n2_dir=/path/to/referenceInfluenza_A_virus_H3N2_Wisconsin/separated_segments
infB_dir=/path/to/referenceInfluenza_B_virus_Yamagata/separated_segments

# Select reference folder:
if [[ ${strain} == "A-H1N1" ]]; then
  reference_dir=${infA_h1n1_dir}
elif [[ ${strain} == "A-H3N2" ]]; then
  reference_dir=${infA_h3n2_dir}
elif [[ ${strain} == "B" ]]; then
  reference_dir=${infB_dir}
else
  echo "LOG: There is no FASTA file for strain ${strain}."
  exit 1
fi
```

<a name="7.2"></a> 
##### 7.2 Align reads by segment using BWA:
```Bash
# In case you install BWA through conda:
conda activate bwa

# Non-host reads from step 4:
r1="sample.unclassified_1.fastq"
r2="sample.unclassified_2.fastq"

for i in $( seq 1 8 ); do
  reference=${reference_dir}/Seg${i}.fasta
  outfile=sample.${strain}.seg${i}.sam

  bwa mem -Y ${reference} ${r1} ${r2} > ${outfile}
done

conda deactivate
```

<a name="7.3"></a> 
##### 7.3 Discard unmapped reads and sort BAMs file using SAMtools:
```Bash
# In case you install SAMtools through conda:
conda activate samtools

for i in $( seq 1 8 ); do
  infile=sample.${strain}.seg${i}.sam
  mapped=sample.${strain}.mapped.seg${i}.bam
  samtools view -F 0x04 -b ${infile} > ${mapped}
  samtools index ${mapped}

  sorted=sample.${strain}.mapped.sorted.seg${i}.bam
  samtools sort ${mapped} > ${sorted}
  samtools index ${sorted}
done

conda deactivate
```

<a name="7.4"></a> 
##### 7.4 Quality Control of sorted BAMs using SAMtools and MosDepth:
```Bash
# In case you install SAMtools through conda:
conda activate samtools

for i in $( seq 1 8 ); do
  infile=sample.${strain}.mapped.sorted.seg${i}.bam
  flagstat=sample.${strain}.mapped.sorted.seg${i}.flagstat
  idxstats=sample.${strain}.mapped.sorted.seg${i}.idxstats

  samtools flagstat ${infile} > ${flagstat}
  samtools idxstats ${infile} > ${idxstats}
done

conda deactivate

# In case you install MosDepth through conda:
conda activate mosdepth

for i in $( seq 1 8 ); do
  infile=sample.${strain}.mapped.sorted.seg${i}.bam
  prefix=sample.${strain}.mapped.sorted.seg${i}

  threads=16

  mosdepth --threads ${threads} --no-per-base ${prefix} ${infile}
done

conda deactivate
```

<a name="8"></a> 
#### 8. Variant Calling with iVar:
```Bash
# In case you install iVar through conda:
conda activate ivar

# Make a pileup and pipe to iVar to call variants:
for i in $( seq 1 8 ); do
  reference=${reference_dir}/Seg${i}.fasta

  infile="sample.${strain}.mapped.sorted.seg${i}.bam"
  out_prefix="sample.${strain}.seg${i}.out_variants"

  samtools mpileup --reference ${reference} ${infile} | ivar variants -r ${reference} -p ${out_prefix}
done

conda deactivate
```

<a name="9"></a> 
#### 9. Consensus FASTA Generation:
```Bash
# In case you install iVar through conda:
conda activate ivar

# Generate consensus FASTA:
# Optionally, you can set different parameters to define minimum thresholds for the consensus 
# (see ivar consensus help).

for i in $( seq 1 8 ); do
  infile="sample.${strain}.mapped.sorted.seg${i}.bam"
  outfile="sample.${strain}.seg${i}.consensus.fasta"

  samtools mpileup -A -Q 0 ${infile} | ivar consensus -p ${outfile} -q 10 -t 0 -m 1
done

conda deactivate
```

<a name="10"></a> 
#### 10. Lineage classification of segment 4 (HA) with NextClade:
```Bash
# Create a multisample FASTA with more sequences separated by strain:
sequences_dir="/path/to/more/sequences"
infile="sample.${strain}.seg4.consensus.fasta"
outfile="multifasta.${strain}.fa"

cat ${sequences_dir}/*.${strain}.seg4.consensus.fa ${infile} > ${outfile}

...
```

<a name="11"></a> 
#### 11. Multisample Alignment and Phylogenetic Analysis of segment 4 (HA):
```Bash
# Align multisample FASTA using MAFFT:
# You should replace the --thread parameter with a proper value that suits your execution environment.
infile="multifasta.${strain}.fa"
outfile="multifasta.${strain}.mafft-aligned.fa"
logfile="multifasta.${strain}.mafft-aligned.log"
mafft --reorder --anysymbol --nomemsave --adjustdirection --thread 48 ${infile} 1> ${outfile} 2> ${logfile}

# Generate phylogenetic tree using IQ-TREE with best-fit model:
# Alternatively, you could choose a model like JC by using the -m parameter.
iqtree -s ${outfile} -nt AUTO

# You will end up with a phylogenetic tree in Newick format (*.treefile extension) which can be represented 
# using your favorite tool such as NCBI Tree Viewer.
```
  
<p align="right">
  <a href="#Influenza" title="Up">
    <img src="https://github.com/genomicsITER/influenza/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
