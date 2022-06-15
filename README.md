<!-- ------------------ HEADER ------------------ -->
<!-- Developed and maintained by Genomics Division
<!-- of the Institute of Technology an Renewable Energy (ITER)
<!-- Tenerife, Canary Islands, SPAIN
<!-- See the "Contact us" section to collaborate with us to growth
<!-- this repository. ;=)

<!-- ------------------ SECTION ------------------ -->
<p align="left">
  <a href="https://www.iter.es" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/ITER_logo.png" width="30%" /> 
    
  </a>
</p>

# Monkeypox
A public repository of resources related to MonkeyPox maintained by ITER.

This is the result of a joint-effort of the following institutions and laboratories:
<ul>
 <li>Servicio de Microbiología, Hospital Universitario Ntra. Sra. de Candelaria, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li>Fundación Canaria Instituto de Investigación Sanitaria de Canarias at the Research Unit, Hospital Universitario Ntra. Sra. de Candelaria, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li>Laboratorio de Inmunología Celular y Viral, Unidad de Farmacología, Facultad de Medicina, Universidad de La Laguna, 38200 San Cristóbal de La Laguna, Spain.</li>
 <li>Genomics Division, Instituto Tecnológico y de Energías Renovables, 38600 Santa Cruz de Tenerife, Spain.</li>
</ul>

<hr>
<!-- ------------------ SECTION ------------------ -->

> June updates. Created the public version of this repository. Enjoy the reading! ;=)

# Table of contents #
<ul>
  <li><a href="#Virological">Virological post: A draft of the first genome sequence of Monkeypox virus associated with the multi-country outbreak in May 2022 from the Canary Islands, Spain</a></li>
  <li><a href="#ilink2">Bioinformatic pipelines</a></li>
  <li><a href="#ilink3">Sequences</a></li>
  <li><a href="#ilink4">References</a></li>
  <li><a href="#ilink5">Acknowledgements</a></li>
  </ul>

<hr>
<!-- ------------------ SECTION 1 ------------------ -->

<a name="ilink1"></a>
## Virological posts ##

[This is an external link](https://github.com/genomicsITER/monkeypox)

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>

<hr>
<!-- ------------------ SECTION 2 ------------------ -->

<a name="ilink2"></a>
## Bioinformatic pipelines ##

<p align="center">
  <a href="https://www.iter.es" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/MPXV_pipeline_v1.png" width="75%" /> 
  </a>
</p>

**Excerpt of Code for Illumina short-reads processing**
```Bash
#!/bin/bash

code here

# End of script
```

**Excerpt of Code for Nanopore long-reads processing**
```Bash
#!/bin/bash

code here

# End of script
```

**List of bioinformatic software used in our pipelines**
<ul>
<li>Reformat FASTQ files to get an interleaved FASTQ file: <a href="https://sourceforge.net/projects/bbmap/">BBMap tools v.38.96</a></li>
<li>Remove Human mapping-reads from your FASTQ files: <a href="https://github.com/ncbi/sra-human-scrubber/">NCBI SRA Human Scrubber v.1.0.2021_05_05</a></li>
<li>Remove Human mapping-reads from your FASTQ files: <a href="https://ccb.jhu.edu/software/kraken2/">Kraken2 v.2.1.2</a></li>
<li>Programming environment of general purpose: <a href="https://www.r-project.org/">R v.4.1.3</a></li>
<li>Compute the depth of coverage and other statistics: <a href="https://github.com/brentp/mosdepth/">Mosdepth v.0.3.3</a></li>
<li>Compute de number of duplicates and other statistics: <a href="https://broadinstitute.github.io/picard/">Picard Tools v.2.18.7</a></li>
<li>Perform the variant calling and consensus: <a href="https://github.com/andersen-lab/ivar/">iVar v.1.3.1</a></li>
<li>Perform the variant calling: <a href="https://csb5.github.io/lofreq/">LoFreq v.2.1.5</a></li>
<li>Get mapping statistics, manipulate BAM files, and generate mpileups for FASTA consensus: <a href="https://github.com/samtools/samtools">SAMtools v.1.6</a></li>
<li>Multiple Sample Alignment: <a href="https://mafft.cbrc.jp/alignment/server/">MAFFT v.7.505]</a></li>
<li>Phylogenomic inference and tree computing: <a href="http://www.iqtree.org/">IQ-TREE v.2.2.0.3</a></li>
<li>Mapping of short-reads: <a href="https://github.com/lh3/minimap2/">Minimap2 v.2.24-r1122</a></li>
<li>Mapping of short-reads: <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2 v.2.4.5</a></li>
<li>Mapping of short-reads: [BWA v.0.7.17-r1188<a href="https://github.com/lh3/bwa/"></a></li>
<li>Framework for analyses and visualization of pathogen genome data (Monkeypox in this case): <a href="https://github.com/nextstrain/monkeypox">Nextstrain</a></li>
<li>Assembler: <a href="https://github.com/rrwick/Unicycler/">Unicycler v.0.5.0</a></li>
<li>Benchmarking and quality control of assemblies: <a href="http://quast.sourceforge.net/">QUAST v.5.0.2</a></li>
<li>Visualization of assemblies: <a href="https://rrwick.github.io/Bandage/">Bandage v.0.9.0</a></li>
<li>Visualization of Kraken 2 reports: <a href="https://ccb.jhu.edu/software/pavian/">Pavian v.1.0</a></li>
<li>Annotation of genomes: <a href="https://pcingola.github.io/SnpEff/">SnpEff v.5.1d</a></li>
<li>Visualization of phylogenetic trees: <a href="http://tree.bio.ed.ac.uk/software/figtree/">Figtree</a></li>
<li>Visualization of phylogenetic trees: <a href="https://bioconductor.org/packages/release/bioc/html/ggtree.html/">ggtree 3.15</a></li>
</ul>

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 3 ------------------ -->

<a name="ilink3"></a>
## Sequences ##

> Consensus FASTA file obtained from Illumina short-reads mapping to [MT903344.1](https://www.ncbi.nlm.nih.gov/nuccore/MT903344.1). See [Virological post for more details](https://virological.org/)

[HUNSC_ITER_MPXV01.fasta](https://github.com/genomicsITER/monkeypox/blob/main/sequences/HUNSC_ITER_MPXV01.fasta)

> Consensus FASTA file obtained from a hybrid de novo Illumina-Nanopore based assembly and [MT903344.1](https://www.ncbi.nlm.nih.gov/nuccore/MT903344.1). See [Virological post for more details](https://virological.org/)

[HUNSC_ITER_MPXV01_ILMN-ONT_hybrid_assembly.fasta](https://github.com/genomicsITER/monkeypox/blob/main/sequences/HUNSC_ITER_MPXV01_ILMN-ONT_hybrid_assembly.fasta)

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>

<hr>
<!-- ------------------ SECTION 4 ------------------ -->

<a name="ilink4"></a>
## References ##

[This is an external link](https://github.com/genomicsITER/monkeypox)

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>

<hr>
<!-- ------------------ SECTION 5 ------------------ -->

<a name="ilink5"></a>
## Acknowledgements ##

This study has been funded by Cabildo Insular de Tenerife (CGIEU0000219140 and "Apuestas científicas del ITER para colaborar en la lucha contra la COVID-19"), Instituto de Salud Carlos III (FI18/00230) cofunded by European Union (ERDF) "A way of making Europe", and and by the agreement with Instituto Tecnológico y de Energías Renovables (ITER) to strengthen scientific and technological education, training, research, development and innovation in Genomics, Personalized Medicine and Biotechnology (OA17/008).

We acknowledge the researchers and their institutions who released the Monkeypox sequences through NCBI GenBank that are being used in 
[Table 1 (EXCEL file)](https://github.com/genomicsITER/monkeypox/blob/main/tables/table1.xlsx).



<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
