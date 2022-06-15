<!-- ------------------ HEADER ------------------ -->
<!-- Developed and maintained by Genomics Division
<!-- of the Institute of Technology an Renewable Energy (ITER)
<!-- Tenerife, Canary Islands, SPAIN
<!-- See the "Contact us" section to collaborate with us to growth
<!-- this repository. ;=)

<!-- ------------------ SECTION ------------------ -->
<p align="left">
  <a href="https://github.com/genomicsITER/monkeypox" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/logos_GH.png" width="75%" /> 
      </a>
</p>

# Monkeypox
A public repository of resources related to MonkeyPox maintained by ITER.

This is the result of an ongoing joint-effort of the following institutions and laboratories:
<ul>
 <li>Servicio de Microbiología, Hospital Universitario Ntra. Sra. de Candelaria, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li>Fundación Canaria Instituto de Investigación Sanitaria de Canarias at the Research Unit, Hospital Universitario Ntra. Sra. de Candelaria, 38010 Santa Cruz de Tenerife, Spain.</li>
 <li>Laboratorio de Inmunología Celular y Viral, Unidad de Farmacología, Facultad de Medicina, Universidad de La Laguna, 38200 San Cristóbal de La Laguna, Spain.</li>
 <li>Genomics Division, Instituto Tecnológico y de Energías Renovables, 38600 Santa Cruz de Tenerife, Spain.</li>
</ul>

<hr>
<!-- ------------------ SECTION ------------------ -->

# Table of contents #
<ul>
  <li><a href="#Virological posts">Virological post: A draft of the first genome sequence of Monkeypox virus associated with the multi-country outbreak in May 2022 from the Canary Islands, Spain</a></li>
  <li><a href="#Bioinformatic pipelines">Bioinformatic pipelines</a></li>
    <ul>
    <li><a href="#Code-Illumina">Code for Illumina short-reads processing</a></li>
    <li><a href="#Code-ONT">Code for Nanopore long-reads processing and hybrid de novo assemby</a></li>
    <li><a href="#List-of-software">List of bioinformatic software used in our pipelines</a></li>
  </ul>
  <li><a href="#Sequences">Sequences</a></li>
  <li><a href="#References">References</a></li>
  <li><a href="#Acknowledgements">Acknowledgements</a></li>
  <li><a href="#License and Attribution">License and Attribution</a></li>
  <li><a href="#Contributing. Contact us">Contributing. Contact us</a></li>
  <li><a href="#Update logs">Update logs</a></li>
</ul>

<hr>
<!-- ------------------ SECTION 1 ------------------ -->

<a name="Virological posts"></a>
## Virological posts ##

[This is an external link](https://github.com/genomicsITER/monkeypox)

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>

<hr>
<!-- ------------------ SECTION 2 ------------------ -->

<a name="Bioinformatic pipelines"></a>
## Bioinformatic pipelines ##

<p align="center">
  <a href="https://www.iter.es" title="Instituto Tecnológico y de Energ&iacute;as Renovables (ITER) / Institute of Technology and Renewable Energy (ITER)">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/MPXV_pipeline_v1.png" width="75%" /> 
  </a>
</p>

<a name="Code-Illumina"></a>
**Code for Illumina short-reads processing**
```Bash
#!/bin/bash

code here

# End of script
```

<a name="Code-ONT"></a>
**Code for Nanopore long-reads processing and hybrid de novo assemby**
```Bash
#!/bin/bash

code here

# End of script
```

<a name="List-of-software"></a>
**List of bioinformatic software used in our pipelines**
<ul>
<li>Conda manual for installation of numerous open-source tools used in these pipelines:<a href="https://docs.conda.io/en/latest/">Conda documentation</a></li>
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

<a name="Sequences"></a>
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

<a name="References"></a>
## References ##

> Published papers

Berthet N, Descorps-Declère S, Besombes C, et al. Genomic history of human monkey pox infections in the Central African Republic between 2001 and 2018. Sci Rep. 2021;11(1):13085. Published 2021 Jun 22. doi: https://doi.org/10.1038/s41598-021-92315-8

Cervantes-Gracia K, Gramalla-Schmitz A, Weischedel J, Chahwan R. APOBECs orchestrate genomic and epigenomic editing across health and disease. Trends Genet. 2021;37(11):1028-1043. doi: https://doi.org/10.1016/j.tig.2021.07.003

Cohen J. Global outbreak puts spotlight on neglected virus. Science. 2022;376(6597):1032-1033. doi:10.1126/science.add2701

Cohen-Gihon I, Israeli O, Shifman O, et al. Identification and Whole-Genome Sequencing of a Monkeypox Virus Strain Isolated in Israel. Microbiol Resour Announc. 2020;9(10):e01524-19. Published 2020 Mar 5. doi: https://doi.org/10.1128/MRA.01524-19

Erez N, Achdout H, Milrot E, et al. Diagnosis of Imported Monkeypox, Israel, 2018. Emerg Infect Dis. 2019;25(5):980-983. doi:10.3201/eid2505.190076

Faye O, Pratt CB, Faye M, et al. Genomic characterisation of human monkeypox virus in Nigeria [published correction appears in Lancet Infect Dis. 2018 Mar;18(3):244]. Lancet Infect Dis. 2018;18(3):246. doi: https://doi.org/10.1016/S1473-3099(18)30043-4

Iizuka I, Saijo M, Shiota T, et al. Loop-mediated isothermal amplification-based diagnostic assay for monkeypox virus infections. J Med Virol. 2009;81(6):1102-1108. doi: https://doi.org/10.1002/jmv.21494

Kraemer MUG, Tegally H, Pigott DM, et al. Tracking the 2022 monkeypox outbreak with epidemiological data in real-time [published online ahead of print, 2022 Jun 8]. Lancet Infect Dis. 2022;S1473-3099(22)00359-0. doi: https://doi.org/10.1016/S1473-3099(22)00359-0

Kugelman JR, Johnston SC, Mulembakani PM, et al. Genomic variability of monkeypox virus among humans, Democratic Republic of the Congo. Emerg Infect Dis. 2014;20(2):232-239. doi: https://doi.org/10.3201/eid2002.130118

Kulesh DA, Loveless BM, Norwood D, et al. Monkeypox virus detection in rodents using real-time 3'-minor groove binder TaqMan assays on the Roche LightCycler. Lab Invest. 2004;84(9):1200-1208. doi: https://doi.org/10.1038/labinvest.3700143

Li D, Wilkins K, McCollum AM, et al. Evaluation of the GeneXpert for Human Monkeypox Diagnosis. Am J Trop Med Hyg. 2017;96(2):405-410. doi:10.4269/ajtmh.16-0567

Li Y, Olson VA, Laue T, Laker MT, Damon IK. Detection of monkeypox virus with real-time PCR assays. J Clin Virol. 2006;36(3):194-203. doi: https://doi.org/10.1016/j.jcv.2006.03.012

Li Y, Zhao H, Wilkins K, Hughes C, Damon IK. Real-time PCR assays for the specific detection of monkeypox virus West African and Congo Basin strain DNA. J Virol Methods. 2010;169(1):223-227. doi: https://doi.org/10.1016/j.jviromet.2010.07.012

Luciani L, Inchauste L, Ferraris O, et al. A novel and sensitive real-time PCR system for universal detection of poxviruses [published correction appears in Sci Rep. 2022 Apr 8;12(1):5961]. Sci Rep. 2021;11(1):1798. Published 2021 Jan 19. doi: https://doi.org/10.1038/s41598-021-81376-4

Maksyutov RA, Gavrilova EV, Shchelkunov SN. Species-specific differentiation of variola, monkeypox, and varicella-zoster viruses by multiplex real-time PCR assay. J Virol Methods. 2016;236:215-220. doi: https://doi.org/10.1016/j.jviromet.2016.07.024

Mucker EM, Hartmann C, Hering D, et al. Validation of a pan-orthopox real-time PCR assay for the detection and quantification of viral genomes from nonhuman primate blood. Virol J. 2017;14(1):210. Published 2017 Nov 3. doi: https://doi.org/10.1186/s12985-017-0880-8

Patrono LV, Pléh K, Samuni L, et al. Monkeypox virus emergence in wild chimpanzees reveals distinct clinical outcomes and viral diversity. Nat Microbiol. 2020;5(7):955-965. doi: https://doi.org/10.1038/s41564-020-0706-0

Pecori R, Di Giorgio S, Paulo Lorenzo J, Nina Papavasiliou F. Functions and consequences of AID/APOBEC-mediated DNA and RNA deamination [published online ahead of print, 2022 Mar 7]. Nat Rev Genet. 2022;1-14. doi: https://doi.org/10.1038/s41576-022-00459-8

Shchelkunov SN, Shcherbakov DN, Maksyutov RA, Gavrilova EV. Species-specific identification of variola, monkeypox, cowpox, and vaccinia viruses by multiplex real-time PCR assay. J Virol Methods. 2011;175(2):163-169. doi: https://doi.org/10.1016/j.jviromet.2011.05.002

Tumewu J, Wardiana M, Ervianty E, et al. An adult patient with suspected of monkeypox infection differential diagnosed to chickenpox. Infect Dis Rep. 2020;12(Suppl 1):8724. Published 2020 Jul 6. doi: https://doi.org/10.4081/idr.2020.8724

Yong SEF, Ng OT, Ho ZJM, et al. Imported Monkeypox, Singapore. Emerg Infect Dis. 2020;26(8):1826-1830. doi:10.3201/eid2608.191387

Zhao K, Wohlhueter RM, Li Y. Finishing monkeypox genomes from short reads: assembly analysis and a neural network method. BMC Genomics. 2016;17 Suppl 5(Suppl 5):497. Published 2016 Aug 31. doi: https://doi.org/10.1186/s12864-016-2826-8



> Preprint papers

Enhanced surveillance of monkeypox in Bas-Uélé, Democratic Republic of Congo: the limitations of symptom-based case definitions
Gaspard Mande, Innocent Akonda, Anja De Weggheleire, Isabel Brosius, Laurens Liesenborghs, Emmanuel Bottieau, Noam Ross, Guy -Crispin Gembu, Robert Colebunders, Erik Verheyen, Ngonda Dauly, Herwig Leirs, Anne Laudisoit
medRxiv 2022.06.03.22275815; doi: https://doi.org/10.1101/2022.06.03.22275815

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 5 ------------------ -->

<a name="Acknowledgements"></a>
## Acknowledgements ##

This study has been funded by Cabildo Insular de Tenerife (CGIEU0000219140 and "Apuestas científicas del ITER para colaborar en la lucha contra la COVID-19"), Instituto de Salud Carlos III (FI18/00230) cofunded by European Union (ERDF) "A way of making Europe", and and by the agreement with Instituto Tecnológico y de Energías Renovables (ITER) to strengthen scientific and technological education, training, research, development and innovation in Genomics, Personalized Medicine and Biotechnology (OA17/008).

We acknowledge in [Table 1 (EXCEL file)](https://github.com/genomicsITER/monkeypox/blob/main/tables/table1.xlsx) the researchers and their institutions who released the Monkeypox sequences through NCBI GenBank that are being used in our studies. 

We also gratefully acknowledge the authors, originating and submitting laboratories of the genetic sequences and metadata for sharing their work as shown in [Nextstrain](https://nextstrain.org/monkeypox) and:
<ul>
<li>Hadfield et al, Nextstrain: real-time tracking of pathogen evolution, Bioinformatics (2018)</li>
  <li>Sagulenko et al, TreeTime: Maximum-likelihood phylodynamic analysis, Virus Evolution (2017)</li>
</ul>

And we also want to thank the contributions of a number of researchers and laboratories sharing their preliminary results through the [Virological](https://virological.org/) website.

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 6 ------------------ -->

<a name="License and Attribution"></a>
## License and Attribution ##

This repository and data exports are published under the CC BY 4.0 license. Please, acknowledge authors, originating and submitting laboratories of the genetic sequences and metadata and open source software used in this work.

Please cite as: "Monkeypox repository of the Refence Laboratory for the Epidemiological Surveillance of Pathogens in the Canary Islands (accessed on YYYY-MM-DD)".

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 7 ------------------ -->

<a name="Contributing. Contact us"></a>
## Contributing. Contact us ##

> Want to share your relevant links? Place a Direct Message to @resocios or @adrmunozb on Twitter (see below).

 <p align="left">
  <a href="https://www.iter.es/areas/area-genomica/" title="Contact us at the Genomics Division of the Institute of Technology and Renewable Energy (ITER), Tenerife, Canary Islands, Spain">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/ITER_logo.png" width="30%" /> 
  </a>
</p>

By AMB <a href="https://twitter.com/adrmunozb" title="Follow to @resocios on Twitter" >@adrmunozb <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/Twitter_Social_Icon_Circle_Color.png" width="32px" /></a> and JMLS <a href="https://twitter.com/resocios" title="Follow to @resocios on Twitter" >@resocios <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/Twitter_Social_Icon_Circle_Color.png" width="32px" /></a>


<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>


<hr>
<!-- ------------------ SECTION 8 ------------------ -->

<a name="Update logs"></a>
## Update logs ##

> June 13, 2022 updates. Created the public version of this repository. Enjoy the reading! ;=)

<p align="right">
  <a href="#Monkeypox" title="Up">
    <img src="https://github.com/genomicsITER/monkeypox/blob/main/images/home-icon.png" style="float: right; margin: 10px; padding: 2px;" />
  </a>
</p>
