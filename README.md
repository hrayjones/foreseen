<p align="center">
  <img src="https://github.com/hrayjones/foreseen/blob/main/foreseen.png?raw=true" width="300">
</p>

<p align="center">
***** UNDER DEVELOPMENT *****
</p>

__TODO__:

Make diagram for the main workflow

Add requirements

Add installation instructions

Add basic usage instructions

Update the main workflow with specifics

</p>

## Foreseen: An automated pipeline for 4C-seq experimental design, with or without genetic variants of interest

4C-seq is a powerful method for profiling the 3D genome at specific genomic regions of interest. Depending on the design, a 4C-seq experiment can be used to identify chromatin interactions in a genotype-specific fashion, which has clear utility in resolving the effect of genotype on gene regulatory mechanisms. However, designing 4C-seq is complicated and time-consuming. Moreover, a manual design process does not necessarily explore all possibilities for experimental design. The foreseen pipeline automates the design process, taking as input a genomic locus and running through all possible combinations of valid restriction enzymes and primers to suggest the most optimal experimenal parameters for the 4C-seq experiment.

### Requirements
- Python version XXX
- Primer 3
- ...

### Installation
*installation instructions*

### Basic usage
*basic usage instructions*

#### The modes are as follows:
1) __Basic__ (classic 4C-seq; no genotype-specificity required)
2) __Allele-specific__ (4C-seq that is focused on a specific allele of a single nucleotide polymorphism [SNP] of interest)
3) __Allele-aware__ (4C-seq that yields separable sequencing results for both alleles of a SNP of interest in a heterozygous sample)

## How does Foreseen work?
Foreseen works off the basic principles, as outlined by prior publications, for the best 4C-seq experimental designs [[1]](#1)[[2]](#2).

The main steps in Foreseen are as follows:
1) __Expand region of interest__

Given a region of interest, which may be in the form of genomic coordinates or sequence (fasta), extend the edges in order to capture the most useful restriction cut sites in relatively close proximity to the region of interest. If the input region is <1000 bp, the default action is to expand by 2000bp each way. If the region is >1000bp, then it is expanded by 1000bp each way.

2) __Explore restriction enzyme combinations with in silico digestion__

Foreseen uses input lists of valid restriction enzymes that can be used for 4C-seq experiments. It checks that the resultant viewpoint fragments are "non blind", i.e. they are flanked by both restriction enzyme 1 (RE1) and restriction enzyme 2 (RE2), and that they are of an optimal length for digestion and ligation in the 4C-seq experiment (200-1500bp). If the _allele-aware_ mode is used, the enzymes will only be selected if the given SNP falls within 100bp of a restriction cut site. If the _allele-specific_ mode is used, the SNP must affect a (non-N) base within one of the restriction cut sites.

3) __Identify primers__

Reading and non-reading primers are first identified separately, using Primer3 software. Primers must fall within close proximity to the cut sites, as recommended in [[1]](#1) and [[2]](#2). This is defined in part by the expected length of sequencing reads, supplied by the user. If the _allele-aware_ mode is used, the SNP must fall between the primary restriction cut site and the reading primer (or the secondary restriction cut site and the non-reading primer, in the case of paired end sequencing). Primers are then paired up based on their molecular properties (GC content, melting temperature and product size).

4) __Generate final design score__

A score is generated for all valid combinations of restriction enzymes and primers. The output is a ranked list of experimental designs. 

 
## References
<a id="1">[1]</a> 
Krijger, Peter H.L., Geert Geeven, Valerio Bianchi, Catharina R.E. Hilvering, and Wouter de Laat (2020).
4C-Seq from Beginning to End: A Detailed Protocol for Sample Preparation and Data Analysis. 
Genome Architecture, 170 (January):17–32.

<a id="2">[2]</a> 
Miranda, Mélanie, Daan Noordermeer, and Benoit Moindrot (2022). 
Detection of Allele-Specific 3D Chromatin Interactions Using High-Resolution In-Nucleus 4C-Seq. 
Methods in Molecular Biology (Clifton, N.J.) 2532:15–33.

