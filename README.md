# tcr-snps

Extract reference and variant sequences for SNPs and other polymorphisms in the human T-cell receptor loci

#### version 3.4, [Jamie Heather](http://jamieheather.github.io), March 2018
<hr>

### About

T-cell receptors (TCRs) are variable antigen receptors produced through somatic DNA recombination in developing T-cells, which fulfil key roles in maintaining a functional adaptive immune system. Recent advances in high-throughput DNA sequencing (HTS) have resulted in large numbers of TCR sequences being produced, with more data, library preparations and analysis techniques being published on a regular basis. 

Polymorphism in the human TCR loci can impact upon the expression levels and antigen binding properties of rearranged receptors, however such information is frequently either not considered or not detected in current analyses. Appreciation of effects of any variation beyond that described across the different TCR gene alleles contained in IMGT/GENE-DB is particularly rare.

This script seeks to produce files containing reference and alternative polymorphic sequences from across the alpha and beta chain loci, in a format convenient for feeding into downstream analyses. This should allow us to test the effect of these variants on other repertoire parameters, so as to better probe and profile high-throughput TCRseq repertoire data.

### Data sources

* Ensembl IDs were found by searching [Ensembl](https://ensembl.org) (only human) for the IMGT-designated gene symbol
    * Note that forward slashes ('/') are omitted in Ensembl and gnomAD gene symbols (e.g. TRAV36/DV7 becomes TRAV36DV7)
* Polymorphism types and coordinates were sourced from [gnomAD](http://gnomad.broadinstitute.org/), searching on Ensembl ID and downloading the csv of variant sites
* Refererence sequences were downloaded using the [UCSC](http://genome.ucsc.edu/) DAS browser (see [this gist](https://gist.github.com/JamieHeather/b03cc8a330a69c622c3e5ffbc8fb7550)); alternate sequences are inferred from these given the description of the variant.
* Polymorphisms that fall into known germline V and J genes were harvested from [IMGT/GENE-DB](https://www.imgt.org/genedb/) by downloading the V/J-REGION fasta sequence for all genes with multiple alleles, determining alternative positions relative to the prototypical allele for each gene (*01) and extracting the polymorphic regions.

Note that it seems that there remain TCR genes which lack both designated Ensembl IDs and gnomAD entries (or at least not ones that are readily found). This was vastly improved in the latest round of data collection (**v3.4 update**), but there remain a number of TRBV genes which lack Ensembl IDs/corresponding SNP data. Those are listed here - and note that all but TRBV3-2 are IMGT-designated functional genes:

* TRBV3-2 (pseudogene)
* TRBV4-3
* TRBV5-8
* TRBV6-3
* TRBV6-9
* TRBV7-8

There are also a number of TRBV genes which have been allocated Ensembl IDs since this project began, but which yet still lack gnomAD entries. Those are the highest numerical Ensembl IDs in this analysis, perhaps suggested that they are just yet to be processed. Again note that these are with one exception annotated as functional genes, and indeed contains a number of frequently used segments (including the entire TRBJ1 cluster):

* TRBV10-3 (ENSG00000275791)
* TRBV11-3 (ENSG00000276597)
* TRBV12-3 (ENSG00000274752)
* TRBV12-4 (ENSG00000276953)
* TRBV12-5 (ENSG00000275158)
* TRBV13 (ENSG00000276405)
* TRBV14 (ENSG00000275743)
* TRBV15 (ENSG00000276819)
* TRBV16 (ENSG00000275243)
* TRBV17 (ENSG00000277880) (ORF)
* TRBV18 (ENSG00000276557)
* TRBV25-1 (ENSG00000282499)
* TRBV6-2 (ENSG00000283063)
* TRBV7-2 (ENSG00000282939)
* TRBV7-9 (ENSG00000278030)
* TRBJ1-1 (ENSG00000282320)
* TRBJ1-2 (ENSG00000282420)
* TRBJ1-3 (ENSG00000282133)
* TRBJ1-4 (ENSG00000281958)
* TRBJ1-5 (ENSG00000282173)
* TRBJ1-6 (ENSG00000282780)

Those genes which lack Ensembl IDs and/or gnomAD entries therefore only have IMGT recorded variants in these data (where they exist).

Also note that while IMGT covers all of the genes, it only stores allelic information within the V and J gene sequences themselves, ignoring any neighbouring sequence which may be functionally relevant. Similarly sometimes the supposed reference sequence at a given position does not match the actual hg19 sequence exactly, and there are also ambiguous polymorphisms, thus not all described variants survive the process (particularly among duplications and deletions).

Should new data come available users may wish to redownload the raw data and re-run the code with it; bear in mind it takes an hour or two to run, due to the need to pull each sequence from UCSC.

### Running the script

This script was written for Python 2.7. Just get the repo from github and run the script, e.g.:

```
git clone https://github.com/JamieHeather/tcr-snps.git
python tcr-snps.py
```

Note this will overwrite any previously SNP files files in the output directory.

This script also requires the **Levenshtein package**, which can be installed via pip:
```
pip2 install python-levenshtein
```

### Changes in v3.4:

* Moved some files around around
    * TCR details csv moved to Raw Data
    * Moved output data to own directory
* Updated the input data, tweaked the output format
    * Replaced older [ExAC](http://exac.broadinstitute.org/) data with the latest release of gnomAD [version 2.0.2](https://storage.googleapis.com/gnomad-public/release/2.0.2/README.txt) (slightly altered filtering)
    * Added some additional TCR genes which were previously absent in the database
    * Added a (commented) heading line per output SNP file
* Updated README/code accordingly






