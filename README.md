<h1> tcr-snps</h1>

<h3>Extract reference and variant sequences for SNPs and other polymorphisms in the human T-cell receptor loci</h3>

<h4>version 3.3, [Jamie Heather](http://jamieheather.github.io), January 2017</h4>
<hr>

T-cell receptors (TCRs) are variable antigen receptors produced through somatic DNA recombination in developing T-cells, which fulfil key roles in maintaining a functional adaptive immune system. Recent advances in high-throughput DNA sequencing (HTS) have resulted in large numbers of TCR sequences being produced, with more data, library preparations and analysis techniques being published on a regular basis. 

Polymorphism in the human TCR loci can impact upon the expression levels and antigen binding properties of rearranged receptors, however such information is frequently either not considered or not detected in current analyses. Appreciation of effects of any variation beyond that described across the different TCR gene alleles contained in IMGT/GENE-DB is particularly rare.

This script seeks to produce files containing reference and alternative polymorphic sequences from across the alpha and beta chain loci, in a format convenient for feeding into downstream analyses. This should allow us to test the effect of these variants on other repertoire parameters, so as to better probe and profile high-throughput TCRseq repertoire data.

SNPs (and other polymorphisms) are sourced using the following two methods:

* Polymorphism types and coordinates were sourced from [ExAC](http://exac.broadinstitute.org/)/[gnomAD](http://gnomad.broadinstitute.org/); the refererence sequence was downloaded using the [UCSC](http://genome.ucsc.edu/) DAS browser and the alternate sequence inferred.
* Polymorphisms that fall into known germline V and J genes were harvested from IMGT/GENE-DB by downloading the V/J-REGION fasta sequence for all genes with multiple alleles, determining alternative positions relative to the prototypical allele for each gene (*01) and extracting the polymorphic regions. 

Note that it seems that not all TCR genes currently have entries in ExAC/gnomAD (or even seemingly Ensembl gene IDs, or at least not ones that are readily found), and so only have IMGT recorded variants (where they exist). Conversely while IMGT covers all of the genes, it only stores allelic information within the V and J gene sequences themselves. 

Similarly sometimes the supposed reference sequence at a given position does not match the actual hg19 exactly, and there are also ambiguous polymorphisms, thus not all described variants surive the process (particularly among duplications and deletions).

Should new data come available users may wish to re-run the code with it - bear in mind it takes an hour or two to run, due to the need to pull each sequence from UCSC.

<h4>New in v3.3:</h4>

* Added IMGT-recorded TCR SNPs
* Split functions out into separate document for readability

Genes containing polymorphisms were taken from [IMGT/GENE-DB](http://www.imgt.org/genedb/) and outputting the V- or J-REGION for those genes with more than 1 allele. All output files were combined into the file 'Raw_Files/imgt_polymorphic_genes.fasta'. Polymorphisms are inferred using the [Levenshtein](https://pypi.python.org/pypi/python-Levenshtein/0.12.0) Python module.







