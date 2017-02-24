<h1> tcr-snps</h1>

<h3>Extract reference and variant sequence for SNPs and other polymorphisms in the human T-cell receptor loci</h3>

<h4>version 3.2, Jamie Heather, January 2017</h4>
<hr>

T-cell receptors (TCRs) are variable antigen receptors produced through somatic DNA recombination in developing T-cells, which fulfil key roles in maintaining a functional adaptive immune system. Recent advances in high-throughput DNA sequencing (HTS) have resulted in large amounts of TCR sequence data produced. 

Polymorphism in the human TCR loci can impact upon the expression levels and antigen binding properties of rearranged receptors, however such information is frequently either not considered or not detected in current analyses, especially any variation beyond that described in the different TCR gene alleles contained in IMGT/GENE-DB.

This script uses SNP (and other polymorphism) data extracted from ExAC/gnomAD and the UCSC DAS browser to generate reference and polymorphism DNA sequences in order to better probe and profile HTS TCR repertoire data.

Note that it currently seems that not all TCR genes have entries in these databases (or even seemingly Ensembl gene IDs, or at least not ones that are readily found), and so are missing from this current analysis. 

Similarly sometimes the supposed reference sequence at a given position does not match the actual hg19 exactly, and there are also ambiguous polymorphisms, thus not all described variants surive the process (particularly among duplications and deletions).
 
<h4>New in v3.2:</h4>
* Added IMGT-recorded TCR SNPs
* Split functions out into separate document for readability

Genes containing polymorphisms were taken from [IMGT/GENE-DB](http://www.imgt.org/genedb/) and outputting the V- or J-REGION for those genes with more than 1 allele. All output files were combined into the file 'Raw_Files/imgt_polymorphic_genes.fasta'. Polymorphisms are inferred using the [Levenshtein](https://pypi.python.org/pypi/python-Levenshtein/0.12.0) Python module.







