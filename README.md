# Duplicated genes in Theobroma cacao genome
*Master genomics projet, Paris Saclay University*

**Aim of the project** : find duplicated sequences across a plant genome. Case study : *Theobroma Cacao*, the cacao tree. For the comparison of nucleic acid sequences would be  time- and ressources-consumming, one should better stick with the information conveyed in the amino acid sequences. Also because silencious mutations do not change the local properties of proteins.

## Project summary

### Isoform selection

Each gene display several transcripts. Because two transcripts of the same gene are likely to be detected as a gene duplication, isoform filtering has to be done. 
Here, we chose to keep the longest isoform, assuming that it countains the most information about the gene.
This filtering also allow to significatively reduce the size BLASTp instance.

### BLASTp run

BLASTp is run in order to align every selected peptide with every other. Such alignment can be quite long, even if the use of peptides sequences makes it faster. A tab like file is outputed by the BLASTp. That one is used to cluster protein into duplication families.

### Clustering of BLAST results

The FTAG Finder v3 pipeline on the Galaxy platform is used to cluster sequences into families based on the BLAST results


### Tools and packages

This project uses Python3, BioPython and pandas packages.
BLAST are run locally with command lines.
