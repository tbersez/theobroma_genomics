# Duplicated genes in Theobroma cacao genome
*Master genomics projet, Paris Saclay University*

**Aim of the project** : find duplicated sequences across a plant genome. Case study : *Theobroma Cacao*, the cacao tree. For the comparison of nucleic acid sequences would be  time- and ressources-consumming, one should better stick with the information conveyed in the amino acid sequences. Also because some mutation do not really change the local properties of the protein.

### Isoform selection

Each gene display several transcripts. Because two transcripts of the same gene are really likely to be detected as a gene duplication, some isoform filtering has to be done. 
Here, we chose to keep the longest isoform, assuming that it countains the most information about the gene

### BLASTp run

BLASTp is run in order to align every selected peptide with every other. Such alignment can be quite long, even if the use of peptides sequences makes it faster.

### Clustering of BLAST results

Then, the FTAG Finder v3 pipeline on the Galaxy platform is used to cluster sequences into families based on the BLAST results


### Tools and packages

This project uses Python3, the BioPython and the pandas packages.
BLAST are run locally with command lines
