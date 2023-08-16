
The pipelines takes as an input pair-end raw reads from illumina Mi-seq, and a sample barcode file. It outputs a folder with the R1 and R2 files for each sample, the R1 and R2 are sorted by the primers of M1S and Z.

Input files:

1. The read 
2. The primers sequences available online at the table below.
3. The sample barcode file

Output files:

* folder with the R1 and R2 files for each sample


Pipeline container:

* Docker: immcantation/suite:4.3.0


Sequence processing steps:

1. generate barcode fasta file
2. MaskPrimer the R1 and R1 with the score method and position 0
5. MaskPrimer the R1 and R1 with align
6. Combine the output files and create the samples files. 
