Example for GARM

Files included:
- 454 assembly done by Newbler, using a 3kb and a 20kb PE library (~30x)
Scaffolds are used in this example

- Illumina assembly done by Velvet, using PE 76bp PCR free Illumina reads (~42x)
Kmer=31, insert size 150, coverage cutoff 10 expected coverage 42, min_pair_count 80
Contigs are used in this example.

You can download the reference (all chromosomes) from the Sanger Institute FTP site:
ftp://ftp.sanger.ac.uk/pub/pathogens/P_berghei/June_2011/

1.- Create a text file (eg. genomes.txt) with the location of the files provided in this example
in one column and in a second column assign a "tag" to each file. For example:

/path/454Scaffolds.fna 	454
/path/contigs.fa	ILLUM

2.- After setting up GARM as described in the README.txt inside the GARM tarball,
run GARM like this:

GARM.pl -g genomes.txt -o my_merged_assembly

3.- The program should start running, and a directory named with the selected prefix
will be created (my_merged_assembly) and a log file.

4.- in the directory you will find your results:
*final_merged_contigs.fasta.- Fasta file with your merged contigs
*final_merged.read.placed .- Location of your original reads in the merged contigs
   This file has the WashU read.placed file format, (check the README.txt)
*final_bin.fasta:
	Contigs with no overlaps
	Contigs excluded by the minimum length filter (200 bases default)
	Contigs with overlaps but not merged for some reason... (AMOS problem)
*final.stats .- genome stats for the results

If you have scaffolds (like here) you also get:

*final_merged_scaffolds.fasta .- Fasta file with the reconstructed scaffolds
*final_merged_scaffolds.agp .- AGP file for the scaffolding information

5.- The results can be compared using ABACAS (http://abacas.sourceforge.net/)
using the reference. (Some snapshots will be available at the GARM website soon!).
http://garm-meta-assem.sourceforge.net/
