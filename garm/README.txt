GARM README

INTRODUCTION
GARM (Genome Assembly, Reconciliation and Merging) is a pipeline of several scripts to combine results from different genome assemblies. The results (contigs and/or scaffolds) can come from different sequencing technologies (capillary, 454, Illumina, SOLiD, etc) libraries or just different programs using the same sequencing technology. 

REQUIREMENTS:

Since version 0.7.3, the necessary third-party software* is already included. Please, configure your enviromental variables (described below) pointing to the locations included in the package. If you try a different version and works, please let us know!. 

-Perl v5.8 or above
-Perl Modules: 
	Parallel::ForkManager
	List::MoreUtils
The best way to install a perl module would be CPAN
-AMOS (tested on v3.0.0) software (http://sourceforge.net/apps/mediawiki/amos/index.php?title=AMOS)
-MUMmer (tested on v3.22) software (nucmer) (http://mummer.sourceforge.net/)
-R language (for GC plots only)

Please, make sure that everything is installed and running properly.

Depending on which shell you have, please add these environmental variables:

If you have tcsh then in your .cshrc add:

setenv GARMBIN /path/where/you/installed/GARM_version/bin
setenv GARMLIB /path/where/you/installed/GARM_version/lib
setenv MUMBIN  /path/where/you/installed/GARM_version/MUMmer3.22
setenv AMOSBIN /path/where/you/installed/GARM_version/amos-3.0.0/bin
setenv AMOSLIB /path/where/you/installed/GARM_version/amos-3.0.0/lib
setenv PATH {$PATH}:/path/where/you/installed/GARM_version	# add GARM.pl to PATH

In bash shell, could be .bashrc or .bash_profile add:


GARMBIN=/path/where/you/installed/GARM_version/bin
GARMLIB=/path/where/you/installed/GARM_version/lib
MUMBIN=/path/where/you/installed/GARM_version/MUMmer3.22
AMOSBIN=/path/where/you/installed/GARM_version/amos-3.0.0/bin
AMOSLIB=/path/where/you/installed/GARM_version/amos-3.0.0/lib
PATH=$PATH:/path/where/you/installed/GARM_version	# add GARM.pl to PATH

export $GARMBIN:$GARMLIB:$MUMBIN:$AMOSBIN:$AMOSLIB

Once you defined the environmental variables you can run the script "/path/where/you/installed/GARM_version/config.pl" to chech if everything it's okay.

INPUT:
- 2 fasta files with contigs or scaffolds from your assemblies.
- A text file (info_genomes_file) with the full path to those genome assembly results. The format is a 2 column file:
	First column: Full path file location of the assembly result
	Second column: A tag for the sequencing technology (Sanger, 454, Illumina, SOLiD, etc).
- For scaffolds, the program will take the number of N's between contigs as the estimated gap size (10 Ns is the minimum). The program will automatically select which result is better and will select it as a reference. In case you want to specifically use one scaffolds, write in the  info_genomes_file a third column with a "*" symbol next to the file you want to use.  

GARM USAGE:

Usage: GARM.pl -g info_genomes_file -o <prefix>
Options:
   -g: Text file with the path to two fasta genome assemblies to merge
   -o: Prefix for output directory and files
   -h: Complete help menu with advanced options

ADVANCE OPTIONS (Do not change unless you know what you are doing):
   -c <n>: Run with multi-threads with n cpus (default n=2)

EXAMPLE (HOW TO RUN)

- The info_genomes_file could have any name you like (eg. genomes.txt or config.file, etc). It should look like this:

/home/assembly_res/ABYSS/k31/kmer31-contigs.fasta ABYSS
/home/assembly_res/Newbler/454Scaffolds.fna 454

- The second column could be any tag that you will like to use to distinguish your assemblies (eg. ILLUM, SOLiD, 454, SOAP, Velvet, etc.). Just be sure that they are not the same. This tag would be used in case that you have contig names that are the same in the 2 results.

- Also check that your fasta files are in the location you wrote and that they are in a fasta format with no empty lines between sequences. 

- After setting up properly all the requirements and writing your info_genomes_file you can run GARM like this:

/path/where/GARM/is/installed/GARM.pl -g genomes.txt -o my_genome

OUTPUT:

- Inside the directory that you specified as output, you will find:

final_merged_contigs.fasta - This are the things in common to both results. Fasta file with your merged contigs.
final_merged.read.placed.filtered - Position of each original contig in the merged contigs. The file is in the WashU read placed format.
final_bin.fasta - Contigs with no overlaps
	            - Contigs excluded by the minimum length filter (200 bases default)
	            - Contigs with overlaps but not merged for some reason... (AMOS problem)
final.stats - genome stats for the results


  Column   Description
  ------   ----------------------------------------------------------
    1      NCBI ti number for read (or *, if none known)
    2      read name
    3      left trimmed position on the original read
    4      number of bases in trimmed read
    5      orientation on contig (0 = forward, 1 = reverse)
    6      contig name
    7      supercontig name (or *, if none known)
    8      approximate start position of the trimmed read in the contig
    9      approximate start position of the trimmed read on supercontig (or *, if none known)

- If there were scaffolds:

final_merged_scaffolds.fasta - Fasta file with the reconstructed scaffolds
final_merged_scaffolds.agp - AGP file for the scaffolding information

Any questions:

Alejandro Sanchez alexsf[at]ibt.unam.mx
Karel Estrada karel[at]ibt.unam.mx




