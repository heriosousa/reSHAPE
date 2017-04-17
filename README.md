# reSHAPE

reSHAPE (Sousa's Hybrid Assembly Pipeline) is an assembly pipeline which combines different genome assembly tools using different sequencing technologies.

## Copyright and license

reSHAPE is free, open source software. Source code from third-part software used in this pipeline can be found at:

- A5-miseq: https://sourceforge.net/projects/ngopt/
- GARM: https://sourceforge.net/projects/garm-meta-assem/
- Canu: https://github.com/marbl/canu/
- QUAST: https://sourceforge.net/projects/quast/

Current used versions:

- A5-miseq (2016-08-25 version)
- GARM (version 0.7.5)
- Canu (version 1.3)
- QUAST (version 4.5)

## Requirements

- 64-bit Linux
- Perl v5.8 or above
- Perl Modules: 
	Parallel::ForkManager
	List::MoreUtils
- Python 2.7+

## Installation 

- After having Perl and required modules installed you can install each third-part software one-by-one or run the following command:

source install.sh

## Usage examples

reshape.py <MiSeq 1 FASTQ> <MiSeq 2 FASTQ> <PacBio FASTQ> <genome_size> <Preassembly FASTA>

reshape.py <MiSeq 1 FASTQ> <MiSeq 2 FASTQ> <PacBio FASTQ> <genome_size>

Where <genome_size> should be provided according to http://canu.readthedocs.io/en/latest/parameter-reference.html#parameter-reference

## Quick start

Go to samples directory and decompress both MiSeq files. Right after, go back to reSHAPE base directory and execute:

./reshape.py samples/miseq1.fastq samples/miseq2.fastq samples/pacbio.fastq 2m samples/asm.fasta

Final assembly will be located on: results_datetime/garm-final/final_merged_contigs.fasta

You can also find assembly statistics of each stage on: results_datetime/quast/report.txt
