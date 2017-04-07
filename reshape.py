#!/usr/bin/env python

import sys, subprocess, os, shutil
from datetime import datetime

#Checking arguments

if len(sys.argv) < 5 or len(sys.argv) > 6:
	print("Usage: reshape.py <MiSeq 1 FASTQ> <MiSeq 2 FASTQ> <PacBio FASTQ> <genome_size> <Preassembly FASTA>")
	print("Or: reshape.py <MiSeq 1 FASTQ> <MiSeq 2 FASTQ> <PacBio FASTQ> <genome_size>")
	exit(-1)

#Creating output directory
output_dir = "results_" + datetime.now().strftime('%d%m%Y_%H%M%S')
subprocess.call(["mkdir", "-p", output_dir])

#Executing a5-miseq (1st stage)
subprocess.call(["./a5-miseq/bin/a5_pipeline.pl", sys.argv[1], sys.argv[2], "assembly.out"])
subprocess.call(["mkdir", "-p", output_dir + "/a5-miseq"])
files = os.listdir(os.getcwd())
for f in files:
	if f.startswith("assembly.out."):
		shutil.move(f, output_dir + "/a5-miseq")
print("---------------------------------------------------------------------------------------------------------------------------------------------------------------")

#Executing GARM (2nd stage - OPTIONAL)
if len(sys.argv) == 6:
	genomes_file = open('genomes.txt', 'w')
	genomes_file.write(os.path.abspath(output_dir + "/a5-miseq/assembly.out.contigs.fasta A5") + '\n')
	genomes_file.write(os.path.abspath(sys.argv[5]) + " INPUT")
	genomes_file.close()
	subprocess.call(["./garm/GARM.pl", "-g", "genomes.txt", "-o", "garm-temp"])
	open('genomes.txt', 'w').close()
	shutil.move("GARM_garm-temp.log", "garm-temp")
	shutil.move("garm-temp", output_dir)

print("---------------------------------------------------------------------------------------------------------------------------------------------------------------")

#Executing Canu (3rd stage)
subprocess.call(["./canu/Linux-amd64/bin/canu", "-p", "asm", "-d", output_dir + "/canu", "genomeSize=" + sys.argv[4], "-pacbio-raw", sys.argv[3]])

print("---------------------------------------------------------------------------------------------------------------------------------------------------------------")

#Executing GARM (4th stage)
genomes_file = open('genomes.txt', 'w')
if len(sys.argv) == 5:
	genomes_file.write(os.path.abspath(output_dir) + "/a5-miseq/assembly.out.contigs.fasta A5" + '\n')
else:
	genomes_file.write(os.path.abspath(output_dir) + "/garm-temp/final_merged_contigs.fasta GARM" + '\n')
genomes_file.write(os.path.abspath(output_dir) + "/canu/asm.contigs.fasta CANU")
genomes_file.close()
subprocess.call(["./garm/GARM.pl", "-g", "genomes.txt", "-o", "garm-final"])
subprocess.call(["rm", "genomes.txt"])
shutil.move("GARM_garm-final.log", "garm-final")
shutil.move("garm-final", output_dir)

print("---------------------------------------------------------------------------------------------------------------------------------------------------------------")

#Storing output name files for each stage
a5_output = output_dir + "/a5-miseq/assembly.out.contigs.fasta"
garm_temp_output = output_dir + "/garm-temp/final_merged_contigs.fasta"
canu_output = output_dir + "/canu/asm.contigs.fasta"
garm_final_output = output_dir + "/garm-final/final_merged_contigs.fasta"

#Executing QUAST
if len(sys.argv) == 6:
	subprocess.call(["./quast/quast.py", a5_output, garm_temp_output, canu_output, garm_final_output, "-o", output_dir + "/quast", "--labels" , "a5-miseq,garm-temp,canu,garm-final", "--gene-finding"])
else:
	subprocess.call(["./quast/quast.py", a5_output, canu_output, garm_final_output, "-o", output_dir + "/quast", "--labels" , "a5-miseq,canu,garm-final", "--gene-finding"])
print("---------------------------------------------------------------------------------------------------------------------------------------------------------------")
print("Final assembly is in: " + garm_final_output)
print("Assembly statistics in each stage is in: " + output_dir + "/quast/report.txt")
