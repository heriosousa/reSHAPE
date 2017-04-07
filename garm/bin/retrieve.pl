#!/usr/bin/env perl
use strict;

die "Usage: retrieve.pl PWD PREFIX\n" unless (@ARGV);

my $pwd = shift;
my $prefix = shift;

my $path = "$pwd/$prefix";

if (-s "$path/scaffolding/new_scaffolds.txt.fasta.agp.fixed.fasta.all_split.fasta") {
	system ("cp $path/scaffolding/new_scaffolds.txt.fasta.agp.fixed.fasta.all_split.fasta $pwd/$prefix/final_merged_with_reovl_contigs.fasta");
	system ("cp $path/scaffolding/new_scaffolds.txt.fasta.agp.fixed.fasta $pwd/$prefix/final_merged_with_reovl_scaffolds.fasta");
	system ("cp $path/scaffolding/new_scaffolds.txt.fasta.agp.fixed.fasta.agp $pwd/$prefix/final_merged_with_reovl_scaffolds.agp");
	system ("cp $path/scaffolding/$prefix.read.placed $pwd/$prefix/final_merged.read.placed.filtered");
} 
elsif (-s "$path/scaffolding/new_scaffolds.txt.fasta.all_split.fasta") {
	system ("cp $path/scaffolding/new_scaffolds.txt.fasta.all_split.fasta $pwd/$prefix/final_merged_contigs.fasta");
	system ("cp $path/scaffolding/new_scaffolds.txt.fasta $pwd/$prefix/final_merged_scaffolds.fasta");
	system ("cp $path/scaffolding/new_scaffolds.txt.fasta.agp $pwd/$prefix/final_merged_scaffolds.agp");
	system ("cp $path/scaffolding/$prefix.read.placed $pwd/$prefix/final_merged.read.placed.filtered");
}
elsif (-s "$path/merging/$prefix.merged.fasta") {
	system ("cp $path/merging/$prefix.merged.fasta $pwd/$prefix/final_merged_contigs.fasta");
	system ("cp $path/merging/$prefix.read.placed $pwd/$prefix/final_merged.read.placed");
}
else {
	die "Not found: $path/new_scaffolds.txt.fasta.agp.fixed.fasta.all_split.fasta or $path/new_scaffolds.txt.fasta.all_split.fasta\n";
}
