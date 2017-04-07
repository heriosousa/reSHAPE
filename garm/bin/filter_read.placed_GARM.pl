#!/usr/bin/env perl

use strict;

die "Usage: filter_read.placed_GARM.pl <read.placed>\n" unless @ARGV;

my $rpf = $ARGV[0];
my %scaffs;
my %contigs;

open IN, "$rpf";
open OUT, ">$rpf.filtered";
open LOG, ">filter_read.placed.log";

while(<IN>) {
	my @line = split (/\s+/);
	my $pos = $line[7];
	my $length = $line[3] - $line[2];
	push (@{$scaffs{$line[5]}{$pos}}, $_);
	push (@{$contigs{$line[5]}{$pos}}, "$_");
}

foreach my $scf (sort { (split /contig/, $a)[1] <=> (split /contig/, $b)[1] } keys %scaffs) {
	#print "I'm in $scf"; <STDIN>;
	foreach my $pos (sort {$a <=> $b} keys %{$scaffs{$scf}}) {
		my $count = scalar (@{$scaffs{$scf}{$pos}});
		#print $count; <STDIN>;
		if ($count > 1 and $contigs{$scf}{$pos}) {
			my %length;
			my $index = 0;
			foreach my $line (@{$contigs{$scf}{$pos}}) {
				my $l = (split /\s+/, $line)[3] - (split /\s+/, $line)[2];
				push (@{$length{$index}{$l}}, $line);
			}
			foreach my $index (sort {$a <=> $b} keys %length) {
				foreach my $l (keys %{$length{$index}}) {
					my $c = scalar @{$length{$index}{$l}};
					if ($c > 1) {
						foreach my $hap (@{$length{$index}{$l}}) {
							chomp $hap;
							print LOG "$hap\tPROBABLE HAPLOTYPE\n";
						}
					} else {
						my $newline = ${$length{$index}{$l}}[0];
						print OUT "$newline";
					}
				}
			}
		} else {
			my $newline = ${$scaffs{$scf}{$pos}}[0];
			print OUT "$newline"; #<STDIN>;
		}
	}
}

close OUT;
close LOG;
