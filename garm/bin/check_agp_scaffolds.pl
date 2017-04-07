#!/usr/bin/env perl

use strict;
use lib  "$ENV{GARMLIB}";
use GetData;

die "Usage: check_agp_scaffolds.pl <AGP_file>\n" unless scalar @ARGV == 1;

my $agp = shift;
open AGP, "$agp";
my %contigs;
my $previous = "";
while (<AGP>) {
	next if (/^\#/);
	unless (/fragment/) {
		my @line = split (/\s+/,$_);
		if ($contigs{$line[5]}) {
			my $index = (scalar @{$contigs{$line[5]}}) - 1;
			if ($line[0] eq ${$contigs{$line[5]}}[$index]) {
				#print;
				print "$previous"; #<STDIN>;
				$previous = $_;
			} else {	
				push (@{$contigs{$line[5]}}, $line[0]); 
			}
		} else {
			push (@{$contigs{$line[5]}}, $line[0]);
			$previous = $_;
		}
	}
}
foreach my $contig (keys %contigs) {
	my $total = scalar @{$contigs{$contig}};
	if ($total == 2) {
		print "join ${$contigs{$contig}}[0] and ${$contigs{$contig}}[1] by $contig\n";
	} 
}
