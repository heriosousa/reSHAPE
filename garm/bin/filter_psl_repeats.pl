#!/usr/bin/env perl 

use strict;

unless (scalar @ARGV == 2) {
    die "Usage: filter_psl_repeats.pl <prefix.coords> <minlen>\n"; 
}

my $psl = $ARGV[0];
my $minlen = $ARGV[1];
my %end;
my %begin;

open IN, "$psl" || die "Can't open $psl\n";
while (<IN>) {
	if (/^\d+/) {
		my @line = (split);
		my $qry = $line[9];
		my $ref = $line[13];
		my $rlen = $line[14];
		my $rstart = $line[15];
		my $rend = $line[16];
		if (($rend == $rlen) && ($rend-$rstart >= $minlen)) {
			$end{$ref} = $qry;
		}
		elsif (($rstart <= 20) && ($rend-$rstart >= $minlen)) {
			$begin{$ref} = $qry;
		}
	}
}
close IN;

open REPEATS, ">$ARGV[0].repeats" || die "Can't create or write to $ARGV[0].repeats\n";

foreach my $contig (sort keys %begin) {
	if ($end{$contig}) {
		print REPEATS "$contig\tBOTH\t$begin{$contig}\n";
	}
	else {
		print REPEATS "$contig\tBEGIN\t$begin{$contig}\n";
	}
}

foreach my $contig (sort keys %end) {
	next if $begin{$contig};
	print REPEATS "$contig\tEND\t$end{$contig}\n";
}

close REPEATS;


