#!/usr/bin/env perl

use strict;
use lib "$ENV{GARMLIB}";
use GetData;

my $scaffs = $ARGV[0];

open SCAFF, ">$scaffs.contigs_in_scaffs.fasta" || die "Can't open $scaffs.contigs_in_scaffs.fasta\n: $!\n";;
open SINGLE, ">$scaffs.contigs_singletons.fasta" || die "Can't open $scaffs.contigs_singletons.fasta\n: $!\n";;
open AGP, ">$scaffs.agp";

my %seqs = GetData::loadfasta($scaffs);

foreach my $key (sort {(split /\_/, $a)[$#_] <=> (split /\_/, $b)[$#_]} keys %seqs) {
	my $conta = 1;
	my $agpcont = 1;
	my $pos = 1;
	if ($seqs{$key} =~ /N{6,}/) {
		my @gaps = $seqs{$key} =~ m/N{6,}/g;
		my @contigs = split (/N{6,}/, $seqs{$key});
		foreach my $contig (@contigs) {
			my $l = length $contig;
			my $end = $pos+$l-1;
			$key = (split /\s+/, $key)[0]; ### WILL TAKE ONLY THE FIRST FIELD
			print SCAFF ">$key\_$conta\n$contig\n"; #<STDIN>;
			print AGP "$key\t$pos\t$end\t$agpcont\tW\t$key\_$conta\t1\t$l\t+\n";
			$agpcont++;
			$pos += $l;
			my $gapend = 0;
			if ($conta-1 < $#contigs) {
				my $lgap = length $gaps[$conta-1];
				$gapend = $pos+$lgap-1;
				print AGP "$key\t$pos\t$gapend\t$agpcont\tN\t$lgap\tfragment\tyes\n"; #<STDIN>;
				$pos = $gapend+1;
			$conta++;
			$agpcont++;
			} 
		}
	} else {
		my $seq = $seqs{$key};
		$key = (split /\s+/, $key)[0]; ### WILL TAKE ONLY THE FIRST FIELD
		print SINGLE ">$key\n$seq\n"; #<STDIN>;
	}
}
close AGP;
close SCAFF;
close SINGLE;
system("cat $scaffs.contigs_in_scaffs.fasta $scaffs.contigs_singletons.fasta > $scaffs.all_split.fasta");
