#!/usr/bin/env perl

use strict;
use lib "$ENV{GARMLIB}";
use GetData;

unless (scalar @ARGV == 4) {
	die "Usage: generate_fasta_scf.pl <merged_contigs.fasta> <singletons.fasta> <extras.fasta> <AGP file>\n";
}

my $mergedfile = shift;
my $singlets = shift;
my $extras = shift;
my $agp = shift;
#my $rpf = shift;

my %mergedfasta = GetData::loadfasta($mergedfile);
my %singlefasta = GetData::loadfasta($singlets);
my %extrasfasta = GetData::loadfasta($extras);

open AGP, "$agp";
my %newscaff;
while (<AGP>) {
	next if (/^\#/);
	my @line = split (/\s+/,$_);
	if (/fragment/) {
		if ($line[5] >= 0) {
			my $gap = "N" x $line[5];
			$newscaff{$line[0]} .= $gap;
			#print $gap; <STDIN>;
		} else {
			my $gap = "N" x 10;
			$newscaff{$line[0]} .= $gap;
		}
	} else {
		my $strand = $line[8];
		my $seq = "";
		if ($strand eq "-") {
			#print "minus"; <STDIN>;
			if ($mergedfasta{$line[5]}) {
				my $revcomp = GetData::revcomp($mergedfasta{$line[5]});
				$seq = $revcomp;
			} elsif ($singlefasta{$line[5]}) {
				my $revcomp = GetData::revcomp($singlefasta{$line[5]});
				$seq = $revcomp;
				delete $singlefasta{$line[5]};
			} elsif ($extrasfasta{$line[5]}) {
				my $revcomp = GetData::revcomp($extrasfasta{$line[5]});
				$seq = $revcomp;
				delete $extrasfasta{$line[5]};
			} else {
				print "$line[5]\tI can't find this sequence... ARGH!!!!\n";
			}
		} else {
			if ($mergedfasta{$line[5]}) {
				$seq = $mergedfasta{$line[5]};
			} elsif ($singlefasta{$line[5]}) {
				$seq = $singlefasta{$line[5]};
				delete $singlefasta{$line[5]};
			} elsif ($extrasfasta{$line[5]}) {
				$seq = $extrasfasta{$line[5]};
				delete $extrasfasta{$line[5]};
			} else {
				print "$line[5]\tI can't find this sequence... ARGH!!!!\n";
			}
		}
		$newscaff{$line[0]} .= $seq;
		#print "$seq"; <STDIN>;
	}
}
close AGP;


#open RPF, "$rpf";
#my %merged;
#while (<RPF>) {
#	my $scaff = (split /\s+/, $_)[1];
	#print $scaff; <STDIN>;
#	my @names = split (/\_/, $scaff);
#	if (scalar @names == 3) {
#		my $header = $names[1];
#		$merged{$header} = (split /\s+/, $_)[5];
#	}
#}
#close RPF;

open OUT, ">$agp.fasta";
foreach my $scf (sort { (split /\_/, $a)[$#_] <=> (split /\_/, $b)[$#_] } keys %newscaff) {
	#print $scf; <STDIN>;
	#my $header = (split /\_/, $scf)[$#_];
	#print "$header"; <STDIN>;
	#if ($merged{$header} and ($merged{$header} =~ /merged/) and $newscaff{$scf}) {
	print OUT ">$scf\n$newscaff{$scf}\n";
	#} #else {
#		print "$header doesn't exist in RPF or just point to individual contigs or empty\n";
#	}
	#print OUT ">$scf\n$newscaff{$scf}\n";
}

open BIN, ">final_leftovers_bin.fasta";
foreach my $leftover (keys %singlefasta) {
	if ($singlefasta{$leftover}) {
		print BIN ">$leftover\n$singlefasta{$leftover}\n";
	}
}
foreach my $leftover (keys %extrasfasta) {
	if ($extrasfasta{$leftover}) {
		my @fields = split (/\_/, $leftover);
		if (scalar @fields == 2) {
			print BIN ">$leftover\n$extrasfasta{$leftover}\n";
		} 
		elsif (scalar @fields > 3) {
			print BIN ">$leftover\n$extrasfasta{$leftover}\n";
		}
	}
}


close BIN,
close OUT;

