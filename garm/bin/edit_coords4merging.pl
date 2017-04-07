#!/usr/bin/env perl

use strict;
use lib "$ENV{GARMLIB}";
use GetData;

unless (@ARGV) {
	die "Usage: edit_coords4merging.pl <genome1.fasta> <genome2.fasta> <file.coords>\n";
}

my $file1 = "$ARGV[0]";
my $file2 = "$ARGV[1]";
my $coords = "$ARGV[2]";

my $list1 = `grep ">" $file1 | tr -d ">"`;
my $list2 = `grep ">" $file2 | tr -d ">"`;
#system ("cat $file1 $file2 > for_merging.fasta");

die "fasta files empty or wrong location...\n" unless ($list1 and $list2);

my @headers1 = split (/\n/, $list1);
my @headers2 = split (/\n/, $list2);
my $refcount = scalar @headers1;
#print $refcount; <STDIN>;

open ONE, ">$file1.index";
my %allheaders;
my $count = 1;
foreach my $name (@headers1) {
	chomp $name;
	$allheaders{$name} = $count;
	$count++;
	print ONE "$name\t$allheaders{$name}\n"; #<STDIN>;
}

open TWO, ">$file2.index";
foreach my $name (@headers2) {
	chomp $name;
	$refcount++;
	$allheaders{$name} = $refcount;
	print TWO "$name\t$allheaders{$name}\n";
}
close ONE,
close TWO,


open IN, "$coords";
open OUT, ">edited.coords";

while (<IN>) {
	chomp;
	my @line = split (/\t/, $_);
	if ($allheaders{$line[17]} and $allheaders{$line[18]}) {
		$line[17] = "$allheaders{$line[17]}";
		$line[18] = "$allheaders{$line[18]}";
		my $edited = join("\t", @line);
		#print $edited; <STDIN>;
		print OUT "$edited\n";
	}
#	elsif ($allheaders{$line[17]}) {
#		print "genome2 $line[18] was excluded?"; <STDIN>;
#	}
#	elsif ($allheaders{$line[18]}) {
#		print "genome1 $line[17] was excluded?"; <STDIN>;
#	} else {
#		print "genome1 $line[17] and genome2 $line[18] were excluded?"; <STDIN>;
#	}
}

close IN;
close OUT;
