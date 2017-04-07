#!/usr/bin/env perl

use strict;
use lib "$ENV{GARMLIB}";
use GetData;
my $mumbin = $ENV{MUMBIN};
my $amosbin = $ENV{AMOSBIN};

die "Usage: check_coords_annotated.pl <AGP file> <splimer.coords> <splimer.notAnnotated.coords> <param> <maxtrim> <conide> <ovlide> <genome1> <genome2>\n" unless (scalar @ARGV == 9);

my $AGP = shift;
my $annotated = shift;
my $notannotated = shift;
my $param = shift;
my $maxtrim = shift;
my $conide = shift;
my $ovlide = shift;
my $genome1 = shift;
my $genome2 = shift;

my %contigs = GetData::loadAGP($AGP);
my %coords = GetData::loadcoords($annotated);
my %notannot = GetData::loadcoords($notannotated);
my %refs = GetData::loadfasta($genome1);
my %qrys = GetData::loadfasta($genome2);

mkdir ("rerun");
open TRIM, ">trimmed_contigs.list";
open TAILS, ">contained_but_overhangs.list";
open RERUN, ">rerun_nucmer_ref.sh";
foreach my $contig (keys %contigs) {
	unless ($coords{$contig}) {
		if ($notannot{$contig}) {
			#print "$contig chucked out...\n"; <STDIN>;
			my $howmany = scalar @{$notannot{$contig}};
			#print "$howmany options for contig $contig\n"; <STDIN>;
			if ($howmany > 1) {
				foreach my $line (@{$notannot{$contig}}) {
					#print "check $line"; <STDIN>;
				}
			} else {
				processline(${$notannot{$contig}}[0]);
			}
		}
	}
}

sub processline {
####

my $line = shift;
#print "$line"; <STDIN>;
my $begin = (split /\s+/, $line)[0];
my $end = (split /\s+/, $line)[1];
my $left = $begin - 1;
my $length = (split /\s+/, $line)[11];
my $right = $length - $end;
my $tag = (split /\s+/, $line)[19];
my $qry = (split /\s+/, $line)[18];
my $ref = (split /\s+/, $line)[17];
my $pct1 = 1;
my $pct2 = 1;
if ($left and $left <= 200) {
	$pct1 = $left/$length;
}
if ($right and $left <= 200) {
	$pct2 = $right/$length;
}
if ($tag =~ /CONTAINED/) {
	if ($pct1 <= 0.25 and $pct2 <= 0.25) {
		open REF, ">rerun/$ref.fasta";
		print REF ">$ref\n";
		my $trimref = substr($refs{$ref}, $left, $length-$right);		
		print REF "$trimref\n";
		close REF;
		open QRY, ">rerun/$qry.fasta";
		print QRY ">$qry\n$qrys{$qry}\n";
		close QRY;
		print RERUN "$mumbin/nucmer --maxmatch --minmatch $param --mincluster $param --prefix=$qry.vs.$ref rerun/$ref.fasta rerun/$qry.fasta && $mumbin/show-coords -H -c -l -o -r -I $conide $qry.vs.$ref.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'CONTAIN|IDENTITY' > clipped.coords && $mumbin/show-coords -H -c -l -o -r -I $ovlide $qry.vs.$ref.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'BEGIN|END' >> clipped.coords\n";
		print TRIM "$ref $left $right\n";
	}
	elsif ($pct1 and $pct1 <= 0.25) {
		my $trimref = substr($refs{$ref}, $left);	
		open REF, ">rerun/$ref.fasta";
		print REF ">$ref\n";
		print REF "$trimref\n";
		close REF;
		open QRY, ">rerun/$qry.fasta";
		print QRY ">$qry\n$qrys{$qry}\n";
		close QRY;
		print RERUN "$mumbin/nucmer --maxmatch --minmatch $param --mincluster $param --prefix=$qry.vs.$ref rerun/$ref.fasta rerun/$qry.fasta && $mumbin/show-coords -H -c -l -o -r -I $conide $qry.vs.$ref.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'CONTAIN|IDENTITY' > clipped.coords && $mumbin/show-coords -H -c -l -o -r -I $ovlide $qry.vs.$ref.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'BEGIN|END' >> clipped.coords\n";
		print TRIM "$ref $left 0\n";
	} 
	elsif ($pct2 and $pct2 <= 0.25) {
		open REF, ">rerun/$ref.fasta";
		my $trimref = substr($refs{$ref}, 0, $length-$right);
		print REF ">$ref\n";
		print REF "$trimref\n";
		close REF;
		open QRY, ">rerun/$qry.fasta";
		print QRY ">$qry\n$qrys{$qry}\n";
		close QRY;
		print RERUN "$mumbin/nucmer --maxmatch --minmatch $param --mincluster $param --prefix=$qry.vs.$ref rerun/$ref.fasta rerun/$qry.fasta && $mumbin/show-coords -H -c -l -o -r -I $conide $qry.vs.$ref.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'CONTAIN|IDENTITY' > clipped.coords && $mumbin/show-coords -H -c -l -o -r -I $ovlide $qry.vs.$ref.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'BEGIN|END' >> clipped.coords\n";
		print TRIM "$ref 0 $right\n";
	} else {
		print TAILS "$ref\tOVERHANG > 200 bases or more than 25% of the sequence\n";
	}
} else {
	#print "check tag\n$line\n"; <STDIN>;
}

####
}
