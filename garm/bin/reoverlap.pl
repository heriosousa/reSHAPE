#!/usr/bin/env perl

use strict;
use lib "$ENV{GARMLIB}";
use GetData;

unless (scalar @ARGV == 2) {
	die "Usage: reoverlap.pl <merged_contigs.fasta> <AGP file>\n";
}
my $bin = "~as9/bin";
my $mumbin = "$ENV{MUMBIN}";
my $amosbin = "$ENV{AMOSBIN}";
my $garmbin = "$ENV{GARMBIN}";

mkdir ("reoverlap");
chdir ("reoverlap");

my $merged = shift;
my $agp = shift;

system("ln -s ../$merged");
system("ln -s ../$agp");

my $minovl = 50;
my $ignore = 10;
my $ide = 98;

my %mergedfasta = GetData::loadfasta($merged);

open AGP, "$agp";

my %newscaff;
my %order;
while (<AGP>) {
	my @line = split (/\s+/,$_);
	push (@{$order{$line[0]}}, $line[5]);
	unless (/fragment/) {
		my $strand = $line[8];
		my $seq = "";
		if ($strand eq "-") {
			if ($mergedfasta{$line[5]}) {
				my $revcomp = GetData::revcomp($mergedfasta{$line[5]});
				$seq = $revcomp;
					
			} else {
				print "$line[5] not found in merged sequences...\n";
			}
		} else {
			if ($mergedfasta{$line[5]}) {
				$seq = $mergedfasta{$line[5]};
			} else {
				print "$line[5] not found in merged sequences...\n";
			}
		}
		$newscaff{$line[5]} = $seq;
	}
}
close AGP;

my %reovl;
my %pairs;
open REOVL, ">reoverlap.fasta";
foreach my $scaff (sort {(split /\_/, $a)[$#_] <=> (split /\_/, $b)[$#_] } keys %order) {
	my $total = (scalar @{$order{$scaff}});
	for (my $index = 1; $index < $total-2; $index += 2) {
		if (${$order{$scaff}}[$index] <= 10) {
			my $ref = ${$order{$scaff}}[$index-1];
			my $qry = ${$order{$scaff}}[$index+1];	
			$reovl{$ref}++;
			$reovl{$qry}++;
			$pairs{$ref}{$qry} = ${$order{$scaff}}[$index];
		}
	}
}

my $count = 1;
my %newindex;
foreach my $seq (keys %reovl) {
	if ($newscaff{$seq}) {
		print REOVL ">$seq\n$newscaff{$seq}\n";
		$newindex{$seq} = $count;
		$count++;
	} else {
		print "$seq not found in hash newscaffs\n"; #<STDIN>;
	}
}
if (-s "reoverlap.fasta") {
	system ("$mumbin/nucmer --maxmatch -minmatch $minovl --prefix=reoverlap reoverlap.fasta reoverlap.fasta && $mumbin/show-coords -H -c -l -o -r -I $ide reoverlap.delta | $amosbin/nucmerAnnotate -noid -ignore $ignore | egrep 'END' > reoverlap.coords && $garmbin/filter_GARM_coords_v0.3.pl reoverlap.coords && mv filtered.coords reoverlap.filtered.coords"); 

	open COORDS, "reoverlap.filtered.coords";
	open NEWCOORDS, ">reoverlap.filtered.edited.coords";
	open REJECT, ">reoverlap.rejected.coords";
	while (<COORDS>) {
		chomp;
		my @line = split (/\t/,);
		if ($pairs{$line[17]}{$line[18]} and ($line[19] eq "[END]") and ($line[3] < $line[4])) {
			$line[17] = $newindex{$line[17]};
			$line[18] = $newindex{$line[18]};
			my $newline = join("\t", @line);
			print NEWCOORDS "$newline\n"; #<STDIN>;
		} else {
			print REJECT "$_\n"; #<STDIN>;
			
		}
	}

	open RESH, ">GARM.merge.sh";
	print RESH "$amosbin/toAmos -s reoverlap.fasta -o reoverlap.afg && \ \n ";
	print RESH "$amosbin/bank-transact -c -z -b reoverlap.bnk -m reoverlap.afg && \ \n";
	print RESH "cp -R reoverlap.bnk bkp.reoverlap.bnk && \ \n";
	print RESH "$amosbin/nucmer2ovl -ignore $ignore -tab reoverlap.filtered.edited.coords | $amosbin/sort2 > reoverlap.ovl && \ \n";
	print RESH "$amosbin/ovl2OVL reoverlap.ovl > reoverlap.OVL && \ \n";
	print RESH "$amosbin/bank-transact -z -b reoverlap.bnk -m reoverlap.OVL && \ \n";
	print RESH "$amosbin/tigger -b reoverlap.bnk && \ \n";
	print RESH "$amosbin/make-consensus -B -e 0.01 -b reoverlap.bnk -w 15 && \ \n";
	print RESH "$amosbin/bank2contig reoverlap.bnk > reoverlap.contig && \ \n";
	print RESH "$garmbin/contig2rpf_minimus2.l reoverlap.contig > reoverlap.read.placed && \ \n";
	print RESH "$amosbin/bank2fasta -b reoverlap.bnk > reoverlap.merged.fasta && \ \n";
	print RESH "$amosbin/listReadPlacedStatus -S -E reoverlap.bnk > reoverlap.singletons && \ \n";
	print RESH "$amosbin/dumpreads -e -E reoverlap.singletons reoverlap.bnk > reoverlap.singletons.fasta && \ \n";
	print RESH "perl -e 'while (<>) { my \@line = split; my \$newname = \"reoverlap_\".(split /\_/, \$line[5])[1]; \$line[5] = \$newname; my \$newline = join(\" \", \@line); print \"\$newline\\n\"; }' reoverlap.read.placed > reoverlap.edited.read.placed && \ \n";
	print RESH "sed 's/>/>reoverlap_contig/' reoverlap.merged.fasta > reoverlap.merged.fasta2 && \ \n";
	print RESH "mv reoverlap.merged.fasta2 reoverlap.merged.fasta && \ \n";
	print RESH "cp reoverlap.merged.fasta reoverlap.edited.read.placed ../; cd ../ && \ \n";
	print RESH "$garmbin/fix_scaffs_with_reoverlap.pl reoverlap.merged.fasta reoverlap.edited.read.placed new_scaffolds.txt.fasta.agp new_scaffolds.txt.fasta.all_split.fasta && \ \n";
	print RESH "$garmbin/scaff_split.pl new_scaffolds.txt.fasta.agp.fixed.fasta\n";	
	
	close RESH;

	system("sh GARM.merge.sh; wait");
} else {
# 	system ("cd ..; $garmbin/scaff_split.pl new_scaffolds.txt.fasta.agp.fixed.fasta\n");	$merged
    system ("cd ..; $garmbin/scaff_split.pl $merged\n");

}
