#!/usr/bin/env perl 

use strict;
use lib "$ENV{GARMLIB}";
use GetData;

die "Usage: fix_scaffs_with_reoverlap.pl <reoverlap.merged.fasta> <reoverlap.edited.read.placed> <new_scaffolds.txt.fasta.agp> <new_scaffolds.txt.fasta.all_split.fasta>\n" unless scalar @ARGV == 4;

my $reovlfile = shift;
my $rpffile = shift;
my $agpfile = shift;
my $newscafffile = shift;

die "Files must contain info:\nreoverlap.merged.fasta\nreoverlap.edited.read.placed\n" unless (-s "reoverlap.merged.fasta" and -s "reoverlap.edited.read.placed");

my %reovl = GetData::loadfasta($reovlfile);
my %oriseqs = GetData::loadfasta($newscafffile);
my %rpf = GetData::loadrpf($rpffile);

open AGP, "$agpfile";

my %scaff;
while (<AGP>) {
	my @line = split(/\t/);
	push(@{$scaff{$line[0]}}, $_);
}

my %seqs;
my %dirs;
$dirs{'FWD'} = "+";
$dirs{'REV'} = "-";

foreach my $scf (sort {(split /\_/, $a)[$#_] <=> (split /\_/, $b)[$#_] } keys %scaff) {
	#print "scaff $scf"; <STDIN>;
	my $total = scalar @{$scaff{$scf}};
	for (my $index = 1; $index <= $total-2; $index+=2) {
		#print "index $index I'm checking ${$scaff{$scf}}[$index]"; <STDIN>;
		my $left = (split /\s+/, ${$scaff{$scf}}[$index-1])[5];
		my $right = (split /\s+/, ${$scaff{$scf}}[$index+1])[5];
		#print "$left $rigth"; <STDIN>;
		if ($rpf{$left} and $rpf{$right}) { ##GAP CLOSED!!!
			my $dir1 = (split /\s+/, $rpf{$left})[1];
			my $dir2 = (split /\s+/, $rpf{$right})[1];
			my $mer = (split /\s+/, $rpf{$right})[2];
			if ($dir1 eq $dir2) {
				my $seq = $reovl{$mer};
				if ($dir1 eq "REV") {
					$seq = reverse $reovl{$mer};
				}
				$seqs{$mer} = $seq;
				
			} else {
				print "$left $dir1\t$right $dir2\tProblems with orientation...\n"; <STDIN>;
			}
			#print "GAP CLOSED HERE!"; <STDIN>;
			#print "disappear index $index I'm checking ${$scaff{$scf}}[$index]"; <STDIN>;
			my $closedlen = length($reovl{$mer});
			my $posini = (split /\s+/, ${$scaff{$scf}}[$index-1])[1];
			my $posfinal = $closedlen+$posini; 
			my $piece = (split /\s+/, ${$scaff{$scf}}[$index-1])[3];
			my $newline = "$scf\t$posini\t$posfinal\t$piece\tW\t$mer\t1\t$closedlen\t$dirs{$dir1}\n"; #<STDIN>;
			${$scaff{$scf}}[$index-1] = $newline;
			my $oldgap = ${$scaff{$scf}}[$index+2];
			$posini += $closedlen+1;
			my $gap = (split /\s+/, $oldgap)[5];
			$posfinal += $gap;
			my $newgap = "$scf\t$posini\t$posfinal\t0\tN\t$gap\tfragment\tyes\n";
			${$scaff{$scf}}[$index+2] = $newgap;
			splice (@{$scaff{$scf}},$index, 2);
			$index = $index-2;
			$total = scalar @{$scaff{$scf}};
		} 
	}

}

open FIXED, ">$agpfile.fixed";
open FASTA, ">$agpfile.fixed.fasta";
foreach my $scf (sort {(split /\_/, $a)[$#_] <=> (split /\_/, $b)[$#_] } keys %scaff) {
	my $conta = 1;
	my $previous = "";
	print FASTA ">$scf\n";
	foreach my $line (@{$scaff{$scf}}) {
		my @oldline = split (/\s+/, $line);
		if ($oldline[6] eq "fragment") {
			my $gap = "N" x $oldline[5];
			print FASTA "$gap";
		} else {
			if ($oriseqs{$oldline[5]}) {
				print FASTA "$oriseqs{$oldline[5]}";
			} 
			elsif ($seqs{$oldline[5]}) {
				print FASTA "$seqs{$oldline[5]}";
			}
		}
		if ($previous) {
			$oldline[1] = $previous+1; 	 
			if ($oldline[4] eq "N") {
				$oldline[2] = $oldline[1]+$oldline[5];
			} else { 
				$oldline[2] = $oldline[1]+$oldline[7];
			}
			$previous = $oldline[2];
		} else {
			$previous = $oldline[2];
		}
		$oldline[3] = $conta;
		my $join = join("\t", @oldline);
		print FIXED "$join\n"; #<STDIN>;
	}
	print FASTA "\n";
}
