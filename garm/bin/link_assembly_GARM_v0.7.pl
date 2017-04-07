#!/usr/bin/perl -w

use strict;
#use Statistics::Descriptive;
#use List::MoreUtils qw(uniq);
use lib "$ENV{GARMLIB}";
use MoreUtils;
use GetData;
#use Statistics::Descriptive;
my $mumbin = "$ENV{MUMBIN}";
my $amosbin = "$ENV{AMOSBIN}";

die "Usage: link_assembly_GARM.pl <merged.fasta> <read.placed.filtered> <scaffolds.agp> <extras.fasta> <singletons.fasta> <overlap ide>\n" unless (scalar @ARGV) == 6; 
open COSA, ">cosastderr"; # karel

my $mergedfile = shift;
my $rpffile = shift;
my $agpfile = shift;
my $extrasfile = shift;
my $singletfile = shift;
my $overlapide = shift;

my %extras = GetData::loadfasta($extrasfile);
my $extraseqs = scalar keys %extras;

my %singletons = GetData::loadfasta($singletfile);
my $singletseqs = scalar keys %singletons;

my %mergedseqs = GetData::loadfasta($mergedfile);
my $merseqs = scalar keys %mergedseqs;

my %scaffolds = SCAFFS($agpfile);
my $scaffsinagp = scalar keys %scaffolds;

my %reoverlapseqs;

open LOG, ">link_assembly_GARM.log";

print LOG "### INPUT ###\n";
print LOG "$merseqs sequence(s) in merged contigs file $mergedfile\n"; #<STDIN>;
print LOG "$extraseqs sequence(s) in extras file $extrasfile\n"; #<STDIN>;
print LOG "$singletseqs sequence(s) in singletons file $singletfile\n"; #<STDIN>;
print LOG "$scaffsinagp scaffolds in the AGP file\n\n"; #<STDIN>;

my %orient;

$orient{"0"} = "FWD";
$orient{"1"} = "REV";

### CHECK IF I HAVE WHOLE SCAFFOLDS WITH NO OVL ###

my %scaffsrecover;
foreach my $extra (keys %extras) {
	my $scaff = GETNAMESCAFF($extra);
	if ($scaff)	{
		$scaffsrecover{$scaff}++;
	}
}
foreach my $singlet (keys %singletons) {
	my $scaff = GETNAMESCAFF($singlet);
	if ($scaff)	{
		$scaffsrecover{$scaff}++;
	}
}


#open RECOVER, ">recover.agp" || die "$0 : Can't write to recover.agp\n";
if (-s "recovered_scaffolds_from_extras.agp") {
	system("rm -rf recovered_scaffolds_from_extras.agp");
}

my $scaffrecovextras = 0;
my $contigrecovextras = 0;
open RECOV, ">recovered_scaffolds_from_extras.fasta" || die "$0 : Can't write to recovered_scaffolds_from_extras.fasta\n";
foreach my $scaff (keys %scaffsrecover) {
	if ($scaffolds{$scaff}) {
		my $howmany = scalar @{$scaffolds{$scaff}};
		#print "$scaff - $scaffsextras{$scaff} - $howmany"; <STDIN>;
		if ($scaffsrecover{$scaff} eq $howmany) {
			#print "$scaff has to be recovered from extras file\n"; #<STDIN>;
			$scaffrecovextras++;
			$contigrecovextras += $howmany;
			system("grep -w $scaff $agpfile >> recovered_scaffolds_from_extras.agp");		
			print RECOV ">$scaff\n"; #<STDIN>;
			foreach my $index (0..$howmany-1) {
				my $contig = (split /\s+/, $scaffolds{$scaff}[$index])[0];
				my $seq = "";
				$seq = $extras{$contig} if $extras{$contig};
				$seq = $singletons{$contig} if $singletons{$contig};
				print RECOV "$seq"; #<STDIN>;
				delete $scaffsrecover{$contig};
				delete $extras{$contig} if $extras{$contig};
				delete $singletons{$contig} if $singletons{$contig};
				if ((split /\s+/, $scaffolds{$scaff}[$index])[2]) {
					my $gap = "N" x (split /\s+/, $scaffolds{$scaff}[$index])[2];
					print RECOV "$gap";
				}
			}
			print RECOV "\n"; #<STDIN>;
		}
	}
}
close RECOV;

print LOG "### RECOVERED ###\n";
print LOG "$scaffrecovextras scaffolds recovered from extras file\n"; #<STDIN>;
print LOG "$contigrecovextras contigs in those scaffolds recovered\n\n"; #<STDIN>;

open RPF, "$rpffile";
my %scaffcontigs;
my %merged;
while (<RPF>) {
	chomp;
	my @line = split(/\s+/, $_);
	my $contig = $line[1];
	my $mergedcontig = $line[5];
	my $scaff = GETNAMESCAFF($contig);
	if ($scaff and $scaffolds{$scaff}) {
		push(@{$scaffcontigs{$scaff}}, $_);
	}
	$merged{$mergedcontig} = 1;
}

#my $inlist = scalar keys %merged;
#print "$inlist merged contigs in the RPF list"; <STDIN>;

my %usage; ### KEEP RECORD OF THE USED SCAFFOLDS
my %newscaffs; ### STORE ORIENTED SEQS FOR A NEW SCAFFOLD
my %gapdist; ### USE THIS FOR COMPRESSION OR EXTENSION DETECTION

print LOG "### REBUILDING SCAFFOLDING INFORMATION ###\n"; #<STDIN>; 
foreach my $name (sort keys %scaffcontigs) {
	#print "$name\n"; #<STDIN>;
	my $inrpf = scalar @{$scaffcontigs{$name}};
	my $scaffcount = scalar @{$scaffolds{$name}};
	my %dirs;
	#my $mergedcount = 1;
	
	for (my $index = 0; $index < $scaffcount-1; $index++) {		
		#print $scaffolds{$name}[$index]; <STDIN>;
		if ($index == $scaffcount-1) { ### END OF THE SCAFFOLD
			#print "last scaffold\n$scaffolds{$name}[$index]\n"; <STDIN>;
			my $leftcontig = $scaffolds{$name}[$index];
			my $pattern = (split /\s+/, $leftcontig)[0];
			my $line = (grep(/$pattern\s+/, @{$scaffcontigs{$name}}))[0];
			my $mergedcontig = "";
			my $mergedlen = "";
			my $orient = "";
			my $newinfo = "$scaffolds{$name}[$index]";
			if ($line) {
				$mergedcontig = (split /\s+/, $line)[5];
				$mergedlen = length $mergedseqs{$mergedcontig};
				if (%dirs) {
					#print "there is the dirs hash"; <STDIN>;
					$orient = CHECKORIENT(\%dirs) if %dirs;
					#print $orient; <STDIN>;
				} else {
					$orient = (split /\s+/, $line)[4];
				}
				$newinfo = "$mergedcontig $mergedlen 0 $orient";
			}
			else {
				$newinfo = $leftcontig;
			}
			
			undef %dirs;
						
			push(@{$usage{$mergedcontig}}, $name) if $mergedcontig;
			push(@{$newscaffs{$name}}, $newinfo);
			
		} else {
			my $leftcontig = $scaffolds{$name}[$index];
			my $pattern1 = (split /\s+/, $leftcontig)[0];
			my $line1 = (grep(/$pattern1\s+/, @{$scaffcontigs{$name}}))[0];
			my $gap = (split /\s+/, $leftcontig)[2];
			my $rightcontig = $scaffolds{$name}[$index+1];
			my $pattern2 = (split /\s+/, $rightcontig)[0];
			my $line2 = (grep(/$pattern2\s+/, @{$scaffcontigs{$name}}))[0];
			#print "$line1\t$line2"; <STDIN>;
			if ($line1 and $line2) {
				#print "both exist\n"; #<STDIN>;
				my $mergedcontig1 = (split /\s+/, $line1)[5];
				my $pos1 = (split /\s+/, $line1)[7];
				my $dir1 = (split /\s+/, $line1)[4];
				my $len1 = (split /\s+/, $line1)[3];
				my $mergedcontig2 = (split /\s+/, $line2)[5];
				my $pos2 = (split /\s+/, $line2)[7];
				my $dir2 = (split /\s+/, $line2)[4];
				my $len2 = (split /\s+/, $line2)[3];
				#print "$mergedcontig1 - $mergedcontig2"; <STDIN>;					
				if ($mergedcontig1 eq $mergedcontig2) { ### CHECK ORIENT AND DISTANCE
					#print "$mergedcontig1 same merged contig for $pattern1 $pattern2\n"; #<STDIN>;
					$dirs{$orient{$dir1}}++;
					$dirs{$orient{$dir2}}++;
					my($distance,$fragment) = CHECKDIST($pos1, $len1, $pos2, $len2);
					#print "fragment=$fragment distance=$distance gap=$gap"; <STDIN>;
					$gapdist{$mergedcontig1}{$fragment} = "distance=$distance gap=$gap";
				} else {
											
					#print "change of contig from $mergedcontig1 to $mergedcontig2 for $pattern1 $pattern2"; <STDIN>;
					#print "recalculate gap"; <STDIN>;
					$dirs{$orient{$dir1}}++;
					my $orient = CHECKORIENT(\%dirs);						
					undef %dirs;
					$dirs{$orient{$dir2}}++;
					my $mergedlen1 = length $mergedseqs{$mergedcontig1};
					my $mergedlen2 = length $mergedseqs{$mergedcontig2};
					#print "$pos1, $len1, $dir1, $mergedlen1, $pos2, $len2, $dir2, $mergedlen2, $gap"; <STDIN>;
					
					#print "both exist $leftcontig\t$rightcontig"; <STDIN>;	
					
					my $newgap = NEWGAP($pos1, $len1, $dir1, $mergedlen1, $pos2, $len2, $dir2, $mergedlen2, $gap);
					
					my $newinfo = "$mergedcontig1 $mergedlen1 $newgap $orient";	
					#print "-> change $name $newinfo"; <STDIN>;
						
					push (@{$usage{$mergedcontig1}}, $name);
					push (@{$newscaffs{$name}}, $newinfo);
				}
							
			} elsif ($line1) {
				#print "line1 yes\n$line1\nline2 no\n-"; <STDIN>;
				my $mergedcontig1 = (split /\s+/, $line1)[5];
				my $pos1 = (split /\s+/, $line1)[7];
				my $dir1 = (split /\s+/, $line1)[4];
				my $len1 = (split /\s+/, $line1)[3];
				my $orient = $dir1;
				if (%dirs) {
					#print "line1 dirs exist...".(keys %dirs)[0]; <STDIN>;
					$orient = CHECKORIENT(\%dirs);
					undef %dirs;
				}
				my $mergedlen1 = length $mergedseqs{$mergedcontig1};
				
				my $recovered = (split /\s+/, $scaffolds{$name}[$index+1])[0];
				my $pos2 = 1;
				my $dir2 = 0;
				my $len2 = (split /\s+/, $scaffolds{$name}[$index+1])[1];
				my $mergedlen2 = $len2;
				
				#print "left exist right recover $leftcontig\t$rightcontig"; <STDIN>;	
				
				my $newgap = NEWGAP($pos1, $len1, $dir1, $mergedlen1, $pos2, $len2, $dir2, $mergedlen2, $gap);	
				my $newinfo = "$mergedcontig1 $mergedlen1 $newgap $orient";
					
				push (@{$usage{$mergedcontig1}}, $name);
				push (@{$newscaffs{$name}}, $newinfo);
				push (@{$newscaffs{$name}}, "$recovered $len2");
				$index++;
				undef %dirs;
			} elsif ($line2) { ### 
				#print "line1 no\n-\nline2 yes\n$line2"; <STDIN>;
				my $recovered = (split /\s+/, $scaffolds{$name}[$index])[0];
				my $pos1 = 1;
				my $dir1 = 0;
				my $len1 = (split /\s+/, $scaffolds{$name}[$index])[1];
				my $gap1 = (split /\s+/, $scaffolds{$name}[$index])[2];
				my $mergedlen1 = $len1;
				
				my $mergedcontig2 = (split /\s+/, $line2)[5];
				my $pos2 = (split /\s+/, $line2)[7];
				my $dir2 = (split /\s+/, $line2)[4];
				my $len2 = (split /\s+/, $line2)[3];
				my $orient = $dir2;
				if (%dirs) {
					#print "line2 dirs exist...but should it?".(keys %dirs)[0]; <STDIN>;
					$orient = CHECKORIENT(\%dirs);
					#undef %dirs;
				}
				my $mergedlen2 = length $mergedseqs{$mergedcontig2};
				
				#print "left recoverted right exist $leftcontig\t$rightcontig"; <STDIN>;
				
				my $newgap = NEWGAP($pos1, $len1, $dir1, $mergedlen1, $pos2, $len2, $dir2, $mergedlen2, $gap);
				my $newinfo = "$mergedcontig2 $mergedlen2 $newgap $orient";
				
				push (@{$usage{$mergedcontig2}}, $name);
				push (@{$newscaffs{$name}}, "$recovered $len1 $gap1 0");
				push (@{$newscaffs{$name}}, $newinfo);
				$index++;
			} else {
				#print "Two contigs in a row don't exist!!! $pattern1 $pattern2"; <STDIN>;
				push (@{$newscaffs{$name}}, $scaffolds{$name}[$index]);
				push (@{$newscaffs{$name}}, $scaffolds{$name}[$index+1]);
				$index++;
				undef %dirs;
			}
		}	
	}
}

my $mergedused = scalar keys %usage;
print LOG "$mergedused merged contigs with scaffolding information\n"; #<STDIN>;

my $mergedsinglet = 0;
open SINGLET, ">merged_contigs_singletons.fasta" || die "$0 : Can't write to merged_contigs_singletons.fasta\n";
foreach my $check (sort { (split /g/, $a)[2] <=> (split /g/, $b)[2] } keys %merged) {
	unless ($usage{$check}) {
		#print "$check singleton merged contig\n"; #<STDIN>;
		print SINGLET ">$check\n$mergedseqs{$check}\n"; #<STDIN>;
		$mergedsinglet++;
	}
}
close SINGLET;

print LOG "$mergedsinglet merged contigs with no scaffolding information\n"; #<STDIN>;


my @fusedscaffolds;
my $inonecontig = 0;

my %complete;
my $rescaffold = 0;

open NEWAGP, ">rescaffold.agp" || die "$0 : Can't write to rescaffold.agp\n";
open NEWAGPSEQ, ">rescaffold.fasta" || die "$0 : Can't write to rescaffold.fasta\n";
foreach my $old (sort { scalar @{$newscaffs{$b}} <=> scalar @{$newscaffs{$a}} } keys %newscaffs) {
	### CHECK LATER ###	
	my @check = @{$newscaffs{$old}};
	my $howmany = scalar @check;
	if ($howmany == 1) {
		$inonecontig++;
		my $contig = (split /\s+/, $check[0])[0];
		$complete{$old} = $contig; ### REPORT THIS SCAFFOLD AS ALL GAPS CLOSED
		### CHECK IF ITS FUSED ### 
		unless ($usage{$contig}) {
			my $seq = $mergedseqs{$contig};
			my $dir = 0;
			$dir = (split /\s+/, $check[0])[3] if (split /\s+/, $check[0])[3];
			$seq = GetData::revcomp($mergedseqs{$contig}) if $dir;
			print NEWAGPSEQ ">$contig\_$old\_no_gaps\n$seq\n";
			delete $mergedseqs{$contig} if $mergedseqs{$contig};
			delete $newscaffs{$old} if $newscaffs{$old};
		}
	}
}

foreach my $old (sort { scalar @{$newscaffs{$b}} <=> scalar @{$newscaffs{$a}} } keys %newscaffs) {
	
	### RESOLVE FUSIONS BEFORE RESCAFFOLDING ###
	foreach my $used (sort { scalar @{$usage{$b}} <=> scalar @{$usage{$a}} } keys %usage) {
		my $howmany = scalar @{$usage{$used}};
		next if $howmany == 1; ## NO NEED TO CHECK A CONTIG RELATED TO JUST ONE SCAFFOLD
		### HERE I PASS THE ARRAY WITH SCAFFOLDS LINKED BY A CERTAIN CONTIG ###
		SPECIAL(\@{$usage{$used}}, $used);
	}
	############################################
	
}

foreach my $old (sort { scalar @{$newscaffs{$b}} <=> scalar @{$newscaffs{$a}} } keys %newscaffs) {
	unless ($complete{$old}) {
		REESCAFFOLD(\@{$newscaffs{$old}}, $old);
	}

}

my $completecount = scalar keys %complete;

close NEWAGP;
close NEWAGPSEQ;

print LOG "$inonecontig scaffolds merged in one merged contig\n"; #<STDIN>;
print LOG "$completecount scaffolds with more than 1 contig merged in a bigger scaffold\n";
print LOG "$rescaffold scaffolds re-scaffolded\n\n"; #<STDIN>;

my $mergedafter = scalar keys %mergedseqs;
my $singletafter = scalar keys %singletons;
my $extrasafter = scalar keys %extras;

print LOG "### LEFTOVERS ###\n";
#print LOG "$mergedafter sequence(s) in merged contigs after rescaffolding\n"; #<STDIN>;
print LOG "$extrasafter sequence(s) in extras after rescaffolding\n"; #<STDIN>;
print LOG "$singletafter sequence(s) in singletons after rescaffolding\n"; #<STDIN>;

#open NOOVL, ">contigs_no_overlaps.fasta" || die "$0 : Can't write to contigs_no_overlaps.fasta\n"; 
LEFTOVERS(\%extras, "contigs_no_overlaps.fasta");
LEFTOVERS(\%singletons, "contigs_with_ovl_but_failed.fasta");
close LOG;


## RECONCILIATON PART ###
#my $checkgapdist = scalar keys %gapdist;
#print "checkgapdist $checkgapdist"; <STDIN>;
system("rm -f probable_CE_problems.txt") if (-s "probable_CE_problems.txt");
open CE, ">probable_CE_problems.txt" || die "$0 : Can't write to probable_CE_problems.txt\n"; 
foreach my $merged (keys %gapdist) {
	#print CE "$merged\t";
	foreach my $region (sort keys %{$gapdist{$merged}}) {
		#print "region=$region\t$gapdist{$merged}{$region}"; <STDIN>;
		print CE "$merged\tregion=$region\t$gapdist{$merged}{$region}\n";
	}
}
close CE;
system("grep NO_OVERLAP reoverlap/reoverlap.log >> probable_CE_problems.txt");


sub LEFTOVERS {
###
my $hash = shift;
my $file = shift;
my %leftovers = %$hash;

open BIN, ">$file" || die "$0 : Can't write to $file\n"; ;

foreach my $name (keys %leftovers) {
	print BIN ">$name\n$leftovers{$name}\n";
}
close BIN;
###
}

#####
# NEED TO FIX THE FUSED SCAFFOLD REPORT
# SUB SPECIAL
# SUB RESOLVE
# SUB RESCAFFOLD
####

sub SPECIAL {
###
my $array = shift;
my $link = shift;
my @linkedscaffolds = @$array;

foreach my $scaff (sort { scalar @{$newscaffs{$a}} <=> scalar @{$newscaffs{$b}} } @linkedscaffolds) {
#	if ($complete{$scaff}) {
#		print "$scaff has all its contig in $link linking contig... no point on checking\n"; <STDIN>;
#	} else {
	unless ($complete{$scaff}) {
		my $contigs = scalar @{$newscaffs{$scaff}};
		#print "$scaff ($contigs) linked to others by $link"; <STDIN>;	
		my %matches;
		foreach my $elem (@{$newscaffs{$scaff}}) {
			my $check = (split /\s+/, $elem)[0];
			#print "check $check"; <STDIN>;
			foreach my $other (@linkedscaffolds) {
				next if ($other eq $scaff);
				next if ($complete{$other} or $complete{$scaff});
				push(@{$matches{$other}}, grep(/^$check /, @{$newscaffs{$other}}));
			}
		}
		foreach my $other (keys %matches) {
			my $count = scalar @{$matches{$other}};
			#print "matches=$count with $other"; <STDIN>;
			if ($count == $contigs) {
				next if $complete{$scaff};
				$complete{$scaff} = $other;
				#print "$scaff is merged completely in $other"; <STDIN>;
			} else {
				#print "\t$scaff is merged partially in $other\nNeed to create a new array!!!"; <STDIN>;
				#print "original scaffold\n";<STDIN>; 
				#PRINT(\@{$newscaffs{$other}});
				my @newarray = RESOLVE($link, $scaff, $other);
				#print "new array\n"; <STDIN>;
				#PRINT(\@newarray);
				#print "check"; <STDIN>;
				@{$newscaffs{$other}} = @newarray;
			}
		}
	}
}


###
}

sub PRINT {
my $array = shift;
my @scaff = @$array;

my $index = 0;
foreach my $elem (@scaff) {
	print "index $index = $elem\n";
	$index++;
}


}


sub RESOLVE {
###
my $link = shift;
my $scaff = shift;
my $other = shift;
my @newarray;

my($small, $big) = sort { scalar @{$newscaffs{$scaff}} <=> scalar @{$newscaffs{$other}} } ($scaff, $other);

my $smallcount = scalar @{$newscaffs{$small}};
my $bigcount = scalar @{$newscaffs{$big}};

#print "link $link"; <STDIN>;
#print "small $small $smallcount"; <STDIN>;
#PRINT(\@{$newscaffs{$small}});
#print "big $big $bigcount"; <STDIN>;
#PRINT(\@{$newscaffs{$big}});

my $fit = "";

my $goback = 0;
my $add = 0;
my $found = 0;
my $smallindex = 0;
for (my $index = 0; $index <= $smallcount-1; $index++) {
	#print "$small $index $newscaffs{$small}[$index]"; <STDIN>;
	my $check = (split /\s+/, $newscaffs{$small}[$index])[0];
	if ($found) {
		$add++;
	} else {
		if ($check eq $link) {
			$found = 1;
			$smallindex = $index; 
		} else {
			$goback++;
		}
	}
}
#print "link $link is the $smallindex element of the array and has $goback elements before and $add elements after\n"; #<STDIN>;
my $stop = 0;

for (my $bigindex = 0; $bigindex <= $bigcount-1; $bigindex++) {
	last if $stop;
	#print "in big index $bigindex or $bigcount-1"; <STDIN>;
	my $check = (split /\s+/, $newscaffs{$big}[$bigindex])[0];
	if ($check eq $link) {
		#print "$check index $index goback=$goback\n"; <STDIN>;
		#push(@newarray, $newscaffs{$big}[$bigindex]);
		if ($goback) {
			#print "I have to go back only"; <STDIN>;
			if ($bigindex == 0 and $smallindex == $smallcount-1) {
				### END OF SMALL AND BEGIN OF BIG ###
				@newarray = @{$newscaffs{$small}};
			} else {
				### INTERLACE CONTIGS OF BOTH SCAFFOLDS
				$bigindex -= $goback;
				for (my $index = ($smallindex-$goback); $index <= $smallindex; $index++) {
					if ($bigindex < 0) {
						#print "out the start of big contig..."; <STDIN>;
						push(@newarray, $newscaffs{$small}[$index]);
						$bigindex++;
					} else {
						my @smallelem = ELEMENTS($small, $index, 0);
						my @bigelem = ELEMENTS($big, $bigindex, 0);
						if ($smallelem[0] eq $bigelem[0]) {
							push(@newarray, $newscaffs{$big}[$bigindex]);
							next;
						}
						my($first, $second) = INTERLACE(\@smallelem, \@bigelem);
						
						$newarray[$bigindex] = $first;
						push(@newarray, $second);
						$bigindex++;
					}
				}
			}
			$goback = 0;
		} 
		if ($add) {

			#print "I have to add only"; <STDIN>;
			if ($bigindex == $bigcount-1 and $smallindex == 0) {
				@newarray = (@newarray, @{$newscaffs{$small}});
			} else {
				for (my $index = ($smallindex+1); $index <= $smallcount-1; $index++) {
					if ($bigindex > $bigcount-1) { ## END OF BIG SCAFFOLD
						#print "reach end of big contig..."; <STDIN>;
						push(@newarray, $newscaffs{$small}[$index]);
					} else {
						my @smallelem = ELEMENTS($small, $index, 0);
						my @bigelem = ELEMENTS($big, $bigindex); ### THIS IS THE LINK
						my($first, $second, $mod) = INTERLACE(\@smallelem, \@bigelem);
						#print "mod $mod"; <STDIN>;
						#if ($mod) { ## MODIFY BIG
						#	$newscaffs{$big}[$bigindex] = $first;
						#} else { ## MOD SMALL
						#	$newscaffs{$small}[$index] = $first;
						#}
						#$newarray[$bigindex] = $first;
						push(@newarray, $first);
						push(@newarray, $second);
						$bigindex++;
					}
				}
				$bigindex--;
			}
			$add = 0;
		}
	} else {
		#print "pushing from bigscaff $newscaffs{$big}[$bigindex]"; <STDIN>;
		push(@newarray, $newscaffs{$big}[$bigindex]);
	}
}
$complete{$small} = $big;
#print "resolved fuck yeah!"; <STDIN>;

return (@newarray);
###
}

sub INTERLACE {
###
my $array1 = shift;
my $array2 = shift;
my @smallelem = @$array1;
my @bigelem = @$array2;

my $smallend = $smallelem[1]+$smallelem[2]; #LENGTH+GAP
my $bigend = $bigelem[1]+$bigelem[2]; #LENGTH+GAP

my $newgap = $bigelem[2]-$smallend;

if (($newgap >= 0) or ($newgap < 0 and $smallend <= $bigend)) { ### SMALL CONTIG GOES INSIDE THE BIG GAP
															    ### OR BIG CONTIG OVERLAPS WITH SMALL 
	$bigelem[2] = $newgap;
	my $newbig = join(" ", @bigelem);
	my $samesmall = join(" ", @smallelem);
	
	return($newbig, $samesmall);
	#push(@newarray, $newbig);
	#push(@newarray, $samesmall);

}
elsif ($smallend > $bigend) { 	## SMALL CONTIG OVERLAPS WITH BIG OR
								## BIG CONTIG GOES INSIDE SMALL GAP
	$newgap = $smallelem[2]-$bigend;
	$smallelem[2] = $newgap;
	my $newsmall = 	join(" ", @smallelem);
	my $samebig = join(" ", @bigelem);
	return($newsmall, $samebig);
}

###
}



sub ELEMENTS {
###
my $scaff = shift;
my $index = shift;
my $array = shift;

#my @newarray = @$array if $array;

my $contig = "";
my $gap = 0;
my $len = 0;
my $dir = 0;

#if ($array) {
#	$contig = (split /\s+/, $newarray[$index])[0];
#	$len = (split /\s+/, $newarray[$index])[1];
#	$gap = (split /\s+/, $newarray[$index])[2] if (split /\s+/, $newarray[$index])[2];
#	$dir = (split /\s+/, $newarray[$index])[3] if (split /\s+/, $newarray[$index])[3];
#} else {
	$contig = (split /\s+/, $newscaffs{$scaff}[$index])[0];
	$len = (split /\s+/, $newscaffs{$scaff}[$index])[1];
	$gap = (split /\s+/, $newscaffs{$scaff}[$index])[2] if (split /\s+/, $newscaffs{$scaff}[$index])[2];
	$dir = (split /\s+/, $newscaffs{$scaff}[$index])[3] if (split /\s+/, $newscaffs{$scaff}[$index])[3];
#}

my @elem = ($contig, $len, $gap, $dir); 				
return (@elem);
###
}


sub REESCAFFOLD {
###
my $array = shift;
my $old = shift;
my @check = @$array;
my @scaff;
	if (scalar @check == 1) {
		@scaff = @check;
	} else {
		@scaff = CLOSEGAPS(\@check);
	}
	my $howmany = scalar @scaff;
	if (scalar @scaff == 1) {
		$inonecontig++;
		#print "check if empty $scaff[0]"; <STDIN>;
		my $contig = (split /\s+/, $scaff[0])[0];
		#print "$mergedseqs{$contig}"; <STDIN>;
		my $seq = $mergedseqs{$contig};
		my $dir = (split /\s+/, $scaff[0])[3];
		$seq = GetData::revcomp($mergedseqs{$contig}) if $dir;
		print NEWAGPSEQ ">$contig\_no_gaps\n$seq\n";
		
	} else {
		$rescaffold++;
		#print "scaffold $old"; <STDIN>;
		my $conta = 1;
		my $begin = 1;
		my $end = 0;
		my $seq = "";
		print NEWAGPSEQ ">$old\_reconstructed\n";
		for (my $index = 0; $index <= $howmany-1; $index++) {
			my @line = (split /\s+/, $scaff[$index]); 
			my $sign = "+";
			$sign = "-" if $line[3];
			if ($line[2]) {
				my $gap = $line[2];
				$end += $line[1];
				#print "@line"; <STDIN>;
				print NEWAGP "$old\t$begin\t$end\t$conta\tW\t$line[0]\t1\t$line[1]\t$sign\n"; #<STDIN>;
				$begin = $end+1;
				$conta++;		
				if ($sign eq "-") {
					$seq .= GetData::revcomp($mergedseqs{$line[0]}) if $mergedseqs{$line[0]};
					$seq .= GetData::revcomp($extras{$line[0]}) if $extras{$line[0]};
					$seq .= GetData::revcomp($singletons{$line[0]}) if $singletons{$line[0]};
				} else {
					$seq .= $mergedseqs{$line[0]} if $mergedseqs{$line[0]};
					$seq .= $extras{$line[0]} if $extras{$line[0]};
					$seq .= $singletons{$line[0]} if $singletons{$line[0]};
				}
				delete $mergedseqs{$line[0]} if $mergedseqs{$line[0]};
				delete $extras{$line[0]} if $extras{$line[0]};
				delete $singletons{$line[0]} if $singletons{$line[0]};
				
				if ($gap < 0) {
					$end = $begin+1;
					print NEWAGP "$old\t$begin\t$end\t$conta\tU\t$gap\tfragment\tyes\n"; #<STDIN>;
					$seq .= "n";
				} else {
					$end += $gap-1;
					print NEWAGP "$old\t$begin\t$end\t$conta\tN\t$gap\tfragment\tyes\n"; #<STDIN>;
					$seq .= "N" x $gap;
				}
				$conta++;
				$begin = $end+1;
			} else {
				$begin = $end+1;
				$end += $line[1]-1;
				print NEWAGP "$old\t$begin\t$end\t$conta\tW\t$line[0]\t1\t$line[1]\t$sign\n"; #<STDIN>;
				#print "last contig, no gap"; <STDIN>;
				if ($sign eq "-") {
					$seq .= GetData::revcomp($mergedseqs{$line[0]}) if $mergedseqs{$line[0]};
					$seq .= GetData::revcomp($extras{$line[0]}) if $extras{$line[0]};
					$seq .= GetData::revcomp($singletons{$line[0]}) if $singletons{$line[0]};
				} else {
					$seq .= $mergedseqs{$line[0]} if $mergedseqs{$line[0]};
					$seq .= $extras{$line[0]} if $extras{$line[0]};
					$seq .= $singletons{$line[0]} if $singletons{$line[0]};
				}
				delete $mergedseqs{$line[0]} if $mergedseqs{$line[0]};
				delete $extras{$line[0]} if $extras{$line[0]};
				delete $singletons{$line[0]} if $singletons{$line[0]};
			}
		}
		print NEWAGPSEQ "$seq\n";
	}

 
###
}

sub CLOSEGAPS {
###
my $array = shift;
my @scaff = @$array;
my $howmany = scalar @scaff;
#die "close gaps for $howmany\n" if $howmany == 1 ;
for (my $index = 0; $index < $howmany-1; $index++) {
	my @line = (split /\s+/, $scaff[$index]); 
	my $sign = "+";
	$sign = "-" if $line[3];
	if ($line[2]) {
		my $gap = $line[2];
		my $reovl = 0;
		my $reovlseq = "";
		if ($gap <= 0) {
			my $sign2 = "+";
			$sign2 = "-" if (split /\s+/, $scaff[$index+1])[3];
			($reovl, $reovlseq) = REOVL((split /\s+/, $scaff[$index])[0], (split /\s+/, $scaff[$index+1])[0], $sign, $sign2, $gap);
		}
		if ($reovl) {
			my $add = (split /\s+/, $scaff[$index+1])[1]-$reovl;
			my $newgap = 0;
			$newgap = (split /\s+/, $scaff[$index+1])[2] if (split /\s+/, $scaff[$index+1])[2];
			my $newlen += ($line[1]+$add);
			my $previous = $line[0];
			$line[0] .= "_".(split /\s+/, $scaff[$index+1])[0];
			#print "reoverlapped $line[0]\n$reovlseq"; <STDIN>;
			$mergedseqs{$line[0]} = $reovlseq;
			delete $mergedseqs{$previous} if $mergedseqs{$previous};
			delete $mergedseqs{(split /\s+/, $scaff[$index+1])[0]} if $mergedseqs{(split /\s+/, $scaff[$index+1])[0]};
			delete $extras{$previous} if $extras{$previous};
			delete $extras{(split /\s+/, $scaff[$index+1])[0]} if $extras{(split /\s+/, $scaff[$index+1])[0]};
			delete $singletons{$previous} if $singletons{$previous};
			delete $singletons{(split /\s+/, $scaff[$index+1])[0]} if $singletons{(split /\s+/, $scaff[$index+1])[0]};
			my $newinfo = "$line[0] $newlen $newgap 0"; 				
			$scaff[$index] = $newinfo;
			$scaff[$index+1] = "remove";
			my @newscaff;
			foreach my $elem (@scaff) {
				next if ($elem eq "remove");
				push(@newscaff, $elem);
			}
			@scaff = @newscaff;
			$howmany = scalar @scaff;
			$index--;
		}
	}

}

return (@scaff);

###
}

sub REOVL {
###

mkdir ("reoverlap") unless (-d "reoverlap");

my $leftname = shift;
my $rightname = shift;
my $sign = shift;
my $sign2 = shift;
my $gap = shift;	

my $left = "";
my $right = "";



open ONE, ">reoverlap/$leftname.fasta" || die "$0 : Can't write to reoverlap/$leftname.fasta\n";
open TWO, ">reoverlap/$rightname.fasta" || die "$0 : Can't write to reoverlap/$rightname.fasta\n";

if ($sign eq "-") {
	$left = GetData::revcomp($mergedseqs{$leftname}) if $mergedseqs{$leftname};
	$left = GetData::revcomp($extras{$leftname}) if $extras{$leftname};
	$left = GetData::revcomp($singletons{$leftname}) if $singletons{$leftname};
} else {
	$left = $mergedseqs{$leftname} if $mergedseqs{$leftname};
	$left = $extras{$leftname} if $extras{$leftname};
	$left = $singletons{$leftname} if $singletons{$leftname};
}

if ($sign2 eq "-") {
	$right = GetData::revcomp($mergedseqs{$rightname}) if $mergedseqs{$rightname};
	$right = GetData::revcomp($extras{$rightname}) if $extras{$rightname};
	$right = GetData::revcomp($singletons{$rightname}) if $singletons{$rightname};
} else {
	$right = $mergedseqs{$rightname} if $mergedseqs{$rightname};
	$right = $extras{$rightname} if $extras{$rightname};
	$right = $singletons{$rightname} if $singletons{$rightname};
}



print ONE ">$leftname\n$left\n";
print TWO ">$rightname\n$right\n";
close ONE;
close TWO;

open OVLLOG, ">>reoverlap/reoverlap.log" || die "$0 : Can't write to reoverlap/reoverlap.log\n";
my $ide = $overlapide;
print OVLLOG "$leftname.vs.$rightname\tgap=$gap"; #<STDIN>;
#print LOG "cd reoverlap; $mumbin/nucmer --maxmatch --prefix=$leftname.vs.$rightname $leftname.fasta $rightname.fasta; wait; $mumbin/show-coords -H -c -l -o -r -I $ide $leftname.vs.$rightname.delta > $leftname.vs.$rightname.coords; cd ..\n";
system ("cd reoverlap; $mumbin/nucmer --maxmatch --prefix=$leftname.vs.$rightname $leftname.fasta $rightname.fasta; wait; $mumbin/show-coords -H -c -l -o -r -I $ide $leftname.vs.$rightname.delta | grep END | sort -k10,10nr > $leftname.vs.$rightname.coords; cd .."); 

my $reovl = 0;
my $reovlseq = "";
if (-s "reoverlap/$leftname.vs.$rightname.coords") {
	my $line = `cat reoverlap/$leftname.vs.$rightname.coords`;
	chomp $line;
	$reovl = (split /\s+/, $line)[7];
	#print $reovl; <STDIN>;
	print OVLLOG "\treoverlap=$reovl\n";
	system("rm reoverlap/$leftname.vs.$rightname.delta");
	my $frag = substr $right, $reovl-1;
	$reovlseq = $left.$frag;
} else {
	print OVLLOG "\tNO_OVERLAP\n";
	system("rm reoverlap/$leftname.vs.$rightname.coords reoverlap/$leftname.vs.$rightname.delta reoverlap/$leftname.fasta reoverlap/$rightname.fasta");
}

return ($reovl, $reovlseq);
###
}

sub NEWGAP {
###
#print "I'm in newgap"; <STDIN>;
my $pos1 = shift;
my $len1 = shift;
my $dir1 = shift;
my $mergedlen1 = shift;
my $pos2 = shift;
my $len2 = shift;
my $dir2 = shift;
my $mergedlen2 = shift;
my $gap = shift;

my $tail1 = 0;
my $tail2 = 0;
my $decision = "";
if ($dir1 == 1 and $dir2 == 1) { ## 1-1
	#$decision = "1-1";
	$tail1 = $pos1-1;
	$tail2 = $mergedlen2-($pos2+$len2);
} elsif ($dir2 and $dir1 == 0) { ## 0-1
	#$decision = "0-1";
	$tail1 = $mergedlen1-($pos1+$len1);
	$tail2 = $mergedlen2-($pos2+$len2);
} elsif ($dir1 and $dir2 = 0) { ## 1-0
	#$decision = "1-0";
	$tail1 = $pos1-1;
	$tail2 = $pos2-1;
} elsif ($dir1 == 0 and $dir2 == 0) { ## 0-0
	#$decision = "0-0";
	$tail1 = $mergedlen1-($pos1+$len1);
	$tail2 = $pos2-1;
}

my $newgap = $gap-($tail1+$tail2);

return ($newgap);

###
}

sub CHECKDIST {
###

my $pos1 = shift;
my $len1 = shift;
my $pos2 = shift;
my $len2 = shift;
#print "$pos1-$pos2"; <STDIN>;
my $distance = 0;
my $fragment = "";
if ($pos1 > $pos2) {
	#print "REV\n";
	$distance = $pos1-($pos2+$len2);
	my $region = $pos1+$len1;
	$fragment = "$pos2<->$region";
} else {
	#print "FWD\n";
	$distance = $pos2-($pos1+$len1);
	my $region = $pos2+$len2;
	$fragment = "$pos1<+>$region"; 
}
#print "distance=$distance"; <STDIN>;
return ($distance, $fragment);

###
}


sub CHECKORIENT {
###
my $dirs = shift;
my $checkdir = scalar keys %$dirs;
#print "scalar keys $checkdir\n";
my $orient = 0;
if ($checkdir == 1) {
	my $dir = (keys %$dirs)[0];
	#print "orientation=$dir"; <STDIN>;
	$orient = 1 if ($dir eq "REV");
} elsif ($checkdir > 1) {
	my($one, $two) = sort { $$dirs{$b} <=> $$dirs{$a} } (keys %$dirs);
	my $more = $$dirs{$one};
	my $less = $$dirs{$two};
	#print "++++orientation disagreement\n$one=$more\n$two=$less\n"; #<STDIN>;
	$orient = 1 if ($one eq "REV");
} else {
	#print "hash empty? $checkdir\n"; #<STDIN>;
}

#print "return $orient"; <STDIN>;
return ($orient);

###
}


sub SCAFFS {
###
my $agpfile = shift;
open SCAFF, "$agpfile";
my %scaff;
my $info = "";
my $index = 0;
while (<SCAFF>) {
	my @line = split (/\s+/, $_);
    print COSA "1 @line\n"; # karel
	my @names = split (/\_/, $line[0]);
    print COSA "2 @names\n"; # karel
	#my $scaffname = $names[$#names];
	#print "$scaffname"; <STDIN>;
	#print "routing scaffs $scaffname"; <STDIN>;
	my $scaffname = $line[0];
	if (/fragment/) {
		${$scaff{$scaffname}}[$index] .= " $line[5]";
			
	} else {
		$info = "$line[5] $line[7]";
		push (@{$scaff{$scaffname}}, $info);
		$index = (scalar (@{$scaff{$scaffname}}))-1;
	}
}
close SCAFF;
return %scaff;
###
}

sub GETNAMESCAFF {
###

#print "I'm in getname"; <STDIN>;
my $contig = shift;
#print $contig; <STDIN>;
my $index = shift;
my @scaff = split(/\_/, $contig);
print COSA "3 @scaff\n";
#my $scaffindex = pop @scaff;
my $scaffindex = $scaff[$#scaff] if $scaff[$#scaff] =~ /\d+/;
print COSA "4 @scaff algo esta mal\n" unless $scaffindex;
#print COSA "4 $scaffindex\n";
my $scaffname = join("_", @scaff);
 print COSA "5 $scaffname\n";
#print $scaffname; <STDIN>;

if ($scaffolds{$scaffname}) {
	if ($index) {
		return ($scaffname, $scaffindex)
	} else {
		return ($scaffname);
	}
} else {
	$scaffname = 0;
	if ($index) {
		$scaffindex = 0;
		return ($scaffname, $scaffindex)
	} else {
		return ($scaffname);
	}
}
###
}

sub mysort {
	my @linea = split(/\_/, $a);
    print COSA "6 @linea\n";
	my @lineb = split(/\_/, $b); 
	my $index1 = $linea[$#linea];
	my $index2 = $lineb[$#lineb];
	my $scaff1 = $linea[$#linea-1];
	my $scaff2 = $lineb[$#lineb-1];
	#print "$scaff1 cmp $scaff2 and $index1 <=> $index2"; <STDIN>;
	($scaff1 cmp $scaff2) or ($index1 <=> $index2); #$index1 <=> $index2;
}
