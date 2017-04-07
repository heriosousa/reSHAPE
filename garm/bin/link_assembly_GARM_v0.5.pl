#!/usr/bin/perl -w

use strict;
#use Statistics::Descriptive;
#use List::MoreUtils qw(uniq);
use lib "$ENV{GARMLIB}";
use MoreUtils;
use GetData;
#use Statistics::Descriptive;

die "Usage: link_assembly_GARM.pl <merged.fasta> <read.placed.filtered> <scaffolds.agp> <extras.fasta>\n" unless (scalar @ARGV) == 4; 

my $mergedfile = shift;
my $rpffile = shift;
my $agpfile = shift;
my $extrasfile = shift;

my %extras = GetData::loadfasta($extrasfile);
my %mergedseqs = GetData::loadfasta($mergedfile);

my %AGP = GetData::loadAGP($agpfile);
my %origaps = SCAFFS($agpfile);
my %orient;

$orient{"0"} = "FWD";
$orient{"1"} = "REV";

open LOG, ">link_assembly_GARM.log";

open MISSING, "$extrasfile";
my %missing;
while (<MISSING>) {
	if (/^>/) {
		tr/>//d;
		my $name = (split /\s+/, $_)[0];
		$missing{$name} = 1; 
	}	
}

open RPF, "$rpffile";
my %links;
my %rpf;
while (<RPF>) {
	chomp;
	my @line = split(/\s+/, $_);
	my $contig = $line[1];
	my $orient = $line[4];
	my $merged = $line[5];
	my $pos = $line[7];
	$links{$contig} = "$merged $orient";
	$rpf{$contig} = "$pos";
}

my %newscaff;
my %linkcount;
my $scaffcount = 0;
my $previous = "none";
my %oriscaffs;

foreach my $contig (sort mysort keys %AGP) {
	my $scaff = (split /\_/, $contig)[$#_-1];
	$oriscaffs{$scaff}++;
	unless ($scaff eq $previous) {
		$previous = $scaff;
		$scaffcount++;
	}
	if ($missing{$contig}) {
		push (@{$newscaff{$scaffcount}}, $contig);
	} else {
		if ($links{$contig}) {
			push (@{$linkcount{$links{$contig}}{$scaff}}, $contig);
			unless (scalar @{$linkcount{$links{$contig}}{$scaff}} > 1) {
				push (@{$newscaff{$scaffcount}}, $links{$contig});
				
			}
		} else { # singletons maybe
			#print "wierd $contig"; <STDIN>;
			push (@{$newscaff{$scaffcount}}, $contig);
		}
	}
}

my %expected;
my %observed;
foreach my $index (sort {$a <=> $b} keys %newscaff) {
	#print "new_scaff_$index"; <STDIN>;
	my $elems = scalar @{$newscaff{$index}};
	#print $elems; <STDIN>;
	#my $arrows = "1";
	foreach my $newcontig (@{$newscaff{$index}}) {
		if ($linkcount{$newcontig}) {
			my $scaffsin = scalar keys %{$linkcount{$newcontig}};
			#print "This merged contig ($newcontig) has $scaffsin diff scaffs"; <STDIN>;
			#if ($scaffsin > 1) {
				#print "Have to find links between different scaffolds"; <STDIN>;
			#} else {
				foreach my $oldscaff (sort keys %{$linkcount{$newcontig}}) {
					#print "+ $oldscaff"; <STDIN>;
					next unless ($oriscaffs{$oldscaff});
					my $oldcontigs = scalar @{$linkcount{$newcontig}{$oldscaff}};
#					print "There are $oldcontigs contig(s) of $oldscaff ($oriscaffs{$oldscaff}) in this $newcontig\n"; #<STDIN>; 
					if ($oldcontigs == $oriscaffs{$oldscaff}) {
						print LOG "$oldscaff\t$newcontig\tCONTAINED_ALL\n"; #<STDIN>;
						GAPS($oldscaff, $newcontig, 0);
						undef @{$origaps{$oldscaff}};
						delete $oriscaffs{$oldscaff};
					} else {
						if ($oldcontigs > 1) {
							#my $join = join(" -> ", @{$linkcount{$newcontig}{$oldscaff}});
							#print "calculate de gap(s) for $oldcontigs contigs in $newcontig\n$join"; <STDIN>;
							GAPS($oldscaff, $newcontig, 1);
						} else {
							#NEWGAP($oldscaff, $newcontig);
							#print "$newcontig must be linked to other ($oldscaff=$oriscaffs{$oldscaff} here=$oldcontigs)"; <STDIN>;
							my $onlymerged = ${$linkcount{$newcontig}{$oldscaff}}[0];
							my $onlyindex = (split /\_/, $onlymerged)[$#_];
							#print "onlyindex $onlyindex"; <STDIN>;
							my $newgap = 0;
							my $mergedlen = length($mergedseqs{(split /\s+/, $newcontig)[0]});;
							#print "length of $newcontig $mergedlen"; <STDIN>;
							#my $test = scalar $oriscaffs{$oldscaff};
							#print "-> ($onlyindex-1) == $test-1)"; <STDIN>;
							unless (($onlyindex-1) == $oriscaffs{$oldscaff}-1) {
								my $onlypos = $rpf{$onlymerged};
								my $onlylen = (split /\s+/,$origaps{$oldscaff}[$onlyindex-1])[1];
								my $onlygap = (split /\s+/,$origaps{$oldscaff}[$onlyindex-1])[2];
								$newgap = $onlygap-($mergedlen - ($onlypos+$onlylen));
								my $dir = $orient{(split /\s+/, $newcontig)[1]};
								#print "\n$onlymerged onlygap-(mergedlen-(onlypos+onlylen)) = $onlygap-($mergedlen - ($onlypos+$onlylen))\n"; <STDIN>;
								if ($dir eq "REV") {
									$newgap = $onlygap-$onlypos;
									#print "\nonlygap-onlypos = $onlygap-$onlypos\n"; <STDIN>;
								}
								#print "only new gap $newgap"; <STDIN>;
							}
							my $newinfo = "$newcontig $mergedlen $newgap";
							if (($onlyindex-1) == $oriscaffs{$oldscaff}-1) {
								$newinfo = "$newcontig $mergedlen";
							}
							$origaps{$oldscaff}[$onlyindex-1] = $newinfo;
						}
					}
					#print "$newcontig($oldscaff:$oldcontigs)"; <STDIN>;
				}
			#}
		} #else {
			#print "--**+++ $newcontig"; <STDIN>;
		#}
		#unless ($arrows >= $elems) {
			#print "->";
			#$arrows++;
		#}
	}
	#print "\n";

}

open CE, ">CE_tags.txt";
foreach my $merged (keys %observed) {
	foreach my $region (keys %{$observed{$merged}}) {
		my $obsval = $observed{$merged}{$region};
		my $expval = $expected{$merged}{$region};
		my $diff = $expval-$obsval;
		my $abs = abs($diff);
		if ($abs > 500) {
			if ($diff > 0) {
				print CE "PROB_CE\t$merged region $region COMPRESSION\tprevious_gap=$expval filled_gap=$obsval\n";
			} else {
				print CE "PROB_CE\t$merged region $region EXTENSION\tprevious_gap=$expval filled_gap=$obsval\n";
			}
		}
	}


}
#my($obsdev, $obsmean) = stats(\%observed);
#my($expdev, $expmean) = stats(\%expected);

#print "Gap expected mean = $expmean"; <STDIN>;
#print "Gap expected stdev = $expdev\n"; <STDIN>;
#print "Gap observed mean = $obsmean"; <STDIN>;
#print "Gap observed stdev = $obsdev\n"; <STDIN>;

my %newagp;
my %newlinks;
my $newscf = 0;
my %dirs;
$dirs{'0'} = "+";
$dirs{'1'} = "-";
foreach my $old (sort keys %origaps) {
	#print "original $old\n"; <STDIN>; 
	if (scalar @{$origaps{$old}}) {
		$newscf++;
		my $posini = 1;
		my $posfinal = 0;
		my $piece = 1;
		my $tag = "GARMscaffold_$newscf";
		my @unique = MoreUtils::uniq(@{$origaps{$old}});
		foreach my $agp (@unique) {
			my $merged = (split /\s+/, $agp)[0];
			my $orient = "";
			my $mergedlen = 0;
			my $gap = 0;
#			my $chain = "->";
			if ($merged =~ /merged_contig/) {
				$orient = $dirs{(split /\s+/, $agp)[1]};
				$mergedlen = (split /\s+/, $agp)[2];
				if ((split /\s+/, $agp)[3]) {
					$gap = (split /\s+/, $agp)[3];
				}
#				$chain = (split /\s+/, $agp)[4];
			} else {
				$mergedlen = (split /\s+/, $agp)[1];
				$orient = "+";
				if ((split /\s+/, $agp)[2]) {
					$gap = (split /\s+/, $agp)[2];
				}
				unless ($extras{$merged}) {
					print LOG "$merged\tCHECK\texcluded for some reason or is in singletons\n"; #<STDIN>;
				}
			}
			#print "$old - agp $agp -- $merged - $mergedlen"; <STDIN>;
			$posfinal += $mergedlen;
			my $line = "$tag\t$posini\t$posfinal\t$piece\tW\t$merged\t1\t$mergedlen\t$orient\t$old\n";
			$posini += $mergedlen;
			$piece++;
			push (@{$newagp{$tag}}, $line);
			push (@{$newlinks{$merged}}, $tag);
			if ($gap) {
				$posfinal += $gap;
				my $gapline = "$tag\t$posini\t$posfinal\t$piece\tN\t$gap\tfragment\tyes\n";
				$posini += $gap;
				$piece++;
				push (@{$newagp{$tag}}, $gapline);
			}
		}
	}
}

open NEWAGP, ">new_scaffolds.txt";
foreach my $new (sort {(split /\_/, $a)[$#_] <=> (split /\_/, $b)[$#_]} keys %newagp) {
	#print "$new\t"; <STDIN>;
	foreach my $line (sort {(split /\s+/, $a)[3] <=> (split /\t/, $b)[3]} @{$newagp{$new}}) {
		print NEWAGP $line; #<STDIN>;
	}
}
foreach my $merged (keys %newlinks) {
	if (scalar @{$newlinks{$merged}} > 1) {
		my $join = join(" ", @{$newlinks{$merged}});
		print LOG "$merged links scaffods $join\n"; #<STDIN>;
	}
}


sub mysort {
	my @linea = split(/\_/, $a);
	my @lineb = split(/\_/, $b); 
	my $index1 = $linea[$#linea];
	my $index2 = $lineb[$#lineb];
	my $scaff1 = $linea[$#linea-1];
	my $scaff2 = $lineb[$#lineb-1];
	#print "$scaff1 cmp $scaff2 and $index1 <=> $index2"; <STDIN>;
	($scaff1 cmp $scaff2) or ($index1 <=> $index2); #$index1 <=> $index2;
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
	my @names = split (/\_/, $line[0]);
	my $scaffname = $names[$#names];
	#print "routing scaffs $scaffname"; <STDIN>;
	#my $scaffname = $line[0];
	if (/fragment/) {
		${$scaff{$scaffname}}[$index] .= " $line[5]";
			
	} else {
		$info = "$line[5] $line[7]";
		push (@{$scaff{$scaffname}}, $info);
		$index = (scalar (@{$scaff{$scaffname}}))-1;
	}
}
return %scaff;
###
}

sub GAPS {  
####
my $scaff = shift;
my $merged = shift;
#print "merged $merged"; <STDIN>;
my $notall = shift;
my $dir = $orient{(split /\s+/, $merged)[1]};
#print "dir according to this scaff $scaff $dir"; <STDIN>;
my @local = @{$linkcount{$merged}{$scaff}};
my $total = scalar @local;
if ($dir eq "REV") {
	@local = sort { (split /\_/, $a)[$#_] <=> (split /\_/, $b)[$#_] } @local;
}
my $chain = join("->", @local);

## check if is continuous
for (my $index = 0; $index < $#local; $index++) {
	my $current = $local[$index];
	my $next = $local[$index+1];
	#print "local $current -> $next"; <STDIN>;
	
	my $cindex = (split /\_/, $current)[$#_];
	my $nindex = (split /\_/, $next)[$#_];
	
	## check if are continuous ##
	my $continous = $nindex - $cindex;
	if ($continous == 1) {	
		#print "$cindex -> $nindex"; <STDIN>;
		#print "$origaps{$scaff}[$cindex] -> $nindex"; <STDIN>;
		my $cpos = $rpf{$current};
		my $npos = $rpf{$next};
		#print "cpos $cpos npos $npos "; <STDIN>;
		my $clength = (split /\s+/, $origaps{$scaff}[$cindex-1])[1];
		my $nlength = (split /\s+/, $origaps{$scaff}[$nindex-1])[1];
		my $expgap = (split /\s+/, $origaps{$scaff}[$cindex-1])[2];
	
		my $obsgap = $npos-($cpos+$clength);
		my $check = "$npos-($cpos+$clength)";
		my $region = ($cpos+$clength)."-".$npos;
		if ($dir eq "REV") {
			$obsgap = $cpos-($npos+$nlength);
			$check = "$cpos-($npos+$nlength)";
			$region = ($npos+$nlength)."-".$cpos;
		}
		#print "expected=$expgap\tobserved=$obsgap"; <STDIN>;
		$expected{$merged}{$region} = $expgap;
		$observed{$merged}{$region} = $obsgap;
	} else {
		print LOG "contig missing between $current and $next\n"; #<STDIN>;
	}
}

if ($notall) {
	my $firstmerged = $local[0];
	#print "first $firstmerged"; <STDIN>;
	my $firstindex = (split /\_/, $firstmerged)[$#_];
	my $lastmerged = $local[$#local];
	#print "last $lastmerged"; <STDIN>;
	my $lastindex = (split /\_/, $lastmerged)[$#_];
	my $newgap = 0;
	my $mergedlen = 0;
	unless (($lastindex-1) == scalar $oriscaffs{$scaff}-1) {
		my $lastpos = $rpf{$lastmerged};
		my $lastlen = (split /\s+/,$origaps{$scaff}[$lastindex-1])[1];
		my $lastgap = (split /\s+/,$origaps{$scaff}[$lastindex-1])[2];
		$mergedlen = length($mergedseqs{(split /\s+/, $merged)[0]});
		$newgap = $lastgap-($mergedlen - ($lastpos+$lastlen));
		#print "lastgap-(mergedlen-(lastpos+lastlen)) = $lastgap-($mergedlen - ($lastpos+$lastlen))";
		if ($dir eq "REV") {
			$newgap = $lastgap-$lastpos;
			#print "lastgap-lastpos = $lastgap-$lastpos"; <STDIN>;
		}
	}
	my $newinfo = "$merged $mergedlen $newgap";
#	my $newinfo = "$merged $mergedlen $newgap $chain";
	if (($lastindex-1) == $oriscaffs{$scaff}-1) {
		$newinfo = "$merged $mergedlen";
#		$newinfo = "$merged $mergedlen 0 $chain";
	}
	#print "HERE!"; <STDIN>;
	for (my $newindex=$firstindex-1; $newindex <= $lastindex-1; $newindex++) {
		#print "I will modify $origaps{$scaff}[$newindex]"; <STDIN>;
		$origaps{$scaff}[$newindex] = $newinfo;
	}
	#splice(@{$origaps{$scaff}}, $firstindex, $total-1);
}

}

sub stats {
####

my $ref = shift;
my %gaps = %$ref;
my $stat = Statistics::Descriptive::Sparse->new();
foreach my $contig (keys %gaps) {
	#print $contig; <STDIN>;
	foreach my $region (keys %{$gaps{$contig}}) {
		#print $region; <STDIN>;
		my $gap = $gaps{$contig}{$region};
		$stat->add_data($gap);
	}
}
my $stdev = $stat->standard_deviation();
my $mean = $stat->mean();
return ($stdev, $mean);
###
}
