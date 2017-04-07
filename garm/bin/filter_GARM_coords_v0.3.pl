#!/usr/bin/env perl

use strict;

unless (@ARGV) {
    die "Usage: filter_GARM_coords.pl <splimer.coords>\n"; 
}
my $file = $ARGV[0];
#my $name = (split /\./, $file)[0];
my $name = $file;
open COORDS, "$file" || die "Can't open $file\n";

my @contain;
my @overlaps;
my @identical;
my %exclude;

open OVL, ">$name.filtered.overlap.coords" || die "Can't open $name.filtered.overlap.coords\n";
open REPEATS, ">$name.repeats.coords" || die "Can't open $name.repeats.coords\n";
open CON, ">$name.contain.coords" || die "Can't open $name.contain.coords\n";
open IDE, ">$name.ide.coords" || die "Can't open $name.ide.coords\n";
open IDELOG, ">$name.identical.log" || die "Can't open $name.identical.log\n";
open LOG, ">$name.filter.log" || die "Can't open $name.filter.log\n";
open EXCLUDE, ">$name.excluded" || die "Can't open $name.excluded\n";
open RESCUE, ">$name.reevaluated.overlaps.coords" || die "Can't open $name.reevaluated.overlaps.coords\n";

while (<COORDS>) {
	#print; <STDIN>;
	my @line = split (/\s+/, $_);
	my $ide = $line[9];
	my $ref = $line[17];
	my $qry = $line[18];
	my $tag = $line[19];
	next if ($ref eq $qry);
	if ($tag eq "[IDENTITY]") {
		push (@identical, $_);
	}
	elsif ($tag eq "[CONTAINS]" || $tag eq "[CONTAINED]") {
		#print; <STDIN>;
		push (@contain, $_);
	} else {
		my $qstart = $line[3];
		my $qend = $line[4];
		if ($tag eq "[BEGIN]") {
			push (@overlaps, $_);
		} elsif ($tag eq "[END]") {
			push (@overlaps, $_);
		}	
	}
}

#print "Printing contigs contained in larger ones...\n";
CONTAIN(\@contain);

#print "Printing identical contigs ...\n";
foreach my $line (@identical) {
	print IDE "$line";
	my $qry = (split /\s+/, $line)[18];
	my $ref = (split /\s+/, $line)[17];
	my $ide = (split /\s+/, $line)[9];
	print IDELOG "$qry\t$ref\t$ide% IDE\n";
}
#print "DONE!\n";
close IDELOG;
close CON;


my %checkreps;
foreach my $line (@overlaps) {
	#print "$line"; <STDIN>;
	my $ref = (split /\s+/, $line)[17];
	my $qry = (split /\s+/, $line)[18];
	my $tag = (split /\s+/, $line)[19];
	$checkreps{$ref}{$tag}++;
}

my %probreps;
foreach my $line (@overlaps) {
	my $ref = (split /\s+/, $line)[17];
	my $qry = (split /\s+/, $line)[18];
	my $tag = (split /\s+/, $line)[19];
	if ($checkreps{$ref}{$tag} > 1 ) {
		if ($checkreps{$ref}{$tag} > 1) {
			print LOG "PROBREPS\t$ref\t$tag\t$checkreps{$ref}{$tag}\n";
			#$exclude{$ref}++;
		}
		push (@{$probreps{$ref}{$tag}}, $line);	
		#print REPEATS "$line";
	} else {
		print OVL "$line";
	}
}
PROBREP(\%probreps);


sub PROBREP {
###

my $hash = shift;
my %probreps = %$hash;

foreach my $ref (keys %probreps) {
	foreach my $tag (keys %{$probreps{$ref}}) {
		my $count = 1;
		foreach my $line (sort { (split /\s+/, $b)[7] <=> (split /\s+/,$a)[7] } @{$probreps{$ref}{$tag}}) {
			if ($count == 1) {
				print RESCUE "$line";
				$count++;
			} else {
				print REPEATS "$line";
			}
		}
	}
}
 

###
}


sub CONTAIN {
########
my $contain = shift;
#my $count = scalar @$contain;
#print $count; <STDIN>;

my %multiple;
foreach my $line (@$contain) {
	my $ref = (split /\s+/, $line)[17];
	my $qry = (split /\s+/, $line)[18];
	my $tag = (split /\s+/, $line)[19];
	push (@{$multiple{$ref}}, $line);
}

foreach my $line (@$contain) {
	my $ref = (split /\s+/, $line)[17];
	my $qry = (split /\s+/, $line)[18];
	my $tag = (split /\s+/, $line)[19];
	if ($tag eq "[CONTAINS]") {
		print CON "$line";
	} else {
		if (scalar @{$multiple{$ref}} > 1) {
			print LOG "PROBHAPLOTYPE\t$ref\n";
			my $join = join("", @{$multiple{$ref}});
			print REPEATS "$join";
			$exclude{$ref}++;
		} else {
			print CON "${$multiple{$ref}}[0]";
		}
	}
}

#print "DONE!\n";
#######
}

foreach my $exclude (sort keys %exclude) {
	print EXCLUDE "$exclude\n";
}

close LOG;
close OVL;
close REPEATS;
close EXCLUDE;

if (-s "clipped.coords") {
	system ("cat $name.filtered.overlap.coords $name.contain.coords $name.reevaluated.overlaps.coords $name.ide.coords clipped.coords > filtered.coords; wait");
} else {
system ("cat $name.filtered.overlap.coords $name.contain.coords $name.reevaluated.overlaps.coords $name.ide.coords > filtered.coords; wait");
}
