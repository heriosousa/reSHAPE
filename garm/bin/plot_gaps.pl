#!/usr/bin/env perl
use strict;

unless (@ARGV) {
	die "Usage: plot_gaps.pl <AGP file>\n";
}

my $agp = shift;
my $garmbin = $ENV{GARMBIN};
open R, ">plot.R";
print R "pdf(file=\"$agp.gapdist\")\n";
print R "x <- scan(pipe(\"$garmbin/gap_dist.pl $agp\"))\n";
print R "hist(x,breaks=200,main=\"$agp GAP DISTRIBUTION\" col=\"gray\")\n";
