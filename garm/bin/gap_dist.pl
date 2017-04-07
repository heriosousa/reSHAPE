#!/usr/bin/env perl

use strict;
use lib "$ENV{GARMLIB}";
use GetData;

die "Usage:gap_dist.pl <AGP file>\n" unless (scalar @ARGV == 1);
my $file = shift;
my @gaps = GetData::gapsAGP($file);

#my $join = join(",", @gaps);
#print $join;
foreach my $gap (@gaps) {
	print "$gap\n";
}
