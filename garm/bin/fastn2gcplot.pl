#!/usr/bin/env perl

=pod

=head1 fastn2length.pl

fastn2gcplot.pl - Perl script to make GC plot of from a fasta or fastq file

=head1 SYNOPSIS

Usage: fastn2gcplot.pl <fast{a,q} file (can be gzipped)> <output_file{pdf,png}> <window length> [-debug]

=head1 DESCRIPTION

Makes a GC plot of a fasta or fastq file (possibly gzipped).
Ouptput is in pdf or png format.  Assumes no extra linebreaks in fasta sequences:
if this is not the case, run fasta2singleLine.pl first.
GC plot is made by calculating the GC content of each window in each sequence
and plotting the histogram.  Ns in the sequences are ignored before calculating the GC
(e.g. ACGTNNNNNNNN would have a GC content of 50%)


=head1 ARGUMENTS

Arg [1] : fasta or fastq filename (can be gzipped with .gz extension)

Arg [2] : output filename: extension decides the ouptut format (pdf or png)

Arg [3] : window size

Arg [optional] : use -debug flag to not delete all temporary files made (eg R script)

=head1 SEE ALSO

See the team wiki page http://scratchy.internal.sanger.ac.uk/wiki/index.php/Team_133

=head1 AUTHOR

Martin Hunt, <mh12@sanger.ac.uk>

=head1 TIME AND DATE

2009-06-29

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 Genome Research Limited. All Rights Reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=cut

use strict;
use warnings;

use File::Temp "tempfile";
use Getopt::Long;

my $infile  = $ARGV[0];
my $outfile = $ARGV[1];
my $window  = $ARGV[2];

my %hist;
my %options = (debug => 0);
my @gc;
my $type = "";

# get the command line options
my $options_ok = GetOptions(\%options, 'debug');

# stop if bad options
if ($#ARGV != 2 or !($options_ok)) {
    print "usage: fastn2gcplot.pl <fast{a,q} file (can be gzipped)> <output_file{pdf,png}> <window length> [-debug]\n",
      "Makes a gc plot of all reads in the file. Output format depends on whether the outfile is\n",
      "called whatever.pdf or whatever.png\n",
      "Use -debug to not clean up:\n",
      "This will leave the tmp R script and .Rout files used to make plot\n\n",
      "IMPORTANT: assumes no extra line breaks in fasta files, if there are, run fasta2singleLine.pl\n\n";
    exit 1;
}

my $in_fh = open_file_reading($infile) or die "error opening $infile\n";


while (<$in_fh>) {
    # determine the file type: fasta or fastq
    unless ($type) {
        if ($_ =~ /^>/) {
            $type = "fasta";
        }
        elsif (/^@/) {
            $type = "fastq";
        }
        else {
            die "file type not recognised\n";
        }
    }

    # this is an id line, so skip it
    $_ = <$in_fh>;

    chomp;
    s/N//gi;    # chuck out the Ns

    my $i = 0;

    # get the gc content of each window in the current sequence
    while ($i < length($_)) {
        my $seq    = substr $_, $i, $window;
        my $length = length $seq;
        my $count  = $seq =~ tr/[c,g,C,G]/a/;
        push @gc, (1.0 * $count / $length);
        $i += $window;
    }

    # skip the header/qual lines if the file is a fastq
    if ($type eq "fastq") {
        $_ = <$in_fh>;
        $_ = <$in_fh>;
    }

}
close $in_fh;

# initialise histogram
foreach (0 .. 100) {
    $hist{$_} = 0;
}

# fill histogram
foreach (@gc) {
    my @a = split /\./, ($_ * 100);
    $hist{ $a[0] }++;
}

# make strings for easy R file writing
my $hist_bins   = "";
my $hist_values = "";

foreach (sort {$a <=> $b} keys %hist) {
    $hist_bins   .= "$_,";
    $hist_values .= "$hist{$_},";
}

$hist_bins   = substr $hist_bins,   0, -1;
$hist_values = substr $hist_values, 0, -1;

# write R script
my ($fh, $fname) = tempfile("fastn2gcplot-tmp-XXXX");

print $fh
  "m=max(c($hist_values))\n",
  "s=sum(c($hist_values))\n",
  "pc_pos = c(0, m/4, m/2, 3*m/4, m)\n",
  "pc_labels= round(100 * pc_pos / s, digits=1)\n";

#pdf or png?
if ($outfile =~ /\.pdf/) {
    print $fh "pdf(file=\"$outfile\")\n";
}
elsif ($outfile =~ /\.png/) {
    print $fh "png(file=\"$outfile\")\n";
}

print $fh
  "par(mar=c(5, 4, 4, 5) + 0.1)\n",
  "barplot(c($hist_values), names.arg=c($hist_bins), xlab=\"GC content (\%)\", ylab=\"Frequency\")\n",
  "axis(4, at=pc_pos, labels=pc_labels)\n",
  "mtext(\"Percentage\", side=4, las=0, line=3)\n",
  "dev.off()\n";

close $fh;

# run R script
system("R CMD BATCH $fname");

# clean up
unless ($options{debug}) {
    unlink $fname;
    unlink "$fname.Rout";
}

sub open_file_reading {
    my $filename = shift;
    my $handle;

    print "Can't find file $filename\n" unless (-e $filename);

    if ($filename =~ /\.gz$/) {
        open $handle, "gunzip -c $filename |" or print "Error opening zipped file $filename";
    }
    else {
        open $handle, $filename or print "Error opening file $filename\n";
    }
    return $handle;
}
