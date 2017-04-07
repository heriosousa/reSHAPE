#!/usr/bin/perl

# split_fasta.pl - split a FASTA file into multiple files with no more than 1000 sequences per file.
#
# Some analysis scripts have a limit on the number of sequences they will accept in a single file,
# and this script is an easy way to split up a large file into multiple smaller ones.

use warnings;
use strict;
#use DaemonsLog::Log;

if (@ARGV == 0) {
  die "Usage: $0 file.fasta\n";
}

for my $filename (@ARGV) {
  split_fasta($filename);
}


sub split_fasta {
  my ($filename) = @_;

  my $chunk = 1;
  my $seq_in_chunk = 0;
  open (my $fh, "<", $filename) or die "Failed to open $filename: $!\n";
  my $chunk_filename = sprintf("%s.%03d", $filename, $chunk);
  if (-e $chunk_filename) {
    die "File '$chunk_filename' already exists.\n";
  }
  open (my $out, ">", $chunk_filename) or die "Failed to open $chunk_filename for output: $!\n";
  print "Creating chunk file $chunk_filename\n";

  while (<$fh>) {
    if (/^>/) {
      if (++$seq_in_chunk > 1000) {
        close $out;
        $chunk++;
        open ($out, ">", "$filename.$chunk") or die "Failed to open $filename.$chunk: $!\n";
        $seq_in_chunk = 1;
      }
    }
    print $out $_;
  }
  close $fh;
}
