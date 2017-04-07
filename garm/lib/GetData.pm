package GetData;

use strict;
use warnings;
use File::Path;
require Exporter;


our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(loadcoords loadAGP filtercontigs checkheaders loadfasta loadrpf gapsAGP splitscaff headlength revcomp);

sub loadcoords {

my $file = shift;
open COORDS, "$file";
my %coords;
while (<COORDS>) {
	my @line = (split /\s+/, $_);
	unless ($line[0]) {
		shift @line;
	}
	my $ref = $line[17];
	my $join = join("\t", @line);
	push (@{$coords{$ref}}, $join); 
}
return %coords;	
		
###	
}

sub loadAGP {
###
my $AGP = shift;
open AGP, "$AGP";
my %contigs;
while (<AGP>) {
	unless (/fragment/) {
		my @line = (split /\s+/, $_);
		my $scaff = $line[0];
		my $contig = $line[5];
		$contigs{$contig} = $scaff;
	}
}
return %contigs;

###
}


sub gapsAGP {
###

my $AGP = shift;
open AGP, "$AGP";
my @gaps;
while (<AGP>) {
	if (/fragment/) {
		my $gap = (split /\t/)[5];
		print "$gap\n";
		push(@gaps, $gap);
	}
}
return @gaps;
###
}


sub filtercontigs {
##

my $file = shift;
my $bplimit = shift;
my $shortest = shift;
my %contigs;
my $id = "";
my $filename = "";
if ($file =~ /\//) {
	my @strip = split (/\//, $file);
	$filename = $strip[$#strip];
} else {
	$filename = $file;
}

my $modfile = "";
if ($shortest >= $bplimit) {
	$bplimit = $shortest;
	$modfile = $file;
} else {
	open FILE, "$file" || die "can't open $file\n";
	open FILTER, ">$filename.min_$bplimit" || die "can't open $filename.min_$bplimit";
	open LEFTOVER, ">$filename.leftovers" || die "$filename.leftovers";

	while (<FILE>) {
   	     chomp;
   	     if (/^>/) {
   	             $id = $_;
   	     } else {
   	             tr/a-z/A-Z/;
   	             $contigs{$id} .= $_;
   	     }
	}
	foreach my $contig (sort {length $contigs{$b} <=> length $contigs{$a}} keys %contigs) {
   	     my $long = length $contigs{$contig};
   	     if ($long >= $bplimit) {
   	             print FILTER "$contig\n$contigs{$contig}\n";
   	     } else {
				print LEFTOVER "$contig\n$contigs{$contig}\n";
			}

	}

$modfile = "$filename.min_$bplimit";
}
close FILTER;
close LEFTOVER;

return $modfile;
##
}

sub loadrpf {
###

my $file = shift;
open RPF, "$file";
my %rpf;
my %orient;
$orient{"0"} = "FWD";
$orient{"1"} = "REV";

while (<RPF>) {
	my $contig = (split /\s+/, $_)[5];
	my $read = (split /\s+/, $_)[1];
	my $pos = (split /\s+/, $_)[7];
	my $dir = (split /\s+/, $_)[4];	
	$rpf{$read} = "$pos $orient{$dir} $contig";
}

return %rpf;
###
}


sub checkheaders {
##

my $file = shift @_;
my $num = shift @_;
my $type = shift @_;
my $filename = "";
if ($file =~ /\//) {
	$filename = (split /\//, $file)[$#_]; 
} else {
	$filename = $file;
}

open FILE, "$file";
open OUT, ">$filename.modname.fasta" || die "can't open $filename.modname.fasta\n";

while (<FILE>) {
        if (/^>/) {
               	tr/>//d;
				my $name = (split /\s+/, $_)[0]; 
				print OUT ">$type.genome$num\_$name\n";
        } else {
                print OUT;
        }
}
my $result = "$filename.modname.fasta";
close OUT;
return $result;
##
}

sub loadfasta {

my $file = shift @_;
open FASTA, "$file";

my %hash;
my $name = "";
while (<FASTA>) {
	chomp;
	if (/^>/) {
		tr/>//d;
		$name = $_; 
	} else {
		$hash{$name} .= $_;
	}
}
close FASTA;

return %hash;

}

sub revcomp {
###
my $dna = shift;
my $revcomp = reverse($dna);
$revcomp =~ tr/ACGTacgt/TGCAtgca/;

return $revcomp;
###
}

sub headlength {
###
my $file = shift @_;
my %seqs = loadfasta($file);

foreach my $seq (keys %seqs) {
	my $size = length $seqs{$seq};
	$seqs{$seq} = $size;
}

return %seqs;
###
}


sub splitscaff {
####

my $scaffs = shift;

open SCAFF, ">$scaffs.contigs_in_scaffs.fasta" || die "Can't open $scaffs.contigs_in_scaffs.fasta\n: $!\n";;
open SINGLE, ">$scaffs.contigs_singletons.fasta" || die "Can't open $scaffs.contigs_singletons.fasta\n: $!\n";;
open AGP, ">$scaffs.agp";

my %seqs = loadfasta($scaffs);

foreach my $key (sort keys %seqs) {
   my $conta = 1;
   my $agpcont = 1;
   my $pos = 1;
   if ($seqs{$key} =~ /N{5,}/) {
	   my @gaps = $seqs{$key} =~ m/N{5,}/g;
	   my @contigs = split (/N{5,}/, $seqs{$key});
	   foreach my $contig (@contigs) {
			my $l = length $contig;
			my $end = $pos+$l-1;
			print SCAFF ">$key\_$conta\n$contig\n"; #<STDIN>;
			print AGP "$key\t$pos\t$end\t$agpcont\tW\t$key\_$conta\t1\t$l\t+\n";
			$agpcont++;
			$pos += $l;
			my $gapend = 0;
			if ($conta-1 < $#contigs) {
				
				my $lgap = length $gaps[$conta-1];
				$gapend = $pos+$lgap-1;
				print AGP "$key\t$pos\t$gapend\t$agpcont\tN\t$lgap\tfragment\tyes\n"; #<STDIN>;
				$pos = $gapend+1;
			    $conta++;
			    $agpcont++;
			} 
	   }
   } else {
	   #my $l = length $seqs{$key};
	   print SINGLE ">$key\n$seqs{$key}\n"; #<STDIN>;
   }
}
close AGP;
close SCAFF;
close SINGLE;

my $singlet = "$scaffs.contigs_singletons.fasta";
my $scf = "$scaffs.contigs_in_scaffs.fasta";
return ($singlet, $scf);

###
}



1;
__END__

