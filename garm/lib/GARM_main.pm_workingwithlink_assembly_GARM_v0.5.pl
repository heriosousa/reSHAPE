package GARM_main;

use strict;
use warnings;
use File::Path;
use lib "$ENV{GARMLIB}";
use GetData;
require Exporter;
my $garmbin = $ENV{GARMBIN};
my $mumbin = $ENV{MUMBIN};
my $amospath = $ENV{AMOSBIN};

our @ISA = qw(Exporter);

our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(usage loadinput checkhost input genomestats splimer merging scaffolding retrievemerged retrievescaff);

sub input {
####

my $genomes = shift;
my $pwd = shift;
my $prefix = shift;
my $mlength = shift;
my $repeatend = shift;
my $repeatsfile = shift;
my $binsize = shift;
my $gcplot = shift;

my @scffiles;
my @qcgen;
my $log = "";
my $filterlog = "";
unless (scalar @$genomes == 2) {
	die "Two genomes are required for the merging. Please, check your input info_genomes_file again\n";
}

my $count = 1;
my $override = "";
foreach my $elem (@$genomes) {
	my $scf = 0;
	my $scaffname = "";
	chomp $elem;
	$log .= "## GENOME$count\nFASTA FILE:\t$elem ##\n";
	my @line = split (/\s+/, $elem);
	my $origenome = $line[0];
	my $genome = $origenome;
	my $type = "NOSPEC";
	if ($line[1]) {
		$type = $line[1];
	}
	$log .= "SEQ TECH SPECIFIED:\t$type\n\n";
	#print $log; <STDIN>;
	if (-z $genome) { ## FILES MUST NOT BE EMPTY... ##
		die "$genome doesn't exist or is empty! Check the $genome file, please.\n";
	} else {
		($genome, $scf, $scaffname, $filterlog) = FILTER($genome, $count, $type, $pwd, $prefix, $mlength, $repeatend, $repeatsfile, $binsize, $gcplot);
		$log .= $filterlog;
		push (@qcgen, $genome);
		if ($line[2] and $line[2] eq "*") {
			$override = $scaffname;
		
		} else {
			if ($scf and $scaffname) {
				push (@scffiles, $scaffname);
			}
		}
	}
	$count++;
}
if ($override) {
	undef @scffiles;
	push (@scffiles, $override);
}
return ($log, \@qcgen, \@scffiles);

###
}


sub usage {

my $version = shift;
my $help = shift;

print "GARM version $version\n";
print "Usage: GARM.pl -g info_genomes_file -o <prefix>\n";
print "Options:\n";	
print "   -g: Text file with the path to two fasta genome assemblies to merge\n";
print "   -o: Prefix for output directory and files\n";
print "   -h: Complete help menu with advanced options\n\n";

if ($help) {
print "PLEASE, CHECK THE DOCUMENTATION FIRST!!!\n";
print "ADVANCE OPTIONS (Do not change unless you know what you are doing):\n";
print "   -c <n>: Run with multi-threads with n cpus (default n=2)\n";
print "   -r: fasta file with repeats\n";
print "   -t: Max length of the end sequence unaligned (defautl 50 bp)\n";
print "   -l: Minimun contig length on each assembly[Defaul 200]\n";
print "   -m: Minimum match (anchor) [Default 200]\n";
print "   -x: Generate GC histograms output [REQUIRE R] \n";
print "   -w: Bin length for GC histograms (default: ave contig size)\n\n";
}
print "For more details, check the DOCUMENTATION\n\n";
print "The info_genomes_file should be a plain text file with 3 colums:\n";
print "first column:\tthe path to the fasta file of one assembly\n";
print "second column:\ta tag (has to be different for each file)\n";
print "You need two files to procced with the merging.\n\n";
exit;
}

sub loadinput {
###

my $input = shift;
open FILE, "$input";
my @genomes;
while (<FILE>) {
	unless (/^$/) {
		chomp;
		push (@genomes, $_);
	}
}

return @genomes;



###
}

sub checkhost {
###

my $host = `hostname`;
chomp $host;
my $q = "normal";
if ($host =~ /seq1/) {
	$q = "phrap";
}

return $q;

###
}

sub FILTER {

my $filterpsl = "$garmbin/filter_psl_repeats.pl";
my $genome = shift;
my $count = shift;
my $type = shift;
my $pwd = shift;
my $prefix = shift;
my $mlength = shift;
my $repeatend = shift;
my $repeatsfile = shift;
my $binsize = shift;
my $gcplot = shift;

my $scfc = 0;
my $scaffname = "";
my $log = "";

$log .= "# GENOME$count FILTERING #\n";

chdir ("$pwd/$prefix/input");

$log .= "# MODIFYING HEADERS #\n";
$genome = GetData::checkheaders($genome, $count, $type);
my %headlength = GetData::headlength($genome);
my $shortest = "";
foreach my $head (sort { $headlength{$a} <=> $headlength{$b} } keys %headlength) {
	$shortest = $headlength{$head};
	last;
}
chomp $shortest;
$log .= "# FILTERING CONTIGS/SCAFFOLDS BY MIN LENGTH $mlength #\n";
$genome = GetData::filtercontigs($genome, $mlength, $shortest, $pwd);

$log .= "# SPLITTING SCAFFOLDS WITH MORE THAN ONE CONTIG #\n";
my ($singlet, $scf) = GetData::splitscaff("$genome");

my $allsplit = "$genome\_split_all.fasta";
if (-s "$scf" and -s "$singlet") {
	$log .= "# SCAFFOLDS AND SINGLETON CONTIGS #\n";
	system ("cat $scf $singlet > $allsplit");
	mkdir("$pwd/$prefix/scaffolding") unless (-d "$pwd/$prefix/scaffolding");
	system("cp $genome.agp $pwd/$prefix/scaffolding");
	$scaffname = "$genome";
	$scfc++;
}

elsif (-s "$scf") {
	$log .= "# ALL SCAFFOLDS HAD 2 OR MORE CONTIGS #\n";
	system ("cat $scf > $genome\_split_all.fasta");
	system ("rm $singlet");
	mkdir("$pwd/$prefix/scaffolding") unless (-d "$pwd/$prefix/scaffolding");
	system("cp $genome.agp $pwd/$prefix/scaffolding");
	$scaffname = "$genome";
	$scfc++;
}
elsif (-s "$singlet") {
	$log .= "# NO SCAFFOLDS... ALL SINGLE CONTIGS #\n";
	system ("cat $singlet > $genome\_split_all.fasta");
	system ("rm $scf");
}

if ($gcplot) {
	if ($binsize) {
		system ("$garmbin/fastn2gcplot.pl $allsplit $allsplit.pdf $binsize");
	} else {
		my $contigave = `$garmbin/garm_stats $allsplit | head -n 1 | tr -d "," | awk '{print \$9}'`;
		system ("$garmbin/fastn2gcplot.pl $allsplit $allsplit.pdf $contigave");
	}
}

mkdir ("$pwd/$prefix/merging/") unless (-d "$pwd/$prefix/merging");
system ("cp $allsplit $pwd/$prefix/merging/");

if ($repeatsfile and -e "$repeatsfile") {
	$log .= "# REPEAT INFORMATION PROVIDED #\n"; #<STDIN>;
	$log .= "# RUNNING BLAT TO FIND REPEATS IN $genome #\n"; #<STDIN>;
	system ("blat $genome\_split_all.fasta $repeatsfile -minScore=300 genome$count.psl && $filterpsl genome$count.psl $repeatend");
	my $reps = "genome$count.psl.repeats";
	system ("cat $reps >> $pwd/$prefix/merging/repeats.txt");
}

$log .= "INPUT PART - DONE\n";

if ($scfc) {
	return ($allsplit, $scfc, $scaffname, $log);
} else {
	return ($allsplit, $scfc, 0, $log);
}

####
}

sub genomestats {
###
my $genome1 = shift;
my $genome2 = shift;

my $stats = "";
my $count1 = `grep -c ">" $genome1`; 
my $count2 = `grep -c ">" $genome2`;
my $size1 = `$garmbin/garm_stats $genome1 | head -n 1 | tr -d "," | awk '{print \$3}'`;
my $size2 = `$garmbin/garm_stats $genome2 | head -n 1 | tr -d "," | awk '{print \$3}'`;
my $ave1 = `$garmbin/garm_stats $genome1 | head -n 1 | tr -d "," | awk '{print \$9}'`;
my $ave2 = `$garmbin/garm_stats $genome2 | head -n 1 | tr -d "," | awk '{print \$9}'`;
chomp ($count1, $count2, $size1, $size2, $ave1, $ave2);
$stats .= "\tContigs\tTotal bases\tAve.Contig size\n";
$stats .= "$genome1\t$count1\t$size1\t$ave1 <--- Best stats!\n";
$stats .= "$genome2\t$count2\t$size2\t$ave2\n";

							
if ($count1 == 1 or $count2 == 1) {
	print "Some of your genome files only have one sequence... Perhaps running AMOScmp or ABACAS will be better?\n";
}
## USE ALWAYS THE LESS FRAGMENTED BIGGER SET AS A REFERENCE ##
if ($count2 >= $count1 and $size2 <= $size1 and $ave2 <= $ave1) { 
	return ($genome1, $genome2, $stats);
} elsif ($size2 <= $size1 and $ave2 <= $ave1){
	return ($genome1, $genome2, $stats);
} elsif ($ave2 <= $ave1) {
	return ($genome1, $genome2, $stats);
} else {
	$stats = "";
	$stats .= "\tContigs\tTotal bases\tAve.Contig size\n";
	$stats .= "$genome1\t$count1\t$size1\t$ave1\n";
	$stats .= "$genome2\t$count2\t$size2\t$ave2 <--- Best stats!\n";
	return ($genome2, $genome1, $stats);
}
###
}


sub splimer {
###
	my $filter = "$garmbin/filter_GARM_coords_v0.2.pl";
	my $ref = shift;
	my $qry = shift;
	my $match = shift; 
	my $filteride = shift;
	my $overlapide = shift;
	my $splimer = shift;
	my $cpus = shift;
	my $maxtrim = shift;
	my $random = shift;
	if ($cpus) {
		system ("$splimer $ref $qry $cpus $match $filteride $overlapide $maxtrim");
	}
	else {
		system("$splimer $ref $qry $match $filteride $overlapide $maxtrim $random"); ###
	}
	
###
}

sub merging {
###

my $genome1 = shift;
my $genome2 = shift;
my $maxtrim = shift;
my $prefix = shift;

open MERGE, ">GARM.merging.sh";
print MERGE "awk '{print \$18}' filtered.coords | sort | uniq > ref.include && \ \n";
print MERGE "awk '{print \$19}' filtered.coords | sort | uniq > qry.include && \ \n";
print MERGE "$garmbin/grep_fasta.l -l ref.include $genome1 > $genome1.filtered && \ \n";
print MERGE "$garmbin/grep_fasta.l -l qry.include $genome2 > $genome2.filtered && \ \n";
print MERGE "$garmbin/edit_coords4merging.pl $genome1.filtered $genome2.filtered filtered.coords && \ \n";
print MERGE "cat $genome1.filtered $genome2.filtered > $prefix.2merge.fasta && \ \n";
print MERGE "$amospath/toAmos -s $prefix.2merge.fasta -o $prefix.afg && \ \n";
print MERGE "$amospath/bank-transact -c -z -b $prefix.bnk -m $prefix.afg && \ \n";
print MERGE "cp -R $prefix.bnk bkp.$prefix.bnk && \ \n";
print MERGE "$amospath/nucmer2ovl -ignore $maxtrim -tab edited.coords | $amospath/sort2 > $prefix.ovl && \ \n";
print MERGE "$amospath/ovl2OVL $prefix.ovl > $prefix.OVL && \ \n";
print MERGE "$amospath/bank-transact -z -b $prefix.bnk -m $prefix.OVL && \ \n";
print MERGE "$amospath/tigger -b $prefix.bnk && \ \n";
print MERGE "$amospath/make-consensus -B -e 0.03 -b $prefix.bnk -w 15 && \ \n";
print MERGE "$amospath/bank2contig $prefix.bnk > $prefix.contig && \ \n";
print MERGE "$garmbin/contig2rpf_minimus2.l $prefix.contig > $prefix.read.placed && \ \n";
print MERGE "$amospath/bank2fasta -b $prefix.bnk > $prefix.merged.fasta \ && \n";
print MERGE "$amospath/listReadPlacedStatus -S -E $prefix.bnk > $prefix.singletons && \ \n";
print MERGE "$amospath/dumpreads -e -E $prefix.singletons $prefix.bnk > $prefix.singletons.fasta && \ \n";
print MERGE "sed 's/>/>merged_contig/' $prefix.merged.fasta > $prefix.merged.fasta2 && \ \n";
print MERGE "mv $prefix.merged.fasta2 $prefix.merged.fasta && \ \n";
print MERGE "grep \">\" $genome1 | tr -d \">\" > $genome1.headers && \ \n";
print MERGE "grep \">\" $genome2 | tr -d \">\" > $genome2.headers && \ \n";
print MERGE "fgrep -f ref.include -w -v $genome1.headers > $genome1.no_overlaps && \ \n";
print MERGE "fgrep -f qry.include -w -v  $genome2.headers > $genome2.no_overlaps && \ \n";
print MERGE "$garmbin/grep_fasta.l -l $genome1.no_overlaps $genome1 > extras.fasta && \ \n";
print MERGE "$garmbin/grep_fasta.l -l $genome2.no_overlaps $genome2 >> extras.fasta\n";
close MERGE;

###
}

sub scaffolding {

my $scf1 = shift;
$scf1 .= ".agp";
my $pwd = shift;
my $prefix = shift;
chdir ("$pwd/$prefix/scaffolding");
open LINK, ">GARM.linkfiles.sh";
open SCF, ">GARM.scaffolding.sh";
print LINK "ln -s $pwd/$prefix/merging/ref.include . && \ \n";
print LINK "ln -s $pwd/$prefix/merging/qry.include . && \ \n";
print LINK "ln -s $pwd/$prefix/merging/$prefix.merged.fasta . && \ \n";
print LINK "ln -s $pwd/$prefix/merging/$prefix.singletons.fasta . && \ \n";
print LINK "ln -s $pwd/$prefix/merging/$prefix.read.placed . && \ \n";
print LINK "ln -s $pwd/$prefix/merging/extras.fasta\n";
print SCF "$garmbin/filter_read.placed_GARM.pl $prefix.read.placed && \ \n";
print SCF "$garmbin/link_assembly_GARM_v0.5.pl $prefix.merged.fasta $prefix.read.placed.filtered $scf1 extras.fasta && \ \n";
print SCF "$garmbin/generate_fasta_scf_v0.3.pl $prefix.merged.fasta $prefix.singletons.fasta extras.fasta new_scaffolds.txt && \ \n";
print SCF "$garmbin/scaff_split.pl new_scaffolds.txt.fasta && \ \n";
print SCF "$garmbin/reoverlap.pl new_scaffolds.txt.fasta.contigs_in_scaffs.fasta new_scaffolds.txt.fasta.agp\n";


close LINK;
close SCF;
	
}

sub retrievemerged {

my $pwd = shift;
my $prefix = shift;
open RET, ">GARM.retrieve.merged.sh";
print RET "$garmbin/retrieve.pl $pwd $prefix && \ \n";
print RET "cat $pwd/$prefix/merging/extras.fasta $pwd/$prefix/merging/$prefix.singletons.fasta > $pwd/$prefix/final_bin.fasta\n";
print RET "echo -e \"GARM merged contigs:\\n\" > $pwd/$prefix/final.stats; $garmbin/garm_stats $pwd/$prefix/final_merged*_contigs.fasta >> $pwd/$prefix/final.stats\n";
print RET "echo -e \"GARM bin:\\n\" >> $pwd/$prefix/final.stats; $garmbin/garm_stats $pwd/$prefix/final_bin.fasta >> $pwd/$prefix/final.stats\n";

close RET;
}

sub retrievescaff {
my $pwd = shift;
my $prefix = shift;
my $path= "*.leftovers";
if( grep {-s $_} glob($path) ){
    open RET, ">GARM.retrieve.scaff.sh";
    print RET "$garmbin/retrieve.pl $pwd $prefix && \ \n";
    print RET "cat $pwd/$prefix/scaffolding/final_leftovers_bin.fasta $pwd/$prefix/input/*.leftovers > $pwd/$prefix/final_bin.fasta\n";
    print RET "echo -e \"GARM reconstructed scaffolds:\\n\" > $pwd/$prefix/final.stats; $garmbin/garm_stats $pwd/$prefix/final_merged*_scaffolds.fasta >> $pwd/$prefix/final.stats\n";
    print RET "echo -e \"GARM merged contigs:\\n\" >> $pwd/$prefix/final.stats; $garmbin/garm_stats $pwd/$prefix/final_merged*_contigs.fasta >> $pwd/$prefix/final.stats\n";
    print RET "echo -e \"GARM bin:\\n\" >> $pwd/$prefix/final.stats; $garmbin/garm_stats $pwd/$prefix/final_bin.fasta >> $pwd/$prefix/final.stats\n";

    close RET;
}else{
    open RET, ">GARM.retrieve.scaff.sh";
    print RET "$garmbin/retrieve.pl $pwd $prefix && \ \n";
    print RET "cat $pwd/$prefix/scaffolding/final_leftovers_bin.fasta > $pwd/$prefix/final_bin.fasta\n";
    print RET "echo -e \"GARM reconstructed scaffolds:\\n\" > $pwd/$prefix/final.stats; $garmbin/garm_stats $pwd/$prefix/final_merged*_scaffolds.fasta >> $pwd/$prefix/final.stats\n";
    print RET "echo -e \"GARM merged contigs:\\n\" >> $pwd/$prefix/final.stats; $garmbin/garm_stats $pwd/$prefix/final_merged*_contigs.fasta >> $pwd/$prefix/final.stats\n";
    print RET "echo -e \"GARM bin:\\n\" >> $pwd/$prefix/final.stats; $garmbin/garm_stats $pwd/$prefix/final_bin.fasta >> $pwd/$prefix/final.stats\n";

    close RET;
    
}

}

1;
__END__

