#!/usr/bin/env perl

###################
# Genome Assembly #
# Reconciliation  #
# and Merging     #
###################

############## 
use strict;
use Getopt::Std;
use lib "$ENV{GARMLIB}";
use GARM_main;
##############

my $version = "0.7.5";
my $random = rand(1000);

### ENV PATHS ###
my $pwd = $ENV{PWD};
my $garmbin = $ENV{GARMBIN};
my $mumbin = $ENV{MUMBIN};
#################
## DETECT IF RUNNING IN SANGER CLUSTER
my $host = `hostname`;
chomp $host;
my $sanger = "";
#################

my $splimer = "$garmbin/splimer_fork.pl";
my $cpus = 2;
my $stdevfactor =1;

#############

### SOME GLOBAL DEFAULT VALUES ###
my $filteride = 97; #contain identity value
my $overlapide = 90; #identity for the overlaps at the ends of contigs (BEGIN | END) # 99 karel
my $conerror = 0.03;
my $maxtrim = 100;
my $prefix = "";
my $mlength = 200; #default
my $match = 200;
my $repeatend = $match;
my $repeatsfile = "";
my $binsize = 0;
my $extrahelp = 0;
my $gcplot = 0;
my $q = ""; ## ONLY FOR LFS QUEUE SYSTEMS...
### GET OPTIONS ###
my %opts;
getopts('hsxm:c:g:o:l:t:w:', \%opts);
###################

GARM_main::usage($version, 1) if $opts{h};

unless ($opts{g} and $opts{o}) {
	print "No info_genomes_file or not output prefix were specified...\n\n";
	GARM_main::usage($version, 0);
}

$prefix = $opts{o};
if ($opts{l}) {
	$mlength = $opts{l};
}

open LOG, ">GARM_$prefix.log";
if ($opts{p} and $opts{f}) {
	warn "Options -p and -f are not compatible!!! Choose one or the other\n";
	GARM_main::usage($version, 1);
}

print LOG "### OPTIONS ###\n\n";
print LOG "Prefix:\t$prefix\n";

$prefix = $opts{o};
if ($opts{l}) {
	$mlength = $opts{l};
}
print LOG "Minimal contig length:\t$mlength\n";

if ($opts{m}) {
	$match = $opts{m};
}
if ($match > $mlength) {
	print LOG "The length of the contig is shorter than the min match parameter\n";
	print LOG "Adjustting match parameter to $mlength\n"; #<STDIN>;
	$match = $mlength;
}
print LOG "Minimal match length:\t$match\n";


if ($opts{c}) {
	$cpus = $opts{c};
	print LOG "Using $garmbin/splimer_fork.pl with $cpus threads\n";
}

mkdir ("$pwd/$prefix") || die "Failed to create $pwd/$prefix: $!\n" unless -d "$pwd/$prefix";
mkdir ("$pwd/$prefix/input") || die "Failed to create $pwd/$prefix/input: $!\n" unless -d "$pwd/$prefix/input";

if ($opts{s}) {
	$stdevfactor = $opts{s};
}
if ($opts{t}) {
	$maxtrim = $opts{t};
}
if ($opts{x}) {
	$gcplot = 1;
	if ($opts{w}) {
		$binsize = $opts{w};
	}
}

my $start = `date`;
my $time = "";
chomp $start;
my $genomesfile = $opts{g};
my $stats = "";

print LOG "\nGARM version $version started: $start\n";
print LOG "1.- PROCESSING INPUT FROM $genomesfile ###\n"; #<STDIN>;
print "GARM version $version started: $start\n";
print "1.- PROCESSING INPUT FROM $genomesfile ###\n"; #<STDIN>;
my @inputgenomes = GARM_main::loadinput("$genomesfile");
my ($inputlog, $getqcgenomes, $getscffiles) = GARM_main::input(\@inputgenomes, $pwd, $prefix, $mlength, $repeatend, $repeatsfile, $binsize, $gcplot);
my @qcgenomes = @$getqcgenomes;
my @scffiles = @$getscffiles;

$time = `date`;
if ($inputlog =~ /INPUT PART - DONE/) {
	print "PROCESSING INPUT - DONE! $time\n";
} else {
	die "PROCESSING INPUT - FAILED...Please check the $pwd/$prefix/input directory\n";
}

print LOG "$inputlog\n";
my $scf1 = "";
my $scf2 = "";
print LOG "2.- CHECKING IF SCAFFOLDING INFORMATION IS AVAILABLE ###\n";
print "2.- CHECKING IF SCAFFOLDING INFORMATION IS AVAILABLE ###\n";
if (scalar @scffiles == 2) {
	print LOG "### TWO SCAFFOLDING FILES... SELECTING THE BEST ###\n";
	print "TWO SCAFFOLDING FILES... SELECTING THE BEST\n";
	($scf1, $scf2, $stats) = GARM_main::genomestats($scffiles[0], $scffiles[1]);
	print LOG "$stats\n";
} elsif (scalar @scffiles ==1) {
	print LOG "### SCAFFOLDING FILE AVAILABLE ###\n";
	print "SCAFFOLDING FILE AVAILABLE\n";
	$scf1 = $scffiles[0];
# 	my $stats = "$garmbin/garm_stats $scf1"; karel
# 	print LOG "Using scaffolds from as ref\n$stats\n"; karel
}

$time = `date`;
chomp $time;
print LOG "3.- RUNNING MERGING PART $time ###\n"; #<STDIN>;
print "3.- RUNNING MERGING PART $time ###\n"; #<STDIN>;
if (scalar @qcgenomes == 2) {
	my $genome1 = $qcgenomes[0];
	my $genome2 = $qcgenomes[1];
	print LOG "## CALCULATING STATS FOR THE GENOMES ##\n"; #<STDIN>;
	if ($scf1) {
		if ($genome2 =~ /$scf1/) {
			my $fake1 = $genome2;
			my $fake2 = $genome1;
			($fake1, $fake2, $stats) = GARM_main::genomestats($fake1, $fake2);
			print LOG "Stats of contigs before merging:\n$stats\n";
		} else {
			($genome1, $genome2, $stats) = GARM_main::genomestats($genome1, $genome2);
			print LOG "Stats of contigs before merging:\n$stats\n";
		}
	} else {
		($genome1, $genome2, $stats) = GARM_main::genomestats($genome1, $genome2);
		print LOG "Stats of contigs before merging:\n$stats\n";
	}
	chdir("$pwd/$prefix/merging/");
	print LOG "## RUNNING SPLIMER ##\n"; #<STDIN>;
	print "RUNNING SPLIMER ##\n"; #<STDIN>;
	GARM_main::splimer($genome1, $genome2, $match, $filteride, $overlapide, $splimer, $cpus, $maxtrim, $random);
	#print LOG "## REMOVING REPEATS AND PROBABLE REPEATS ##\n";
	#print "REMOVING REPEATS AND PROBABLE REPEATS ##\n";cat /fre
	#my $agp = (split /\_/, $genome1)[0].".agp";
	my @agpname = split (/\_/, $genome1);
	splice(@agpname, -2, 2);
	my $agp = join("_", @agpname).".agp";
	my $forcheck = "$pwd/$prefix/input/$agp splimer.coords splimer.notAnnotated.coords $match $maxtrim $filteride $overlapide $genome1 $genome2";
	print LOG "## RUNNING O-L-C PIPELINE ##\n";
	print "RUNNING O-L-C PIPELINE ##\n";
	chdir ("$pwd/$prefix/merging");
	GARM_main::merging($genome1, $genome2, $maxtrim, $prefix);
	# GARM_main::retrievemerged($pwd, $prefix);
	### SCAFFOLDING ###
	if ($scf1) {
		GARM_main::scaffolding($scf1, $pwd, $prefix, $stdevfactor);
		GARM_main::retrievescaff($pwd, $prefix);
		print LOG "## RUNNING GARM.merging.sh script and GARM.scaffolding.sh...\n";
		print LOG "sh GARM.merging.sh; wait; cd $pwd/$prefix/scaffolding; sh GARM.linkfiles.sh && sh GARM.scaffolding.sh\n";
		chdir ("$pwd/$prefix/merging");
		system("sh GARM.merging.sh; wait");
		system("$garmbin/check_coords_annotated.pl $forcheck; wait; sh rerun_nucmer_ref.sh; wait; $garmbin/filter_GARM_coords_v0.3.pl splimer.coords; wait");
		system("cd $pwd/$prefix/scaffolding; sh GARM.linkfiles.sh; wait; sh GARM.scaffolding.sh; wait; sh GARM.retrieve.scaff.sh; wait");
		$time = `date`;
		chomp $time;
		print LOG "\nGARM version $version finished: $time\n";
	} else {
        GARM_main::retrievemerged($pwd, $prefix);
		print LOG "## RUNNING GARM.merging.sh script...\n";
		print LOG "sh GARM.merging.sh && sh GARM.retrieve.merged.sh\n";
		chdir ("$pwd/$prefix/merging");
		system("sh GARM.merging.sh; wait; sh GARM.retrieve.merged.sh");
		system("$garmbin/check_coords_annotated.pl $forcheck; wait; sh rerun_nucmer_ref.sh; wait; $garmbin/filter_GARM_coords_v0.3.pl splimer.coords; wait");
		system("sh GARM.retrieve.merged.sh");
		$time = `date`;
		chomp $time;
		print LOG "\nGARM version $version finished: $time\n";
	}
	
}
