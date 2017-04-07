#!/usr/bin/env perl

use strict;
use Parallel::ForkManager;
my $garmbin = $ENV{GARMBIN};
my $mumbin = $ENV{MUMBIN};
my $amosbin = $ENV{AMOSBIN};
unless (@ARGV) {
	die "Usage: splimer_fork.pl <db.fasta> <qry.fasta> <num CPUS to use> <minmatch> <filter_identity> <ovl_identity> <maxtrim>\n";
}

my $qry = $ARGV[1];
my $db = $ARGV[0];
my $cpus =$ARGV[2]; ;
my $param = $ARGV[3];
my $conide = $ARGV[4];
my $ovlide = $ARGV[5];
my $maxtrim = $ARGV[6];

# my $agp = (split /\_/, $db)[0].".agp"; # karel
my @split= split('\.', $db);
pop(@split) for 1..2;
my $agp =join('.',@split).'.fasta.agp';

my $forcheck = "../input/$agp splimer.coords splimer.notAnnotated.coords $param $maxtrim $conide $ovlide $db $qry";

my $chunk = 10000000; ## split size in bases (10Mb default)
my $prefix = "splimer";
my $bin = $garmbin;
#my $splitter = "$garmbin/split_fasta.pl";
my $splitter = "$garmbin/fastn_split.py";

### CHECK THE STATS OF THE DB ###
my $total = `$garmbin/garm_stats $db | head -1 | awk '{print \$3}' | tr -d ,`;
if ($total < 80000000) {
	$chunk = int($total/10);
} 

mkdir ("splits") || die "Failed to create splits directory: $!\n" unless -d "splits";
#system ("$splitter $db; wait; mv $db.001 $db.1; mv $db.* splits/");
system ("$splitter $db split. $chunk; wait; rename 's/.fasta//' split.*.fasta; mv split.* splits/");

my $ls = `ls splits/split.*`;
die "The split directory is empty... something wrong\n" unless $ls;
my @list = split (/\n/, $ls);

my $split = "";
if (scalar @list > 1) {
	my $pm = new Parallel::ForkManager($cpus);
	foreach my $file (@list) {
		chomp $file;
		$split = (split /\./, $file)[$#_];
		$pm->start and next;
			system("$mumbin/nucmer --maxmatch --minmatch $param --mincluster $param --prefix=$qry.vs.$split $file $qry && $mumbin/show-coords -H -c -l -o -r -I $conide $qry.vs.$split.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'CONTAIN|IDENTITY' > $qry.vs.$split.coords && $mumbin/show-coords -H -c -l -o -r -I $ovlide $qry.vs.$split.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'END|BEGIN' >> $qry.vs.$split.coords && $mumbin/show-coords -H -c -l -o -r -I $conide $qry.vs.$split.delta | egrep 'CONTAIN|IDENTITY' > $qry.vs.$split.notAnnotated.coords && $mumbin/show-coords -H -c -l -o -r -I $ovlide $qry.vs.$split.delta | egrep 'BEGIN|END' >> $qry.vs.$split.notAnnotated.coords; wait");
		$pm->finish;
	}
	my $last = (split /\./, $list[$#list])[$#_];	
#	system("sleep 5m; cat $qry.vs.split*.coords > $prefix.coords && cat $qry.vs.split*$prefix.notAnnotated.coords > $prefix.notAnnotated.coords && $garmbin/check_coords_annotated.pl $forcheck && sh rerun_nucmer_ref.sh; wait; $garmbin/filter_GARM_coords_v0.3.pl splimer.coords && rm $qry.vs.*.coords $qry.vs.*.delta; rm -rf splits");
	system("sleep 2m; cat $qry.vs.[0-9]*.coords > $prefix.coords && cat $qry.vs.*.notAnnotated.coords > $prefix.notAnnotated.coords && $garmbin/check_coords_annotated.pl $forcheck && sh rerun_nucmer_ref.sh; wait; $garmbin/filter_GARM_coords_v0.3.pl splimer.coords; wait; rm $qry.vs.*.coords $qry.vs.*.delta; rm -rf splits ");
} else {
	print STDERR "Only one split... Not running with fork...\n";
	chomp $list[0];
	my $file = $list[0];
	$split = (split /\./, $file)[$#_];
	system("$mumbin/nucmer --maxmatch --minmatch $param --mincluster $param --prefix=$qry.vs.$split $file $qry && $mumbin/show-coords -H -c -l -o -r -I $conide $qry.vs.$split.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'CONTAIN|IDENTITY' > $prefix.coords && $mumbin/show-coords -H -c -l -o -r -I $ovlide $qry.vs.$split.delta | $amosbin/nucmerAnnotate -ignore $maxtrim | egrep 'END|BEGIN' >> $prefix.coords && $garmbin/filter_GARM_coords_v0.3.pl splimer.coords && $mumbin/show-coords -H -c -l -o -r -I $conide $qry.vs.split.\${LSB_JOBINDEX}.delta | egrep 'CONTAIN|IDENTITY' > $qry.vs.split.\${LSB_JOBINDEX}.$prefix.notAnnotated.coords && $mumbin/show-coords -H -c -l -o -r -I $ovlide $qry.vs.split.\${LSB_JOBINDEX}.delta | egrep 'BEGIN|END' >> $qry.vs.split.\${LSB_JOBINDEX}.$prefix.notAnnotated.coords; wait");
	system("sleep 2m; cat $qry.vs.split*.coords > $prefix.coords && cat $qry.vs.split*$prefix.notAnnotated.coords > $prefix.notAnnotated.coords && $garmbin/check_coords_annotated.pl $forcheck && sh rerun_nucmer_ref.sh; wait; $garmbin/filter_GARM_coords_v0.3.pl splimer.coords; wait; rm $qry.vs.*.coords $qry.vs.*.delta; rm -rf splits ");
}



