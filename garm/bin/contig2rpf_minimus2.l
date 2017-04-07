#!/usr/bin/env perl
# takes a '.contig' file and creates an 'rpf' file with some fields missing
#use lib "/home/alejandr/perl"; 
#use yo qw(:DEFAULT); 
use strict; 
my ($contig); 
@ARGV == 1 || die "contig2rpf.l <algo.contig>\n"; 

#CTGSCF();

while(<>){
    if( /^##(contig\d+)/i || /^##(\d+)/ ){
        $contig="$1";
		if ($contig !~ /contig/) {
			$contig = "merged_contig$contig";
		} 
        next; 
    }
    next unless /^#([^(]+).*checksum.\s+\{(\d+)\s+(\d+)\}\s+\<(\d+)\s+(\d+)\>/; 
    my $name=$1; 
    my $read_from=$2; 
    my $read_to=$3; 
    my $contig_from=$4; 
    my $contig_to=$5; 
    my $strand; 
    if( $read_to<$read_from ){
        $strand=1; 
        my $tmp=$read_from; 
        $read_from=$read_to; 
        $read_to=$tmp; 
    }else{
        $strand=0; 
    }
    my $len= $read_to -$read_from; 
    
    print "* $name $read_from $len $strand $contig * $contig_from *\n"; 
}


#sub CTGSCF {

#my $path = "/home/alexsf/CA_results/taenia_assemble_sanger050209_bicholike/9-terminator/";




#}
