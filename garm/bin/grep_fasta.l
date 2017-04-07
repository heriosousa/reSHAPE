#!/usr/bin/env perl
use lib "$ENV{GARMLIB}"; 
use bio; 

use Getopt::Std;
my(%opts); 

getopts("hosql:", \%opts); 
help() if $opts{'h'}; 
my $regex=shift unless $opts{'l'}; 

die "grep_fasta.l: las opciones -q y -s no son compatibles\n" if $opts{'s'} && $opts{'q'}; 
my @files=@ARGV; 
@files >=1 or help(); 

my $cnt=0; 

if( $opts{'l'} ){
    die "grep_fasta.l: las opciones -l y -s no son compatibles\n" if $opts{'s'}; 
    w_list(); 
}else{
    w_regex(); 
}
print STDERR "$cnt sequence found\n"; 

    
sub w_regex{
    foreach my $file ( @files ){
        if( $opts{q} ){
            tie *FASTA, 'qual', "$file";  
        }else{
            tie *FASTA, 'fasta', "$file";  
        }
        
        
        if( $opts{'s'} ){ # busca en la secuencia
            while(my($name, $seq)= <FASTA>){
                next unless $$seq=~ m{$regex}; 
                print ">$name\n$$seq\n"; 
                exit if $opts{'o'}; 
                $cnt++; 
            }
        }else{
            while(my($name, $seq)= <FASTA>){
                next unless $name=~ m{$regex}; 
                print ">$name\n$$seq\n"; 
                exit if $opts{'o'}; 
                $cnt++; 
            }
        }
    }
}

sub w_list{
    my $IN = sopen("< $opts{'l'}", 0);  
    my %wanted; 
    while(<$IN>){
        next if /^#/; 
        next if/^\s*$/; 
        my($name)=(split); 
        $wanted{$name}++; 
    }
    foreach my $file ( @files ){
        if( $opts{q} ){
            tie *FASTA, 'qual', "$file";  
        }else{
            tie *FASTA, 'fasta', "$file";  
        }
        while(my($name, $seq)= <FASTA>){
            next unless $wanted{$name}; 
            print ">$name\n$$seq\n"; 
            exit if $opts{'o'}; 
            $cnt++; 
        }
    }
}
sub sopen{
# $fh=sopen($filename, [$verbose]) returns a filehandle_rference
    my $filename = shift;
    my $in = do { local *IN };
    open $in, $filename or die "$filename: $!";
    print STDERR "myopen $filename\n" if $_[0]; 
    return $in;
}


sub help{
die <<"aqui"; 
grep_fasta.l [options] <regex> <fasta_files>
options:
-h         This help
-s         Search in sequence, not in seq_names
-o         Exit after first match
-q         Fasta are quals;  not compatible with -s
-l <file>  Take sequence names from 'file'; names should be exact; 
           incompatible with -s; 
           
Note: seq_names are as this (\$seq_name)=~ m/^>(\\S+)/
aqui
; 

}
