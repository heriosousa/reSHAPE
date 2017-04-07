# Mis funciones mas comunes
package bio; 
use Exporter (); 
@ISA = qw (Exporter); 
@EXPORT= qw (); 
@EXPORT_OK = qw (reverse_seq make_meres balance_idx frec_meres); 

############## bio #################

push  @EXPORT_OK, 'revierte_seq'; 
sub revierte_seq{
    my $nucs= ref($_[0]) ? scalar(reverse(${$_[0]} ) ) : scalar(reverse($_[0])); 
    $nucs=~ tr/ACGTacgt/TGCAtgca/;
    return \$nucs; 
}

push  @EXPORT_OK, 'sixway'; 
sub sixway{
# $prots=sixway($dna) returns a ref to an array with the 6 way
# translations of dna. dna can be scalar or ref_to_scalar. STOP codons
# will translate to 'X'
    my @prot; 
    if( ref $_[0] ){
        for(my $i=0; $i< 3; $i++){
            push @prot, code2aa( $_[0], $i); 
        }
        my $rev=reverse(${$_[0]}); 
        $rev=~ tr/ACGTacgt/TGCAtgca/; 
        for(my $i=0; $i< 3; $i++){
            push @prot, code2aa( \$rev, $i); 
        }
        
    }else{
        for(my $i=0; $i< 3; $i++){
            push @prot, code2aa( \$_[0], $i);
        }
        my $rev=reverse($_[0]); 
        $rev=~ tr/ACGTacgt/TGCAtgca/; 
        for(my $i=0; $i< 3; $i++){
            push @prot, code2aa( \$rev, $i); 
        }
    }
    for (@prot){ tr/*/X/ }; 
    return \@prot; 
    
}

push @EXPORT_OK, 'indel_dup_DNA'; 
sub indel_dup_DNA ($$$) {
# $newseq=indel_dup_DNA($seq, $num, $max_size) return a mutated seq with
# num macro_mutations,  whose average size is $max_size/2 and where frec
# of del ~= frec (ins + dupl). if $max_size is <1,  $max_size is set to
# $len * $max_size. $max_size can never be larger than $len/3. To
# prevent destroying the seq, some bundaries apply

    my($seq, $num, $max_size)=@_; 
    my $len0 = length($seq); 
    $max_size= $max_size <1 ? $max_size * $len0 : $max_size; 
    $max_size= $max_size> $len0 /3 ? int($len0 /3) : $max_size; 
    
    for(my $i=0; $i< $num; $i++){
        my $len = length($seq); 
        
        my $pos = int(rand($len)); 
        my $range=int(rand($max_size)); 
        my $case; 
        if( $len >$len0 ){        # forza delete
            $case=int(rand(2)); 
        }elsif( $len<$len0 ){     # impide delete
            $case=int(rand(2)) +2; 
        }else{                    # el mismo chance de crecer o disminuir
            $case=int(rand(4)); 
        }
        if( $case <2){                          # delete
            substr($seq, $pos, $range)=''; 
        }elsif( $case <3 ){                     # insert
            substr($seq, $pos, 0)=rand_DNA($range); 
        }else{                                  # duplicate
            substr($seq, $pos, 0)=substr($seq, $pos, $range); 
        }
    }
    return $seq; 
}

push @EXPORT_OK, 'muta_DNA_pam'; 
sub muta_DNA_pam ($$){
# $new_seq=  muta_DNA_pam($seq, $pam) returns a muated vesrion of seq
# which distance is aprox PAM
    my($seq, $pam)=@_; 
    my $len = length($seq); 
    my @nucs=qw( A C G T); 
    for(my $i=0; $i< $len/100 * $pam; $i++){
        my $pos = rand($len); 
        my $nuc= $nucs[rand(4)]; 
        substr($seq, $pos, 1)=$nuc; 
    }
    return $seq; 
}

push @EXPORT_OK, 'rand_DNA'; 
sub rand_DNA ($){
# $seq=rand_DNA($len) return a random DNA sequence of desired len
    my($len)=shift; 
    my @nucs=qw( A C G T); 
    my $seq=''; 
    for(my $i=0; $i< $len; $i++){
        $seq.= $nucs[rand(4)]; 
    }
    return $seq; 
}

push @EXPORT_OK, 'one_let_codon'; 
sub one_let_codon{
# %hash=one_let_codon() returns a hash with one letter_names for codons
    return qw(
AAA K AAC n AAG k AAT N ACA 2 ACC t ACG 3 ACT T AGA 4 AGC 7 AGG 5 AGT 6 ATA ! 
ATC i ATG M ATT I CAA Q CAC h CAG q CAT H CCA O CCC p CCG o CCT P CGA 8 CGC r 
CGG 9 CGT R CTA J CTC b CTG j CTT B GAA E GAC d GAG e GAT D GCA @ GCC a GCG & 
GCT A GGA : GGC g GGG = GGT G GTA U GTC v GTG u GTT V TAA % TAC y TAG $ TAT Y 
TCA Z TCC s TCG z TCT S TGA * TGC c TGG W TGT C TTA L TTC f TTG l TTT F --- -     
    ); 
}

push @EXPORT_OK, 'hori2vert_aln'; 
sub hori2vert_aln{
# \@cols= hori2vert_aln(\@rows) verticalizes a block of text switching
# rows to cols. @rows is a unidimensional arrar of strings of equal
# lengths. It returns a ref to an array of newline terminated strings
# which are the cols of the original block

    my $lines=shift; 
    my(@rows, @cols); 
    foreach my $line ( @$lines ){
        chomp $line; 
        push @rows, [split('', $line)]; 
    }
    for(my $i=0; $i< @{$rows[0]}; $i++){
        for(my $j=0; $j< @rows; $j++){
            $cols[$i].=$rows[$j][$i]; 
        }
         $cols[$i].="\n"; 
    }
    return \@cols; 
}

push @EXPORT_OK, 'nuc_coords'; 
sub nuc_coords{
# my($gene, $strand, $start, $end)=nuc_coords($query, $from, $to, $nuc_len)
# Given the query from a blat hit at the protein level,
# returns the coordinates in the nucleotide CDS. It's merit is that it
# takes into acount the phase in which the CDS was read to produce the
# query protein. Here $query has to make a match like this:
# $query=~/_[fr][012]$/; 
    my ($query, $from, $to, $len)=@_; 
    my($gene, $start, $end, $strand, $fase); 
    $query=~/^(\S+)_([fr])([012])$/ or die "Bad query name:$query\n"; 
    $gene=$1; 
    $strand=$2; 
    $fase=$3; 
    my $s= ($from *3) +$fase; 
    my $e=  ($to *3) +$fase; 
    if( $2 eq 'r' ){
        $start=$len- $e ; 
        $end=  $len -$s; 
    }else{
        $start= $s; 
        $end=   $e; 
    }
    return $gene, $strand, $start, $end; 
}

push @EXPORT_OK, 'contigs'; 
sub contigs{
# my($left, $right, $members)=contigs(\@start, \@end [, $gap]) joins
# fragments in contigs if their distance is <= $gap. If you want true
# contig, make $gap negative. It return refs to 3 parallel arrays
# describing the contigs: left of the contig,  right of the contig, and
# the members in the contig. Members are shown as string of ':'
# separated numbers. The numbers are the indices of the fragments as
# present in the input arrays.
    my($start, $end, $gap)=@_; 
    my($l, $r, $mbs, @l, @r, @mbs, @pos); 
    
    $gap||=0; 
    my $mi_sort= sub {
        return $start->[$a] <=> $start->[$b]
        ||
        $end->[$a] <=> $end->[$b]; 
    }; 
    
    @pos= sort $mi_sort  (0..$#{$start}); 
    
    ($l, $r)=($start->[$pos[0]], $end->[$pos[0]]); 
    $mbs="0"; 
    for(my $i=1; $i< @pos; $i++){
        my $pos=$pos[$i]; 
        if( $r+$gap >= $start->[$pos] ){
            $mbs.=":$pos"; 
            $r=$end->[$pos] if  $end->[$pos] >$r; 
        }else{ 
            push @l, $l; 
            push @r, $r; 
            push @mbs, $mbs; 
            $l=$start->[$pos]; 
            $r=$end->[$pos]; 
            $mbs="$pos"; 
        }
    }
    push @l, $l; 
    push @r, $r; 
    push @mbs, $mbs; 
    return \@l, \@r, \@mbs; 
    
}


push @EXPORT_OK,  'read_fasta'; 
sub read_fasta{
# ($seq, $names)=read_fasta(<fasta_file>, [<regex_to_catch_name>]) reads
# all sequences in afasta file and returns a ref to hash with the
# sequences,  and a ref to an array with the names. Note that the regex
# must match in the heading line and must contain parentesis to catch
# the name
# Nota: <fasta_file> puede ser un refa un array de nombres; 
    my ($file, $rgx)=@_;  
    $rgx||='^>(.+)';    # default: take all line
    my(@names, %seq, $name, %ya, @files); 
    
    @files= ref($file) ? @$file : $file; 
    foreach my $file_to_open ( @files ){
        open IN, $file_to_open or die "Cant read $file_to_open\n"; 
        while(<IN>){
            if(/$rgx/){ 
                $name=$1;  
                $ya{$name}++ && die "read_fasta: repeated name: $name\n"; 
                push @names, $name;  
            }else{ 
                chomp;  
                $seq{$name}.=$_; 
            }
        }
    }
    
    return \%seq, \@names; 
}

push @EXPORT_OK,  'overlap'; 
sub overlap{
# ($size, $from, $to)=overlap($f1, $t1, $f2, $t2) returns the size and
# the ends of the overlap of both ranges. If there is no overlap,  both
# ends are -1
    my($a, $b, $c, $d)=@_; 
    if($b<$a ||$d<$c){
        warn "ranges are wrong: $a $b $c $d\n"; 
        return (-1, -1, -1); 
    }
    my $from= $a>$c ? $a : $c; 
    my $to=   $b<$d ? $b : $d; 
    my $len=  $to-$from +1; 
    return ($len<=0) ? ($len, -1, -1) : ($len, $from, $to); 
}


push @EXPORT_OK,  'cnt_codons'; 
{
my @codons; 
sub cnt_codons{
# @codon_cnt= cnt_codons($orf) return in an array with the count of the
# 64 codons in the $orf sequence. The order of codons is alfabetical.
# You can generate the list with: @codons =bio::make_meres(3);


    my ($orf)=lc(shift); 
    my %cnt; 
    @codons=make_meres(3) unless @codons; 
    @cnt{@codons}=(0) x 64; 
    while($orf=~/(...)/g){
        $cnt{$1}++; 
    }
    return @cnt{@codons}; 
}
}

push @EXPORT_OK,  'indentity'; 
sub identity{
# $ident=identity($seq1, $seq2) compares $seq1 and $seq2 and returns the
# percentage of identities. Initial and final gaps are ignored.
# Positions where one seq has gap '-' but not the other are counted as
# negative cases. When both have '-',  position is ignored.
    my($seq1, $seq2)=@_; 
    my($i, $j, $k)=(0) x3; 
    $seq1=~s/-+$//; 
    $seq2=~s/-+$//; 
    # skip initial gaps
    while(my $c1= substr($seq1, $i, 1) and my $c2= substr($seq2, $i, 1) ){
        $i++; 
        next if $c1 eq '-' or $c2 eq '-'; 
        last; 
    }
    
    $i--; 
    while(my $c1= substr($seq1, $i, 1) and my $c2= substr($seq2, $i, 1) ){
        $i++; 
        next if $c1 eq '-' && $c2 eq '-'; 
        $j++; 
        $k++ if $c1 eq $c2; 
    }
    return $j?100*$k/$j:undef; 
}

push @EXPORT_OK, 'rand_peptide';
{
my($order, $pie); 
sub rand_peptide{
# $prot=rand_peptide($len) returns a random peptide of length $len, 
# generated using the average aa. frecuencies from all genomes in the
# database
    my($len)=shift; 
    unless(defined $order){
        my(@valid)=qw( A C D E F G H I K L M N P Q R S T V W Y); 
        my %u_frec; 
        eval q{ use store_db; } ; 
        @u_frec{@valid}=store_db::universal_frec(); 
        ($order, $pie)=yo::pie(\%u_frec)
    }
    return join('', yo::pie_hit($order, $pie, $len))
    
}
eval q{ no store_db; } ;
}
push @EXPORT_OK, 'switch_pamize'; 
# $seq=switch_pamisize($seq, $pam) switches 2 randomly chosen residues
# from $seq. The process is repeated to achive $pam point accepted
# mutations. Note: since each cicle changes 2 residues, the resulting pam
# is the best aproximation smaller than the desired pam
sub switch_pamize{
    my ($seq, $pam)=@_; 
    my @array=split('', $seq); 
    my $len=@array; 
    my $stop=$pam*$len/200; 
    return $seq if $len <2; 
    for(my $i=0; $i<$stop; ){
        my $k=int rand($len); 
        my $l=int rand($len); 
        next if $l==$k; 
        @array[$k,$l]=@array[$l,$k]; 
        $i++; 
    }    
    return join('', @array); 

}
push @EXPORT_OK, 'shuffle_seq'; 
sub shuffle_seq{
# $randseq=shuffle_seq($seq|\$seq) returns a shuffled version of $seq.
# Note: you must take care not to include a "\n" at the end of seq or it
# will get randomized too
    my $seq=shift; 
    my @array; 
    if(ref $seq){
        @array=split('', $$seq); 
    }else{
        @array=split('', $seq); 
    }
    yo::shuffle_array(\@array); 
    return join('', @array); 
}

push  @EXPORT_OK, 'faa2hash'; 
sub faa2hash{
# $seqs=faa2hash($faa_file) returns a hash ref containing all the seqs in a
# '.faa' file;

    my $faa_file=shift; 
    my(%seqs); 
    open FAA, $faa_file or die "Cant read $faa_file at bio::faa2vect()\n"; 
    $_=<FAA>; 
    /^>(gi\|)?([A-Za-z0-9._+-]+)/i or die "First line in $faa_file is no fasta-like\n"; 
    my $name=$2; 
    while (<FAA>){
        if( /^>(gi\|)?([A-Za-z0-9._+-]+)/i){
            $name=$2; 
            next; 
        }
        $seqs{$name}.=$_; 
    }
    foreach my $seq (values %seqs){
        $seq=~s/\s+//g; 
    }
    return(\%seqs); 
}

push @EXPORT_OK, 'pdb2seq'; 
sub pdb2seq{
# %seqs=pdb2seq(*FILEHANDLE) return a hash where keys are the names
# of the peptide chains in a "pdb" file, and the values are the
# corresponding sequences. If there is only one chain, the only
# key might be a space character: ' '. 
# Presumably FILEHANDLE is associated to a real pdb file, but may
# also be a pipe from something producing a pdb file.
# Non standard aminoacids are represented by 'X'
    my($fh)=shift; 
    my %three2one=three2one(); 
    my($chain, @aa, $tri, %seq); 
    while(<$fh>){
        last if /^ATOM/; 
        next unless /^SEQRES/; 
        if( /^SEQRES.....(.).......(([A-Z]{3} )+)/ ){
            $chain=$1; @aa=split(" ", $2);
            for $tri (@aa){
                $seq{$chain}.=($three2one{$tri}?$three2one{$tri}:'X'); 
            }
        }else{
#             warn "bad line:$_"
        }
    }
    return %seq; 
}


push @EXPORT_OK, 'align'; 
sub align{
# align(\$seq1, \$seq2) corre bl2seq y regresa un array con los resultados.
# Los 2 ultimos datos del array son referncias a strings que representan
# la primera y segunda secuencias ya alineadas.
# El uso recomendado es asi:

# ($score, $expect, $ident, $positive, $seq1_start, $seq2_start, 
#     $align1, $align2)= align(\$seq1, \$seq2); 

# Score, Expect, Identities y Positive son lo que son
# $seq1_start, $seq2_start son las primeros residuos del las proteinas
# que aparecen en los alineamientos,  ya que bl2seq no pone gaps iniciales

    my ($seq1, $seq2)=@_; 
    open S1, ">/var/tmp/align_seq1"||die "Bio.pm: Cant write to align_seq1\n"; 
    print S1 $$seq1, "\n"; 
    close S1, 
    open S1, ">/var/tmp/align_seq2"||die "Bio.pm: Cant write to align_seq2\n"; 
    print S1 $$seq2, "\n"; 
    close S1, 
    open ALIGN, "bl2seq -p blastp -i /var/tmp/align_seq1 -j /var/tmp/align_seq2|"
        or die "cant run bl2seq\n"; 
    my($score, $expect, $ident, $positive); 
    while(<ALIGN>){
        last if /^Query:/; 
        if(/^ Score =\s+(\d+).* Expect =\s+(.*)/){
            $score=$1, $expect=$2; 
            next; 
        }elsif(/^ Identities =[^(]+\((\d+)\%\), Positives =[^(]+\((\d+)/){
            $ident=$1, $positive=$2; 
        }
    }
    my($seq1_start, $seq2_start, $align1, $align2, $ya1, $ya2); 
    if(/^Query:\s+(\d+)\s+(\S+)\s+\d+/){
        $seq1_start=$1 unless $ya1++; 
        $align1.=$2; 
    }; 
    while (<ALIGN>){
        last if /^(\s+Score =|Lambda\s+K\s+H)/; 
        if(/^Query:\s+(\d+)\s+(\S+)\s+\d+/){
            $seq1_start=$1 unless $ya1++; 
            $align1.=$2; 
            next; 
        }elsif(/^Sbjct:\s+(\d+)\s+(\S+)\s+\d+/){
            $seq2_start=$1 unless $ya2++; 
            $align2.=$2; 
            next; 
        }    
    }
    return ($score, $expect, $ident, $positive, $seq1_start, $seq2_start, \$align1, \$align2); 



}

push @EXPORT_OK, 'codon_list'; 
sub codon_list{
# codon_list() returns a list of codons->aa in alphabetical order. 
# Ambiguous nucleotides are not considered
# Use (codon_list())[0..127] to get uppercase codons, 
# or  (codon_list())[128..255] to get lowercase codons
return qw(
AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T AGA R AGC S AGG R AGT S 
ATA I ATC I ATG M ATT I CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L GAA E GAC D GAG E GAT D 
GCA A GCC A GCG A GCT A GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S TGA * TGC C TGG W TGT C 
TTA L TTC F TTG L TTT F aaa K aac N aag K aat N aca T acc T acg T act T 
aga R agc S agg R agt S ata I atc I atg M att I caa Q cac H cag Q cat H 
cca P ccc P ccg P cct P cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
gaa E gac D gag E gat D gca A gcc A gcg A gct A gga G ggc G ggg G ggt G 
gta V gtc V gtg V gtt V taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
tga * tgc C tgg W tgt C tta L ttc F ttg L ttt F 
); 
}

push @EXPORT_OK, 'three2one'; 
sub three2one{
# %hash=three2one(); 
# returns a hash list where key are 3 upercase-letter AA symbols and values
# are one uppercase letter AA symbols.
    return qw(
    ALA A ARG R ASN N ASP D CYS C GLN Q GLU E GLY G HIS H ILE I
    LEU L LYS K MET M PHE F PRO P SER S THR T TRP W TYR Y VAL V 
    )  ;   
}
push @EXPORT_OK, 'code2aa'; 
{
my %aa; 
my $table='';   
sub code2aa{
# code2aa($seq, [$fase], [$table]) or code2aa(\$seq, [$fase]) translates $seq 
# $seq can be passed as value or reference (faster for long strings); 
# fase can be specified (0|1|2), default=0
# posible incomplete codons at begining or end are ignored; 
# codons with strange symbols translate to 'X', <STOP> to '*'
    my (%aa, $seq, $prot, $fase); 
    
    $_[2]||=1; 
    unless( keys %aa && $_[2] ==$table ){
        $table=($_[2] <23) ? $table :1 ; 
        %aa= ncbi_table($_[2]); 
    }
    
    $seq   = ref($_[0])?$_[0]:\$_[0]; 
    $fase  = defined ($_[1])? '.' x $_[1]:''; 
    $$seq=~/($fase)/g; 
    while ($$seq=~/(...)/g){
        $prot.= defined($aa{$1}) ? $aa{$1} : 'X'; 
    }
    return $prot; 
}
}
push @EXPORT_OK, 'cnt_ntuples'; 
sub cnt_ntuples{
# cnt_ntuples($n, \$string) return a %hash containing the count of all
# n-long words in string
    my($n, $string)=@_; 
    my %cnt; 
    my $len =length($$string) - $n; 
    for( my $i=0; $i<$len; $i++){
        $cnt{substr($$string, $i, $n)}++; 
    }
    return %cnt; 
}

push @EXPORT_OK, 'frec_meres'; 
sub frec_meres{
# frec_meres(merlen, seq) return the frequency of all merlen-long
# nucleotide meres frequncies are sorted alphabeticaly and noramlized to
# the expected random  values (a frec of 1 implies equilibrium) seq must
# be in lower case letters .Non-standard meres are excluded
    my($merelen, $seq)=@_; 
    my ($mere, $i, $seqlen, $total, @meres, @out); 
    @meres=make_meres($merelen); 
    $seqlen=length($seq)-$merelen+1; 
    my %cnt=(); 
    for my $mere (@meres){
        $cnt{$mere}=0; 
    }
    for $i(0.. $seqlen){
        $mere=substr($seq, $i, $merelen);  
        $cnt{$mere}++; 
    }
    grep {push @out, $cnt{$_}; $total+=$cnt{$_}} @meres; 
    @out=map {$_ * @meres /$total} @out; 
    
    return @out; 
}

sub balance_idx{
# balance_idx(@meres) returns the balance index.
    my(@meres)=@_; 
    my ($factor, $balance)=(0, 1); 
    grep {$factor+=$_} @meres; 
    $factor=($#meres+1)/$factor; 
    for (@meres){
        $balance*=$_*$factor; 
    }
    return $balance
}

sub reverse_seq{
# reverse_seq($seq) returns the complementary sequence of the DNA $seq
    my($seq)=shift; 
    $seq=~tr/ACGTacgt/TGCAtgca/; 
    return scalar reverse $seq; 
}


push @EXPORT_OK, 'make_all_words'; 
sub make_all_words{
# @words = make_all_words(\@alphabet, $size) returns an alfabeticaly
# sorted list of all the words of length $size, built with a given
# alphabet. The alphabet should be in the form of an anonymous array
    my($alphabet, $size)=@_; 
    my($len, $symb, $word, @prev, @new); 
    $prev[0]=''; 
    for($len=0; $len<$size; $len++){
        @new=(); 
        for $symb(@{$alphabet}){
            for $word (@prev){
                push(@new, "$symb$word"); 
            }
        }
        @prev=@new; 
    }
    return @prev; 
}

sub make_meres{
# make_meres($len)
# Produce recursivamente la lista de los n-meres de nucleotidos
    my($len)=shift; 
    my(@list, @newlist, $mere, $nuc); 
    $len--; 
    if($len<0){
        return '';  
    }elsif($len==0){
        return qw(a c g t)
    }else{
        @list=make_meres($len); 
    }
    for $mere (@list){
        for $nuc (qw (a c g t)){
            push (@newlist, $mere.$nuc); 
        }
    }
    return @newlist; 
}

push @EXPORT_OK, "make_markov"; 
sub make_markov{
# make_markov($size, $merlen, @frecs) regresa un genoma Markoviano 
# Si falla regresa undef
    my($size, $merlen, @frecs)=@_; 
    my(%ac, %ct)=(); 
    {
# llena el hash %ct con los valores de @frec        
        my(@meres)=make_meres($merlen); 
        return undef if $#meres != $#frecs; 
        for my $mere (@meres){
            $ct{$mere}=shift @frecs; 
        }
    }
    {
# equivale a acum_mat
        my(@meres)=make_meres($merlen-1); 
        my($val, $total, $fst, $snd);
        for $fst (@meres){
            $total=0; 
            for $snd (qw(a c g t)){
                $val=$ct{$fst.$snd};
                $total+=$val;
                $ac{$fst.$snd}=$total;
            }
        }
    }
    {
#genera el markoviano        
        my($res, $rnd, $i, $prev, $genome); 
        my(@meres)=make_meres($merlen-1); 
        $genome='';
        $prev=$meres[0]; 
        for($i=0;$i<$size;$i++){
            $rnd=rand($ac{$prev.'t'});
            if($rnd<$ac{$prev.'a'}){
                $res='a';
            }elsif($rnd<$ac{$prev.'c'}){
                $res='c';
            }elsif($rnd<$ac{$prev.'g'}){
                $res='g';
            }else{
                $res='t';
            }
            $genome.=$res;
            $prev=substr($prev.$res, 1); 
        }
        return $genome; 
    }
}

sub pam_score{
# pam_score($seq1, $seq2) califica el parcido de dos seqs del mismo
# tamanio alineadas, sin gaps. Regresa el score
    my $pam=pam_250(); 
    my(@sec1)=split('', $_[0]); 
    my(@sec2)=split('', $_[1]); 
    my($score, $i, $len)=(0, 0, 0); 
    $len=@sec1<@sec2?@sec1:@sec2; 
    for $i(0..$len-1){
        $score+=$$pam{$sec1[$i]}{$sec2[$i]}; 
    }
    return $score; 
}


sub pam_250{
# genera un hash llamado %pam250 con scope a este package que contiene
# la tabla de pam_250 pa los siguientes simbolos:
# A R N D C Q E G H I L K M F P S T W Y V B Z X
# es util a otras sub de este mismo package,  pero tambien puede ser usada
# como %bio::pam250
    my %pam250=(); 
@{$pam250{A}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(2 -2 0 0 -2 0 0 1 -1 -1 -2 -1 -1 -4 1 1 1 -6 -3 0 0 0 0);
@{$pam250{R}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-2 6 0 -1 -4 1 -1 -3 2 -2 -3 3 0 -4 0 0 -1 2 -4 -2 -1 0 0);
@{$pam250{N}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(0 0 2 2 -4 1 1 0 2 -2 -3 1 -2 -4 -1 1 0 -4 -2 -2 2 1 0);
@{$pam250{D}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(0 -1 2 4 -5 2 3 1 1 -2 -4 0 -3 -6 -1 0 0 -7 -4 -2 3 3 0);
@{$pam250{C}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-2 -4 -4 -5 12 -5 -5 -3 -3 -2 -6 -5 -5 -4 -3 0 -2 -8 0 -2 -4 -5 0);
@{$pam250{Q}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(0 1 1 2 -5 4 2 -1 3 -2 -2 1 -1 -5 0 -1 -1 -5 -4 -2 1 3 0);
@{$pam250{E}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(0 -1 1 3 -5 2 4 0 1 -2 -3 0 -2 -5 -1 0 0 -7 -4 -2 2 3 0);
@{$pam250{G}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(1 -3 0 1 -3 -1 0 5 -2 -3 -4 -2 -3 -5 -1 1 0 -7 -5 -1 0 -1 0);
@{$pam250{H}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-1 2 2 1 -3 3 1 -2 6 -2 -2 0 -2 -2 0 -1 -1 -3 0 -2 1 2 0);
@{$pam250{I}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-1 -2 -2 -2 -2 -2 -2 -3 -2 5 2 -2 2 1 -2 -1 0 -5 -1 4 -2 -2 0);
@{$pam250{L}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-2 -3 -3 -4 -6 -2 -3 -4 -2 2 6 -3 4 2 -3 -3 -2 -2 -1 2 -3 -3 0);
@{$pam250{K}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-1 3 1 0 -5 1 0 -2 0 -2 -3 5 0 -5 -1 0 0 -3 -4 -2 1 0 0);
@{$pam250{M}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-1 0 -2 -3 -5 -1 -2 -3 -2 2 4 0 6 0 -2 -2 -1 -4 -2 2 -2 -2 0);
@{$pam250{F}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-4 -4 -4 -6 -4 -5 -5 -5 -2 1 2 -5 0 9 -5 -3 -3 0 7 -1 -5 -5 0);
@{$pam250{P}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(1 0 -1 -1 -3 0 -1 -1 0 -2 -3 -1 -2 -5 6 1 0 -6 -5 -1 -1 0 0);
@{$pam250{S}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(1 0 1 0 0 -1 0 1 -1 -1 -3 0 -2 -3 1 2 1 -2 -3 -1 0 0 0);
@{$pam250{T}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(1 -1 0 0 -2 -1 0 0 -1 0 -2 0 -1 -3 0 1 3 -5 -3 0 0 -1 0);
@{$pam250{W}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-6 2 -4 -7 -8 -5 -7 -7 -3 -5 -2 -3 -4 0 -6 -2 -5 17 0 -6 -5 -6 0);
@{$pam250{Y}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(-3 -4 -2 -4 0 -4 -4 -5 0 -1 -1 -4 -2 7 -5 -3 -3 0 10 -2 -3 -4 0);
@{$pam250{V}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(0 -2 -2 -2 -2 -2 -2 -1 -2 4 2 -2 2 -1 -1 -1 0 -6 -2 4 -2 -2 0);
@{$pam250{B}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(0 -1 2 3 -4 1 2 0 1 -2 -3 1 -2 -5 -1 0 0 -5 -3 -2 2 2 0);
@{$pam250{Z}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(0 0 1 3 -5 3 3 -1 2 -2 -3 0 -2 -5 0 0 -1 -6 -4 -2 2 3 0);
@{$pam250{X}}{qw(A R N D C Q E G H I L K M F P S T W Y V B Z X)}=
	qw(0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0);
# print Dumper(%pam250); 
    return \%pam250; 
}

push @EXPORT_OK, 'testseq'; 
sub testseq{
# $seq=testseq([$index]) returns  one of several protein seqs,  usefull
# for testing. If index is omited seq[0] is used; 
    my $seq; 
    my $index=defined $_[0]?$_[0]:0; 
    if($index==0){
        $seq="MSGKMTGIVKWFNADKGFGFITPDDGSKDVFVHFSAIQNDGYKS
             LDEGQKVSFTIESGAKGPAAGNVTSL"; 
    }elsif($index==1){
        $seq="MKPNIHPEYRTVVFHDTSVDEYFKIGSTIKTDREIELDGVTYPY
             VTIDVSSKSHPFYTGKLRTVASEGNVARFTQRFGRFVSTKKGA"; 
    }elsif($index==2){
        $seq="MSSDYAGELMIWIMLATLAVVFVVGFRVLTSGARKAIRRLSDRL
             NIDVVPVESMVDQMGKSAGDEFLRYLHRPDESHLQNAAQVLLIWQIVIVDGSEQNLLQ
             WHRILQKARLAAPITDAQVRLALGFLRETEPEMQDINAFQMRYNAFFQPAEGVHWLH"; 
    }elsif($index==3){
        $seq="MDKNELVQKAKLAEQAERYDDMAACMKSVTEQGAELSNEERNLLSVAYKNVVGARRSSWR
            VVSSIEQKTEGAEKKQQMAREYREKIETELRDICNDVLSLLEKFLIPNASQAESKVFYLK
            MKGDYYRYLAEVAAGDDKKGIVDQSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYE
            ILNSPEKACSLAKTAFDEAIAELDTLSEESYKDSTLIMQLLRDNLTLWTSDTQGDEAEAG
            EGGEN"; 
    }elsif($index==4){
        $seq="MKQGLQLRLSQQLAMTPQLQQAIRLLQLSTLELQQELQQALESN
             PLLEQIDTHEEIDTRETQDSETLDTADALEQKEMPEELPLDASWDTIYTAGTPSGTSG
             DYIDDELPVYQGETTQTLQDYLMWQVELTPFSDTDRAIATSIVDAVDETGYLTVPLED
             ILESIGDEEIDIDEVEAVLKRIQRFDPVGVAAKDLRDCLLIQLSQFDKTTPWLEEARL
             IISDHLDLLANHDFRTLMRVTRLKEDVLKEAVNLIQSLDPRPGQSIQTGEPEYVIPDV
             LVRKHNGHWTVELNSDSIPRLQINQHYASMCNNARNDGDSQFIRSNLQDAKWLIKSLE
             SRNDTLLRVSRCIVEQQQAFFEQGEEYMKPMVLADIAQAVEMHESTISRVTTQKYLHS
             PRGIFELKYFFSSHVNTEGGGEASSTAIRALVKKLIAAENPAKPLSDSKLTSLLSEQG
             IMVARRTVAKYRESLSIPPSNQRKQLV"; 
   }elsif($index==5){
        $seq="MNVIAILNHMGVYFKEEPIRELHRALERLNFQIVYPNDRDDLLK
             LIENNARLCGVIFDWDKYNLELCEEISKMNENLPLYAFANTYSTLDVSLNDLRLQISF
             FEYALGAAEDIANKIKQTTDEYINTILPPLTKALFKYVREGKYTFCTPGHMGGTAFQK
             SPVGSLFYDFFGPNTMKSDISISVSELGSLLDHSGPHKEAEQYIARVFNADRSYMVTN
             GTSTANKIVGMYSAPAGSTILIDRNCHKSLTHLMMMSDVTPIYFRPTRNAYGILGGIP
             QSEFQHATIAKRVKETPNATWPVHAVITNSTYDGLLYNTDFIKKTLDVKSIHFDSAWV
             PYTNFSPIYEGKCGMSGGRVEGKVIYETQSTHKLLAAFSQASMIHVKGDVNEETFNEA
             YMMHTTTSPHYGIVASTETAAAMMKGNAGKRLINGSIERAIKFRKEIKRLRTESDGWF
             FDVWQPDHIDTTECWPLRSDSTWHGFKNIDNEHMYLDPIKVTLLTPGMEKDGTMSDFG
             IPASIVAKYLDEHGIVVEKTGPYNLLFLFSIGIDKTKALSLLRALTDFKRAFDLNLRV
             KNMLPSLYREDPEFYENMRIQELAQNIHKLIVHHNLPDLMYRAFEVLPTMVMTPYAAF
             QKELHGMTEEVYLDEMVGRINANMILPYPPGVPLVMPGEMITEESRPVLEFLQMLCEI
             GAHYPGFETDIHGAYRQADGRYTVKVLKEESKK"; 
    }else{
        $seq="MSTVDKEELVQKAKLAEQSERYDDMAQAMKSVTETGVELSNEERNLLSVAYKNVVGARRS
            SWRVISSIEQKTEASARKQQLAREYRERVEKELREICYEVLGLLDKYLIPKASNPESKVF
            YLKMKGDYYRYLAEVATGDARNTVVDDSQTAYQDAFDISKGKMQPTHPIRLGLALNFSVF
            YYEILNSPDKACQLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDTQGDEA
            EPQEGGDN"; 
    }
    $seq=~tr/ \t\n//d; 
    return $seq; 
}

sub ncbi_table{
# %hash=ncbi_table[$table] return a hash with the corresponding
# translation table
    my($table)=shift; 
    my @ncbi_table; 
$ncbi_table[1]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA * TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga * tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[2]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA * AGC S AGG * AGT S ATA M ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA W TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga * agc S agg * agt S ata M atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga W tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[3]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA M ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA T CTC T CTG T CTT T 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA W TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata M atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta T ctc T ctg T ctt T 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga W tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[4]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA W TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga W tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[5]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA S AGC S AGG S AGT S ATA M ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA W TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga S agc S agg S agt S ata M atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga W tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[6]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA Q TAC Y TAG Q TAT Y TCA S TCC S TCG S TCT S 
	TGA * TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa Q tac Y tag Q tat Y tca S tcc S tcg S tct S 
	tga * tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[9]= {qw(
	AAA N AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA S AGC S AGG S AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA W TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa N aac N aag K aat N aca T acc T acg T act T 
	aga S agc S agg S agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga W tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[10]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA C TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga C tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[11]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA * TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga * tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[12]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG S CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA * TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg S ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga * tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[13]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA G AGC S AGG G AGT S ATA M ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA W TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga G agc S agg G agt S ata M atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga W tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[14]= {qw(
	AAA N AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA S AGC S AGG S AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA Y TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA W TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa N aac N aag K aat N aca T acc T acg T act T 
	aga S agc S agg S agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa Y tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga W tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[15]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG Q TAT Y TCA S TCC S TCG S TCT S 
	TGA * TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag Q tat Y tca S tcc S tcg S tct S 
	tga * tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[16]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG L TAT Y TCA S TCC S TCG S TCT S 
	TGA * TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag L tat Y tca S tcc S tcg S tct S 
	tga * tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[21]= {qw(
	AAA N AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA S AGC S AGG S AGT S ATA M ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG * TAT Y TCA S TCC S TCG S TCT S 
	TGA W TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa N aac N aag K aat N aca T acc T acg T act T 
	aga S agc S agg S agt S ata M atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag * tat Y tca S tcc S tcg S tct S 
	tga W tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
$ncbi_table[22]= {qw(
	AAA K AAC N AAG K AAT N ACA T ACC T ACG T ACT T 
	AGA R AGC S AGG R AGT S ATA I ATC I ATG M ATT I 
	CAA Q CAC H CAG Q CAT H CCA P CCC P CCG P CCT P 
	CGA R CGC R CGG R CGT R CTA L CTC L CTG L CTT L 
	GAA E GAC D GAG E GAT D GCA A GCC A GCG A GCT A 
	GGA G GGC G GGG G GGT G GTA V GTC V GTG V GTT V 
	TAA * TAC Y TAG L TAT Y TCA * TCC S TCG S TCT S 
	TGA * TGC C TGG W TGT C TTA L TTC F TTG L TTT F 
	aaa K aac N aag K aat N aca T acc T acg T act T 
	aga R agc S agg R agt S ata I atc I atg M att I 
	caa Q cac H cag Q cat H cca P ccc P ccg P cct P 
	cga R cgc R cgg R cgt R cta L ctc L ctg L ctt L 
	gaa E gac D gag E gat D gca A gcc A gcg A gct A 
	gga G ggc G ggg G ggt G gta V gtc V gtg V gtt V 
	taa * tac Y tag L tat Y tca * tcc S tcg S tct S 
	tga * tgc C tgg W tgt C tta L ttc F ttg L ttt F 
	)}; 
    return %{$ncbi_table[$table]}; 
}

package fasta;
use Carp;
sub TIEHANDLE {
# USAGE: while( (\$name, \$seqref) = <FASTA> ){etc}
    my $class=shift; 
    my $fasta=shift; 
    open (my $self, "<$fasta") || croak "can't open $fasta: $!";
    while(<$self>){
        next unless /^>(\S+)\s*(.*)$/; 
        $$self->{nexthead}=$2; 
        $$self->{nextname}=$1; 
        last; 
    }
    $$self->{eof}= defined  $$self->{nextname} ? 0 : 1;
    return bless $self, $class; 
}
sub READLINE {
    my $self = shift;
    return () if $$self->{eof}; 
    my $name=$$self->{nextname}; 
    my $head=$$self->{nexthead}; 
    my $seq=''; 
    while(<$self>){
        if( /^>(\S+)\s*(.*)$/ ){
            $$self->{nextname}=$1; 
            $$self->{nexthead}=$2; 
            return $name, \$seq, $head; 
        }
        chomp; 
        $seq.=$_; 
    }
    $$self->{eof}=1; 
    return ($name, \$seq, $head) if $seq; 
    
    return <$self>;
}

sub CLOSE{
    my $self = shift;
    return close $self;
}
#==========
package qual;
use Carp;
sub TIEHANDLE {
# USAGE: while( (\$name, \$seqref) = <FASTA> ){etc}
    my $class=shift; 
    my $fasta=shift; 
    open (my $self, "<$fasta") || croak "can't open $fasta: $!";
    while(<$self>){
        next unless /^>(\S+)\s*(.*)$/; 
        $$self->{nextname}=$1; 
        $$self->{nexthead}=$2; 
        last; 
    }
    $$self->{eof}= defined  $$self->{nextname} ? 0 : 1;
    return bless $self, $class; 
}
sub READLINE {
    my $self = shift;
    return () if $$self->{eof}; 
    my $name=$$self->{nextname}; 
    my $head=$$self->{nexthead}; 
    my $seq=''; 
    while(<$self>){
        if( /^>(\S+)\s*(.*)/ ){
            $$self->{nextname}=$1; 
            $$self->{nexthead}=$2; 
            return $name, \$seq, $head; 
        }
        chomp; 
        $seq.="$_ "; 
    }
    $$self->{eof}=1; 
        return ($name, \$seq, $head) if $seq; 
    
    return <$self>;
}

sub CLOSE{
    my $self = shift;
    return close $self;
}

1; 

__END__

Regex que encuentra el maximo segmento sin codones de stop
$max_orf_rgx='(...(?!(taa|tag|tga)))+';

Regex que encuentra el maximo CDS
$cds_rgx='atg(...(?!(taa|tag|tga)))+...(taa|tag|tga)'; 

sub generate_mere_list{
# generate_mere_list($merlen) returns an alfabeticaly sorted list of
# all the DNA multimers of length $merlen
    my($merlen)=shift; 
    my($len, $nuc, $mere, @prev, @new); 
    $prev[0]=''; 
    for($len=0; $len<$merelen; $len++){
        @new=(); 
        for $nuc(qw( a c g t)){
            for $mere (@prev){
                push(@new, "$mere$nuc"); 
            }
        }
        @prev=@new; 
    }
    return @prev; 
}

sub overlap{
# ($size, $from, $to)=overlap($f1, $t1, $f2, $t2) returns the size and
# the ends of the overlap of both ranges. If there is no overlap,  both
# ends are -1
    my($a, $b, $c, $d)=@_; 
    if($b<$a ||$d<$c){
        warn "ranges are wrong: $a $b $c $d\n"; 
        return (-1, -1, -1); 
    }
    if($a>$d){
        return (0, -1, -1); 
    }elsif($a>$c){
        if($b<$d){
            return ($b-$a+1, $a, $b)
        }else{
            return ($d-$a+1, $a, $d)
        }
    }else{
        if($b<$c){
            return (0, -1, -1)
        }elsif($b<$d){
            return ($b-$c+1, $c, $b)
        }else{
            return ($d-$c+1, $c, $d); 
        }
    }
}
