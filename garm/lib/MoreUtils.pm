package MoreUtils;

use strict;
require Exporter;

our @ISA = qw(Exporter);
our %EXPORT_TAGS = ( 'all' => [ qw(
	
) ] );
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(uniq);



sub uniq (@) {
    my %seen = ();
    grep { not $seen{$_}++ } @_;
}



# Function aliases
*first_index = \&firstidx;
*last_index  = \&lastidx;
*first_value = \&firstval;
*last_value  = \&lastval;
*zip         = \&mesh;
*distinct    = \&uniq;

1;

__END__

