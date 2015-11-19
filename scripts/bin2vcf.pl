#!/usr/bin/perl
# bin2vcf.pl     wkretzsch@gmail.com     2015/11/18 09:38:24
use warnings;
use strictures;
$| = 1;

our $VERSION = "0.001";
$VERSION = eval $VERSION;

use autodie;
use Carp;
use List::Util qw/max/;
use POSIX qw/log10/;

use Getopt::Std;
my %args;
getopts( 'o:', \%args );

my $glTag = "GL";
$glTag = "PL" if $args{P};

# print header
print <<_EOF_;
##fileformat=VCFv4.2
##source=bin2vcf.plv$VERSION
##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihood comprised of comma separated floating point log10-scaled likelihoods">
_EOF_

my $header = <>;
chomp $header;

my @header = split "\t", $header;
my @samples = @header[ 3 .. $#header ];
print join( "\t",
    ( '#CHROM', qw/POS ID REF ALT QUAL FILTER INFO FORMAT/, @samples ) )
  . "\n";

# parse site lines
while (<>) {
    chomp;
    my @line = split "\t";
    my ( $ref, $alt ) = get_alleles( $line[2] );
    print join( "\t", ( @line[ 0, 1 ], '.', $ref, $alt, qw/. . ./, $glTag ) );

    # parse GLs
    my @gls = @line[ 3 .. $#line ];
    croak "Wrong number of samples in GLs" if @gls != @samples;
    for (@gls) {
        my ( $het, $ahom ) = split ' ', $_;
        my $rhom = max( 1 - $het - $ahom, 0 );
        print "\t"
          . join( ',',
            map { $_ == 0 ? -999 : log10($_) } ( $rhom, $het, $ahom ) );
    }
    print "\n";
}

sub get_alleles {
    my $alls = shift;
    my ( $ref, $alt );
    if ( $alls =~ m/\s/ ) {
        ( $ref, $alt ) = split /\s/, $alls;
    }
    elsif ( length $alls == 2 ) {
        $ref = substr $alls, 0, 1;
        $alt = substr $alls, 1, 1;
    }
    else {
        croak "Could not parse alleles [$alls]";
    }
    return ( $ref, $alt );
}

__END__

=head1 NAME

bin2vcf.pl

=head1 SYNOPSIS

gzip -dc input.bin | bin2vcf.pl | bgzip -c > input.vcf.gz

=head1 DESCRIPTION

Converts text in SNPTools .bin format to a VCF.

=head1 AUTHOR

Warren Winfried Kretzschmar, E<lt>winni@winnis-MacBook-Pro2.localE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Warren Winfried Kretzschmar

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
