#!/usr/bin/perl
# hapLeg2haps.pl                   wkretzsch@gmail.com
#                                  28 Nov 2013
our $VERSION = '0.002';
$VERSION = eval $VERSION;
print STDERR "hapLeg2haps.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

=head1 Changes

0.002 Only first four cols of legend file are used now
      Added -t option

0.001 First version

=cut

use warnings;
use strict;
$| = 1;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;

use Getopt::Std;
my %args;
getopts( 'l:h:c:t:', \%args );

my $chrom = $args{c} || croak "need to define chrom with -c";

my $cmd = openCmd( $args{l} );
open( my $legFH, $cmd );
my $header     = <$legFH>;                # get rid of header
my @header     = split( ' ', $header );
my $typeFilter = undef;
if ( @header > 4 && $header[4] eq 'type' ) {
    $typeFilter = $args{t};
}

$cmd = openCmd( $args{h} );
open( my $hapFH, $cmd );

my $sitePrinted = 0;
my $siteSeen    = 0;
while ( my $hap = <$hapFH> ) {
    $siteSeen++;
    my $leg = <$legFH>;
    chomp $leg;
    my @leg = split( ' ', $leg );
    croak "ERROR: $args{h} is longer than $args{l}" unless defined $leg;
    if ( !defined $typeFilter || $leg[4] eq $typeFilter ) {
        print $chrom. ' ' . join( ' ', @leg[ 0 .. 3 ] ) . ' ' . $hap;
        $sitePrinted++;
    }
    if ( $siteSeen % 1000 == 0 ) {
        print STDERR "$sitePrinted / $siteSeen\r";
    }
}
print STDERR "\n";
croak "ERROR: $args{h} is shorter than $args{l}" if defined <$legFH>;

sub openCmd {
    my $file = shift;
    my $cmd;
    if ( $file =~ m/\.gz/ ) {
        $cmd = "gzip -dc $file |";
    }
    else {
        $cmd = " < $file";
    }
    return $cmd;
}

__END__

=head1 NAME

hapLeg2haps.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION

merges hap and legend file to haps file

=head1 AUTHOR

Warren Winfried Kretzschmar, E<lt>wkretzsch@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Warren Winfried Kretzschmar

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
