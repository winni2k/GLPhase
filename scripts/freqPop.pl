#!/usr/bin/perl
# freqPop.pl                   wkretzsch@gmail.com
#                              26 Nov 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "freqPop.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

use warnings;
use strict;
$| = 1;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;
use List::MoreUtils qw/firstidx/;

use Getopt::Std;
my %args;
getopts( 'df:l:p:L:', \%args );
my $DEBUG = $args{d} || 1;

my $freqFile = $args{f} || $args{L};
croak "need to define a pop"        unless defined $args{p};
croak "need to pass in freq file"   unless -e $freqFile;
croak "need to pass in legend file" unless -e $args{l};

my $openCmd;
if ( $args{l} =~ m/\.gz$/ ) {
    $openCmd = " gzip -dc $args{l} |";
}
else {
    $openCmd = " < $args{l}";
}
open( my $legFH, $openCmd );
<$legFH>;    #remove header

# open freq file
if ( $freqFile =~ m/\.gz$/ ) {
    $openCmd = " gzip -dc $freqFile |";
}
else {
    $openCmd = " < $freqFile";
}
open( my $freqFH, $openCmd );
my $head = <$freqFH>;    # remove header
chomp $head;
my @head   = split( ' ', $head );
my $str    = $args{p} . '\.aaf';
my $popAAF = firstidx { $_ =~ m/($str)/i } @head;
$popAAF = 4 if $popAAF < 0;

# print header
print "$args{p}\n";

my $hit = 0;
LEGLINE: while (<$legFH>) {
    my $found = 0;
    chomp;
    my @L = split(q/ /);
  FREQLINE: while ( my $f = <$freqFH> ) {
        chomp $f;
        my @F = split( q/\s/, $f );
        if ( $L[1] eq $F[1] ) {
            $hit = 1;
            if ( $L[2] eq $F[2] && $L[3] eq $F[3] ) {
                print $F[$popAAF] . "\n";
                $found = 1;
                $hit   = 0;
                last FREQLINE;
            }
        }
    }
    if ( $hit == 1 ) {
        croak "Found position (@L) but it did not have correct alleles.";
    }
    die "$L[1] could not be found" unless $found;

}

__END__

=head1 NAME

freqPop.pl

=head1 SYNOPSIS
   
freqPop.pl -f freqFile -l legendFile -p EUR

=head1 DESCRIPTION

Pulls the frequencies from a freq file (-f) for all positions in the legend file (-l) and creates a population frequency file.

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

