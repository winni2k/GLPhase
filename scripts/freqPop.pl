#!/usr/bin/perl
# freqPop.pl                   wkretzsch@gmail.com
#                              26 Nov 2013
our $VERSION = '0.002';
$VERSION = eval $VERSION;
print STDERR "freqPop.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

=head1 CHANGES

0.2  Tue Dec 17 10:25:35 GMT 2013
     Added support for filtering on -g instead of -l

=cut

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
getopts( 'g:df:l:p:L:', \%args );
my $DEBUG = $args{d} || 1;

my $freqFile = $args{f} || $args{L};
croak "need to define a pop" unless defined $args{p};
croak "need to pass in freq file" unless -e $freqFile;
croak "need to pass in legend or gen file"
  unless ( defined $args{l} xor defined $args{g} );

my $openCmd;
my $siteFile = defined $args{l} ? $args{l} : $args{g};
croak "file does not exist: $siteFile" unless -e $siteFile;
if ( $siteFile =~ m/\.gz$/ ) {
    $openCmd = " gzip -dc $siteFile |";
}
else {
    $openCmd = " < $siteFile";
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

    # unshift chromosome id if legFH is a gen file
    shift @L if defined $args{g}; 
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

Pulls the frequencies from a freq file (-f) for all positions in the legend (-l) or gen (-g) file and creates a population frequency file.

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

