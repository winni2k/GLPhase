#!/usr/bin/perl
# hapsFilter.pl                   wkretzsch@gmail.com
#                                 14 Oct 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "hapsFilter.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

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
getopts( 'dp:s:', \%args );
my $DEBUG = $args{d} || 1;

if ( $args{p} ) {
    print STDERR "Filtering on sites in $args{p}\n";
    print STDERR "Using exact matching\n";
}

open( my $fh, "<", $args{p} );
my @positions;
my ( $lastChr, $lastSite );
while (<$fh>) {
    chomp;
    my @line = split(/\t/);
    croak "positions file needs to have exactly two tab separated entries"
      unless @line == 2;
    unless ( defined $lastSite ) {
        ( $lastChr, $lastSite ) = @line;
    }
    else {
        croak "positions file needs to be sorted"
          if ( $line[0] eq $lastChr && $line[1] < $lastSite );
    }
    push @positions, @line;
}

print STDERR "Filtering on " . ( @positions / 2 ) . " sites.\n";
print STDERR "Writing haps file.\n";
my $searchChrom = shift @positions;
my $searchSite  = shift @positions;
while (<>) {
    chomp;
    m/^(\S+) (\S+) (\S+)/;
    my $chrom = $1;
    my $site = $2 eq '.' ? $3 : $2;
    if ( $chrom eq $searchChrom ) {
        if ( $site eq $searchSite ) {
            print $_."\n";
            ( $searchChrom, $searchSite ) = splice( @positions, 0, 2 );
        }
        elsif ( $site > $searchSite ) {
            croak
              "Could not find position $searchChrom:$searchSite in haps input";
        }
    }
}

if (@positions) {
    croak "The following positions are not in the haps input:\n"
      . join( ' ', @positions );
}

__END__

=head1 NAME

hapsFilter.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for hapsFilter.pl, 
created by template.el.

It looks like the author of this script was negligent
enough to leave the stub unedited.


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
