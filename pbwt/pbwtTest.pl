#!/usr/bin/perl
# pbwtTest.pl                   wkretzsch@gmail.com
#                               06 Nov 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "pbwtTest.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

use warnings;
use strict;
$| = 1;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;
use Term::ANSIColor;

use Getopt::Std;
my %args;
getopts( 'drj:', \%args );
my $DEBUG = $args{d} || 1;

srand(2);

my @aaHaps = (
    [qw/0 1 0 0 0/], [qw/0 0 1 0 0/], [qw/0 0 0 0 0/], [qw/0 0 1 0 0/],
    [qw/1 0 0 0 0/], [qw/0 0 0 1 0/],
);
@aaHaps = RandHaps( 5, 6 );

my $numHaps  = @aaHaps;
my $numSites = @{ $aaHaps[0] };

VisualizeMat( \@aaHaps, "Haps", 0 );

# this is k = 1
my @aaDk = ( [ 0, (0) x ( $numHaps - 1 ) ] );
my @aaAk = ( [ 0 .. $numHaps - 1 ] );

my $maxK = $numSites;
for my $k ( 0 .. $maxK - 1 ) {
    BuildPrefixAndDivergenceArrays( \@aaHaps, \@aaDk, \@aaAk, $k );
}

#add sentinel array to Dk
# push @aaDk, [($maxK+1) x $numHaps];

VisualizeMat( \@aaDk, "Dk", 0 );
VisualizeMat( \@aaAk, "Ak", 0 );

# visualize Ak with differing k
for my $sortRow ( 1 .. $maxK ) {

    # find longest match
    my $matchLength = 2;
    my @matches;
    if ( $sortRow > $matchLength ) {
        @matches =
          ReportLongMatches( \@aaDk, \@aaHaps, \@aaAk,
            $numHaps, $sortRow, $matchLength );
        print "Matches greater length $matchLength: ";
        for my $set (@matches) {
            print " " . join( " ", @{$set} ) . ";";
        }
        print "\n";

    }
    print "Dk[$sortRow]: " . join( ' ', @{ $aaDk[$sortRow] } ) . "\n";
    print "Ak[$sortRow]: " . join( ' ', @{ $aaAk[$sortRow] } ) . "\n";
    my $raaSortedHaps = SortHapsByAk( $aaAk[$sortRow], \@aaHaps );
    VisualizeMat( $raaSortedHaps, "Haps sorted by Ak row $sortRow", 0, \@matches );
}

sub ReportLongMatches {

    my ( $raaDk, $raaHaps, $raaAk, $M, $k, $L ) = @_;
    my ( $v, $u, @a, @b ) = qw/0 0/;

    croak "k needs to be greater 0" unless $k > 0;

    my $hapLength = @{$raaAk} - 1;
    my @matches;

    for my $i ( 0 .. $M - 1 ) {
        if ( $raaDk->[$k][$i] > $k - $L - 1 ) {
            if ( @a > 0 && @b > 0 ) {
                push @matches, [ @a, @b ];
            }
            @a = ();
            @b = ();
        }
        if ( $k == $hapLength ) {
            if ( $i % 2 == 0 ) {
                push @a, $raaAk->[$k][$i];
            }
            else {
                push @b, $raaAk->[$k][$i];
            }
        }
        elsif ( $raaHaps->[ $raaAk->[$k][$i] ][$k] == 0 ) {
            push @a, $raaAk->[$k][$i];
        }
        else {
            push @b, $raaAk->[$k][$i];
        }
    }
    push @matches, [ @a, @b ] if @a && @b;

    return @matches;

}

sub RandHaps {

    my ( $rows, $cols ) = @_;

    my @aaHaps;
    for my $row ( 1 .. $rows ) {
        my @hap;
        for my $col ( 1 .. $cols ) {
            push @hap, int( rand(2) );
        }
        push @aaHaps, \@hap;
    }
    return @aaHaps;
}

sub SortHapsByAk {

    my $raAk    = shift;
    my $raaHaps = shift;

    my @retHaps;
    for my $idx ( @{$raAk} ) {
        push @retHaps, $raaHaps->[$idx];
    }
    return \@retHaps;
}

sub VisualizeMat {

    my $raa   = shift;
    my $title = shift;
    my $base  = shift;
    my $raMatches = shift;

    my @matches;
    for my $set (@{$raMatches}){
        push @matches, @{$set};
    }
    my %matches = map { $_ => 1 } @matches;

    print $title. "\n";
    my $rowNum = $base;
    for my $row ( @{$raa} ) {
        my $allNum = 0;
        for my $all (@{$row}){
            $allNum ++;
            print ' ' if $allNum > 1;
            if(defined $matches{$rowNum}){
                print colored [ 'yellow' ], $all;
            }
            else{
                print $all;
            }
        print join( ' ', ( "$rowNum: ", @{$row} ) ) . "\n";
        $rowNum++;
    }
    print "\n";
}

sub BuildPrefixAndDivergenceArrays {

    my ( $raaHaps, $raaDk, $raaAk, $k ) = @_;

    # k is 0 based
    croak "k plus 1 needs to be > 0"
      unless ( $k >= 0 && int($k) == $k );

    my ( $p, $q ) = ( $k + 1, $k + 1 );
    my ( @a, @b, @d, @e );
    my $M = @{$raaHaps};    # M is number of haps

    for my $i ( 0 .. $M - 1 ) {
        $p = $raaDk->[$k][$i] if $raaDk->[$k][$i] > $p;
        $q = $raaDk->[$k][$i] if $raaDk->[$k][$i] > $q;
        if ( $raaHaps->[ $raaAk->[$k][$i] ][$k] == 0 ) {
            push @a, $raaAk->[$k][$i];
            push @d, $p;
            $p = 0;
        }
        else {
            push @b, $raaAk->[$k][$i];
            push @e, $q;
            $q = 0;
        }
    }

    push @{$raaAk}, [ @a, @b ];    # adding sentinel value at end
    push @{$raaDk}, [ @d, @e ];
}

__END__

=head1 NAME

pbwtTest.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION


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
