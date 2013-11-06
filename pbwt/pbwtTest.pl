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

use Getopt::Std;
my %args;
getopts( 'drj:', \%args );
my $DEBUG = $args{d} || 1;

my $numSites = 5;
my @aaHaps   = (
    [qw/0 1 0 0 0/], [qw/0 0 1 0 0/], [qw/0 0 0 0 0/], [qw/0 0 1 0 0/],
    [qw/0 0 0 1 0/]
);

VisualizeMat(\@aaHaps, "Haps");

# this is k = 1
my @aaDk = ( [qw/2 0 0 0 0/] );
my @aaAk = ( [ 0 .. $numSites - 1 ] );

BuildPrefixAndDivergenceArrays( \@aaHaps, \@aaDk, \@aaAk, 1 );
VisualizeMat(\@aaDk, "Dk");
VisualizeMat(\@aaAk, "Ak");


sub VisualizeMat {

    my $raa = shift;
    my $title = shift;
    
    print $title."\n";
    for my $row (@{$raa}){
        print join(' ', @{$row})."\n";
    }
    print "\n";
}

sub BuildPrefixAndDivergenceArrays {

    my ( $raaHaps, $raaDk, $raaAk, $k ) = @_;

    # k is 1 based
    croak "k needs to be >0" unless ( $k > 0 && int($k) == $k );

    my ( $u, $v, $p, $q ) = ( 0, 0, $k + 1, $k + 1 );
    my ( @a, @b, @d, @e );
    my $M = @{$raaHaps};    # M is number of haps

    for my $i ( 0 .. $M - 1 ) {
        $p = $raaDk->[$k-1][$i] if $raaDk->[$k-1][$i] > $p;
        $q = $raaDk->[$k-1][$i] if $raaDk->[$k-1][$i] > $q;
        if($raaHaps->[$i][$k-1] == 0){
            push @a, $raaAk->[$k-1][$i];
            push @d, $p;
            $p = 0;
            $u = $u+1;
        }
        else{
            push @b, $raaAk->[$k-1][$i];
            push @e, $q;
            $q = 0;
            $v = $v+1;
        }
    }

    push @{$raaAk}, [@a, @b];
    push @{$raaDk}, [@d, @e];
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
