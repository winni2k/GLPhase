#!/usr/bin/perl
# vcf2marker.pl                   wkretzsch@gmail.com
#                                 17 Sep 2013

use warnings;
use strict;
$| = 1;
use Data::Dumper;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use Carp;

use Getopt::Std;
my %args;
getopts( 'd:v:', \%args );
my $DEBUG = $args{d} || 1;

croak "File $args{v} does not exist" unless -f $args{v};

my ( $name, $path, $suffix ) = fileparse( $args{v}, ".gz" );

my $useGz = $suffix ? 1 : 0;

my $cmd = '';
$cmd .= "gzip -dc "    if $useGz;
$cmd .= "dd bs=1M if=" if !$useGz;

$cmd .= "$args{v}";
$cmd .= q( | grep -v '^#');
$cmd .= q( | awk '{print $1, $2,$3,$4,$5}');
$cmd .= " | ";

open( my $fh, $cmd );
while (<$fh>) {
    chomp;
    my @line = split(/\s/);
    my @out = ( $line[2], $line[1], @line[ 3 .. 4 ] );
    $out[0] = join( ':', @line[ 0 .. 1 ] ) if $out[0] eq '.';
    print join( "\t", @out ) . "\n";
}

__END__

=head1 NAME

vcf2marker.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION

This script takes a vcf through stdin and turns it in to a marker file.

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
