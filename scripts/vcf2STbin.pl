#!/usr/bin/perl
# vcf2STbin.pl                   wkretzsch@gmail.com
#                                18 Sep 2013

use warnings;
use strict;
$| = 1;
use Data::Dumper;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use Croak;
use ListUtil qw/sum/;

use Getopt::Std;
my %args;
getopts( 'd', \%args );
my $DEBUG = $args{d} || 1;

my $inFile = shift;

croak "Input VCF file was not given or does not exist" unless -f $inFile;

my ( $name, $path, $suffix ) = fileparse( $inFile, ".gz" );

my $useGz = $suffix ? 1 : 0;

( $name, $path, $suffix ) = fileparse( $inFile, ".vcf.gz" );

my $cmd = '';
$cmd .= "gzip -dc "    if $useGz;
$cmd .= "dd bs=1M if=" if !$useGz;

$cmd .= "$inFile";
$cmd .= " | ";

# write the rest of the lines
my $lineNum = 0;
open( my $fh,    $cmd );
open( my $fhOut, " | gzip -c > $name.bin" );
my $headerWritten = 0;
LINE: while (<$fh>) {
    $lineNum++;
    chomp;
    next if m/^##/;
    my @line = split(/\t/);

    # write header
    if ( $headerWritten == 0 ) {

        # first three header columns
        print $fhOut join( "\t", qw/chr pos allele/ );

        # remove invariant columns
        splice( @line, 0, 9 );

        # print samples to header
        map { print $fhOut "\t$_" } @line;
        print $fhOut "\n";

        $headerWritten++;
        next LINE;
    }

    # first 3 cols
    my $alleles = join( '', @line[ 3 .. 4 ] );
    print $fhOut join( "\t", ( @line[ 0 .. 1 ], $alleles ) );

    # figure out where PLs are
    my @fields = split( /:/, $line[8] );
    my $PLFieldNum;
    for my $field ( 0 .. $#fields ) {
        if ( $fields[$field] eq 'PL' ) {
            $PLFieldNum = $field;
            last;
        }
    }
    croak "PL field does not exist in line $lineNum" unless defined $PLFieldNum;

    # print renormed likelihoods
    splice( @line, 0, 9 );
    for my $sampField (@line) {
        my @field = split( /:/, $sampField );
        my @pls   = split( /,/, $field[$PLFieldNum] );
        my @gls = map { 10 ^ ( ( -$_ ) / 10 ) } @pls;
        my $sum = sum @gls;
        @gls = map { $_ / $sum } @gls;

        # only keep het and homo alt
        print "\t" . join( ' ', @gls[ 1 .. 2 ] );
    }
    print "\n";
}

__END__

=head1 NAME

vcf2STbin.pl

=head1 SYNOPSIS

# convert $MYFILE.vcf.gz to bin file named $MYFILE.bin
vcf2STbin.pl $MYFILE.vcf.gz   

=head1 DESCRIPTION

Converts vcf to SNPTools bin file.
VCF needs to have PL fields

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
