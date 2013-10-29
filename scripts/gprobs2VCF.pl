#!/usr/bin/perl
# gprobs2VCF.pl                   wkretzsch@gmail.com
#                                 21 Oct 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "gprobs2VCF.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

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
getopts( 's:g:o:', \%args );

our @gGTTable = qw(0/0 1/0 1/1);

croak "Need to defined input gprobs file" unless $args{g};
my $inCMD = $args{g} =~ m/\.gz/ ? "gzip -dc $args{g} |" : "< $args{g}";
open( my $inFH, $inCMD );

my @samples = readHeader($inFH);
print STDERR "Samples: " . @samples . "\n";

croak "Need to define output file with .vcf or .vcf.gz file ending"
  unless $args{o} =~ m/\.vcf(|\.gz)$/;
my $outIsGz = $args{o} =~ m/\.gz$/;
my $openCMD = $outIsGz ? " | bgzip -c > $args{o} " : "> $args{o}";
open( my $outFH, $openCMD );

# print VCF header
print STDERR "printing VCF header\n";
printHeader( $outFH, \@samples, $args{s} );
print STDERR "printing VCF body\n";
printBody( $inFH, $outFH, \@samples );

sub printBody {

    my ( $inFH, $outFH, $raSamples ) = @_;

    my $numSamples = @{$raSamples};
    my $siteNum = 0;
    while (<$inFH>) {
        $siteNum ++;
        print STDERR "Processing site $siteNum\n" if $siteNum % 1000 == 0;
        chomp;
        my @line = split(/ /);
        my ( $chrom, $pos, $id ) = qw/. . ./;
        if ( $line[0] =~ m/:/ ) {
            ( $chrom, $pos ) = split( /:/, $line[0] );
        }
        else {
            $id = $line[0];
        }

        my @out = ( $chrom, $pos, $id, @line[ 1 .. 2 ], qw/. . . GT/ );
        splice( @line, 0, 3 );
        croak "incorrect number of lines in gprobs file"
          if $numSamples * 3 != @line;
        for my $sampNum ( 0 .. ( $numSamples - 1 ) ) {
            my @samp = splice( @line, 0, 3 );
            my $maxIdx = 0;
            map { $maxIdx = $_ if $samp[$_] > $samp[$maxIdx] } 1 .. 2;
            push @out, $gGTTable[$maxIdx];
        }
        print $outFH join( "\t", @out ) . "\n";
    }
}

sub readHeader {

    my $inFH = shift;

    my $line = <$inFH>;
    chomp $line;
    my $firstCols = 'marker alleleA alleleB';
    croak "malformed header: first three columns are not '$firstCols'"
      unless $line =~ m/^${firstCols}/;
    my @line = split( / /, $line );
    splice( @line, 0, 3 );
    croak "number of samples in header is not evenly divisible by three"
      unless @line % 3 == 0;
    my @return;

    for my $sampNum ( 0 .. $#line ) {
        push @return, $line[$sampNum] if $sampNum % 3 == 0;
    }
    return @return;

}

sub printHeader {

    my $outFH     = shift;
    my $raSamples = shift;
    my $source = shift;

    croak "source (-s) can not have any white space" if $source =~ m/\s/;

    my @out;
    push @out, q(##fileformat=VCFv4.2);
    push @out, q(##source=).$source if defined $source;
    push @out, q(##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">);
    push @out,
      join( "\t",
        ( qw/#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT/, @{$raSamples} ) );

    map { print $outFH $_ . "\n" } @out;

}

__END__

=head1 NAME

gprobs2VCF.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for gprobs2VCF.pl, 
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
