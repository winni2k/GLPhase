#!/usr/bin/perl
# vcfFuse.pl                   wkretzsch@gmail.com
#                              16 Oct 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "vcfFuse.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

use warnings;
use strict;
$| = 1;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;
use POSIX qw/ceil floor/;

use Getopt::Std;
my %args;
getopts( 'd:v:', \%args );

# read input files
opendir( my $dh, $args{d} ) || die "can't opendir $args{d}: $!";
my @vcfs = sort grep { m/\.vcf(\.gz|)$/ && -f "$args{d}/$_" } readdir($dh);
closedir $dh;

print STDERR "Found " . @vcfs . " vcf files to fuse\n";

my ( $firstVcf, $secVcf );
my $firstVcfName = shift @vcfs;
$secVcf = loadVcf( $firstVcfName, $args{d} );

my $vcfCmd = $args{v} =~ m/\.gz/ ? " | gzip -c > $args{v}" : " > $args{v}";
open( my $outFH, $vcfCmd );

map { print $outFH $_ . "\n" } @{ $secVcf->{headers} };
print $outFH "##vcfFuse.pl.version=$VERSION\n";
print $outFH "##NumberOfFusedChunks=" . ( @vcfs + 1 ) . "\n";
print $outFH $secVcf->{chromHeader} . "\n";

# load two files at a time
while (@vcfs) {
    my $vcf = shift @vcfs;
    print STDERR "Loading $vcf\n";
    $firstVcf = $secVcf;
    $secVcf = loadVcf( $vcf, $args{d} );

    # make sure both vcfs are on same chrom
    croak "all input vcfs must be of same chromosome"
      if $firstVcf->{chrom} ne $secVcf->{chrom};

    # make sure all loaded vcfs have same number of individuals
    my $indIdx  = 0;
    my @secInds = @{ $secVcf->{individuals} };
    for ( @{ $firstVcf->{individuals} } ) {
        croak
          "individuals in $secVcf->{file} don't match previous individual(s)"
          unless $_ eq $secInds[$indIdx];
        $indIdx++;
    }

    # find position where overlap starts
    my $overlapStartIdx = 0;
    my $secVcfFirstPos  = $secVcf->{positions}->[0];
    for ( @{ $firstVcf->{positions} } ) {
        last if $_ eq $secVcfFirstPos;
        $overlapStartIdx++;
    }

    # find index where overlap ends
    # using BED basing for the overlap region
    my $overlapEndIdx    = $overlapStartIdx;
    my $firstVcfNumSites = $#{ $firstVcf->{positions} } + 1;

    # splice away overlapping sites if we have an overlap
    unless ( $overlapStartIdx == $firstVcfNumSites - 1 ) {

        my $firstVcfLastPos = $firstVcf->{positions}->[ $firstVcfNumSites - 1 ];
        for ( @{ $secVcf->{positions} } ) {
            $overlapEndIdx++;
            last if $_ eq $firstVcfLastPos;
        }

        my $overlapLength = $overlapEndIdx - $overlapStartIdx;
        my $firstVcfLastKeepIdx =
          ceil( $overlapLength / 2 ) + $overlapStartIdx - 1;
        printLines( $firstVcf, $outFH, $firstVcfLastKeepIdx );
        my $secVcfFirstKeepIdx = ceil( $overlapLength / 2 );
        splice( @{ $secVcf->{positions} }, 0, $secVcfFirstKeepIdx );
        splice( @{ $secVcf->{siteLines} }, 0, $secVcfFirstKeepIdx );

    }
    else {
        printLines( $firstVcf, $outFH, $firstVcfNumSites - 1 );
    }

    # print out new vcf
}

#print out remaining lines
printLines( $secVcf, $outFH, $#{$secVcf->{siteLines}});

sub printLines {

    my $vcf     = shift;
    my $fh      = shift;
    my $lastIdx = shift;

    for ( 0 .. $lastIdx ) {
        print $fh $vcf->{siteLines}->[$_] . "\n";
    }
}

sub loadVcf {

    my $vcfName = shift || confess;
    my $dir     = shift || confess;

    my $vcfObject = {
        headers     => [],
        chromHeader => undef,
        individuals => [],
        file        => "$dir/$vcfName",
        chrom       => undef,
        positions   => [],
        siteLines   => [],
    };
    my $openCmd =
      $vcfName =~ m/\.gz$/ ? " gzip -dc $dir/$vcfName |" : " < $dir/$vcfName";
    open( my $vcfFH, $openCmd );

    my $chrom          = undef;
    my $headerComplete = 0;
    while (<$vcfFH>) {
        chomp;
        if (m/^#/) {
            croak "found header line after #CHROM line" if $headerComplete;
            if (m/^##/) {
                push @{ $vcfObject->{headers} }, $_;
            }
            elsif (m/^#CHROM/) {
                $vcfObject->{chromHeader} = $_;
                my @line = split(/\t/);
                splice( @line, 0, 9 );
                push @{ $vcfObject->{individuals} }, @line;
                $headerComplete = 1;
            }
        }
        else {
            push @{ $vcfObject->{siteLines} }, $_;
            m/^(\S+)\t(\S+)\t/;
            if ( defined $chrom ) {
                croak "input vcfs may only contain sites from one chromosome"
                  if $1 ne $chrom;
            }
            else {
                $chrom = $1;
            }

            push @{ $vcfObject->{positions} }, $2;

        }
    }
    $vcfObject->{chrom} = $chrom;
    return $vcfObject;
}

__END__

=head1 NAME

vcfFuse.pl

=head1 SYNOPSIS
   
vcfFuse.pl -d vcfDir -v fusedVCF.gz

=head1 DESCRIPTION

Takes all vcf (or vcf.gz) files in dir, sorts them by !!!basename!!! and fuses them.

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
## Please see file perltidy.ERR
