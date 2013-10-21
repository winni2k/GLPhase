#!/usr/bin/perl
# hapLegBin.pl                 wkretzsch@gmail.com
#                              14 Oct 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "hapLegBin.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

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
getopts( 'h:l:n:s:', \%args );
my $numSites = $args{n} || 1024;
my $stride   = $args{s} || 512;
croak "stride (-s) must be smaller than chunk size (-n)"
  unless $stride < $numSites;
my $hap = $args{h} || croak "Need to specify hap input file with -h";
my $leg = $args{l} || croak "Need to specify legend input file with -l";

my $lastLine = qx(tail -1 $leg);
chomp $lastLine;
$lastLine =~ m/^(\S+) (\S+)/;
my $lastChrom      = $1;
my $lastSite       = $2;
my $lastSiteLength = 9;    # nine works... #length $lastSite;

my $openCmd = $args{h} =~ m/\.gz$/ ? "gzip -dc $args{h} |" : "< $args{h}";
open( my $hapFH, $openCmd );

open( my $legFH, "<", $leg );
my $header = <$legFH>;     # discard header
chomp $header;
croak "header does not match 'id position a0 a1'"
  unless $header eq 'id position a0 a1';
my @outHapBuffer;
my @outLegBuffer;
my @outSites;

while (<$legFH>) {
    m/^(\S+) (\S+)/;
    my $chrom = $1;
    my $site  = $2;

    push @outLegBuffer, $_;
    my $hapLine = <$hapFH> || croak "error: hap file shorter than legend file";
    push @outHapBuffer, $hapLine;
    push @outSites,     $site;
    if ( @outLegBuffer == $numSites + $stride ) {
        print STDERR "$outSites[0]\t$outSites[ $numSites - 1 ]\n";
        my $outHapCmd = sprintf(
"| gzip -c > ${chrom}_%0${lastSiteLength}u_%0${lastSiteLength}u.hap.gz",
            $outSites[0], $outSites[ $numSites - 1 ] );
        open( my $outHapFH, $outHapCmd );
        open(
            my $outLegFH,
            ">",
            sprintf(
                "${chrom}_%0${lastSiteLength}u_%0${lastSiteLength}u.legend",
                $outSites[0], $outSites[ $numSites - 1 ]
            )
        );
        print $outLegFH $header . "\n";
        for my $outIdx ( 0 .. ( $numSites - 1 ) ) {
            print $outLegFH $outLegBuffer[$outIdx];
            print $outHapFH $outHapBuffer[$outIdx];
        }
        close($outHapFH);
        close($outLegFH);
        splice( @outHapBuffer, 0, $stride );
        splice( @outLegBuffer, 0, $stride );
        splice( @outSites,     0, $stride );
    }
}

croak "error: hap file longer than legend file" if defined <$hapFH>;

# print the remaining lines to an output file
my $outHapCmd = sprintf(
    "| gzip -c > ${lastChrom}_%0${lastSiteLength}u_%0${lastSiteLength}u.hap.gz",
    $outSites[0], $outSites[$#outSites] );
open( my $outHapFH, $outHapCmd );
open(
    my $outLegFH,
    ">",
    sprintf(
        "${lastChrom}_%0${lastSiteLength}u_%0${lastSiteLength}u.legend",
        $outSites[0], $outSites[$#outSites]
    )
);
print $outLegFH $header . "\n";
for my $outIdx ( 0 .. $#outSites ) {
    print $outLegFH $outLegBuffer[$outIdx];
    print $outHapFH $outHapBuffer[$outIdx];
}

__END__

=head1 NAME

hapsBin.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for hapsBin.pl, 
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
