#!/usr/bin/perl
# hapsBin.pl                   wkretzsch@gmail.com
#                              14 Oct 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "hapsBin.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

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
getopts( 'di:n:s:', \%args );
my $DEBUG    = $args{d} || 1;
my $numSites = $args{n} || 1024;
my $stride   = $args{s} || 512;
croak "stride (-s) must be smaller than chunk size (-n)"
  unless $stride < $numSites;
my $haps = $args{i} || croak "Need to specify haps input file with -i";

my $openCmd = $args{i} =~ m/\.gz$/ ? "gzip -dc $args{i} |" : "< $args{i}";
my $lastLine =
  $args{i} =~ m/\.gz$/ ? qx(gzip -dc $args{i} | tail -1) : qx(tail -1 $args{i});
chomp $lastLine;
$lastLine =~ m/^(\S+) (\S+) (\S+)/;
my $lastChrom      = $1;
my $lastSite       = $2 eq '.' ? $3 : $2;
my $lastSiteLength = length $lastSite;

open( my $fh, $openCmd );
my @outBuffer;
my @outSites;
while (<$fh>) {
    chomp;
    m/^(\S+) (\S+) (\S+)/;
    my $chrom = $1;
    my $site = $2 eq '.' ? $3 : $2;

    push @outBuffer, $_;
    push @outSites,  $site;
    if ( @outBuffer == 2 * $numSites ) {
        my $outCmd = sprintf(
"| gzip -c > ${chrom}_%0${lastSiteLength}u_%0${lastSiteLength}u.haps.gz",
            $outSites[0], $outSites[ $numSites - 1 ] );
        open( my $outFH, $outCmd );
        for my $outLine ( @outBuffer[ 0 .. ( $numSites - 1 ) ] ) {
            print $outLine. "\n";
        }
        close($outFH);
        splice( @outBuffer, 0, ( $numSites - $stride ) );
        splice( @outSites,  0, ( $numSites - $stride ) );
    }
}

# print the remaining lines to an output file
my $outCmd = sprintf(
    "| gzip -c > ${lastChrom}_%0${lastSiteLength}u_%0${lastSiteLength}u.haps.gz",
    $outSites[0], $outSites[$#outSites] );
open( my $outFH, $outCmd );
for my $outLine (@outBuffer) {
    print $outLine. "\n";
}
close($outFH);

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
