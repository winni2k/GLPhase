#!/usr/bin/perl
# doseKeepSnps.pl                   wkretzsch@gmail.com
#                                   25 Sep 2013

use warnings;
use strict;
$| = 1;
use Data::Dumper;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);

use Getopt::Std;
my %args;
getopts( 'l:', \%args );
my $siteList = $args{l} || 0;

my %sites;
if ($siteList) {
    print STDERR "Filtering on Site list $siteList\n";
    my $cmd;
    if ( $siteList =~ m/\.gz$/ ) {
        $cmd = "gzip -dc $siteList |";
    }
    else {
        $cmd = " < $siteList";
    }
    open my $fhSL, $cmd
      or die "could not open $siteList with command $cmd: $!, $?";
    while (<$fhSL>) {
        chomp;
        m/([a-zA-Z:0-9\-]+)/;
        my $site = $1;
        if ( exists $sites{$site} ) {
            $sites{$site} = 0;
        }
        else {
            $sites{$site} = 1;
        }
    }
}

my $header = <>;
die "header malformed" unless $header =~ m/^marker/;
print $header;
while (<>) {
    chomp;
    my @line    = split(/ /);
    my $printOK = 0;
    for my $allele ( @line[ 1 .. 2 ] ) {
        $printOK++ if $allele =~ m/^[ATGCatgc]$/;
    }
    if ($siteList) {
        $printOK = 0 unless ( exists $sites{ $line[0] } && $sites{ $line[0] } );
    }
    print $_ . "\n" if $printOK == 2;
}

__END__

=head1 NAME

doseKeepSnps.pl

=head1 SYNOPSIS
   
doseKeepSnps.pl < in.dose > out.dose

=head1 DESCRIPTION

Removes any sites that are not SNPs.

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
