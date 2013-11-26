#!/usr/bin/perl
# gprobs2hapGen.pl                   wkretzsch@gmail.com
#                                    03 Oct 2013

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
getopts( 'o:', \%args );

confess "-o file is required" unless $args{o};

#confess "chromosome (-c) is required" unless $args{c};

open my $fh, " | gzip -c > $args{o}.gen.gz";
my @samples;
while (<>) {
    chomp;
    my @line = split(/ /);
    unless (@samples) {
        @samples = @line[ 3 .. $#line ];
        printSamples( $args{o}, \@samples );
        next;
    }

## Please see file perltidy.ERR
    my ( $chrom, $loc ) = split( /:/, $line[0] );
    confess "need to implement a markers option"
      unless defined $line[0] =~ m/^[0-9]+:[0-9]+$/;

    print $fh "$chrom $line[0] $loc $line[1] $line[2] "
      . join( ' ', @line[ 3 .. $#line ] ) . "\n";
}
close($fh);

sub printSamples {
    my $base      = shift;
    my $raSamples = shift;
    my @samples   = @{$raSamples};

    my $sampFile = $base . ".sample";
    open my $fh, '>', $sampFile;
    print $fh "ID_1 ID_2 missing\n0 0 0\n";
    my $sampNum = 0;
    for my $samp (@samples) {
        print $fh "$samp $samp 0\n" if $sampNum % 3 == 0;
        $sampNum++;
    }
    close($fh);
}

__END__

=head1 NAME

gprobs2hapGen.pl

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
