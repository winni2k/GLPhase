#!/usr/bin/perl
# doseKeepSamp.pl                   wkretzsch@gmail.com
#                                   25 Sep 2013

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
getopts( 'ds:', \%args );
my $DEBUG = $args{d} || 1;

my $samp = $args{s} || croak "need sample list passed in by -s";

# load samples
open my $fh, '<', $samp;
my %samples;
while (<$fh>) {
    chomp;
    $samples{$_} = 1;
}
close($fh);

my @keepCols;
while (<>) {
    chomp;
    my @line = split(/ /);

    # this is the header
    unless (@keepCols) {
        @keepCols = 0 .. 2;
        for my $headIdx ( 3 .. $#line ) {
            push @keepCols, $headIdx if exists $samples{ $line[$headIdx] };
        }
    }

    my @printLine = map { $line[$_] } @keepCols;
    print join( ' ', @printLine ) . "\n";
}

__END__

=head1 NAME

doseKeepSamp.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION

Stub documentation for doseKeepSamp.pl, 
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
