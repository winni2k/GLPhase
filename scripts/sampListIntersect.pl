#!/usr/bin/perl
# sampListIntersect.pl                   wkretzsch@gmail.com
#                                        03 Oct 2013

use warnings;
use strict;
$| = 1;
use Data::Dumper;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);

use Getopt::Std;
my %args;
getopts( 'iu', \%args );

die "Please only specify -u or -i" if $args{i} && $args{u};
die "Need to specify -u or -i" unless $args{i} || $args{u};

my @sampLists = @ARGV;
my $sListNum  = 0;
my @ahSamples;
foreach my $sampList (@sampLists) {
    $sListNum++;
    open my $fh, '<', $sampList;
    my %samples;
    while (<$fh>) {
        chomp;
        $samples{$_} = 1;
    }
    push @ahSamples, \%samples;
}

my %gSamples;
for my $rhSamples (@ahSamples) {
    map { $gSamples{$_}++ } keys %{$rhSamples};
}

my $numSampLists = @sampLists;
if ( $args{i} ) {
    map { print STDOUT "$_\n" if $gSamples{$_} == $numSampLists }
      sort keys %gSamples;
}
elsif ( $args{u} ) {
    map { print STDOUT "$_\n" } sort keys %gSamples;
}
else {
    die "unexpected error";
}

__END__

=head1 NAME

sampListIntersect.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for sampListIntersect.pl, 
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
