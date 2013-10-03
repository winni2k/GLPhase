#!/usr/bin/perl
# vcfReplaceHeaders.pl                   wkretzsch@gmail.com
#                                        01 Oct 2013

use warnings;
use strict;
$| = 1;
use Data::Dumper;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;

=head1 DESCRIPTION

  -c <file> Input conversion table (required).
  -p        allow no matches to conversion table.

=cut

use Getopt::Std;
my %args;
getopts( 'c:p', \%args );
my $convTable = $args{c} || confess "conversion table needed (-c)";

my @search;
my @replace;
open my $fh, '<', $convTable;
while (<$fh>) {
    chomp;
    confess "Unexpected line format in conversion table in line:\n$_"
      unless m/^[0-9A-Za-z\-\_\.]+\s[0-9A-Za-z\-\_\.]+\s*$/;
    my @line = split(/\s/);
    push @search,  $line[0];
    push @replace, $line[1];
}

for my $ra ( \@search, \@replace ) {
    my %seen = ();
    $seen{$_}++ for @{$ra};
    my @broken = grep { $seen{$_} > 1 } @{$ra};
    confess "these samples exist more than once:\n"
      . join( "\n", @broken ) . "\n"
      if @broken;
}

while (<>) {
    print $_ and next unless m/^#CHROM/;
    chomp;
    my $line = $_;
    for my $idx ( 0 .. $#search ) {
        my $search = $search[$idx];
        my $rep    = $replace[$idx];
        my $subNum = $line =~ s/\t${search}\t/\t${rep}\t/g;
        confess "multiple substitutions made on search string $search"
          if $subNum > 1;
        confess "no substitutions made on search string $search"
          if $subNum == 0 && !$args{p};
    }
    print $line. "\n";
}

__END__

=head1 NAME

vcfReplaceHeaders.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for vcfReplaceHeaders.pl, 
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
