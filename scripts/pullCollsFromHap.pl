#!/usr/bin/perl
# pullCollsFromHap.pl                   wkretzsch@gmail.com
#                                       18 Oct 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "pullCollsFromHap.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

use warnings;
use strict;
$| = 1;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;
use autodie;

use Getopt::Std;
my %args;
getopts( 'drj:', \%args );
my $DEBUG = $args{d} || 1;

my $sampCols = shift;
open( my $fh, "<", $sampCols );
my @cols;
while (<$fh>) {
    chomp;
    push @cols, $_;
}

@cols = map { $_ - 3 } @cols;

while (<>) {
    chomp;
    my @line = split(' ');
    my @out = splice( @line, 0, 5 );
    map { push @out, @line[ 2 * $_ .. ( 2 * $_ + 1 ) ] } @cols;
    print join( ' ', @out ) . "\n";
}

__END__

=head1 NAME

pullCollsFromHap.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for pullCollsFromHap.pl, 
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
