#!/usr/bin/perl
# regression_tests.pl                   wkretzsch@gmail.com
#                                       01 Dec 2014
use warnings;
use strict;
use TAP::Harness;
use FindBin qw($Bin);
use File::Path qw(make_path remove_tree);
use File::Slurp;

remove_tree("$Bin/results") if -d "$Bin/results";
make_path("$Bin/results");

# pass along insti version
my %defines;
map { my @l = split( /\s+/, $_ ); $defines{ $l[1] } = $l[2] }
  grep { m/\#define/ } read_file("$Bin/../src/version.hpp");

my $insti = shift @ARGV;
$insti = "$Bin/$insti";

my $harness = TAP::Harness->new(
    {
        test_args => [
            $defines{VERSION_MAJOR}, $defines{VERSION_MINOR},
            $defines{VERSION_REVISION}, $insti,
        ]
    }
);

# find all files in this script's dir with a .t ending
opendir( my $dh, $Bin ) || die "can't opendir $Bin: $!";
my @tests = map { "$Bin/$_" } grep { /\.t$/ && -f "$Bin/$_" } readdir($dh);
closedir $dh;

$harness->runtests(@tests);

__END__

=head1 NAME

regression_tests.pl

=head1 SYNOPSIS
   
t/regression_tests.pl

=head1 DESCRIPTION

Runs all regression tests in t

=head1 AUTHOR

Warren Winfried Kretzschmar, E<lt>wkretzsch@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Warren Winfried Kretzschmar

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
