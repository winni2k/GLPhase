#!/usr/bin/perl
# bamTracker.pl                   wkretzsch@gmail.com
#                                 19 Jul 2013

use warnings;
use strict;
$| = 1;
use Data::Dumper;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use ParseArgs 0.000002 qw/getCommandArguments/;
use POSIX;
use DBI;

my %args = getCommandArguments(
    requiredArgs => {

        #        WD                    => undef,
        DB => "/Net/sparse/data/BGI/bamDB/",
    }
);

print "Available DBI drivers:\n";
my @drivers = DBI->available_drivers('quiet');
my @sources;

foreach my $driver (@drivers) {
    print "$driver\n";
    @sources = eval { DBI->data_sources($driver) };
    if ($@) {
        print "\tError: $@\n";
    }
    elsif (@sources) {
        foreach my $source (@sources) {
            print "\t$source\n";
        }
    }
    else {
        print "\tno known data sources\n";
    }
}

__END__

=head1 NAME

bamTracker.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION

This script shall be used to keep track of BAM files, 
and perform operations on BAM files while noting down what has been done on them.

Uses may include:
- backup
- tracking of broken files
- creating lists of bam files matching certain criteria

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
