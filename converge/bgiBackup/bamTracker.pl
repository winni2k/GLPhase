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
use bamTrackLib;

my %args = getCommandArguments(
    requiredArgs => {

        WD       => undef,

                     # database options
        host     => 'mus.well.ox.ac.uk',
        db       => "bgibams",
        user     => 'winni',
        
                    # diagnosis
                     listDrivers=> 0,

                     # commands
                     addBams=>0,

                     # input/output files
                     bamList=>0,
    }
);

die "Please export your password to DBI_PASS environment var before running this script" unless defined $ENV{DBI_PASS};

my $bto = bamTrackLib->new(operations=>[], inputFile=>'testing.list', db=>$args{db}, host=>$args{host}, user=>$args{user});

$bto->listDrivers() if $args{listDrivers};

my $dbh=DBI->connect('dbi:mysql:'.$args{db}.':'.$args{host}, $args{user}) || die "Error opening database: $DBI::errstr";

$bto->dbHandle($dbh);exit;

# load in bam files from 
#addBams($args{bamList}) if $args{addBams};



my $sth = $dbh->prepare("select * from pet where species = 'bird' or species = 'cat';") || die "Prepare failed: $DBI::errstr\n";

$sth->execute() || die "Could not execute query: $DBI::errstr\n";

while(my ( $id, $name) = $sth->fetchrow_array){
    print "$name has $id\n";
}

$sth->finish();

# disconnect
$dbh->disconnect() || die "Failed to disconnect";

#sub addBam

sub listDrivers {
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
