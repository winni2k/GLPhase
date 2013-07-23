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
use Carp;
use DM;

my %args = getCommandArguments(
    requiredArgs => {

        WD => undef,

        # database options
        host => 'mus.well.ox.ac.uk',
        db   => "bgibams",
        user => 'winni',

        # diagnosis
        listDrivers => 0,

        # commands
        addBamList       => 0,
        removeBamList    => 0,
        validateBamList  => 0,
        backupBamList    => 0,
        backupTargetDir  => 0,
        backupDeviceName => 0,

        # validation jar file
        validateSamJar => 0,

        # DM related commands
        dryRun  => 1,
        numJobs => 1,
    }
);

my $dm = DM->new( dryRun => $args{dryRun}, numJobs => $args{numJobs} );

my $bto =
  bamTrackLib->new( db => $args{db}, host => $args{host}, user => $args{user} );

$bto->registerBams( fileList => $args{addBamList} ) if $args{addBamList};
$bto->dropBams( fileList => $args{removeBamList} ) if $args{removeBamList};
$bto->validateBams( fileList => $args{validateBamList} )
  if $args{validateBamList};
$bto->listDrivers() if $args{listDrivers};

if ( $args{backupBamList} ) {
    croak "need to specify backupTargetDir if backing up"
      unless $args{backupTargetDir};
    croak "need to specify backupDeviceName if backing up"
      unless $args{backupDeviceName};
    croak "need to specify validateSamJar if backing up"
      unless $args{validateSamJar};

    # register bams to be backed up
    $bto->backupBams(
        fileList         => $args{backupBamList},
        backupTargetDir  => $args{backupTargetDir},
        backupDeviceName => $args{backupDeviceName},
        validateSamJar   => $args{validateSamJar},
    );

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
