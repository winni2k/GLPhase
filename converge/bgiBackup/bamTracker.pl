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
                     addBamList=>0,

                     # input/output files
                     bamList=>0,
    }
);

my $bto = bamTrackLib->new( db=>$args{db}, host=>$args{host}, user=>$args{user});

$bto->registerBams(fileList => $args{addBamList}) if $args{addBamList};

$bto->listDrivers() if $args{listDrivers};




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
