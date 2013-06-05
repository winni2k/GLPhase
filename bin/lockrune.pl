#!/usr/bin/perl -w

use strict;
use Fcntl qw(:flock);
use Getopt::Long;

my $help      = 0;
my $lockFile  = undef;
my $sleepTime = 60;
my $result    = GetOptions(
    'help'        => \$help,
    'lockFile=s'  => \$lockFile,
    'sleepTime=i' => \$sleepTime,
);

if ( !$result || !defined $lockFile || $help || @ARGV == 0 ) {
    die "Error: Missing Arguments\n\n"
      . "\tSynopsis:\n"
      . "\t\tlockrune.pl [ --sleepTime integer ] --lockFile file -- command\n"
      . "\t\tlockrune.pl --help\n\n"
      . "\tArguments:\n"
      . "\t\t--lockFile [ file ]\n\t\t--sleepTime [ integer ]\n"
      . "\t\t-- [ command ] the command to be run after --sleepTime seconds.\n"
      . "\t\t--help prints this message.";
}

my @command = @ARGV;
open( my $fh_lockfile, ">", $lockFile )
  or die "Can't open lockfile $lockFile for reading";
flock( $fh_lockfile, LOCK_EX ) or die "Cannot lock $lockFile";

# fork it
my $pid = fork();
if ( $pid == 0 ) {

    # child, run command
    system @command;
    my $rc = $? >> 8;
    exit $rc;
}
elsif ($pid) {

    # parent, sleep $sleepTime seconds, then release file lok
    sleep $sleepTime;
    flock( $fh_lockfile, LOCK_UN ) or die "Cannot unlock $lockFile";

    # return exit status of child
    waitpid( $pid, 0 );
    my $rc = $? >> 8;
    exit $rc;
}
else {
    # something went wrong in fork
    die "could not fork: $!\n";
}
