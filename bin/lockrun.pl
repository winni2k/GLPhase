#!/usr/bin/perl -w

use strict;
use Fcntl qw(:flock);

my ($lockfile, @command) = @ARGV;

open(my $fh_lockfile, ">", $lockfile) or die "Can't open lockfile $lockfile for reading";
flock($fh_lockfile, LOCK_EX) or die "Cannot lock $lockfile";
system @command;
flock($fh_lockfile, LOCK_UN) or die "Cannot unlock $lockfile";
