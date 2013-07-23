# Test file created outside of h2xs framework.
# Run this like so: `perl 01-load.t'
#   wkretzsch@gmail.com     2013/07/19 14:58:07

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More;
use FindBin;
use lib "$FindBin::Bin/..";

BEGIN { plan tests => 10 }

use warnings;
use strict;
$| = 1;
use Data::Dumper;
use bamTrackLib;

use Test::Files;
use Test::Differences;

use File::Path qw(make_path remove_tree);

my $bto = bamTrackLib->new(
    host => 'mus.well.ox.ac.uk',
    db   => "testbgibams",
    user => 'winni'
);

# remove test table if it exists
my $dbh = $bto->dbHandle;
$dbh->do("DROP TABLE IF EXISTS md5sums");
$dbh->do("DROP TABLE IF EXISTS bamNames");
$dbh->commit;

$bto = bamTrackLib->new(
    host               => 'mus.well.ox.ac.uk',
    db                 => "testbgibams",
    user               => 'winni',
    allowRelativePaths => 1,
);

isa_ok( $bto, 'bamTrackLib' );

ok( $bto->connect(), "connection ok" );

### load a bam file list
# load expected bam list
my $bamList = './samples/bam.list';

open( my $fh, '<', $bamList ) or die "could not open $bamList";
my %bamList;
my @bamList;
while (<$fh>) {
    chomp;
    $bamList{$_} = 1;
    push @bamList, $_;

    # remove any md5sums that should not exist
    if ( -e $_ . ".md5" ) { remove_tree( $_ . ".md5" ) }
}

# touch an out of date .ok file
system("touch $bamList[1].ok");
system("touch $bamList[1]");

# touch an in date .ok file
system("touch $bamList[0].ok");

# now start the api
ok( $bto->registerBams( fileList => $bamList ), "registration completed" );

my @bams = sort keys %{ $bto->inputBams() };
@bamList = sort @bamList;
eq_or_diff \@bams, \@bamList, "bamlist parsed correctly";

my %md5sums = %{ $bto->inputBamMD5sums };

#my @md5Files = sort map { $md5sums{$_} } sort keys %md5sums;
#my @expected_md5Files = sort map { $_.".md5.expected"} @bamList;

# check if md5sums match
for my $bam (@bamList) {
    compare_ok(
        $bam . '.md5',
        $bam . '.md5.expected',
        "$bam md5sum is correct"
    );
}

# validate md5sums ... maybe later
# $bto->validateInputBamMD5sums;

# make sure bams were registered correctly
my @bamSampleNames = map { m|^.*/([A-Za-z0-9_]+)|; $1; } @bamList;
@bams =
  $bto->retrieveBams( filterColumns => { sampleName => \@bamSampleNames } );
eq_or_diff \@bams, \@bamList, "bamlist saved in db and retrieved correctly";

@bams = $bto->retrieveBams(
    filterColumns => { sampleName => \@bamSampleNames, passedValidateSam => 1 }
);
eq_or_diff \@bams, $bamList[0],
  "bamlist saved in db and only validated bam retrieved correctly";

## Please see file perltidy.ERR
my $wd = "t/02-workdir";
remove_tree($wd) if -d $wd;
 make_path($wd);

my $backupBam = "$wd/MD_CHW_AAS_10179.head500.bam";
system "cp samples/bams/MD_CHW_AAS_10179.head500.bam $backupBam";

$bto->registerBams(
    file         => $backupBam,
    backup       => 1,
    backupDevice => 'externalHD1'
);
@bams = $bto->retrieveBams( filterColumns => { backup => 1 } );
eq_or_diff \@bams, $backupBam, "backup saved in db retrieved correctly";

@bams = $bto->retrieveBams( filterColumns => { backupDevice => 'externalHD1' } );
eq_or_diff \@bams, $backupBam, "backup saved in db retrieved correctly";


#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

