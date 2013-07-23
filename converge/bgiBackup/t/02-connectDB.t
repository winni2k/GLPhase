# Test file created outside of h2xs framework.
# Run this like so: `perl 01-load.t'
#   wkretzsch@gmail.com     2013/07/19 14:58:07

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More;
use FindBin;
use lib "$FindBin::Bin/..";

BEGIN { plan tests => 14 }

use warnings;
use strict;
$| = 1;
use Data::Dumper;
use bamTrackLib;

use Test::Files;
use Test::Differences;

use File::Path qw(make_path remove_tree);
use File::Copy;

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
my $wd1 = "t/01-workdir";
remove_tree($wd1) if -d $wd1;
make_path($wd1);

my $bamList = $wd1 . '/bam.list';
open( my $fh, '>', $bamList ) or die "could not open $bamList for writing";
for my $bam ( "MD_CHW_AAS_10011.head500.bam", "MD_CHW_AAS_10179.head500.bam" ) {
    copy( "samples/bams/" . $bam,                   $wd1 );
    copy( "samples/bams/" . $bam . ".md5.expected", $wd1 );
    print $fh $wd1 . '/' . $bam . "\n";
}
close($fh);

open( $fh, '<', $bamList ) or die "could not open $bamList";
my %bamList;
my @bamList;
while (<$fh>) {
    chomp;
    $bamList{$_} = 1;
    push @bamList, $_;

    # remove any md5sums that should not exist
    die "no md5 files should exist at $_.md5" if ( -e $_ . ".md5" );
}

# touch an out of date .ok file
system("touch $bamList[1].ok");
system("sleep 1 && touch $bamList[1]");

# touch an in date .ok file
system("touch $bamList[0].ok");

# now start the api
ok( $bto->registerBams( fileList => $bamList ), "registration completed" );

my @bams = sort keys %{ $bto->inputBams() };
@bamList = sort @bamList;
eq_or_diff \@bams, \@bamList, "bamlist parsed correctly";

#my %md5sums = %{ $bto->inputBamMD5sums };

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
  "bamlist saved in db and only validateSam validated bam retrieved correctly";

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

@bams =
  $bto->retrieveBams( filterColumns => { backupDevice => 'externalHD1' } );
eq_or_diff \@bams, $backupBam, "backup saved in db retrieved correctly";

### testing dropping of bams
@bams =
  $bto->retrieveBams( filterColumns => { sampleName => \@bamSampleNames } );
eq_or_diff \@bams, [ @bamList, $backupBam ],
  "backup + bamList saved in db retrieved correctly";

is( $bto->dropBams( file => $bamList[1] ), 1, "one row dropped successfully" );

@bams =
  $bto->retrieveBams( filterColumns => { sampleName => \@bamSampleNames } );
eq_or_diff \@bams, [ $bamList[0], $backupBam ],
  "backup + bamlist - dropped bam retrieved from db correctly";

## now try some validation
# create backup bam and phony md5 sum
my $backupBam2 = "$wd/MD_CHW_AAS_10011.head500.bam";
system "cp samples/bams/MD_CHW_AAS_10011.head500.bam $backupBam2";
system "cp samples/bams/MD_CHW_AAS_10011.head500.bam.md5 $backupBam2.md5";
system
"sed 's/2bea/ffff/' < $backupBam2.md5 >$backupBam2.tmp && mv $backupBam2.tmp $backupBam2.md5";
$bto->registerBams(
    file         => "$backupBam2",
    backup       => 1,
    backupDevice => 'externalHD1'
);

my $bamList2 = "$wd/bam.list";
open( $fh, '>', $bamList2 ) or die "could not open $bamList2 for writing";
map { print $fh "$_\n" } ( @bams, $backupBam, $backupBam2 );
close($fh);

my @brokenBams =
  $bto->validateBams( fileList => $bamList2, validationType => 'md5sum' );
eq_or_diff( \@brokenBams, [$backupBam2],
    "broken bam was successfully identified as broken" )

#########################

  # Insert your test code below, the Test::More module is used here so read
  # its man page ( perldoc Test::More ) for help writing this test script.

