# Test file created outside of h2xs framework.
# Run this like so: `perl 01-load.t'
#   wkretzsch@gmail.com     2013/07/19 14:58:07

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More;
use FindBin;
use lib "$FindBin::Bin/..";

BEGIN { plan tests => 13 }

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

# now start the api
ok( $bto->registerBams($bamList), "registration completed" );

my @bams = sort keys %{ $bto->inputBams() };
@bamList = sort @bamList;
eq_or_diff \@bams, \@bamList, "bamlist parsed correctly";

my %md5sums = %{ $bto->inputBamMD5sums };

#my @md5Files = sort map { $md5sums{$_} } sort keys %md5sums;
#my @expected_md5Files = sort map { $_.".md5.expected"} @bamList;

# check if md5sums match
for my $bam (@bamList) {
    compare_ok( $bam . '.md5', $bam . '.md5.expected', "$bam md5sum is correct" );
}

TODO: {

    todo_skip "not implemented yet", 6;

    # make sure bams were registered correctly
    @bams = $bto->retrieveAllBams;

    my %bams;
    for my $bam (@bams) {
        $bams{$bam} = 1;
    }

    for my $bam ( sort keys %bamList ) {
        ok( exists $bams{$bam}, "$bam was registered in db" );
    }

    # now create any missing md5sums
    ok( $bto->updateMD5sums(), "md5sums created successfully" );

    # check if md5sums match
    for my $bam ( sort keys %bamList ) {
        compare_ok(
            $bam . '.md5',
            $bam . '.md5.expected',
            "$bam md5sum is correct"
        );
    }

    # now validate md5sums
    ok( $bto->validateMD5sums(), "md5sums validated" );
}

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

