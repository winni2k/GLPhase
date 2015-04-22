# Test file created outside of h2xs framework.
# Run this like so: `perl 20-phase_simulated_gls.t'
#   wkretzsch@gmail.com     2015/04/21 15:24:17

#########################

use Test::More;
BEGIN { plan tests => 2 }

use warnings;
use strict;
$| = 1;
use Data::Dumper;
use FindBin qw($Bin);
use File::Path qw(make_path remove_tree);
use File::Copy;
use File::Basename;
use lib $Bin;
use VCFComp ':all';

my $insti  = shift @ARGV;
my $srcDir = "$Bin/../samples/hapGen";
my $resDir = "$Bin/results/" . basename( $0, '.t' );

make_path($resDir);

# copy gls to res dir
my $glBase = "$resDir/simulated_gls";
my $gls    = "$glBase.bin";
copy( "$srcDir/ex.bin", $gls );

my $gMap = "$srcDir/ex.map";

ok( system("$insti -g $gMap -C100 -m 100 -B5 -i5 $gls") == 0,
    "ran insti" );
BGZIPandIndexSTVCFGZ("$gls.vcf.gz");

my $nrd = VCFNRD( "$gls.vcf.gz", "$srcDir/ex.vcf.gz", $resDir );
cmp_ok( $nrd, '<', 5, "simulated hap gen NRD < 5" );

#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.

