# Test file created outside of h2xs framework.
# Run this like so: `perl 10-phase_easy_gls.t'
#   wkretzsch@gmail.com     2014/12/01 14:56:15

#########################

use Test::More;
BEGIN { plan tests => 5 }

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

my $insti        = shift @ARGV;
my $simpleDir    = "$Bin/../samples/simple_gls";
my $resDir       = "$Bin/results/" . basename( $0, '.t' );

make_path($resDir);

# copy gls to res dir
my $simpleGLBase = "$resDir/simple.gls.v1.samePos";
my $simpleGLs    = "$simpleGLBase.bin";
copy( "$simpleDir/simple.gls.v1.samePos.bin", $simpleGLs );

my $gMap = "$Bin/../samples/geneticMap/genetic_map_chr20_combined_b37.txt.gz";

ok( system("$insti -g $gMap -C100 -m 100 -B0 -i10 $simpleGLs") == 0,
    "ran insti" );
BGZIPandIndexSTVCFGZ("$simpleGLs.vcf.gz");

my $code = VCFHapMatch( "$simpleGLs.vcf.gz",
    "$simpleDir/simple.gls.v1.expected.bin.vcf", $resDir );
ok( $code eq 0, "simple haps v1 same pos" ) or diag($code);

###
# do same without genetic map
my $nogmapBase = "$simpleGLBase.noGmap";
ok( system("$insti -C100 -m 100 -B0 -i10 -o $nogmapBase $simpleGLs") == 0,
    "ran insti without GMap" );
BGZIPandIndexSTVCFGZ("$nogmapBase.vcf.gz");

$code = VCFHapMatch( "$nogmapBase.vcf.gz",
    "$simpleDir/simple.gls.v1.expected.bin.vcf", $resDir );
ok( $code eq 0, "simple haps v1 ngmap exit code OK" ) or diag($code);

###
# try running on a broken input file
ok( system("$insti $simpleDir/simple.gls.v2.samePos.err.bin 2>/dev/null") != 0);
