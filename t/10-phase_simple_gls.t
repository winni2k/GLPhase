# Test file created outside of h2xs framework.
# Run this like so: `perl 10-phase_easy_gls.t'
#   wkretzsch@gmail.com     2014/12/01 14:56:15

#########################

use Test::More;
BEGIN { plan tests => 8 }

use warnings;
use strict;
$| = 1;
use Data::Dumper;
use FindBin qw($Bin);
use File::Path qw(make_path remove_tree);
use File::Copy;
use File::Basename;
use File::Slurp;
use lib $Bin;
use VCFComp ':all';

my $insti     = shift @ARGV;
my $simpleDir = "$Bin/../samples/simple_gls";
my $resDir    = "$Bin/results/" . basename( $0, '.t' );

make_path($resDir);

# copy gls to res dir
my $simpleGLBase = "$resDir/simple.gls.v1";
my $simpleGLs    = "$simpleGLBase.bin";
copy( "$simpleDir/simple.gls.v1.bin", $simpleGLs );

my $gMap = "$Bin/../samples/geneticMap/genetic_map_chr20_combined_b37.txt.gz";

ok( system("$insti -g $gMap -C100 -m 100 -B0 -i10 $simpleGLs") == 0,
    "ran insti" );
BGZIPandIndexSTVCFGZ("$simpleGLs.vcf.gz");

my $code = VCFHapMatch( "$simpleGLs.vcf.gz",
    "$simpleDir/simple.gls.v1.expected.bin.vcf", $resDir );
ok( $code eq 0, "simple haps v1" ) or diag($code);

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
# test --region option
{
    my $base = "$simpleGLBase.region";
    ok(
        system(
            "$insti -C100 -m 5 -B0 -i5 --region 20:101363- -o $base $simpleGLs")
          == 0,
        "ran insti on region 20:101363-"
    );
    BGZIPandIndexSTVCFGZ("$base.vcf.gz");

    $code = VCFHapMatch( "$base.vcf.gz",
        "$simpleDir/simple.gls.v1.region.expected.bin.vcf", $resDir );
    ok( $code eq 0, "simple haps v1 only region exit code OK" ) or diag($code);
}

###
# test --samples-file option
TODO: {
    local $TODO = "can't unambiguosly phase just 2 samples";
    
    my $base = "$simpleGLBase.samplesFile";
    write_file( "$base.samples", "samp2\nsamp3\n" );
    ok(
        system(
            "$insti -C100 -m 5 -B0 -i5 -S $base.samples -o $base $simpleGLs")
          == 0,
        "ran insti on no samp1"
    );
    BGZIPandIndexSTVCFGZ("$base.vcf.gz");

    $code = VCFHapMatch( "$base.vcf.gz",
        "$simpleDir/simple.gls.v1.noSamp1.expected.bin.vcf", $resDir );
    ok( $code eq 0, "simple haps v1 only samp2 exit code OK" ) or diag($code);
}
