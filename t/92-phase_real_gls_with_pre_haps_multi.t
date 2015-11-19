# Test file created outside of h2xs framework.
# Run this like so: `perl 20-phase_simulated_gls.t'
#   wkretzsch@gmail.com     2015/04/21 15:24:17

#########################

use Test::More;

BEGIN { plan tests => 8; }

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
use version 0.77;

my $insti = shift @ARGV;

my $insti_version = qx($insti 2>/dev/null | head -1);
$insti_version =~ m/(v\d+\S+)/;
$insti_version = qv($1);

my $sampDir = "$Bin/../samples";
my $srcDir  = "$sampDir/multi_gls";
my $resDir  = "$Bin/results/" . basename( $0, '.t' );

make_path($resDir);

# copy gls to res dir
my $gls = "$resDir/gls.bcf.gz";
my $srcGLs =
"$srcDir/chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bcf.gz";
my $glsBin = "$resDir/gls.bin";
my $srcGLsBin =
"$srcDir/chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.bin";
copy( $srcGLs,       $gls );
copy( $srcGLsBin,    $glsBin );
copy( "$srcGLs.csi", "$gls.csi" );

my $phBase =
"20.5335724-5861377.1024_site_subset.union.filteredAC5.onlyPhased.NM_HOMMAJORv3.inGLSamples.winni_filt_subset.with_multi.ordered";
my $pre_haps  = "$srcDir/$phBase.bcf.gz";
my $phTabHaps = "$srcDir/$phBase.tabhaps.gz";
my $phSamp    = "$srcDir/$phBase.sample";

my $phLoadCommand = "-h $phTabHaps -s $phSamp";

#    $phLoadCommand = "" if $insti_version >= qv("v1.6.7");

my $gMap = "$sampDir/geneticMap/genetic_map_chr20_combined_b37.txt.gz";
my $omniSamples =
"$srcDir/20.5335724-5861377.chip.omni_broad_sanger_combined.20140818.snps.genotypes.samples";
my $omniGenotypes =
"$srcDir/20.5335724-5861377.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz";

{
    ok(
        system(
"$insti -R20:5335724-5861377 -g $gMap -C100 -m 10 -B0 -i5 $phLoadCommand --gls $gls -Fbcf -o $gls"
          ) == 0,
        "ran insti"
    );
    BGZIPandIndexSTVCFGZ("$gls.vcf.gz");
    my $nrd = VCFNRD( "$gls.vcf.gz", $omniGenotypes, $resDir );
    cmp_ok( $nrd, '<', 15, "NRD is sufficiently small" );
}

# how about gls in bin format
{

    ok(
        system(
"$insti -R20:5335724-5861377 -g $gMap -C100 -m 10 -B0 -i5 $phLoadCommand --gls $glsBin -Fbin -o $glsBin"
          ) == 0,
        "ran insti"
    );
    BGZIPandIndexSTVCFGZ("$glsBin.vcf.gz");
    my $nrd = VCFNRD( "$glsBin.vcf.gz", $omniGenotypes, $resDir );
    cmp_ok( $nrd, '<', 15, "NRD is sufficiently small" );

}

# let's try not using pre-haps
my $outBase = "$gls.no_pre_haps";
{
    ok(
        system(
"$insti -R20:5335724-5861377 -g $gMap -C100 -m 10 -B0 -i5 --gls $gls -Fbcf -o $outBase"
          ) == 0,
        "ran insti"
    );
    BGZIPandIndexSTVCFGZ("$outBase.vcf.gz");
    my $nrd = VCFNRD( "$outBase.vcf.gz", $omniGenotypes, $resDir );
    cmp_ok( $nrd, '<', 15, "NRD is sufficiently small" );
}

# let's try only using biallelics
$outBase = "$gls.no_pre_haps.no_multi";
$gls     = "$resDir/gls.no_multi.bcf.gz";
$srcGLs =
"$srcDir/chr20.5335724-5861377.1024_site_subset.HRC.r1.AC5.TGPP3_samples.likelihoods.winniFilter.no_multi.bcf.gz";
copy( $srcGLs,       $gls );
copy( "$srcGLs.csi", "$gls.csi" );

{
    ok(
        system(
"$insti -R20:5335724-5861377 -g $gMap -C100 -m 10 -B0 -i5 --gls $gls -Fbcf $phLoadCommand -o $outBase"
          ) == 0,
        "ran insti"
    );
    BGZIPandIndexSTVCFGZ("$outBase.vcf.gz");
    my $nrd = VCFNRD( "$outBase.vcf.gz", $omniGenotypes, $resDir );
    cmp_ok( $nrd, '<', 15, "NRD is sufficiently small" );

}

