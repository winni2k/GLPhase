# Test file created outside of h2xs framework.
# Run this like so: `perl 20-phase_simulated_gls.t'
#   wkretzsch@gmail.com     2015/04/21 15:24:17

#########################

use Test::More;
BEGIN { plan tests => 6 }

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

# let's run on the multiallelic version of the VCF file
my $gls = "$glBase.gls.multi.vcf.gz";
copy( "$srcDir/ex.gls.multi.vcf.gz", $gls ) or die "could not copy";
copy( "$srcDir/ex.gls.multi.vcf.gz.csi", "$gls.csi" ) or die "could not copy";

ok( system("$insti -C100 -m 20 -B5 -i5 -o $gls -Fbcf $gls") == 0, "ran insti on multiallelics" );
BGZIPandIndexSTVCFGZ("$gls.vcf.gz");

my $nrd = VCFNRD( "$gls.vcf.gz", "$srcDir/ex.multi.vcf.gz", $resDir );
cmp_ok($nrd , '>', 0, "simulated hap gen NRD ($nrd) > 0 in multiallelics" );
cmp_ok( $nrd, '<', 5.5, "simulated hap gen NRD ($nrd) < 5.5 in multiallelics" );

# let's run the almost multiallelic version of the bin file
$gls = "$glBase.gls.almost_multi.vcf.gz";
copy( "$srcDir/ex.gls.almost_multi.vcf.gz", $gls ) or die "could not copy";
copy( "$srcDir/ex.gls.almost_multi.vcf.gz.csi", "$gls.csi" ) or die "could not copy";

ok( system("$insti -C100 -m 20 -B5 -i5 -o $gls -Fbcf $gls") == 0, "ran insti on almost_multiallelics" );
BGZIPandIndexSTVCFGZ("$gls.vcf.gz");

$nrd = VCFNRD( "$gls.vcf.gz", "$srcDir/ex.almost_multi.vcf.gz", $resDir );
cmp_ok($nrd , '>', 0, "simulated hap gen NRD ($nrd) > 0 in almost_multiallelics" );
cmp_ok( $nrd, '<', 5.5, "simulated hap gen NRD ($nrd) < 5.5 in almost_multiallelics" );


