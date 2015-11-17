package VCFComp;

# wkretzsch@gmail.com          08 Dec 2014

use 5.006;
use strict;
use warnings;
use Carp;
use Data::Dumper;
use File::Basename;
use File::Slurp;
use Carp;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use VCFComp ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = (
    'all' => [
        qw(
          VCFHapMatch VCFNRD BGZIPandIndexSTVCFGZ
          )
    ]
);

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw(  );

our $VERSION = '0.01';

# Preloaded methods go here.

sub GetFH {
    my $file = shift;
    open( my ($fh), $file =~ m/\.gz$/ ? "gzip -dc $file |" : "<$file" )
      or die "could not open file $file";
    return $fh;
}

sub HapEq {
    my $tHap   = shift;
    my $eHap   = shift;
    my $tFH    = GetFH($tHap);
    my @tLines = <$tFH>;
    my $eFH    = GetFH($eHap);
    my @eLines = <$eFH>;

    return "hap files have different lengths: $tHap, $eHap"
      if @tLines != @eLines;

    my @flippedToMatch;
    for my $lNum ( 0 .. $#tLines ) {
        my $lineNum = $lNum + 1;
        my $tLine   = $tLines[$lNum];
        chomp $tLine;
        my @tAll = split( / /, $tLine );

        my $eLine = $eLines[$lNum];
        chomp $eLine;
        my @eAll = split( / /, $eLine );

        return "hap files have differing number of haps at line $lineNum"
          if @tAll != @eAll;
        return "hap files have 0 haplotypes" if @tAll == 0;
        return "hap files have uneven number of haplotypes" if @tAll % 2 != 0;

        @flippedToMatch = map { 0 } 0 .. ( $#tAll / 2 ) if @flippedToMatch == 0;

      ALLELE: for my $tHapNum ( 0 .. $#tAll ) {
            my $sampNum        = int( $tHapNum / 2 );
            my $firstHapOfPair = $sampNum * 2 == $tHapNum;
            my $eHapNum =
                $flippedToMatch[$sampNum] == 0
              ? $tHapNum
              : ( $firstHapOfPair ? $tHapNum + 1 : $tHapNum - 1 );
            if ( $tAll[$tHapNum] ne $eAll[$eHapNum] ) {
                if ( $flippedToMatch[$sampNum] == 0 && $firstHapOfPair ) {
                    $flippedToMatch[$sampNum] = 1;
                    redo ALLELE;
                }
                return
                  "Haplotypes do not match at line $lineNum, at Allele number "
                  . ( $tHapNum + 1 )
                  . ":\nTest line:\n$tLine"
                  . "\nExpected line:\n$eLine";

            }
        }
    }

    return 0;
}

sub VCFHapMatch {
    my $testVcf = shift;
    my $expVcf  = shift;
    my $wd      = shift;

    my $tBase = basename( $testVcf, '.vcf.gz' );
    my $eBase = basename( $expVcf,  '.vcf.gz' );
    system("bcftools convert -h $wd/$tBase $testVcf 2>/dev/null")
      and die "could not convert file $wd/$tBase";
    system("bcftools convert -h $wd/$eBase $expVcf 2>/dev/null")
      and die "could not convert file $wd/$tBase";
    return HapEq( "$wd/$tBase.hap.gz", "$wd/$eBase.hap.gz" );
}

sub VCFNRD {
    my $testVCF = shift;
    my $expVCF  = shift;
    my $wd      = shift;
    my $tBase   = basename( $testVCF, '.vcf.gz' );
    system( "bcftools stats $expVCF $testVCF -s " . '$'
          . "(bcftools query -l $expVCF"
          . q( | perl -0 -ane 'print join(",", @F)')
          . " ) > $wd/$tBase.stats" );
    my @lines = read_file("$wd/$tBase.stats");
    @lines = grep { m/^NRDs/ } @lines;
    my @line = split /\t/, $lines[0];

    $line[2] =~ m/^\d*\.?\d+$/ or die "NRD ($line[2]) is not a float";
    return $line[2];
}

sub BGZIPandIndexSTVCFGZ {
    my $vcf = shift;
    system(
"gzip -dc $vcf | bgzip -c > $vcf.tmp && mv -f $vcf.tmp $vcf && bcftools index $vcf"
    ) and die "could not bgzip and index $vcf";
}

1;
__END__
# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

VCFComp - Perl extension for blah blah blah

=head1 SYNOPSIS

   use VCFComp;
   blah blah blah

=head1 DESCRIPTION

Stub documentation for VCFComp, 
created by template.el.

It looks like the author of the extension was negligent
enough to leave the stub unedited.

Blah blah blah.

=head2 EXPORT

None by default.

=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Warren Winfried Kretzschmar

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut

