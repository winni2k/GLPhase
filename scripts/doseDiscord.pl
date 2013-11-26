#!/usr/bin/perl
# doseDiscord.pl                   wkretzsch@gmail.com
#                                  16 Oct 2013
our $VERSION = '0.002';
$VERSION = eval $VERSION;
print STDERR "doseDiscord.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

use warnings;
use strict;
$| = 1;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;

use Getopt::Std;
my %args;
getopts( 'T:t:f:s:e:', \%args );
my $truthFile = $args{T} || croak "Need to submit truth file -T";
my $testFile  = $args{t} || croak "Need to submit test file -t";
my $frqFile   = $args{f};
my $siteDiscFile = $args{s} =~ m/\.gz$/
  || !defined $args{s} ? $args{s} : croak "-s must have .gz ending";
my $exclusionSiteList = $args{e};

print STDERR "outputting bySite discordance\n";

my $frq;
if ( defined $args{f} ) {
    print STDERR "Loading frequency file\n";
    $frq = load_frq( $args{f} );
}

my $excList;
if ( defined $exclusionSiteList ) {
    print STDERR "Loading site exclusion list\n";
    $excList = load_list($exclusionSiteList);

    die unless defined $frq;

    # remove exclusion list from frq
    my $deleteNum = 0;
    for my $key ( @{$excList} ) {
        if ( exists $frq->{$key} ) {
            delete $frq->{$key};
            $deleteNum++;
        }
    }
    print STDERR "\tremoved $deleteNum sites\n";

}

print STDERR "Loading dose files\n";
my $truth = load_dose($truthFile);
my $test  = load_dose($testFile);

print STDERR "Intersecting sites files\n";
( $truth, $test ) = intersect_sites( $truth, $test, $frq );

print STDERR "Calculating discordance\n";
my ( $discord, $numDoses, $numDiscords, $rhSiteDisc ) =
  measure_discordance( $truth, $test );

print
"Discordance: $discord\nNumber of Doses: $numDoses\nNumber of discordances: $numDiscords\n";

if ( defined $siteDiscFile ) {
    print STDERR "printing sitewise discordance\n";
    open( my $fh, " | gzip -c > $siteDiscFile" );
    print $fh q/POS REF ALT eRR tRR eRA tRA eAA tAA eALT tALT eALL tALL MIS/
      . "\n";
    for my $siteKey (
        sort { $truth->{sites}->{$a}->{pos} <=> $truth->{sites}->{$b}->{pos} }
        keys %{$rhSiteDisc} )
    {
        my $site = $truth->{sites}->{$siteKey};
        print $fh join(
            ' ',
            (
                $site->{pos},
                @{ $site->{alleles} },
                @{ $rhSiteDisc->{$siteKey} }
            )
        ) . "\n";
    }
}

sub load_list {

    my $file = shift;
    my $openCmd = $file =~ m/\.gz$/ ? " gzip -dc $file | " : " < $file ";
    open( my $fh, $openCmd ) or croak "could not open file $file; $!";
    my @list;
    while (<$fh>) {
        chomp;
        push @list, $_;
    }
    close($fh);
    return \@list;
}

sub load_frq {

    my $file = shift;
    my $openCmd = $file =~ m/\.gz$/ ? " gzip -dc $file | " : " < $file ";
    open( my $fh, $openCmd ) or croak "could not open file $file; $!";
    my %frq;
    while (<$fh>) {
        chomp;
        my @line = split(/ /);
        $frq{ $line[0] } = $line[1];
    }
    close($fh);
    return \%frq;
}

sub measure_discordance {

    my ( $truthDose, $testDose ) = @_;
    croak "dose files do not have the same number of sites"
      unless keys %{ $truthDose->{sites} } == keys %{ $testDose->{sites} };
    for my $site ( sort keys %{ $truthDose->{sites} } ) {
        croak "dose files do not have the same sites"
          unless exists $testDose->{sites}->{$site};
    }
    croak "sample numbers do not match"
      unless @{ $truthDose->{samples} } == @{ $testDose->{samples} };

    # testDose to make sure samples match
    for my $sampIdx ( 0 .. $#{ $truthDose->{samples} } ) {
        croak "dose samples do not match at samples:\n"
          . $truthDose->{samples}->[$sampIdx] . " and "
          . $testDose->{samples}->[$sampIdx] . "\n"
          . join( ' ', @{ $truthDose->{samples} } ) . "\n"
          . join( ' ', @{ $testDose->{samples} } ) . "\n"
          unless $testDose->{samples}->[$sampIdx] eq
          $truthDose->{samples}->[$sampIdx];
    }

    # now calculate discordance overall
    my $numDoses = 0;
    my $discords = 0;
    my %siteDisc;
    for my $siteKey ( sort keys %{ $truthDose->{sites} } ) {
        my %siteNumDoses = ( 0 => 0, 1 => 0, 2 => 0, NA => 0 );
        my %siteNumDiscords = ( 0 => 0, 1 => 0, 2 => 0 );
        my @truthDoses = @{ $truthDose->{sites}->{$siteKey}->{doses} };
        my @testDoses  = @{ $testDose->{sites}->{$siteKey}->{doses} };
        croak "number of doses do not match" unless @truthDoses == @testDoses;
        for my $idx ( 0 .. $#truthDoses ) {
            $siteNumDoses{ $truthDoses[$idx] }++;
            next if ( $truthDoses[$idx] eq 'NA' || $testDoses[$idx] eq 'NA' );
            $siteNumDiscords{ $truthDoses[$idx] }++
              if $truthDoses[$idx] != $testDoses[$idx];
        }
        my @return = map { ( $siteNumDiscords{$_}, $siteNumDoses{$_} ) } 0 .. 2;
        my $eALT   = $return[2] + $return[4];
        my $tALT   = $return[3] + $return[5];
        my $eALL   = $return[0] + $eALT;
        my $tALL   = $return[1] + $tALT;
        push @return, ( $eALT, $tALT, $eALL, $tALL, $siteNumDoses{NA} );
        $siteDisc{$siteKey} = \@return;
        $discords += $eALL;
        $numDoses += $tALL;
    }

    return ( $discords / $numDoses, $numDoses, $discords, \%siteDisc );

}

sub intersect_sites {

    my ( $dose1, $dose2, $frq ) = @_;

    my $compNum = defined $frq ? 3 : 2;
    my %sites;
    map { $sites{$_}++ }
      ( keys %{ $dose1->{sites} }, keys %{ $dose2->{sites} } );
    map { $sites{$_}++ } keys %{$frq} if defined $frq;

    my @nonInterSites = grep { $sites{$_} != $compNum } sort keys %sites;
    map { delete $dose1->{sites}->{$_} if exists $dose1->{sites}->{$_} }
      @nonInterSites;
    map { delete $dose2->{sites}->{$_} if exists $dose2->{sites}->{$_} }
      @nonInterSites;

    return ( $dose1, $dose2 );
}

sub load_dose {

    my $file = shift;
    my %object = ( samples => [], sites => {} );

    my $openCmd = $file =~ m/\.gz$/ ? " gzip -dc $file | " : " < $file ";
    open( my $fh, $openCmd ) or croak "could not open file $file; $!";
    my $firstLine = 1;
    my @reorderedIdx;
    while (<$fh>) {
        chomp;
        my @line = split(/ /);
        if ($firstLine) {
            my @samples = @line[ 3 .. $#line ];
            $firstLine = 0;

            # create list of indexes that sort names
            @reorderedIdx =
              sort { $samples[$a] cmp $samples[$b] } 0 .. $#samples;
            @samples = map { $samples[$_] } @reorderedIdx;
            $object{samples} = \@samples;
            next;
        }
        my $id = shift @line;
        my @alleles = splice( @line, 0, 2 );
        croak "$id exists twice in file $file" if exists $object{sites}->{$id};
        my ( $chrom, $position ) = split( /\:/, $id );

        # sort doses according to sample name
        @line = map { $line[$_] } @reorderedIdx;
        $object{sites}->{$id} = {
            chrom   => $chrom,
            pos     => $position,
            alleles => \@alleles,
            doses   => \@line
        };
    }
    return \%object;
}

__END__

=head1 NAME

doseDiscord.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for doseDiscord.pl, 
created by template.el.

It looks like the author of this script was negligent
enough to leave the stub unedited.


=head1 AUTHOR

Warren Winfried Kretzschmar, E<lt>wkretzsch@gmail.comE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Warren Winfried Kretzschmar

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
