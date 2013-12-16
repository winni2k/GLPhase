#!/usr/bin/perl
# hapLegSampFilter.pl                   wkretzsch@gmail.com
#                                       12 Dec 2013
our $VERSION = '0.002';
$VERSION = eval $VERSION;
print STDERR "hapLegSampFilter.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

=head1 CHANGES

0.002  Mon Dec 16 13:30:16 GMT 2013
       Removed unnecessary -r and -j options

=cut

package Runall;
use warnings;
use strict;
use Moose;
use MooseX::StrictConstructor;
use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;
use List::Util qw/sum/;
$| = 1;

#input variables
has 'hap'       => ( is => 'rw', default => q// );
has 'leg'       => ( is => 'rw', default => q// );
has 'samp'      => ( is => 'rw', default => q// );
has 'keepPop'   => ( is => 'rw', default => q// );
has 'keepGroup' => ( is => 'rw', default => q// );
has 'keepSex'   => ( is => 'rw', default => q// );
has 'keepSites' => ( is => 'rw', default => q// );
has base        => ( is => 'rw', default => 'out' );

# constants
has _filterFiles =>
  ( is => 'ro', isa => 'HashRef[ArrayRef]', builder => '_build_filterFiles' );
has _inputFiles => (
    is      => 'ro',
    isa     => 'HashRef[Bool]',
    default => sub { { hap => 1, samp => 1, leg => 1 } }
);

# build map of file dependencies for each filtering option
sub _build_filterFiles {
    my %return;
    for (qw/keepPop keepGroup keepSex/) {
        $return{$_} = [qw/hap samp/];
    }
    for (qw/keepSites/) {
        $return{$_} = [qw/hap leg/];
    }
    return \%return;
}

sub outHap {
    my $self = shift;
    return $self->base . ".hap.gz";
}

sub outSamp {
    my $self = shift;
    return $self->base . ".samp.gz";
}

sub outLeg {
    my $self = shift;
    return $self->base . ".leg.gz";
}

# internal variables
has '_requiredFiles' => (
    is  => 'rw',
    isa => 'HashRef',
);
has keepHapCols  => ( is => 'rw', isa => 'ArrayRef[Bool]' );
has keepHapRows  => ( is => 'rw', isa => 'ArrayRef[Bool]' );
has keepSitesMap => ( is => 'rw', isa => 'HashRef[Int]' );

around keepSitesMap => sub {
    my $orig = shift;
    my $self = shift;

    return $self->$orig() if defined $self->$orig();

    confess
      "programming error: keepSites is undefined eventhough it should exist."
      unless $self->keepSites =~ m/./;

    my %sites;
    open( my $fh, '<', $self->keepSites );
    while (<$fh>) {
        chomp;
        croak "keepSites should have only a single column" if $_ =~ m/\s/;
        $sites{$_} = $_;
    }
    close($fh);

    return $self->$orig( \%sites );
};

# create array of sites to keep
around keepHapRows => sub {
    my $orig = shift;
    my $self = shift;

    return $self->$orig() if defined $self->$orig();

    # use wc to read number of rows
    my @keepSites;
    unless ( $self->leg =~ m/./ ) {
        my $cmd;
        if ( $self->hap =~ m/\.gz$/ ) {
            $cmd = "gzip -dc " . $self->hap . " | wc -l";
        }
        else {
            $cmd = " wc -l " . $self->hap;
        }
        my $numLines = qx/$cmd/;
        my @line = split( /\s/, $numLines );
        $numLines = $line[0];
        croak "weird wc. hap empty?" unless $numLines > 0;
        @keepSites = map { 1 } 0 .. ( $numLines - 1 );
    }
    else {
        my $openCmd = $self->openCmd( $self->leg, '<' );
        open( my $inFH, $openCmd );
        my $head = <$inFH>;
        chomp $head;
        croak "Unexpected leg header '$head' in leg file: " . $self->leg
          unless $head =~ m/^id position a0 a1/;

        my $rowNum = 0;
        while (<$inFH>) {
            chomp;
            my @line = split(' ');
            push @keepSites, $self->keepLegLine( \@line );
        }
    }
    print "Keeping " . sum(@keepSites) . " sites\n";
    return $self->$orig( \@keepSites );
};

around keepHapCols => sub {
    my $orig = shift;
    my $self = shift;

    return $self->$orig() if defined $self->$orig();
    my @keepSamples;
    unless ( $self->samp =~ m/./ ) {
        my $cmd;
        if ( $self->hap =~ m/\.gz$/ ) {
            $cmd = "gzip -dc " . $self->hap . q/ | head -1 | awk '{print NF}'/;
        }
        else {
            $cmd = " head -1 " . $self->hap . q/ | awk '{print NF}'/;
        }
        my $numLines = qx/$cmd/;
        chomp $numLines;
        croak "weird awk nf. hap empty?" unless $numLines > 0;
        @keepSamples = map { 1 } 0 .. ( 2 * $numLines - 1 );
    }
    else {
        my $openCmd = $self->openCmd( $self->samp, '<' );
        open( my $inFH, $openCmd );
        my $head = <$inFH>;
        chomp $head;
        croak "Unexpected samp header '$head' in samp file: " . $self->samp
          unless $head =~ m/^sample population group sex/;

        my $rowNum = 0;
        while (<$inFH>) {
            chomp;
            my @line = split(' ');
            my $keep = 0;
            $keep = 1 if $self->keepSampLine( \@line );
            push @keepSamples, ( $keep, $keep );
        }
    }
    print "Keeping " . sum(@keepSamples) . " haplotypes\n";
    return $self->$orig( \@keepSamples );
};

around _requiredFiles => sub {
    my $orig = shift;
    my $self = shift;

    return $self->$orig() if defined $self->$orig();

    my %files = %{ $self->_inputFiles };
    my %counts;
    map { $counts{$_} = 0 } keys %files;
    my $rhFilters = $self->_filterFiles;
    for my $filter ( sort keys %{$rhFilters} ) {
        next unless $self->$filter =~ m/./;
        map { $counts{$_}++ } @{ $rhFilters->{$filter} };
    }
    croak "error in building filterFiles" if keys %counts != 3;

    my @status =
      ( "Reading files:", grep { $counts{$_} > 0 } sort keys %counts );
    print join( ' ', @status ) . "\n";

    return $self->$orig( \%counts );
};

# no filters implemented yet
sub keepLegLine {
    my $self   = shift;
    my $raLine = shift;

    my $position = $raLine->[1];
    return 1 if exists $self->keepSitesMap->{$position};

    return 0;
}

# filtering sample line according to filters
sub keepSampLine {
    my $self   = shift;
    my $raLine = shift;
    my $keep   = 1;
    if ( $self->keepPop =~ m/./ ) {
        $keep = 0 unless $self->keepPop eq $raLine->[1];
    }
    if ( $self->keepGroup =~ m/./ ) {
        $keep = 0 unless $self->keepGroup eq $raLine->[2];
    }
    if ( $self->keepSex =~ m/./ ) {
        $keep = 0 unless $self->keepPop eq $raLine->[3];
    }
    return $keep;
}

sub PrintFilterHap {
    my $self = shift;

    # don't print anything if the hap file wasn't required!
    return unless $self->_requiredFiles->{hap};

    print "Printing hap output file to: " . $self->outHap . "\n";
    my @keepCols = @{ $self->keepHapCols };
    my @keepRows = @{ $self->keepHapRows };

    my $openCmd = $self->openCmd( $self->hap, '<' );
    open( my $hapFH, $openCmd );
    $openCmd = $self->openCmd( $self->outHap, '>' );
    open( my $outHapFH, $openCmd );
    my $numRowsKept = 0;
  LINE: for my $rowNum ( 0 .. $#keepRows ) {
        print STDERR "$numRowsKept/$rowNum\r" if $rowNum % 1000 == 0;
        my $line = <$hapFH>;
        croak "hap file has less rows than legend says it should have!"
          unless defined $line;
        next LINE unless $keepRows[$rowNum];

        $numRowsKept++;

        chomp $line;
        my @line = split( ' ', $line );
        croak "number of cols in hap ("
          . @line
          . ") and cols according to sample file ("
          . @keepCols
          . ") do not agree!"
          unless @keepCols == @line*2;
        my @outLine;
      COL: for my $colNum ( 0 .. $#line ) {
            next COL unless $keepCols[$colNum];
            push @outLine, $line[$colNum];
        }
        print $outHapFH join( ' ', @outLine ) . "\n";
    }
    print STDERR "$numRowsKept/" . @keepRows . "\n";

    croak "hap file longer than legend!" if <$hapFH>;
    close($hapFH);
    close($outHapFH);
}

sub PrintFilterLeg {
    my $self = shift;

    # don't print anything if the hap file wasn't required!
    return unless $self->_requiredFiles->{leg};

    print "Printing leg output file to: " . $self->outLeg . "\n";
    my @keepRows = @{ $self->keepHapRows };
    my $openCmd = $self->openCmd( $self->leg, '<' );
    open( my $inFH, $openCmd );
    $openCmd = $self->openCmd( $self->outleg, '>' );
    open( my $outFH, $openCmd );

    #print header
    my $head = <$inFH>;
    print $outFH $head;

  LINE: for my $rowNum ( 0 .. $#keepRows ) {
        my $line = <$inFH>;
        croak "legend file has less rows than legend says it should have!"
          unless defined $line;
        next LINE unless $keepRows[$rowNum];
        print $outFH $line;
    }
    confess "legend file is longer than expected!" if <$inFH>;
    close($inFH);
    close($outFH);
}

sub PrintFilterSamp {
    my $self = shift;

    # don't print anything if the samp file wasn't required!
    return unless $self->_requiredFiles->{samp};

    print "Printing samp output file to: " . $self->outSamp . "\n";
    my @keepCols = @{ $self->keepHapCols };
    my $openCmd = $self->openCmd( $self->samp, '<' );
    open( my $inFH, $openCmd );
    $openCmd = $self->openCmd( $self->outSamp, '>' );
    open( my $outFH, $openCmd );

    my $head = <$inFH>;
    print $outFH $head;

  LINE: for my $rowNum ( 0 .. $#keepCols ) {
        next unless $rowNum % 2 == 0;    # only use even rows
        my $line = <$inFH>;
        croak "samp file only has $rowNum rows when it should have "
          . @keepCols . "!\n"
          unless defined $line;
        next LINE unless $keepCols[$rowNum];
        print $outFH $line;
    }
    croak "samp file longer than samp!" if <$inFH>;
    close($inFH);
    close($outFH);
}

sub openCmd {
    my $self   = shift;
    my $file   = shift;
    my $rwFlag = shift;

    my $read  = $rwFlag eq '<' ? 1 : 0;
    my $write = $rwFlag eq '>' ? 1 : 0;
    croak "rwFlag is not '<' or '>': $rwFlag" unless $read xor $write;

    my $cmd;
    if ( $file =~ m/\.gz$/ ) {
        if ($read) {
            $cmd = "gzip -dc $file |";
        }
        else {
            $cmd = "| gzip -c > $file";
        }
    }
    else {
        $cmd = "$rwFlag $file";
    }
    return $cmd;
}

sub run {
    my $self = shift;

    $self->PrintFilterHap;
    $self->PrintFilterLeg;
    $self->PrintFilterSamp;
}

no Moose;

use warnings;
use strict;
use ParseArgs qw/getCommandArguments/;
my %args = getCommandArguments(
    requiredArgs => {
        hap       => q//,
        leg       => q//,
        samp      => q//,
        keepPop   => q//,
        keepGroup => q//,
        keepSex   => q//,
        keepSites => q//,
        base      => q/out/,
    }
);

my $runner = Runall->new(%args);

$runner->run();

__END__

=head1 NAME

hapLegSampFilter.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for hapLegSampFilter.pl, 
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
