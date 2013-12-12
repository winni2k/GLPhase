#!/usr/bin/perl
# hapLegSampFilter.pl                   wkretzsch@gmail.com
#                                       12 Dec 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "hapLegSampFilter.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

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
use DM;
$| = 1;

#input variables
has 'noDryRun' => ( is => 'ro', init_arg => 'r', isa => 'Bool', default => 0 );
has 'numJobs'  => ( is => 'ro', init_arg => 'j', isa => 'Int',  default => 1 );
has 'dm'        => ( is => 'rw', isa     => 'DM' );
has 'hap'       => ( is => 'rw', default => q// );
has 'leg'       => ( is => 'rw', default => q// );
has 'samp'      => ( is => 'rw', default => q// );
has 'keepPop'   => ( is => 'rw', default => q// );
has 'keepGroup' => ( is => 'rw', default => q// );
has 'keepSex'   => ( is => 'rw', default => q// );
has base        => ( is => 'rw', default => 'out' );

# constants
has _filterFiles =>
  ( is => 'ro', isa => 'HashRef[ArrayRef]', builder => '_build_filterFiles' );
has _inputFiles => (
    is      => 'ro',
    isa     => 'HashRef[Bool]',
    default => sub { { hap => 1, samp => 1, leg => 1 } }
);

sub _build_filterFiles {
    my %return;
    for (qw/keepPop keepGroup keepSex/) {
        $return{$_} = [qw/hap samp/];
    }
    return \%return;
}

sub outHap {
    my $self = shift;
    return $self->base . ".hap.gz"
}
sub outSamp {
    my $self = shift;
    return $self->base . ".samp.gz"
}
sub outLeg {
    my $self = shift;
    return $self->base . ".leg.gz"
}



# internal variables
has '_requiredFiles' => (
    is  => 'rw',
    isa => 'HashRef',
);
has keepHapCols => ( is => 'rw', isa => 'ArrayRef[Bool]' );
has keepHapRows => ( is => 'rw', isa => 'ArrayRef[Bool]' );

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
            my $keep = 0;
            $keep = 1 if $self->keepLegLine( \@line );
            push @keepSites, $keep;
        }
    }
    print "Keeping " . @keepSites . " sites\n";
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
        @keepSamples = map { 1 } 0 .. ( 2*$numLines - 1 );
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
            push @keepSamples, ($keep, $keep);
        }
    }
    print "Keeping " . @keepSamples . " haplotypes\n";
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

sub BUILD {

    my $self = shift;
    $self->dm(
        DM->new( dryRun => !$self->noDryRun, numJobs => $self->numJobs ) );

}

# no filters implemented yet
sub keepLegLine {
    my $self   = shift;
    my $raLine = shift;

    return 1;
}

# filtering sample line according to filters
sub keepSampLine {
    my $self   = shift;
    my $raLine = shift;
    my $keep   = 1;
    if ( defined $self->keepPop ) {
        $keep = 0 unless $self->keepPop eq $raLine->[1];
    }
    if ( defined $self->keepGroup ) {
        $keep = 0 unless $self->keepGroup eq $raLine->[2];
    }
    if ( defined $self->keepSex ) {
        $keep = 0 unless $self->keepPop eq $raLine->[3];
    }
    return $keep;
}

sub PrintFilterHap {
    my $self = shift;

    # don't print anything if the hap file wasn't required!
    return unless $self->_requiredFiles->{hap};

    print "Printing hap output file to: ".$self->outHap;
    my @keepCols = @{ $self->keepHapCols };
    my @keepRows = @{ $self->keepHapRows };
    my $openCmd  = $self->openCmd( $self->hap, '<' );
    open( my $hapFH, $openCmd );
    $openCmd = $self->openCmd( $self->outHap, '>' );
    open( my $outHapFH, $openCmd );
  LINE: for my $rowNum ( 0 .. $#keepRows ) {
        my $line = <$hapFH>;
        croak "hap file has less rows than legend says it should have!"
          unless defined $line;
        next LINE unless $keepRows[$rowNum];

        chomp $line;
        my @line = split( ' ', $line );
        croak "number of cols in hap ("
          . @line
          . ") and cols according to sample file ("
          . @keepCols
          . ") do not agree!"
          unless @keepCols == @line;
        my @outLine;
      COL: for my $colNum ( 0 .. $#line ) {
            next COL unless $keepCols[$colNum];
            push @outLine, $line[$colNum];
        }
        print $outHapFH join( ' ', @outLine ) . "\n";
    }
    croak "hap file longer than legend!" if <$hapFH>;
    close($hapFH);
    close($outHapFH);
}

sub PrintFilterLeg {
    my $self = shift;

    # don't print anything if the hap file wasn't required!
    return unless $self->_requiredFiles->{leg};
    print "Printing leg output file to: ".$self->outLeg;

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
    print "Printing samp output file to: ".$self->outSamp;

    my @keepCols = @{ $self->keepHapCols };
    my $openCmd = $self->openCmd( $self->samp, '<' );
    open( my $inFH, $openCmd );
    $openCmd = $self->openCmd( $self->outSamp, '>' );
    open( my $outFH, $openCmd );

    my $head = <$inFH>;
    print $head;

  LINE: for my $rowNum ( 0 .. $#keepCols ) {
        my $line = <$inFH>;
        croak "samp file has less rows than samp says it should have!"
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
        r         => 0,
        j         => 1,
        hap       => q//,
        leg       => q//,
        samp      => q//,
        keepPop   => q//,
        keepGroup => q//,
        keepSex   => q//,
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
