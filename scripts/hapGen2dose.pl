#!/usr/bin/perl
# hap2dose.pl                   wkretzsch@gmail.com
#                               08 Oct 2013

use warnings;
use strict;
$| = 1;

use Env qw(HOME);
use Carp;
use threads;
use Thread::Queue 3.02;
use 5.008;
use Getopt::Std;
my %args;
getopts( 'hgs:t:b:', \%args );

my $numThreads = $args{t} || 1;
my $hapFlag    = $args{h} || 0;
my $genFlag    = $args{g} || 0;
croak "Please only specify hap or gen flag" unless $hapFlag xor $genFlag;
my $sampList = $args{s} || croak "Please specify sample list with -s";
croak "Please make sure sampList $sampList exists" unless -e $sampList;
my $cacheSize = $args{b} || 200;

# print header
open( my $fh, '<', $sampList) or confess "could not open file $sampList";
my @out = qw/marker alleleA alleleB/;
while (<$fh>) {
    chomp;
    push @out, $_;
}
print STDOUT join( ' ', @out ) . "\n";
my $numCols = @out - 3;

my $qLines = Thread::Queue->new();
my $qOut   = Thread::Queue->new();
my @threads;
for ( 1 .. $numThreads ) {
    push @threads,
      threads->create(
        'processLinesWorker',
        isHap   => $hapFlag,
        isGen   => $genFlag,
        qLines  => $qLines,
        qOut    => $qOut,
        numCols => $numCols,
      );
}

# end qOut to tell empty outline worker when to return
my $qInLineNum  = Thread::Queue->new();
my $qOutLineNum = Thread::Queue->new();
my $tEmptyOutLineWorker =
  threads->create( 'emptyOutInOrderWorker', $qOut, $qInLineNum, $qOutLineNum );

my $inLineNum  = 0;
my $outLineNum = 0;
while ( my $new_line = <> ) {
    $inLineNum++;

    # Only read next line in to queue if there are less lines in cache
    # than the cache size parameter
    # Otherwise sleep for a second and give an update of how far we are
    while ( $inLineNum % $cacheSize == 0 && $qLines->pending() > $cacheSize ) {
        sleep 1;

        # check status of jobs if we are waiting anyway
        #defined( my $nextOutLineNum = $qOutLineNum->dequeue() )
        if ( my $numPending = $qOutLineNum->pending() ) {
            my @nextOutLineNums = $qOutLineNum->dequeue($numPending);
            $outLineNum = pop @nextOutLineNums;
            print STDERR "Processing line number $outLineNum/$inLineNum\n";
        }
    }
    $qLines->enqueue( $inLineNum, $new_line );
    $qInLineNum->enqueue($inLineNum);
}

# finish input and print queues
$qLines->end();
$qInLineNum->end();

# continue to check on status
while ( $outLineNum < $inLineNum ) {

    sleep 1;
    if ( my $numPending = $qOutLineNum->pending() ) {
        my @nextOutLineNums = $qOutLineNum->dequeue($numPending);
        $outLineNum = pop @nextOutLineNums;
        print STDERR "Processing line number $outLineNum/$inLineNum\r";
    }
}

# wait for all threads to finish
print STDERR "Waiting on worker queues to finish...";
map { $_->join() } @threads;
$tEmptyOutLineWorker->join();
print STDERR "done.\n";

exit 0;

sub processLinesWorker {

    my %args    = @_;
    my $isHap   = $args{isHap};
    my $isGen   = $args{isGen};
    my $qLines  = $args{qLines} || confess;
    my $qOut    = $args{qOut} || confess;
    my $numCols = $args{numCols} || confess;
    confess "isHap and isGen are true" unless $isHap xor $isGen;

    my ( $lineNum, $new_line ) = $qLines->dequeue(2);
    while ( defined $lineNum ) {
        chomp($new_line);
        my @line = split( ' ', $new_line );
        my @out = ("$line[0]:$line[2] $line[3] $line[4]");
        if ($isHap) {
            for my $colNum ( 5 .. $#line ) {
                if ( $colNum % 2 == 0 ) {
                    push @out, $line[$colNum] + $line[ $colNum - 1 ];
                }
            }
        }
        elsif ($isGen) {
            for my $colNum ( 5 .. $#line ) {
                if ( $colNum % 3 == 0 ) {
                    push @out, $line[$colNum] + 2 * $line[ $colNum + 1 ];
                }
            }
        }
        else {
            confess "neither hap nor gen are true";
        }

        my $lineNumCols = @out - 1;
        croak
"Input at line $lineNum has data for $lineNumCols instead of $numCols samples"
          if $lineNumCols != $numCols;

        $qOut->enqueue( $lineNum, join( ' ', @out ) . "\n" );

        ( $lineNum, $new_line ) = $qLines->dequeue(2);
    }
}

sub emptyOutInOrderWorker {
    my ( $qOut, $qInLineNum, $qOutLineNum ) = @_;

    my %cache;
    my $outLineNum = 0;
    my $inLineNum  = 1;
    while ( defined( my $nextInLineNum = $qInLineNum->dequeue() ) ) {
        $inLineNum = $nextInLineNum;
        $outLineNum = emptyOutInOrder( $qOut, $outLineNum, \%cache );
        $qOutLineNum->enqueue($outLineNum);
    }
    while ( $outLineNum < $inLineNum ) {
        $outLineNum = emptyOutInOrder( $qOut, $outLineNum, \%cache );
        $qOutLineNum->enqueue($outLineNum);
        sleep 1;
    }
    $qOutLineNum->end();
    croak "Only $outLineNum out of $inLineNum lines printed"
      if $outLineNum != $inLineNum;
}

sub emptyOutInOrder {
    my $qOut        = shift;
    my $outLineNum  = shift;
    my $rhLineCache = shift;

    while ( $qOut->pending() > 1 ) {
        my ( $lineNum, $out ) = $qOut->dequeue(2);

        # this is the next line to print
        if ( $lineNum + 1 == $outLineNum ) {
            print STDOUT $out;
            $outLineNum++;
        }

        # otherwise put the line in cache
        else {
            $rhLineCache->{$lineNum} = $out;
        }
    }

    # check cache for line with the correct line number
  KEY: for my $key ( sort keys %{$rhLineCache} ) {
        my $searchLineNum = $outLineNum + 1;

        # if the correct line number exists,
        # then pull it out of cache, print it and go to the next key
        if ( exists $rhLineCache->{$searchLineNum} ) {
            print STDOUT $rhLineCache->{$searchLineNum};
            delete $rhLineCache->{$searchLineNum};
            $outLineNum++;
        }

        # if the current key was not the searched for one
        # then none of the subsequent keys will be either.
        else {
            last KEY;
        }
    }

    return $outLineNum;
}

__END__

=head1 NAME

hap2dose.pl

=head1 SYNOPSIS
   
hapGen2dose.pl -t 10 -h -s sample.list <in.hap >out

hapGen2dose.pl -t 10 -g -s sample.list <in.gen >out

=head1 DESCRIPTION



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
