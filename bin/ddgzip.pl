#!/usr/bin/perl -w
our $VERSION = "1.001";
use strict;

=head1 Name

ddgzip

=head1 Author

Warren Kretzschmar

=head1 Description

Run gzip using dd for in- and output buffering, but preserve creation time.
This script will delete the original file.

=head1 Synopsis

ddgzip.pl *  # gzips all files found by the glob

=head1 caveats

ddgzip.pl takes no options to gzip, and only compresses at level -1. 
Buffer of dd is set to 1GB.
This script is meant to run on an ubuntu system.

=cut

foreach my $file (@ARGV) {
    if ( $file =~ /^--version$/ ) { print "$VERSION\n"; exit; }
    system "dd ibs=1G if=$file | gzip -c -1 | "
      . "dd obs=1G of=$file.tmp && touch -r $file $file.tmp && "
      . "mv $file.tmp $file.gz && rm $file";
}
