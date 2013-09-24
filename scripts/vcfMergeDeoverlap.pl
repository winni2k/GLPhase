#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use autodie;

=head1 Name

beagle genotype probability file merger

=head1 Synopsis

# -b    take a list of chromosome regions and gprobs files and merge them to one
# formatting is in BED format (three columns: chromosome, start, end)
# fourth column is file location (gzipped or not is fine) by default
# -c    (3) column in BED file with file location (0 based)
# -o output file (optional, STDOUT otherwise)

gprobs_merger.pl -b nonoverlapping.bed -o outputFile

=cut

my %args;
getopts( "b:o:", \%args );

die 'need to define -o option (output file)' unless defined $args{o};
open( my $fhOUT, " | bgzip -c > $args{o}" ) or die "Can't open $args{o}; $!";

die 'need to define -b option (BED file)' unless defined $args{b};
open( my $fh_bed, '<', "$args{b}" ) or die "Can't open $args{b}; $!";

my $fileCol = $args{c} || 3;

# parse bed file
my @general_header;
my @general_otherHeaders;
while (<$fh_bed>) {
    chomp;
    my @bed_line = split( /\t/, $_ );
    my $vcf_file = $bed_line[$fileCol];
    my $vcf_open = "dd ibs=1G if=$vcf_file | ";
    if ( $vcf_file =~ m/\.gz$/ ) { $vcf_open .= 'gzip -dc | ' }
    open( my $fh_in, "$vcf_open" )
      or die "Can't open $vcf_file; $!";

    # read in header
    my @header;
    my @otherHeaders;

    # read in vcf file and keep only lines specified by bed position
    my $line = 0;    # line counter
    while (<$fh_in>) {
        $line++;
        my $printLine = 0;    # true if we wish to print this line
        chomp;
        my $in_line = $_;
        my @in_line = split( /\t/, $_ );

        unless (@header) {
            croak "no header in vcf file $vcf_file" unless $in_line =~ m/^#/;
            if ( $in_line =~ m/^##/ ) {
                push @otherHeaders, $in_line;
            }
            elsif ( $in_line =~ m/^#CHROM/ ) {
                @header = @in_line;
                unless (@general_header) {
                    @general_header       = @header;
                    @general_otherHeaders = @otherHeaders;
                }
                for my $index ( 0 .. $#general_header ) {
                    croak
"$vcf_file: headers don't match at column $index (0-based)."
                      if $header[$index] != $general_header[$index];
                }
                croak "$vcf_file: headers don't match: "
                  . @otherHeaders . " "
                  . @general_otherHeaders
                  if @otherHeaders != @general_otherHeaders;

                last;
            }
            else {
                croak "Unexpected header line:\n$in_line";
            }
            $printLine = 1 unless @general_header;
        }

        # otherwise decide if we want to keep this line
        else {
            # parse first column
            my $position = $in_line[1]
              || croak "position not available Unexpectedly";

            # print line if it is in the bed interval
            if ( $position > $bed_line[0] && $position <= $bed_line[1] ) {
                $printLine = 1;
            }
        }

        # print the line if we set printLine
        if ($printLine) {
            print $fhOUT $in_line . "\n";
        }
    }
    close($fh_in) or die "could not close $vcf_open; $?, $!";
}
close($fh_bed);
