#!/usr/bin/perl
# haps2hapLegend.pl                   wkretzsch@gmail.com
#                                     12 Oct 2013

use warnings;
use strict;
$|=1;
use Data::Dumper;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use Carp;
use autodie;

use Getopt::Std;
my %args;
getopts( 'i:o:', \%args );
my $DEBUG = $args{d} || 1;

croak "need to specify haps input file" unless $args{i};
croak "input (-i) does not exist" unless -e $args{i};

my $openCmd;
if($args{i} =~ m/\.gz$/){
    $openCmd = "gzip -dc $args{i} |";
}
else{
    $openCmd = " < $args{i}";
}

open (my $fhIN, $openCmd);
open (my $fhHAP, " | gzip -c > $args{o}.hap.gz");
open (my $fhLEGEND, '>', $args{o}.".legend");

# print header
print $fhLEGEND join(q/ /,qw/id position a0 a1/)."\n";

while(<$fhIN>){
    chomp;
    my @line =split(/ /);

    # print legend line
    my $chrom = shift @line;
    my $pos = $line[0] eq '.' ? $line[1] : $line[0];
    my (undef,undef, @alleles) = splice(@line, 0,4);
    print $fhLEGEND join(q/ /, ($chrom,$pos,@alleles))."\n";

    # print hap line
    print $fhHAP join(q/ /, @line)."\n";
}


__END__

=head1 NAME

haps2hapLegend.pl

=head1 SYNOPSIS
   
    # splits a haps.gz file in to out.hap.gz and out.legend files
    haps2hapLegend.pl -i haps.gz -o out

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
