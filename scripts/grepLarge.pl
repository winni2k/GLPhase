#!/usr/bin/perl
# grepLarge.pl                   wkretzsch@gmail.com
#                                09 Oct 2013

use warnings;
use strict;
$|=1;
use Data::Dumper;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Getopt::Std;

my %args;
getopts( 'f:vc:d:', \%args );

# 0 = match whole line
# col is 1 based
my $splitCol = $args{c} || 0; 
my $del = $args{d} || '\s';

open (my $fh, '<', $args{f});
my %strings;
while(<$fh>){
    chomp;
    $strings{$_} = 1;
}

while(<>){
    chomp;
    my $print = 0;
    my @line;
    my $matchString = $_;
    if($splitCol){
        @line = split($del, $_);
        $matchString = $line[$splitCol-1];
    }

    if($strings{$matchString}){
        $print = 1;
    }
    $print = !$print if $args{v};
    print $_ ."\n" if $print;
}


__END__

=head1 NAME

grepLarge.pl

=head1 SYNOPSIS
   
# only looks for exact matches
grepLarge.pl -f searchStringList < searchfile >matches

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
