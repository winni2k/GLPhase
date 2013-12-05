#!/usr/bin/perl
# hapLeg2haps.pl                   wkretzsch@gmail.com
#                                  05 Dec 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "hapLeg2haps.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

use warnings;
use strict;
$|=1;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;

use Getopt::Std;
my %args;
getopts( 'h:l:c:', \%args );

my $chrom = $args{c} || croak "need to specify chromosome";

open(my $hapFH, openCmd($args{h}));
open(my $legFH, openCmd($args{l}));

# discard legend header
<$legFH>;

while(1){
    my $raLeg = getTokens($legFH);
    my $raHap = getTokens($hapFH);
    if(defined $raLeg xor defined $raHap){
        croak "legend and hap files are different lengths";
    }
    last if !defined $raLeg && !defined $raHap;

    #  looks like we have two good lines from here down
    print join(q/ /, ($chrom, @{$raLeg}, @{$raHap}))."\n";
}
      

sub getTokens{
    my $fh = shift;
    my $line = <$fh> || return undef;
    chomp $line;
    my @line = split(' ', $line);
    return \@line;
}

sub openCmd{

    my $file = shift;
    my $cmd;
    croak "need to define input file" unless defined $file;
    if($file =~ m/\.gz$/){
        return "gzip -dc $file |";
    }
    else{
        return "<$file";
    }
}


__END__

=head1 NAME

hapLeg2haps.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for hapLeg2haps.pl, 
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
