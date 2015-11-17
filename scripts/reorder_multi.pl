#!/usr/bin/perl
# reorder_multi.pl     wkretzsch@gmail.com     2015/10/25 08:37:04

use warnings;
use strictures;
$| = 1;

use autodie;
use Carp;

use Getopt::Std;
my %args;
getopts( "r:", \%args );

my %sites;

croak "Need to specify reference site list with -r" unless $args{r};
croak "Reference site list does not exist.\nCould not find $args{r}"
  unless -f $args{r};

my $fh;
if ( $args{r} =~ m/\.gz$/ ) {
    open( $fh, "gzip -dc $args{r} |" );
}
else {
    open( $fh, "<", $args{r} );
}

# storing allele order in hash
while (<$fh>) {
    next if m/^#/;
    chomp;
    my @line = get_sites($_);
    my $key  = make_key( \@line );
    if ( exists $sites{$key} ) {
        push @{ $sites{$key} }, $line[4];
    }
    else {
        $sites{$key} = [ $line[4] ];
    }
}
close $fh;

# reorder the VCF
my %buffer = ( key => "", buf => [] );   # place to store lines to be rearranged
LINE: while (<>) {

    # skip comment lines
    if (m/^#/) {
        print $_;
        next LINE;
    }

    # process the line
    chomp;
    my @line = get_sites($_);
    my $key  = make_key( \@line );

    croak "Could not find site in reference site list: $key:$line[4]"
      unless exists $sites{$key};
    my $ref = $sites{$key};
    croak
"Could not find site alternate allele in reference site list: $key:$line[4]"
      unless defined get_alt_idx( $ref, $line[4] );

    # check if buffer needs to be emptied
    if ( @{ $buffer{buf} } ) {

        # new site encountered, empty buffer before pushing new line on
        if ( $buffer{key} ne $key ) {
            print_lines_in_order( $buffer{buf}, $sites{ $buffer{key} } );
            $buffer{buf} = [];
            $buffer{key} = "";
        }
    }

    # push line into buffer
    $buffer{key} = $key;
    push @{ $buffer{buf} }, \@line;

}

# and last print
print_lines_in_order( $buffer{buf}, $sites{ $buffer{key} } );

exit 0;

sub print_lines_in_order {
    my $lines     = shift;
    my $alt_order = shift;

    map { print join( "\t", @{$_} ) . "\n" }
      reorder_lines( $lines, $alt_order );
}

sub reorder_lines {
    my $lines     = shift;
    my $alt_order = shift;

    return @{$lines} if (@{$lines} < 2);
    my @lines = sort {
        get_alt_idx( $alt_order, $a->[4] )
          <=> get_alt_idx( $alt_order, $b->[4] )
      } @{$lines};
    return @lines
}

sub get_alt_idx {
    my $arr = shift;
    my $alt = shift;

    my @idxs = grep { $alt eq $arr->[$_] } 0 .. $#{$arr};
    croak "Found too many indices" unless @idxs < 2;
    return @idxs ? $idxs[0] : undef;
}

sub get_sites {
    my $line = shift;
    return ( split "\t", $_, 6 );
}

sub make_key {
    my $line = shift;
    return join( ':', @{$line}[ 0, 1, 3 ] );
}

__END__

=head1 NAME

reorder_multi.pl

=head1 SYNOPSIS

reorder_multi.pl -r sites.vcf.gz < out_of_order.vcf > in_order.vcf

=head1 DESCRIPTION

This scripts loads up all sites in sites.vcf.gz.  It then scans through the VCF input through STDIN and rearranges lines such that the alternate allele ordering in STDIN matches sites.vcf.gz.  The result is printed to STDOUT.  

sites.vcf.gz has to be a superset of the sites in STDIN.

Multiallelics need to have been broken as bi-allelics using for example `bcftools norm -m-`.

=head1 AUTHOR

Warren Winfried Kretzschmar, E<lt>winni@winnis-MacBook-Pro2.localE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Warren Winfried Kretzschmar

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
Couldn't open file:perltidy.ERR in mode:w : Permission denied
couldn't open perltidy.ERR Permission denied
