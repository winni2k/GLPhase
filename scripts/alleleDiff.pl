#!/usr/bin/perl
# alleleDiff.pl                   wkretzsch@gmail.com
#                                 19 Oct 2013
our $VERSION = '0.001';
$VERSION = eval $VERSION;
print STDERR "alleleDiff.pl -- $VERSION\nBy\twkretzsch@gmail.com\n\n";

use warnings;
use strict;
$| = 1;

use File::Path qw(make_path);
use File::Basename;
use Env qw(HOME);
use autodie;
use Carp;

my @files = @ARGV;

my %diffAllSites;
my @siteLists;
print STDERR "Loading files\n";
for my $file (@files) {
    push @siteLists, load_expSiteList($file);
}

print STDERR "Comparing files\n";
my %diffSites;
for my $sLNum ( 0 .. $#siteLists ) {
    my @diffSites = findDiffAllSites( $sLNum, @siteLists );
    map { $diffSites{$_}++ } @diffSites;
}

print STDERR "Found ".(keys %diffSites). " differing sites\n";
print STDERR "Printing differing sites\n";
map {print $_ ."\n"} sort keys %diffSites;

# find differing alleles between the file stated by $num
# and any of the other files
sub findDiffAllSites {

    my $num = shift;
    my @SLs = @_;
    croak "too few files" unless @SLs > 1;
    my $compSL = $SLs[$num];
    splice( @SLs, $num, 1 );
    my @diffAllSiteKeys;
    for my $secSL (@SLs) {
        for my $siteKey ( sort keys %{$compSL} ) {
            if ( exists $secSL->{$siteKey} ) {
                my $match = 1;
                map {
                    $match = 0
                      unless $compSL->{$siteKey}->[$_] eq
                      $secSL->{$siteKey}->[$_]
                } 0 .. 1;
                push @diffAllSiteKeys, $siteKey unless $match;
            }
        }
    }
    return (@diffAllSiteKeys);
}

sub load_expSiteList {

    my $file = shift;
    my $openCmd = $file =~ m/\.gz$/ ? " gzip -dc $file | " : " < $file ";
    open( my $fh, $openCmd ) or croak "could not open file $file; $!";
    my %siteList;
    my %blackList;
    while (<$fh>) {
        chomp;
        my @line = split(/ /);

        # throw out multiallelic sites
        my $siteKey = "$line[0]:$line[1]";
        $blackList{$siteKey}++ if exists $siteList{$siteKey};
        $siteList{"$line[0]:$line[1]"} = [ @line[ 2 .. 3 ] ];
    }

    # remove any multiallelic snps
    for my $siteKey ( sort keys %blackList ) {
        delete $siteList{$siteKey};
    }
    close($fh);
    return \%siteList;
}

__END__

=head1 NAME

alleleDiff.pl

=head1 SYNOPSIS
   

=head1 DESCRIPTION
Stub documentation for alleleDiff.pl, 
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
