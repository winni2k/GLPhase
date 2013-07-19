package bamTrackLib;
our $VERSION = '0.001';

# wkretzsch@gmail.com          19 Jul 2013

use 5.006;
use strict;
use warnings;
use Carp;
use Data::Dumper;
use DBI;
use File::Basename;
use Digest::MD5;

use Moose;

has operations => ( is => 'ro', isa => 'ArrayRef[Str]' );
has inputFile  => ( is => 'ro', isa => 'Str', );
has dbHandle   => ( is => 'rw', isa => 'DBI::db' );

has db   => ( is => 'ro', isa => 'Str' );
has host => ( is => 'ro', isa => 'Str' );
has user => ( is => 'ro', isa => 'Str' );

has inputBams => ( is => 'rw', isa => 'HashRef[Str]' );

around inputBams => sub {
    my $orig = shift;
    my $self = shift;

    return $self->$orig()
      unless @_;

    my $input = shift;

    # parse bam list
    my @inputBams;
    unless ( $input eq 'ARRAYREF' ) {
        @inputBams = @{ $self->_parseBamListFile($input) };
    }
    else {
        @inputBams = @{$input};
    }
    my %inputBams = map { $_ => 1 } @inputBams;

    return $self->$orig( \%inputBams );
};

has inputBamMD5sums => ( is => 'rw', isa => 'HashRef[Str]' );

around inputBamMD5sums => sub {
    my $orig = shift;
    my $self = shift;

    unless (@_) {

        # normal behaviour if data has been initialized
        my $retval = $self->$orig();
        return $retval if defined $retval;

        # pull md5s from bams if not defined yet
        my %inputBams = %{ $self->inputBams };

        ### look for md5sums based on bamList names
        # find md5sums not existent
        # run md5sums on any files that still need it
        my %inputBamMD5sums = %inputBams;
        print STDERR "creating missing md5sums...\n";

        my $ctx = Digest::MD5->new;
        for my $bam ( sort keys %inputBams ) {
            my $md5 = $bam . '.md5';
            unless ( -e $md5 ) {
                open( my $fh, '<', $bam ) or croak "Could not open $bam";

                $ctx->addfile($fh);
                my $digest = $ctx->hexdigest;
                close($fh);

                # print out md5sum
                my ( $bamName,$path ) = fileparse($bam);
                print $bamName."\n";
                open( $fh, '>', $md5 )
                  or croak "could not open $md5 for writing";
                print $fh $digest . "  $bamName\n";
                close($fh);

              #                my (undef, $md5Name)=fileparse($md5);
              #                system("cd $path && md5sum $bamName > $md5Name");
            }
            $inputBamMD5sums{$bam} = $md5;
        }

        $self->$orig( \%inputBamMD5sums );
        return $self->$orig();
    }
    return $self->$orig(@_);

};

sub BUILD {
    croak
"Please export your password to DBI_PASS environment var before running this script"
      unless defined $ENV{DBI_PASS};
}

sub connect {
    my $self = shift;
    unless ( defined $self->dbHandle ) {
        my $DBI =
          DBI->connect( 'dbi:mysql:' . $self->db() . ':' . $self->host(),
            $self->user() )
          || croak "Error opening database: $DBI::errstr";
        $self->dbHandle($DBI);
    }
    return 1;
}

# takes a list of bams (file or arrayref) and saves them into the DB
sub registerBams {
    my $self  = shift;
    my $input = shift;

    $self->inputBams($input);
    $self->inputBamMD5sums();

    #
    return 1;
}

sub _parseBamListFile {
    my $self = shift;
    my $file = shift;
    open( my $fh, '<', $file ) or croak "Could not open $file";
    my @bams;
    my $line;
    while (<$fh>) {
        chomp;
        $line++;
        croak
          "File on line $line of bamlist $file does not look like a bam file"
          unless m/\.bam$/;
        push @bams, $_;
    }
    return \@bams;
}

sub listDrivers {
    my $self = shift;

    print "Available DBI drivers:\n";
    my @drivers = DBI->available_drivers('quiet');
    my @sources;

    foreach my $driver (@drivers) {
        print "$driver\n";
        @sources = eval { DBI->data_sources($driver) };
        if ($@) {
            print "\tError: $@\n";
        }
        elsif (@sources) {
            foreach my $source (@sources) {
                print "\t$source\n";
            }
        }
        else {
            print "\tno known data sources\n";
        }
    }
}

no Moose;
1;
__END__

=head1 NAME

bamTrackLib - Perl extension for blah blah blah

=head1 SYNOPSIS

   use bamTrackLib;
   blah blah blah

=head1 DESCRIPTION

This class was designed for use with bamTrack.pl

=head2 EXPORT

None by default.

=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Warren Winfried Kretzschmar

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut

