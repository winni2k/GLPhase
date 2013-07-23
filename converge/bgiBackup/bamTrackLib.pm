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
use namespace::autoclean;
use List::MoreUtils qw/any/;

has operations => ( is => 'ro', isa => 'ArrayRef[Str]' );
has inputFile  => ( is => 'ro', isa => 'Str', );
has dbHandle   => ( is => 'ro', isa => 'DBI::db', writer => '_dbHandle' );

has db   => ( is => 'ro', isa => 'Str', required => 1 );
has host => ( is => 'ro', isa => 'Str', required => 1 );
has user => ( is => 'ro', isa => 'Str', required => 1 );
has allowRelativePaths =>
  ( is => 'ro', isa => 'Bool', required => 1, default => 0 );

has _tables  => ( is => 'rw', isa => 'HashRef[Str]' );
has _headers => ( is => 'rw', isa => 'HashRef[Str]' );

has inputBams => ( is => 'ro', isa => 'HashRef[Str]', writer => '_inputBams' );

around _inputBams => sub {
    my $orig = shift;
    my $self = shift;
    my %args = @_;

    # parse bam list
    my @inputBams;
    croak "type must be defined" unless defined $args{type};
    if ( $args{type} eq 'fileList' ) {
        push @inputBams, @{ $self->_parseBamListFile( $args{input} ) };
    }
    if ( $args{type} eq 'file' ) {
        push @inputBams, $args{input};
    }
    croak "no input files given" unless @inputBams;
    croak "Bams with relative path are not allowed"
      if !$self->allowRelativePaths && grep { $_ !~ m|^/| } @inputBams;

    my %inputBams = map { $_ => 1 } @inputBams;

    # store input bams
    $self->$orig( \%inputBams );

    ### update md5sums based on those input bams
    ## look for md5sums based on bamList names
    # find md5sums not existent
    # run md5sums on any files that still need it
    my %inputBamMD5sums = %inputBams;
    print STDERR "creating missing md5sums...\n";

    my $ctx = Digest::MD5->new;
    for my $bam ( sort keys %inputBams ) {
        my $md5 = $bam . '.md5';

        # create new md5sum if md5 file does not exist
        # or it is younger than the bam modification time
        if ( !-e $md5 || ( stat($md5) )[9] < ( stat($bam) )[9] ) {
            open( my $fh, '<', $bam ) or croak "Could not open $bam";

            $ctx->addfile($fh);
            my $digest = $ctx->hexdigest;
            close($fh);

            # print out md5sum
            my ( $bamName, $path ) = fileparse($bam);
            print STDERR "\t" . $bamName . "\n";
            open( $fh, '>', $md5 )
              or croak "could not open $md5 for writing";
            print $fh $digest . "  $bamName\n";
            close($fh);

            #                my (undef, $md5Name)=fileparse($md5);
            #                system("cd $path && md5sum $bamName > $md5Name");
        }

        $inputBamMD5sums{$bam} = $md5;
    }

    $self->_inputBamMD5sums( \%inputBamMD5sums );

    return \%inputBams;
};

has inputBamMD5sums =>
  ( is => 'ro', isa => 'HashRef[Str]', writer => '_inputBamMD5sums' );

sub BUILD {
    my $self = shift;
    croak
"Please export your password to DBI_PASS environment var before running this script"
      unless defined $ENV{DBI_PASS};

    $self->connect;
    my $dbh = $self->dbHandle;

   #    my $sth = $dbh->prepare('SHOW VARIABLES LIKE "%version%"');
   #    $sth->execute();
   #    while ( my @row = $sth->fetchrow_array ) {
   #    print "@row\n";
   #  }
   #exit;
   #    print STDERR "Using mysql server version $dbh->{mysql_serverversion}\n";

    my %tables;

    # create required tables if they do not exist
    my $table        = 'md5sums';
    my $stringLength = 255;
    my $tableCreateCommand =
"CREATE TABLE IF NOT EXISTS $table ( md5sum VARCHAR(22) PRIMARY KEY, sampleName VARCHAR($stringLength), sampleNum INT, passedValidateSam BOOL, recalibrated BOOL);";
    $dbh->do($tableCreateCommand);
    $tables{$table} = 1;

    $table        = 'bamNames';
    $stringLength = 255;
    $tableCreateCommand =
"CREATE TABLE IF NOT EXISTS $table ( md5sum VARCHAR(22),bamName VARCHAR($stringLength), bamDir VARCHAR($stringLength), bamPath VARCHAR($stringLength), backup BOOL, backupDevice VARCHAR($stringLength));";
    $dbh->do($tableCreateCommand);
    $tables{$table} = 1;

    $self->_tables( \%tables );
    my %headers;
    $headers{$_} = 1
      for
      qw/md5sum sampleName sampleNum passedValidateSam recalibrated bamName bamDir bamPath backup backupDevice/;
    $self->_headers( \%headers );

    $dbh->commit;

}

sub DEMOLISH {
    my $self = shift;
    my $dbh  = $self->dbHandle();
    $dbh->disconnect;
}

sub connect {
    my $self = shift;
    unless ( defined $self->dbHandle ) {
        my $DBI =
          DBI->connect( 'dbi:mysql:' . $self->db() . ':' . $self->host(),
            $self->user(), undef, { RaiseError => 1, AutoCommit => 0 } )
          || croak "Error opening database: $DBI::errstr";
        $self->_dbHandle($DBI);
    }
    return 1;
}

sub retrieveBams {
    my $self = shift;

    my %args = @_;

    my %headers       = %{ $self->_headers };
    my %filterColumns = %{ $args{filterColumns} };
    my %searchHeaders = %filterColumns;

    # we are always searching for bamPaths
    $searchHeaders{bamPath} = 1;

    # save index in array of each header in searchHeaders
    my $colNum = 0;
    map { $searchHeaders{$_} = $colNum++ } sort keys %searchHeaders;

    #    print join( ' ', keys %searchHeaders ) . "\n";

    # make sure each search header exists
    for my $key ( keys %searchHeaders ) {
        croak "unknown header supplied to retrieveBams: '$key'"
          if !exists $headers{$key};
    }

    #
    my $dbh = $self->dbHandle;

    # just get all bams and parse in perl
    # lazy, I know but we'll sort it out later
    my $sth =
      $dbh->prepare( "SELECT "
          . join( ',', sort keys %searchHeaders )
          . " FROM md5sums a INNER JOIN bamNames b USING (md5sum);" );
    $sth->execute or croak $sth->errstr;

    my @returnBams;
    while ( my @row = $sth->fetchrow_array ) {

        # decide whether to keep row
        my $keep = 1;
        for my $header ( sort keys %filterColumns ) {
            my $index = $searchHeaders{$header};
            if ( ref( $filterColumns{$header} ) eq 'ARRAY' ) {
                ( $keep = 0 && last )
                  unless any { $row[$index] eq $_ }
                @{ $filterColumns{$header} };
            }
            else {
                ( $keep = 0 && last )
                  unless defined $row[$index]
                  && $row[$index] eq $filterColumns{$header};
            }
        }

        push @returnBams, $row[ $searchHeaders{bamPath} ] if $keep;
    }
    return @returnBams;

}

#sub validateInputBamMD5sums {
#    my $self = shift;

#    my $rhMD5sums = $self->inputBamMD5sums();

#}

# takes a list of bams (file or arrayref) and saves them into the DB
sub registerBams {
    my $self = shift;
    my %args = @_;

    if ( exists $args{fileList} ) {
        $self->_inputBams( type => 'fileList', input => $args{fileList} );
    }
    if ( exists $args{file} ) {
        $self->_inputBams( type => 'file', input => $args{file} );
    }

    my %availCommands;
    $availCommands{$_} = 1 for qw/register drop/;
    my $command =
      defined $args{command}
      ? $args{command}
      : croak "command needs to be defined";
    croak "invalid command: $command" unless defined $availCommands{$command};

    # save input bam list and create missing md5sums
    my %inputBams = %{ $self->inputBams() };
    my %md5sums   = %{ $self->inputBamMD5sums() };

    # make sure same bam files exist in both hashes
    my $areEqual = ( keys %inputBams == keys %md5sums )
      && 0 == grep { !exists $md5sums{$_} } keys %inputBams;

    croak "programming error. list of bams and md5sums don't match."
      unless $areEqual;

    # pull out only name from inputBamFiles
    my @bamFiles = sort keys %md5sums;
    my @bamNames = map { basename($_) } @bamFiles;
    my @bamDirs  = map { dirname($_) } @bamFiles;

    # pull md5sums from files
    my @md5sums;
    for my $md5File ( map { $md5sums{$_} } @bamFiles ) {
        open( my $fh, '<', $md5File )
          or croak "could not read md5 file $md5File";
        my $line = <$fh>;
        chomp $line;
        $line =~ m/^([0-9a-f]+)/;
        push @md5sums, $1;
    }

    # pull sample names from files
    my @sampleNames;
    my @sampleNums;
    for my $bam (@bamFiles) {
        my $rgHeader = `samtools view -H $bam | grep '^\@RG' |head -1`;
        $rgHeader =~ m/\tSM:([A-Za-z0-9_]+)\t/;
        my $sm = $1;
        push @sampleNames, $sm;
        $sm =~ m/(\d+)$/;
        my $sn = $1;
        push @sampleNums, $sn;
    }

    # look for validation files
    my @validated;
    for my $bam (@bamFiles) {
        my $pushVal = 0;
        $pushVal = 1
          if -e "$bam.ok" && ( stat($bam) )[9] < ( stat("$bam.ok") )[9];
        push @validated, $pushVal;
    }

    # keep track of backup flag
    my @backup = map { exists $args{backup} && $args{backup} } @bamFiles;
    croak "need to define backup Device"
      if exists $args{backup} && $args{backup} && !defined $args{backupDevice};

    my @backupDevice =
      map { exists $args{backupDevice} ? $args{backupDevice} : undef }
      @bamFiles;

    my @columns = ( \@md5sums, \@bamNames, \@bamDirs, \@bamFiles, \@backup,
        \@backupDevice );
    my @headers = qw/md5sum bamName bamDir bamPath backup backupDevice/;
    $self->_insertInTable(
        table   => 'bamNames',
        columns => \@columns,
        headers => \@headers,
    );

    $self->_insertInTable(
        table   => 'md5sums',
        columns => [ \@md5sums, \@sampleNums, \@sampleNames, \@validated ],
        headers => [qw/md5sum sampleNum sampleName passedValidateSam/],
    );

    # connect to database
    #    $self->connect();

    # create a hash sorted by md5
    my %bamByMD5;

    #
    return 1;

}

sub _insertInTable {
    my $self   = shift;
    my %args   = @_;
    my %tables = %{ $self->_tables };
    my $table =
      defined $args{table}
      ? $args{table}
      : croak "programming error: need to define input table";

    croak "programming error: table $table is not one of known tables: "
      . join( ' ', sort keys %tables )
      unless exists $tables{ $args{table} };

    my @raColumns =
      ref( $args{columns} ) eq 'ARRAY'
      ? @{ $args{columns} }
      : croak "programming error: columns needs to be an array ref";
    my @headers =
      ref( $args{headers} ) eq 'ARRAY'
      ? @{ $args{headers} }
      : croak "programming error: headers needs to be an array ref";

    scalar @headers == scalar @raColumns
      or croak
      "programming error: number of columns must match number of headers";

    # make sure every column has same length
    my $nRows = undef;
    for my $rCol (@raColumns) {
        unless ( defined $nRows ) {
            $nRows = @{$rCol};
        }
        croak "programming error: not all columns are same length"
          unless $nRows == @{$rCol};
    }

    # get ready to transmit
    my $dbh = $self->dbHandle;

    # assume table is created

    # insert rows
    my $cmd =
        "INSERT INTO $table("
      . join( ',', @headers ) . ")"
      . " VALUES ("
      . join( ',', map { '?' } @headers ) . ")";

    # add in on duplicate update
    $cmd .= " ON DUPLICATE KEY UPDATE $headers[0] = $headers[0]";
    $cmd .= ";";

    my $sth = $dbh->prepare($cmd);

    for my $rowNum ( 0 .. $#{ $raColumns[0] } ) {
        my @row;
        for my $raCol (@raColumns) {
            push @row, $raCol->[$rowNum];
        }
        eval { $sth->execute(@row) };
        croak $sth->errstr . "\n SQL command: $cmd\n" if $@;
    }

    # commit changes
    $dbh->commit;

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

sub dropBams {
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

__PACKAGE__->meta->make_immutable;
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


