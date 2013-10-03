#!/usr/bin/perl -w
our $VERSION = 1.006;
$VERSION = eval $VERSION;

use strict;
use 5.010;
use feature 'switch';
use List::Util qw/sum/;
use Getopt::Std;
use Carp;
use Scalar::Util::Numeric qw(isint);

use Memoize;
#memoize('toDose');
#memoize('get_geno_probs');

=head1 Name

VCF to beagle dose format converter

=head1 Description

Takes a VCF file through stdin. Will rename ID column to
chromosome_number:position format and calculate dose for each
individual. Can handle genotype likelihoods as well.

Use -p to pass in preferred field. (GP, GL, PL, AP or GT)

Use -u to treat GP field as unscaled probabilities

Use -g to output in gprobs instead of dose format

Use -f to specify frequencies file from which to use Hardy Weinberg
priors to convert GLs or PLs to genotype probabilities

Use -h to create hard called genotypes instead of doses

Use -n for no sanity checking on output (speeds things up)

=head2 Changelog

v1.006 - added some optimizations to increase speed (-n most notably)

v1.005 - added hard call (-h) option

v1.002 - Added support for SNPtools AP field

=cut

# skip comments
my $new_line = 1;

my %args;
getopts( 'np:ugf:h', \%args );

our $gHC       = $args{h};
our $gGProbs   = $args{g};
our $gUnscaled = $args{u};

print STDERR "vcf_to_dose.pl v$VERSION\n\n";

# parse freq file
my $hVaraf;
if ( $args{f} ) {
    $hVaraf = getVarAFs( f => $args{f} );
}

=head1 FORMAT Fields Used

included: true if the flag is contained in the bam header
field: location of flag info in FORMAT field

fields are ordered by preference GP > GL > PL > AP > GT

=cut

# determine preferred field
my @ordered_fields  = qw/GP AP GT GL PL/;
my %supported_flags = (
    GP => { included => 0, field => undef }
    ,    # genotype probability (phred scaled)
    GL => { included => 0, field => undef }
    ,    # genotype likelihood (log10 scaled)
    PL => { included => 0, field => undef }
    ,    # genotype likelihood (rounded, phred scaled)
    AP => { included => 0, field => undef },    # allelic probability
    GT => { included => 0, field => undef },    # genotype
);

my $preferredField = $args{p};
if ( defined $preferredField ) {
    confess "need to specify valid preferred field"
      unless defined $supported_flags{$preferredField};
    unshift @ordered_fields, $preferredField;
}

# skipping comments
while ($new_line) {
    $new_line = <>;
    if ( $new_line =~ m/^##/ ) {

        # parse ID tags to search for
        # only register flags we are looking for
        $new_line =~ m/^##FORMAT=\<ID=([a-zA-Z0-9]+)/;
        my $match = $1;
        if ( defined $match && defined $supported_flags{$match} ) {
            $supported_flags{$match}->{included} = 1;
        }
    }
    else {
        last;
    }
}

# determine which information to use for dose estimate
my $preferred_id = undef;
FIELD_ID: for my $fieldID (@ordered_fields) {
    if ( $supported_flags{$fieldID}->{included} ) {
        $preferred_id = $fieldID;
        last FIELD_ID;
    }
}

confess "no field necessary for dose estimation is contained in input file"
  unless defined $preferred_id;

# tell user what we're doing
if ( $args{g} ) {
    print STDERR "Creating GLs in gprobs format from $preferred_id field\n";
}
else {
    print STDERR "Creating dose estimate from $preferred_id field\n";
}
print STDERR "Outputting hard calls\n" if $args{h};
print STDERR "No sanity checks\n" if $args{n};

# print header
chomp($new_line);
my @header_line = split( /\t/, $new_line );
print STDOUT 'marker alleleA alleleB ';

if ( $args{g} ) {
    print STDOUT
      join( ' ', map { ($_) x 3 } @header_line[ 9 .. $#header_line ] ) . "\n";
}
else {
    print STDOUT join( ' ', @header_line[ 9 .. $#header_line ] ) . "\n";
}
$new_line = <>;

# print body
my $lineNum    = 0;
my $first_line = 1;
while ($new_line) {
    $lineNum++;
    print STDERR "Line number $lineNum\r" unless $lineNum % 1000;
    chomp($new_line);
    my @line = split( /\t/, $new_line );
    my $chromLoc = join( q/:/, @line[ 0 .. 1 ] );
    print STDOUT $chromLoc . ' ' . join( q/ /, @line[ 3 .. 4 ] );

    if ($first_line) {
        $first_line = 0;

        # figure out format ordering
        my @format = split( /\:/, $line[8] );
        foreach my $field_num ( 0 .. $#format ) {
            my $fieldID = $format[$field_num];
            if ( defined $supported_flags{$fieldID} ) {

                # check for missing header
                if ( !$supported_flags{$fieldID}->{included} ) {
                    confess "$fieldID does not have a header line";
                }

                # save location of field
                $supported_flags{$fieldID}->{field} = $field_num;
            }
        }
    }

    # parsing individual specific fields
    foreach my $col ( @line[ 9 .. $#line ] ) {

        my @print_val;

        my @col = split( /\:/, $col );
        my $used_col = $col[ $supported_flags{$preferred_id}->{field} ];
        unless ( defined $used_col ) {
            confess '$used_col is not defined';
        }

        # convert input fields to something 0:1 scaled (like, prob, or genotype)
        my @genotype_probs = get_geno_probs( $preferred_id, $used_col );

        # use hardy weinberg prior if freqs file given
        if ( $args{f} ) {
            @genotype_probs = GLToGP(
                prefID   => $preferred_id,
                likes    => \@genotype_probs,
                varafs   => $hVaraf,
                chromLoc => $chromLoc,
            );
        }

        # turn into genotype or dose, or keep probs if -g
        @print_val = toDose( $preferred_id, @genotype_probs );

        unless ( $args{n} ) {

            # sanity check: are all vals between 0 and 2?
            foreach my $val_num ( 0 .. $#print_val ) {
                unless ( $print_val[$val_num] eq '.' ) {
                    my $val = $print_val[$val_num];
                    $print_val[$val_num] = sprintf( "%.4f", $val );
                    if ( $val > 2 || $val < 0 ) {
                        confess "unexpected output: $val at line $lineNum";
                    }
                }
            }

            if ( $args{h} ) {

                # sanity check: there should be only one print val
                confess
"unexpected output: more than one print val ( @print_val ) at line $lineNum"
                  if @print_val != 1;
                confess
"unexpected output: print val is not an integer ( @print_val ) at line $lineNum"
                  if $print_val[0] != sprintf( '%u', $print_val[0] );

                # sanity check: the print val should be integer
            }
        }

        if ( $args{h} ) {
            print STDOUT q/ / . sprintf( '%u', $print_val[0] );
        }
        else {
            print STDOUT q/ / . join( ' ', @print_val );
        }
    }
    print STDOUT "\n";

    $new_line = <>;
}

sub getVarAFs {

    my %args = @_;
    confess "need to define freqs file with -f" unless defined $args{f};

    # read culled freqs file
    my $fh;

    if ( $args{f} =~ m/\.gz$/ ) {
        open( $fh, "gzip -dc $args{f} |" ) or die "could not open $args{f}";
    }
    else {
        open( $fh, "<", $args{f} ) or die "could not open $args{f}";
    }

    # checking header
    my $trash = <$fh>;
    chomp $trash;
    confess "malformed header in $args{f}"
      unless $trash eq "CHROM\tPOS\talleleA\talleleB\tFreqA\tFreqB";

    # parse each line and save to hash
    my %varaf;
    while (<$fh>) {
        chomp;
        my @line = split(/\t/);
        my $key = join( ':', @line[ 0 .. 1 ] );
        $varaf{$key} = $line[5];
    }
    close($fh);

    return \%varaf;
}

# convert GLs to genotype probs if freqs file given
sub GLToGP {
    my %args = @_;

    confess "need to define preferred ID" unless defined $args{prefID};
    confess "need to define chromLoc"     unless defined $args{chromLoc};
    confess "need to define genotype likes"
      unless ref $args{likes} eq 'ARRAY';
    confess "need to pass in VarAFs" unless ref $args{varafs} eq 'HASH';

    # only run if using GLs or PLs
    confess
"multiplying non-genotype likelihood with hardy-weinberg probably makes no sense. Exiting"
      unless $args{prefID} =~ m/GL|PL/;

    # compute HW priors
    my $varaf = $args{varafs}->{ $args{chromLoc} }
      || confess "$args{chromLoc} does not exist in freq file";
    my @hardyWeinbergPriors = (
        ( 1 - $varaf ) * ( 1 - $varaf ),
        2 * ( $varaf * ( 1 - $varaf ) ),
        $varaf * $varaf
    );

    # convert to probs
    my @likes = @{ $args{likes} };
    my @probs = ();
    for my $idx ( 0 .. 2 ) {
        $probs[$idx] = $likes[$idx] * $hardyWeinbergPriors[$idx];
    }
    return @probs;
}

# convert to dose if not GT field used
sub toDose {

#    confess "need to define preferred ID" unless defined $args{prefID};
#    confess "need to define genotype probs" unless ref $args{probs} eq 'ARRAY';

    my $preferred_id = shift;
    my @gprobs       = @_;

    my @print_val;

    # don't change anything if GT field used
    if ($gGProbs) {
        if ( $preferred_id eq 'GT' ) {
            confess "-g option can't work with GT field";
        }
        @print_val = @gprobs;
    }

    # convert from prob/like to dose otherwise
    else {

        given ($preferred_id) {
            when (/^(GP|GL|PL)$/) {
                if ($gHC) {
                    my $idxMax = 0;
                    $gprobs[$idxMax] > $gprobs[$_]
                      or $idxMax = $_
                      for 1 .. $#gprobs;
                    @print_val = ($idxMax);
                }
                else {
                    @print_val =
                      ( ( $gprobs[1] + $gprobs[2] * 2 ) / sum(@gprobs) );
                }
            }
            when ('AP') {
                if ($gHC) {

                    # convert AP to GP
                    my ( $a1, $a2 ) = @gprobs;
                    my @rGProbs = (
                        ( 1 - $a1 ) * ( 1 - $a2 ),
                        ( $a1 * ( 1 - $a2 ) + $a2 * ( 1 - $a1 ) ),
                        ( $a1 * $a2 )
                    );
                    my $idxMax = 0;
                    $rGProbs[$idxMax] > $rGProbs[$_]
                      or $idxMax = $_
                      for 1 .. $#rGProbs;
                    @print_val = ($idxMax);
                }
                else {
                    @print_val = ( sum(@gprobs) );
                }
            }
            when ('GT') {
                @print_val = @gprobs;
            }
            default {
                confess
                  "could not figure out what to do with field $preferred_id"
            }
        }
    }
    return @print_val;

}

# return data as genotype probs/likes scaled [0,1]
sub get_geno_probs {

    #    confess "need to define preferred ID" unless defined $args{prefID};
    #    confess "need to define used col"     unless defined $args{usedCol};

    my $preferred_id = shift;
    my $used_col     = shift;

    my @print_val;
    given ($preferred_id) {
        when ('GP') {
            my @genotype_phred_probability = split( /\,/, $used_col );
            my @genotype_probability;

            if ($gUnscaled) {
                @genotype_probability = @genotype_phred_probability;
            }
            else {
                @genotype_probability =
                  map { 10**( $_ / -10 ) } @genotype_phred_probability;
            }

            @print_val = @genotype_probability;
        }
        when ('GL') {
            my @genotype_log_likelihoods = split( /\,/, $used_col );
            my @genotype_likelihoods =
              map { 10**($_) } @genotype_log_likelihoods;

            @print_val = @genotype_likelihoods;
        }
        when ('PL') {
            my @genotype_phred_likelihoods = split( /\,/, $used_col );
            my @genotype_likelihoods =
              map { 10**( $_ / -10 ) } @genotype_phred_likelihoods;

            @print_val = @genotype_likelihoods;
        }
        when ('AP') {
            @print_val = split( /\,/, $used_col );
        }
        when ('GT') {

            # check for missing data in 'GT' field
            if ( $used_col =~ /\./ ) {
                @print_val = qw'.';
            }
            else {
                #                print STDERR "$lineNum\t$used_col\n";
                my @alleles = split( qr([\|\/]), $used_col );
                unless ( defined $alleles[1] ) {
                    @print_val = qw/./;
                }
                else {
                    @print_val = ( $alleles[0] + $alleles[1] );
                }
            }
        }
        default {
            confess "could not figure out what to do with field $preferred_id"
        }
    }
    return @print_val;
}

