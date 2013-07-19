# Test file created outside of h2xs framework.
# Run this like so: `perl 01-load.t'
#   wkretzsch@gmail.com     2013/07/19 14:58:07

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use Test::More;
use FindBin;
use lib "$FindBin::Bin/..";

BEGIN { plan tests => 2 };

use warnings;
use strict;
$|=1;
use Data::Dumper;

use_ok( 'bamTrackLib' ); 


ok(1); # If we made it this far, we're ok.



#########################

# Insert your test code below, the Test::More module is used here so read
# its man page ( perldoc Test::More ) for help writing this test script.


