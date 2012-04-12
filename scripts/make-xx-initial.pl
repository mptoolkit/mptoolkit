#!/usr/bin/perl

# make the initial state for an xx spin chain, |up ... up up down down ... down>

use strict;
use warnings;

@ARGV == 1 || die "usage: make-xx-initial.pl <number-of-lattice-sites>";

my $Size = $ARGV[0];

for (my $i = 0; $i < $Size/2; ++$i)
{
   print "-0.5\n";
}
for (my $i = 0; $i < $Size/2; ++$i)
{
   print "0.5\n";
}
