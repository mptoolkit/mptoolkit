#!/usr/bin/perl

use strict;
use warnings;

@ARGV == 9 || die "expected: lattice ham energy minfreq maxfreq increment broadening lv config";

(my $Lattice, my $Ham, my $Energy, my $MinFreq, my $MaxFreq, my $Inc, my $Broad, my $Lv, my $Config) = @ARGV;

for (my $f = $MinFreq; $f < $MaxFreq+$Inc; $f = $f+$Inc)
{
   my $Psi = "$Lattice.gmres.cv-$f-$Broad";
   system("cp $Lv $Psi");
   system("mp-gmres-init $Lattice $Ham $Energy $f $Broad $Psi $Lv $Config $Lattice.gmres-$f-$Broad");
}

for (my $f = $MinFreq; $f < $MaxFreq+$Inc; $f = $f+$Inc)
{
   system("mp-gmres-resume $Lattice.gmres-$f-$Broad");
}
