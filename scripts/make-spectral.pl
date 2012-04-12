#!/usr/bin/perl

use strict;
use warnings;

@ARGV == 9 || die "expected: lattice ham groundstate minfreq maxfreq increment broadening config resources";

(my $Lattice, my $Ham, my $Gs, my $MinFreq, my $MaxFreq, my $Inc, my $Broad, my $Config, my $Res) = @ARGV;

(my $Energy) = (`mp-expectation $Lattice $Gs $Ham` =~ /\((.*),/);
print "Energy is $Energy\n";

system("mp-apply $Lattice 'CHup(0)' $Gs $Gs.lvch");
system("mp-apply $Lattice 'Cup(0)' $Gs $Gs.lvc");

my $n = 1;
for (my $f = $MinFreq; $f < $MaxFreq+$Inc; $f = $f+$Inc)
{
  my $ActualBroad = -(abs($f * $Broad) + 1.0e-6);
  my $Psi = "$Lattice.gmres-$f-$ActualBroad";
  system("cp $Gs.lvch $Psi.cvch");
  system("cp $Gs.lvc $Psi.cvc");
  system("init.pl mp-gmres-init $Lattice $Ham $Energy $f $ActualBroad $Psi.cvch $Gs.lvch $Config $Lattice.gmres-ch$f-$ActualBroad sch$n $Res");
  system("init.pl mp-gmres-init $Lattice $Ham $Energy $f $ActualBroad $Psi.cvc $Gs.lvc $Config $Lattice.gmres-c$f-$ActualBroad sc$n $Res");
  $n = $n+1;
}
