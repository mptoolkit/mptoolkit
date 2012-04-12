#!/usr/bin/perl

use strict;
use warnings;

my $Job = $ARGV[0];
my $Resources = $ARGV[1];

defined($Resources) || die("expected: resume.pl <jobid> <resource-list>");

my $ScriptName = $Job;
if ($ScriptName !~ /\.pbs$/)
{
  $ScriptName = $ScriptName . ".pbs";
}

open(PBSFILE, "<$ScriptName");
while (<PBSFILE>)
{
   if (/^wrapper/)
   {
      my @ThisLine = split(" ", $_);
      my $BasePath = $ThisLine[3];
      print "BasePath is $BasePath\n";
      my $Seq = `cat $BasePath.seq`;
      chop $Seq;
      print "Sequence number is $Seq\n";
      my $exestring = "qsub -q default -l $Resources -o $Job.o$Seq -e $Job.e$Seq -N $Job.$Seq -r y $ScriptName";
      print ("executing: $exestring\n");
      system($exestring);
    }
 }
