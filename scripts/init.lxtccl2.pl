#!/usr/bin/perl
#
use strict;
use warnings;
use File::Basename;

sub Panic
{
   die(@_);
}

my $CopyCommand = "cp";

sub BulkCopy
{
   my $StartCopyTime = time;
   my $Filespec = $_[0];
   my @DestDirList = glob $_[1];
   my $DestDirspec = $_[1];
   my $WhichDir = 0;
   my $SuffixToRemove = $_[2];
   my $SuffixToAdd = $_[3];
   my $FilesCopied = 0;
   $SuffixToAdd = "" unless defined($SuffixToAdd);
   $SuffixToRemove = "" unless defined($SuffixToRemove);
   foreach my $File (glob $Filespec)
   {
      $FilesCopied = 1;
      my $BaseFile = $File;
      my $DestDirectory = $DestDirList[$WhichDir];
      $WhichDir = ($WhichDir + 1) % scalar @DestDirList;
      $BaseFile =~ s/.*\///;          # strip the directory
      $BaseFile =~ s/$SuffixToRemove$//;
      $BaseFile = "$DestDirectory/$BaseFile$SuffixToAdd";
      if (-e "$BaseFile") {`rm -f $BaseFile`};
      print "Copying $File to $BaseFile\n";
      if (fork == 0)
      {
	system($CopyCommand, $File, $BaseFile);
	my $Ret = $?;
	print `ls -l $BaseFile`;
	exit $Ret;
      }
   }
   # wait until all copy operations are complete
   my $Failed = 0;
   my $Pid = wait;
   while ($Pid != -1)
   {
      if ($? != 0)
      {
         print "Process $Pid returned exit code $? : $!\n";
	 $Failed = 1;
       }
       $Pid = wait;
   }
   ($Failed) && Panic("Could not copy files to $DestDirspec");
   if ($FilesCopied)
   {
      my $CopyTime = time - $StartCopyTime;
      print "Time to copy files: $CopyTime\n";
   }
}

# parameters are the usual parameters for the init program, then script name then resource list.

my $Resources = pop(@ARGV);
my $ScriptName = pop(@ARGV);
my $BasePathFileName = pop(@ARGV);

defined($BasePathFileName) || die ("expected: parameters script-name resource-list");

my $BasePath = dirname($BasePathFileName);
my $FileName = basename($BasePathFileName);
print "B=$BasePath F=$FileName\n";

my $ExecName = $ARGV[0];
my $ResumeName = $ExecName;
$ResumeName =~ s/-init/-resume/;

push(@ARGV, $BasePathFileName);

$ENV{"TIME_CHECKPOINT"} = 0;
$ENV{"JOBFS_CHECKPOINT"} = 0;
$ENV{"PBS_JOBFS"} = $BasePath;

my $WorkDir = $ENV{"PWD"};

my $Queue;
defined($Queue = $ENV{'INIT_QUEUE'}) || defined($Queue = $ENV{'PBS_QUEUE'}) || ($Queue = "default");

# get the local directory where the data is stored
defined(my $DataDir = $ENV{'MP_DATA_DIRECTORY'}) || die("MP_DATA_DIRECTORY is undefined.");

system(@ARGV);
($? == 0) || die("Execution of " . join(" ",@ARGV) . " Failed - return is $?");

BulkCopy("$BasePath/$FileName.bin*", "$DataDir", "", ".0");
system("rm $BasePath/$FileName.bin*");
system("echo 0 > $BasePath/$FileName.seq");

open(PBSFILE, ">$ScriptName.pbs");
print PBSFILE "#!/bin/bash
trap \"echo job script caught SIGTERM\" SIGTERM
cd $WorkDir
export PBS_JOBFS=/tmp
wrapper.lxtccl2.pl $ScriptName.pbs $ResumeName $BasePathFileName
";
close(PBSFILE);

print "executing:\nqsub -q $Queue -l $Resources -o $ScriptName.o0 -e $ScriptName.e0 -N $ScriptName.0 -r y $ScriptName.pbs\n";
system("qsub -q $Queue -l $Resources -o $ScriptName.o0 -e $ScriptName.e0 -N $ScriptName.0 -r y $ScriptName.pbs");
