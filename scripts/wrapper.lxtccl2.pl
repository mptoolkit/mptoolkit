#!/usr/bin/perl
#
#
# Wrapper script around the checkpoint/restartable DMRG program
#
# Handles copying the .bin.* and .meta checkpoint files between /fast and /jobfs.
#
# bugs: If job gets SIGTERM while copying files, multiple pid's call Panic resulting
# in multiple mailings.
#
# 2003-03-31: Added other=not_bigsmp to the resource list for resubmitted jobs

use strict;
use warnings;
use File::Basename;
use Cwd;

# some commands to stripe/copy
my $StripeCommand = "echo /usr/sbin/stripe";
my $RemoteStripeCommand = "echo rsh sc1 /usr/sbin/stripe";
#my $CopyCommand = "scfs_cp";
my $CopyCommand = "cp";

sub SendMail
{
   my $Subject = "Job $ENV{'PBS_JOBID'} ($ENV{'PBS_JOBNAME'}) has terminated normally";
   my $PanicString = join(' ', @_);
   my $User = defined($ENV{'PBS_O_LOGNAME'}) ? $ENV{'PBS_O_LOGNAME'} : $ENV{'LOGNAME'};
   my $Host = "physik.rwth-aachen.de";
   my $MailTo = $User . "@" . $Host;
   #`echo "$PanicString\" | mailx -s \"$Subject\" $MailTo\n`;
}

sub Panic
{
   my $Subject = "Job $ENV{'PBS_JOBID'} ($ENV{'PBS_JOBNAME'}) abnormal termination";
   my $PanicString = join(' ', @_);
   my $User = defined($ENV{'PBS_O_LOGNAME'}) ? $ENV{'PBS_O_LOGNAME'} : $ENV{'LOGNAME'};
   my $Host = "physik.rwth-aachen.de";
   my $MailTo = $User . "@" . $Host;
   print "job script terminating on: $PanicString\n";
   #`echo \"$PanicString\" | mailx -s \"$Subject\" $MailTo\n`;
   print `ps -o pid,systime,usertime,time,cmd -u $User`;
   print "\n";
   die $PanicString;
}

sub ParseQStat
{
   my $In = join('', @_);
   $In =~ s/\n\t//g;                # put multi-line fields on one line
   $In =~ s/Job Id:/Job_Id =/;      # Make the Job Id: line consistent format 
   $In =~ s/\n    /\n/g;            # get rid of leading whitespace
   $In =~ s/ = /\n/g;               # make ' = ' into a newline separator
   return split(/\n/, $In);         # form a hash out of it
}

sub ConvertToKB
{
   my $Var = $_[0];
   if ($Var =~ s/mb$//) { return $Var * 1024 }
   elsif ($Var =~ s/gb$//) { return $Var * 1024 * 1024 }
   else {$Var =~ s/kb$// || Panic "size must be in kb, mb or gb." };
   return int($Var);
}

sub SecsToTime
{
   (my $hours, my $mins, my $secs) = (int($_[0] / 3600), (int($_[0] / 60) % 60), $_[0] % 60);
   return join(':', $hours < 10 ? '0'.$hours : $hours,
	            $mins < 10 ? '0'.$mins: $mins,
                    $secs < 10 ? '0'.$secs : $secs);
}

sub TimeToSecs
{
   my ($hours, $mins, $secs) = split(':', $_[0]);
   return (($hours * 60) + $mins) * 60 + $secs;
}

#
# BulkCopy
#
# subroutine to copy files given by a glob in $_[0].
# Destination directory in $_[1].  This can be a glob, if multiple directories are requested.
# in this case, the files are written round-robin into each directory.
# If there is a suffix on the source files that should be removed, put it in $_[2].
# if there is a suffix to add to the destination files, put it in $_[3].
# Uses the command specified in the global variable $CopyCommand to do the actual copying.
#
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
	my $Ret = system($CopyCommand, $File, $BaseFile);
	print `ls -l $BaseFile`;
	exit $Ret / 256;
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

# ignore SIGTERM
$SIG{TERM} = sub { syswrite STDERR,"wrapper.pl caught SIGTERM\n" };
# ignore SIGSEGV (assume that its been sent by qsig, for debugging purposes)
$SIG{SEGV} = sub { syswrite STDERR,"wrapper.pl caught SIGSEGV\n" };
$SIG{USR1} = sub { syswrite STDERR,"wrapper.pl caught SIGUSR1\n" };

use POSIX ':signal_h';

sigaction SIGTERM, new POSIX::SigAction sub { syswrite STDERR,"wrapper caught SIGTERM via sigaction.\n" }
        or die "Error setting SIGTERM handler: $!\n";

# load the arguments
(my $ScriptName, my $ExecName, my $BasePathFileName) = @ARGV;

defined($BasePathFileName) || die("Not enough arguments supplied: expected ScriptName, ExecutableName, PathName");

my $BasePath = dirname($BasePathFileName);
my $FileName = basename($BasePathFileName);

# get the local directory where the data is stored
defined(my $DataDir = $ENV{'MP_DATA_DIRECTORY'}) || Panic "MP_DATA_DIRECTORY is undefined.";

# get the preferred queue to use for the restart
my $PBSQueue;
defined($PBSQueue = $ENV{'RESUB_QUEUE'}) || defined($PBSQueue = $ENV{'PBS_QUEUE'}) || Panic "PBS_QUEUE is undefined.";

# get the path to jobfs directory
defined(my $BinPath = $ENV{'PBS_JOBFS'}) || Panic "PBS_JOBFS is undefined.";
defined(my $JobID = $ENV{'PBS_JOBID'}) || Panic "PBS_JOBID is undefined.";

print "PBS_JOBFS is $ENV{'PBS_JOBFS'}\n";

my $ShortScriptName = $ScriptName;
$ShortScriptName =~ s/\.pbs$//;

# hard-coded maximum size of jobfs allocations (in KB).  If this is exceeded,
# the job must be restarted manually.
#my $MaxJobFS = 10485760;          # 10 GB
my $MaxJobFS = 10485760 * 4;          # 40 GB
# hard-coded conservative guess as to how long it takes to checkpoint & copy, in seconds.
my $CheckpointExtra = 3600;       # 1 hour
#my $CheckpointExtra = 180;
# checkpoint interval, in seconds
#my $CheckpointInterval = 18000;   # 5 hours
my $CheckpointInterval = 3600 * 40;
#my $CheckpointInterval = 300;

# parse 'qstat -f' to get our jobfs, vmem and walltime allocations
my $QStatOutput = `qstat -f $JobID`;
my %List = &ParseQStat($QStatOutput);
print "Run started at " . scalar gmtime(time) . "\n";
print "Output of qstat -f is\n$QStatOutput\n";

# see if the output file(s) already exist; if so, then this is
# likely to be a restarted run, and we want to include the
# previous output in the current stdout/stderr.

my $StdoutFile = $List{'Output_Path'};
$StdoutFile =~ s/.*://;  # remove the hostname from the path
my $StderrFile = $List{'Error_Path'};
$StderrFile =~ s/.*://;

if (-e $StdoutFile)
{
   print "\nOutput file $StdoutFile already exists, output is included here:\n";
   print   "===========\n";
   print `cat $StdoutFile`;
   print   "===========End included output file\n";
 }

if (-e $StderrFile && -s $StderrFile)
{
   print STDERR  "Error file $StderrFile already exists, output is included here:\n";
   print STDERR  "==========\n";
   print STDERR  `cat $StderrFile`;
   print STDERR  "==========End included error file\n";
}

# get the jobfs allocation
#my $ReqJobFS = ConvertToKB($List{'Resource_List.jobfs'});

# get the vmem size (in K)
my $VMem = ConvertToKB($List{'Resource_List.mem'});

# get the walltime allocation
my $Walltime = TimeToSecs($List{'Resource_List.walltime'});

# set the JOBFS_CHECKPOINT variable to our allocated jobfs size minus the mem size (to allow for the metadata)
#my $jc = $ReqJobFS - $VMem;
#if ($jc < $VMem) { $jc = $ReqJobFS * 0.5 }; # this is really only for testing
#$ENV{'JOBFS_CHECKPOINT'} = int($jc);
#print "JOBFS_CHECKPOINT is $ENV{'JOBFS_CHECKPOINT'}\n";

# make sure that the sequence number file exists
(-e "$BasePath/$FileName.seq") || 
      Panic "Sequence number file $BasePath/$FileName.seq does not exist!";

chop(my $Suffix = `cat $BasePath/$FileName.seq`);

BulkCopy("$DataDir/$FileName.bin*.$Suffix", $BinPath, ".$Suffix");

# this is where the main loop starts
my $ForceRequeue = 0; # when this is set, a requeue has been forced irrespective of the walltime remaining.

my $HasRun = 0;  # the loop structure below is seriously fscked, this flag is set to make sure we don't
                 # go into an infinite qsub loop if the allocated walltime is less than the minimum.

# touch all files matching $BasePath/$FileName* .  This prevents output files from being deleted
# if they are stored on a tempoary filesystem, ven if they are not used in a given session
system("touch $BasePath/$FileName.*");

while (1)
{
   # calculate how much walltime is remaining
   %List = &ParseQStat(`qstat -f $JobID`);

   (my $hoursUsed, my $minsUsed, my $secsUsed)
     = defined($List{'resources_used.walltime'}) ? split(':', $List{'resources_used.walltime'}) : (0,0,0);
   my $WalltimeUsed = (($hoursUsed * 60) + $minsUsed) * 60 + $secsUsed;
   my $WalltimeRemain = $Walltime - $WalltimeUsed;

   print "Walltime used: $WalltimeUsed, walltime remaining: $WalltimeRemain\n";

   # if we don't have enough walltime left, requeue and exit.
   if ($WalltimeRemain < 2 * $CheckpointExtra || $ForceRequeue)
   {
      # restart

      # paranoid check that the program has run at least once
      $HasRun || Panic("Allocated walltime is insufficient to run program.");

      my $dir = cwd;
      my $ExecString = "ssh master /opt/pbs/bin/qsub -q $PBSQueue -l nodes=1:ppn=1,mem=$VMem"."kb,walltime="
        .SecsToTime($Walltime)." -o $ShortScriptName.o$Suffix -e $ShortScriptName.e$Suffix "
	  ."-N $ShortScriptName.$Suffix -r y $dir/$ScriptName";

      print "executing: " . $ExecString . "\n";

      # submit the job and get the identifier
      my $NewJobID = `$ExecString`;
      chomp($NewJobID);
      # shift it to the top of the queue
      system("qshifttop", $NewJobID);
      exit(0);
   }
   # otherwise calculate how long we should run for

   my $WalltimeReq = $CheckpointInterval;
   # if there a full checkpoint interval left?
   if ($WalltimeReq > $WalltimeRemain - $CheckpointExtra) { $WalltimeReq = $WalltimeRemain - $CheckpointExtra };
   # if we use a full interval + CheckpointExtra, will we leave enough time to restart in this session?
   if ($WalltimeRemain - ($WalltimeReq + $CheckpointExtra) < 2 * $CheckpointExtra)
   {
      # if not, then we might as well keep going until the end of the session
      $WalltimeReq = $WalltimeRemain - $CheckpointExtra;
   }

   # set the TIME_CHECKPOINT variable to our allocated walltime minus the predetermined reserve
   $ENV{'TIME_CHECKPOINT'} = $WalltimeReq;
   print "TIME_CHECKPOINT is $ENV{'TIME_CHECKPOINT'}\n";

   # run the program
   print "executing: $ExecName $BasePathFileName\n";
   my $exit_value = system($ExecName, $BasePathFileName);
   #print "Raw exit code is $exit_value\n";
   ($exit_value != -1) || Panic "Failed to execute $ExecName.  Reason: $!";
   if ($exit_value > 0 && $exit_value < 256)
   {
      Panic("$ExecName aborted on signal $exit_value");
   }

   $exit_value >>= 8;

   if ($exit_value == 0)
   {
      # normal exit, we should delete the persistent files.
      `rm -f $DataDir/$FileName.bin*.$Suffix`;

      SendMail("Terminated normally,run is complete.");

      print "Terminated normally.\n";
      exit(0);
   };

   if ($exit_value == 1) { Panic "Abnormal termination!" };
   if ($exit_value > 8) { Panic "Abnormal termination, return code is $exit_value " };

   # if we got here, we should have checkpointed OK
   # so copy the bin files back over to $BasePath

   # update newsuffix/oldsuffix
   my $OldSuffix = $Suffix;
   $Suffix = $Suffix + 1;

   BulkCopy("$BinPath/$FileName.bin*", "$DataDir", "", ".$Suffix");

   # 'atomically' notify that our checkpoint files are valid
   `echo $Suffix > $BasePath/$FileName.seq`;

   # remove the old checkpoint files
   `rm -f $DataDir/$FileName.bin*.$OldSuffix`;

   # see if we have terminated abnormally
   if ($exit_value == 2) { Panic "Caught SIGTERM - not auto-requeueing." };

   if ($exit_value == 4) { Panic "Synchronous checkpoint - not auto-requeueing." };
   if ($exit_value == 5) { Panic "Huh?  Jobfs?" };
   if ($exit_value == 6) { Panic "Memory limit exceeded, please teach me how to increase it." };
   if ($exit_value == 7) { print "Caught SIGUSR1, forced checkpoint, restarting.\n" };

   if ($exit_value == 8)
   {
     # increase memory allocation and requeue
     my $NewVMem = int($VMem * 1.2);   # increase mem allocation by 20%
     print "Caught SIGUSR2, Increasing mem from " . $VMem . " to " . $NewVMem . " and requeueing.\n";
     $VMem = $NewVMem;
     $ForceRequeue = 1;
   }

   $HasRun = 1;
}
