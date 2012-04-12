#!/usr/bin/perl
#
# migrate.pl - PBS job 'migration' script.
# Created 2002-07-09 Ian McCulloch
#
# The intention is to 'migrate' a suspended PBS job, ie kill it,
# and resubmit it.  We test that the current state really is 'suspended'
# before killing it.
#
# Only a small selection of parameters to qsub are preserved over the restart, namely
# -N (job name),
# -o (path for stdout)
# -e (path for stderr)
# -l (resource list)
#
# The default values of Error_Path and Output_Path include the Job Id,
# so it would be misleading if we always set the new path to be the same
# as the old.  Thus, we only override Error_Path or Output_Path if they
# don't match the default, on the assumption that if the user has
# changed them, they would want the same path name for the requeued job.
# In this case, the existing error/output file is copied with a '.migrated'
# suffix attached.
#
# The resource list is only checked for vmem,jobfs, walltime, other and ncpus.  Other
# options are NOT present in the resubmitted job.  It is apparantly
# not appropriate to just copy the entire Resource_List verbatim, so if we ever
# want to preserve other entries they need to be added explicitly.
# The current script should work for all jobs with 4 or fewer CPUs however.

use strict;
use warnings;

sub ParseQStat
{
   my $In = join('', @_);
   $In =~ s/\n\t//g;                # put multi-line fields on one line
   $In =~ s/Job Id:/Job_Id =/;      # Make the Job Id: line consistent format
   $In =~ s/\n    /\n/g;            # get rid of leading whitespace
   $In =~ s/ = /\n/g;               # make ' = ' into a newline separator
   return split(/\n/, $In);         # form a hash/array out of it
}

my $JobID = $ARGV[0];

defined($JobID) || die "usage: migrate <jobid>";

# parse 'qstat -f'
my $QStatOutput = `qstat -f $JobID`;
($? == 0) || die "Cannot qstat job $JobID!";
my %Stat = &ParseQStat($QStatOutput);

my ($ExecHost, $FullJobID, $Name, $Queue, $JobFS,
    $Vmem, $WallTime, $ResOther, $ResNcpus, $ErrorPath, $OutputPath)
  = ($Stat{'exec_host'},$Stat{'Job_Id'}, $Stat{'Job_Name'},
     $Stat{'queue'}, $Stat{'Resource_List.jobfs'},
     $Stat{'Resource_List.vmem'}, $Stat{'Resource_List.walltime'},
     $Stat{'Resource_List.other'}, $Stat{'Resource_List.ncpus'},
     $Stat{'Error_Path'}, $Stat{'Output_Path'});

# see if output or error paths have been changed from default
my $UseOutputPath = ($ErrorPath !~ /".o$JobID"/);
my $UseErrorPath = ($ErrorPath !~ /".e$JobID"/);

# get rid of the hostname in the ErrorPath and OutputPath
$ErrorPath =~ s/^.*://;
$OutputPath =~ s/^.*://;

# get the exec host name
my $RemoteHost = $ExecHost;
$RemoteHost =~ s/\/.*//; # keep only the first listed host

# make sure the current state really is 'S'
($Stat{'job_state'} eq "S") || die "job $JobID is not suspended!";

# copy the job script from $RemoteHost
system("rcp", "$RemoteHost:/local1/pbs/mom_priv/jobs/$FullJobID.SC", "/tmp/$FullJobID.SC");
($? == 0) || die "Could not copy job script from $RemoteHost:/local1/pbs/mom_priv/jobs/$FullJobID.SC";

# kill the job
system("qsig",  "-s", "KILL", "$JobID");

# (paranoid) make sure its dead
`qstat $JobID`;
while ($? == 0)
{
   sleep 2;
   `qstat $JobID`;
}

# assemble the qsub command
my $QSubCommand = "qsub -q " . $Queue . " -l walltime=" . $WallTime . ",vmem=" . $Vmem . ",jobfs=" . $JobFS .
                  ",other=" . $ResOther . ",ncpus-" . $ResNcpus . " -N " . $Name;

if ($UseOutputPath)
{
   system("cp -f $OutputPath $OutputPath.migrated"); # make a backup of the old stdout file
   $QSubCommand .= " -o $OutputPath";
}
if ($UseErrorPath)
{
   system("cp -f $ErrorPath $ErrorPath.migrated"); # make a backup of the old stderr file
   $QSubCommand .= " -e $ErrorPath";
}

$QSubCommand .= " /tmp/$FullJobID.SC";

# resubmit the job
system($QSubCommand);
($? == 0) || die "failed to resubmit job: leaving job script in /tmp/$FullJobID.SC";

# finally, remove the temp copy of the job script
system("rm -f /tmp/$FullJobID.SC");

