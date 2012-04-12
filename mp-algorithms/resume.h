// -*- C++ -*- $Id$

#if !defined(RESUME_H_DSHF43758U89JP98)
#define RESUME_H_DSHF43758U89JP98

#include "pheap/pheap.h"
#include "mp/copyright.h"
#include <iostream>

template <typename SolverType>
int GenericResume(int argc, char** argv, std::string const& ProgName)
{
   if (argc != 2)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: " << ProgName << " <base-path>\n";
      return 1;
   }

   {
   pvalue_ptr<SolverType> Solver;
   try
   {
      std::string BasePathFull = argv[1];
      std::string BasePath, FileName;
      boost::tie(BasePath, FileName) = pheap::SplitPathFile(BasePathFull);
      if (BasePath.empty()) BasePath = "./";

      ConfList Conf(BasePath+FileName+".conf");

      long PageCacheSize = Conf.GetBytes("PageCacheSize", 16*1024*1024);
      std::string BinPath = Conf.Get("BinPath", std::string("."));

      if (BinPath[BinPath.size()-1] != '/') BinPath += '/';

      MessageLogger::Logger("PHEAPFS").SetThreshold(Conf.Get("PHeapLogLevel", 0));

      Solver = pheap::OpenPersistent(BinPath+FileName+".bin", PageCacheSize);

      ProcControl::Initialize(argv[0], 
                              Solver->PreviousCPUTime(),
                              Solver->PreviousElapsedTime(),
                              true);

      double TimeLimit = Conf.Get("MaxCPUTime", 0.0);
      if (TimeLimit > 0)
         ProcControl::SetCPUTimeLimit(TimeLimit);

      TimeLimit = Conf.Get("MaxWallTime", 0.0);
      if (TimeLimit > 0)
      {
         std::cout << "MaxWallTime = " << TimeLimit << '\n';
         ProcControl::SetWallTimeLimit(TimeLimit);
      }

#if defined(CONFIG_PBS_WALLTIME)
      TimeLimit = Conf.Get("MaxPBSWallTime", 0.0);
      if (TimeLimit > 0)
      {
         std::cout << "MaxPBSWallTime = " << TimeLimit << '\n';
         ProcControl::SetWallTimeLimit(TimeLimit);
      }
#endif

      std::cout << "Restarting...\n";
      Solver.mutate()->Run(BasePathFull, Conf);
      std::cout << "Finished.\n";
      std::cout << "CPU time this run: " << ProcControl::GetCPUTime() << '\n';
      std::cout << "Total CPU time: " << ProcControl::GetCumulativeCPUTime() << '\n';
      std::cout << "Elapsed time this run: " << ProcControl::GetElapsedTime() << '\n';
      std::cout << "Total elapsed time: " << ProcControl::GetCumulativeElapsedTime() << '\n';
   }
   catch (ProcControl::Checkpoint const& c)
   {
      std::cout << "Initiating checkpoint.  Reason: " << c.Reason() << std::endl;
      Solver.mutate()->SetPreviousCPUTime(ProcControl::GetCumulativeCPUTime());
      Solver.mutate()->SetPreviousElapsedTime(ProcControl::GetCumulativeElapsedTime());
      pheap::ShutdownPersistent(Solver);
      std::cout << "CPU time this run: " << ProcControl::GetCPUTime() << '\n';
      std::cout << "Cumulative CPU time: " << ProcControl::GetCumulativeCPUTime() << '\n';
      std::cout << "Elapsed time this run: " << ProcControl::GetElapsedTime() << '\n';
      std::cout << "Cumulative elapsed time: " << ProcControl::GetCumulativeElapsedTime() << '\n';
      ProcControl::Shutdown();
      return c.ReturnCode();
   }
   catch (std::exception& E)
   {
      std::cerr << "Fatal Exception: " << E.what() << std::endl;
      exit(1);
   }
   catch (...)
   {
      std::cerr << "Unknown exception.\n";
      exit(1);
   }
   }

   pheap::Shutdown();
   ProcControl::Shutdown();

   return 0;
}

#endif
