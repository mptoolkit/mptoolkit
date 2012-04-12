
#include "matrixproduct/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      std::cout << "usage: mp-iortho <wavefunction>\n";
      return 1;
   }

   std::string FName = argv[1];

   pvalue_ptr<InfiniteWavefunction> Psi = pheap::OpenPersistent(FName, mp_pheap::CacheSize(), true);

   std::cout.precision(14);
   std::cout << (1.0 - orthogonality_fidelity(*Psi)) << '\n';

   pheap::Shutdown();
}
