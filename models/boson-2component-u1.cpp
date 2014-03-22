-*- C++ -*- $Id$

#include "models/boson-u1.h"
#include "lattice/unitcell.h"
#include "lattice/localoperator.h"
#include "lattice/infinitelattice.h"
#include "common/prog_options.h"
#include "common/environment.h"

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int NMax = 5;
      std::string LatticeName;
      int L = 0;
      double Lambda = 0.5;
      double Xi = 0;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("NumBosons,N", prog_opt::value(&NMax), FormatDefault("Maximum number of bosons per site", NMax).c_str())
         ("out,o", prog_opt::value(&LatticeName), "output filename [required]")
         ;
      
      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).run(), vm);
      prog_opt::notify(vm);    
      
      if (vm.count("help") || !vm.count("LatticeSize") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: boson-2component-u1 [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Defined operators:\n"
                   << "H_Ja         Nearest-neighbor hopping for component a\n"
                   << "H_Jb         Nearest-neighbor hopping for component b\n"
                   << "H_Jab        On-site tunneling between components a and b\n"
                   << "H_Ua         On-site Coulomb repulsion for component a\n"
                   << "H_Ub         On-site Coulomb repulsion for component b\n"
                   << "H_Uab        On-site Coulomb repulsion between components a and b\n"
            ;
         return 1;
      }

      // the basic lattice site
      LatticeSite Site = BosonU1(NMax);

      // our unit cell is two sites
      UnitCell Cell(Site, Site);

      LocalOperator N(Cell, "N");
      LocalOperator BH(Cell, "BH");
      LocalOperator B(Cell, "N");
      LocalOperator B(Cell, "N2");

      // Total particle number
      Cell.Operator("N") = N[0] + N[1];

      // particle number difference
      Cell.Operator("DN") = N[0] - N[1];

      // Tunneling operator
      Cell.Operator("Jab") = BH[0]*B[1] + B[0]*BH[1];

      // Construct the lattice
      InfiniteLattice Lattice(Cell);

      TriangularMPO DN = sum_unit(Cell, Cell.Operator("DN"));

      TriangularMPO H_Ja = sum_unit(Cell, BH(0)[0]*B(1)[0] + B(0)[0]*BH(1)[0]);
      TriangularMPO H_Jb = sum_unit(Cell, BH(0)[1]*B(1)[1] + B(0)[1]*BH(1)[1]);

      TriangularMPO H_Ua = sum_unit(Cell, N2(0)[0]);
      TriangularMPO H_Ub = sum_unit(Cell, N2(0)[1]);

      TriangularMPO H_Uab = sum_unit(Cell, N(0)[0]*N(0)[1]);
      
      TriangularMPO H_Jab = sum_unit(Cell, Jab(0));

      Lattice.add("DN", DN);
      Lattice.add("H_Ja", H_Ja);
      Lattice.add("H_Jb", H_Jb);
      Lattice.add("H_Jab", H_Jab);
      Lattice.add("H_Ua", H_Ua);
      Lattice.add("H_Ub", H_Ub);
      Lattice.add("H_Uab", H_Uab);

      // make a copy of OpList that exists on the persistent heap
      pvalue_ptr<InfiniteLattice> OList = new InfiniteLattice(Lattice);

      pheap::Initialize(LatticeName, 1, 65536, 655360);
      pheap::ShutdownPersistent(OList);
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
