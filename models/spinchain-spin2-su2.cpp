// -*- C++ -*- $Id$

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-su2.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 2;
      std::string LatticeName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("out,o", prog_opt::value(&LatticeName), "output filename [required]")
         ;
      
      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
					  prog_opt::command_line_style::allow_guessing).
		      run(), vm);
      prog_opt::notify(vm);    
      
      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: spinchain-spin2-su2 [options]\n";
         std::cerr << desc << '\n';
	 std::cerr << "Operators:\n"
		   << "H_p0      - projector onto nearest-neighbor spin 0 bond\n"
		   << "H_p1      - projector onto nearest-neighbor spin 1 bond\n"
		   << "H_p2      - projector onto nearest-neighbor spin 2 bond\n"
		   << "H_p3      - projector onto nearest-neighbor spin 3 bond\n"
		   << "H_p4      - projector onto nearest-neighbor spin 4 bond\n"
		   << "\n"
		   << "H_di     - next-nearest neighbor dipole spin exchange     (S.S)\n"
		   << "H_quad   - next-nearest neighbor quadrupole spin exchange (Q.Q)\n"
		   << "H_hex    - next-nearest neighbor hexapole spin exchange   (T.T)\n"
		   << "H_oct    - next-nearest neighbor octapole spin exchange   (F.F)\n"
	    ;
         return 1;
      }

      LatticeSite Site = CreateSU2SpinSite(Spin);
      UnitCell Cell(Site);
      InfiniteLattice Lattice(Cell);
      UnitCellOperator S(Cell, "S"), Q(Cell, "Q"), T(Cell, "T"), F(Cell, "F"), I(Cell, "I");

      // (Q.Q) = -1/3 (S^2 . S^2) + (S.S)^2 + 1/2 S.S

      TriangularMPO& H_di   = Lattice["H_di"];
      TriangularMPO& H_quad = Lattice["H_quad"];
      TriangularMPO& H_hex  = Lattice["H_hex"];
      TriangularMPO& H_oct  = Lattice["H_oct"];

      H_di   = sum_unit(inner(S(0), S(1)));
      H_quad = sum_unit(inner(Q(0), Q(1)));
      H_hex  = sum_unit(inner(T(0), T(1)));
      H_oct  = sum_unit(inner(F(0), F(1)));

      TriangularMPO c = sum_unit(I(0)); // constant energy shift

      // These magic values for the spin 2 model projectors come from solving the equations
      // in misc/spin2.cpp
      Lattice["H_p0"] = (-1/50.0)*H_di + (1/105.0) *H_quad + (-1/180.0)*H_hex + (1/180.0) *H_oct + (1/25.0)*c;
      Lattice["H_p1"] = (-1/20.0)*H_di + (1/70.0)  *H_quad + (0.0)     *H_hex + (-1/90.0) *H_oct + (3/25.0)*c;
      Lattice["H_p2"] = (-1/20.0)*H_di + (-1/98.0) *H_quad + (1/63.0)  *H_hex + (1/126.0) *H_oct + (5/55.0)*c;
      Lattice["H_p3"] = (0.0)    *H_di + (-4/105.0)*H_quad + (-1/72.0) *H_hex + (-1/360.0)*H_oct + (7/25.0)*c;
      Lattice["H_p4"] = (3/25.0) *H_di + (6/245.0) *H_quad + (1/280.0) *H_hex + (1/2520.0)*H_oct + (9/25.0)*c;

      pheap::ExportObject(LatticeName, Lattice);
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
