
#include "lattice/infinitelattice.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"

namespace prog_opt = boost::program_options;

void DescribeLattice(InfiniteLattice const& L)
{
   std::cout << "Description: " << L.description() << '\n';
   std::cout << "Date: " << L.timestamp() << '\n';
   std::cout << "Command line: " << L.command_line() << '\n';
   std::cout << "SymmetryList: " << L.GetSymmetryList() << '\n';
   std::cout << "Unit cell size: " << L.GetUnitCell().size() << '\n';
   for (int i = 0; i < L.GetUnitCell().size(); ++i)
   {
      std::cout << "   site[" << i << "] is: " << L.GetUnitCell()[i].Description() << '\n';
   }

   std::cout << "\nUnit cell operators:\n";
   for (UnitCell::const_operator_iterator I = L.GetUnitCell().begin_operator();
	I != L.GetUnitCell().end_operator(); ++I)
   {
      std::cout << "   " << std::setw(10) << std::left << I->first << " - transforms: "
		<< I->second.TransformsAs() << ", commutes: "
		<< I->second.Commute() << '\n';
   }
 
   std::cout << "\nLattice arguments:\n";
   if (L.arg_empty())
   {
      std::cout << "   (none)\n";
   }
   for (InfiniteLattice::const_argument_iterator I = L.begin_arg(); I != L.end_arg(); ++I)
   {
      std::cout << "   " << std::setw(10) << std::left << I->first << " - " << format_complex(I->second) << '\n';
   }

   std::cout << "\nLattice operators:\n";
   for (InfiniteLattice::const_operator_iterator I = L.begin_operator(); I != L.end_operator(); ++I)
   {
      std::cout << "   " << std::setw(10) << std::left << I->first << " - " 
		<< I->second.description() << '\n';
   }

   std::cout << "\nLattice functions:\n";
   if (L.function_empty())
   {
      std::cout << "   (none)\n";
   }
   for (InfiniteLattice::const_function_iterator I = L.begin_function(); I != L.end_function(); ++I)
   {
      std::cout << "   " << I->first << "=" << I->second << '\n';
   }
}

int main(int argc, char** argv)
{
   try 
   {
      std::string LatticeName;
      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lattice", prog_opt::value<std::string>(&LatticeName), "lattice")
         ;
      prog_opt::positional_options_description p;
      p.add("lattice", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);
      if (vm.count("help")) 
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-lattice-info [options] <lattice>\n";
         std::cerr << desc << "\n";
         return 1;
      }
      
      pvalue_ptr<InfiniteLattice> Lattice = pheap::OpenPersistent(LatticeName, 
								  mp_pheap::CacheSize(), true);

      DescribeLattice(*Lattice);

   }
   catch(std::exception& e) 
   {
      std::cerr << "error: " << e.what() << "\n";
      return 1;
   }
   catch(...) 
   {
      std::cerr << "Exception of unknown type!\n";
   }
   return 0;
}
