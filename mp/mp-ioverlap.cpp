// -*- C++ -*- $Id$

#include "mps/infinitewavefunction.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "lattice/unitcell.h"
#include "lattice/product-parser.h"
#include "common/statistics.h"

namespace prog_opt = boost::program_options;

struct TransEigenInfo
{
   TransEigenInfo() {}
   TransEigenInfo(QuantumNumber const& q_, std::complex<double> const& x_) : q(q_), x(x_) {}

   QuantumNumber q;
   std::complex<double> x;
};

// sort the magnitude of the eigenvalue in reverse order
bool operator<(TransEigenInfo const& x, TransEigenInfo const& y)
{
   return LinearAlgebra::norm_frob_sq(x.x) > LinearAlgebra::norm_frob_sq(y.x);
}

void PrintFormat(TransEigenInfo const& x, bool ShowRealPart, bool ShowImagPart, 
		 bool ShowCorrLength, bool ShowMagnitude, bool ShowArgument,
		 bool ShowRadians, double ScaleFactor)
{
   std::string SectorStr = boost::lexical_cast<std::string>(x.q);
   std::complex<double> Value = std::pow(x.x, ScaleFactor);
   std::cout << std::setw(11) << SectorStr << ' ';
   if (ShowRealPart)
   {
      std::cout << std::setw(20) << Value.real() << "    ";
   }
   if (ShowImagPart)
   {
      std::cout << std::setw(20) << Value.imag() << "    ";
   }
   if (ShowCorrLength)
   {
      std::cout << std::setw(20) << (-1.0/std::log(LinearAlgebra::norm_frob(Value)))
		<< "    ";
   }
   if (ShowMagnitude)
   {
      std::cout << std::setw(20) << LinearAlgebra::norm_frob(x.x) << "    ";
   }
   if (ShowArgument)
   {
      double Arg =  std::atan2(Value.imag(), Value.real());
      if (!ShowRadians)
	 Arg *= 180.0 / math_const::pi;
      std::cout << std::setw(20) << Arg << "    ";
   }
   std::cout << std::endl;
}


InfiniteWavefunction conj(InfiniteWavefunction const& Psi)
{
   InfiniteWavefunction Ret;
   Ret.C_old = conj(Psi.C_old);
   Ret.C_right = conj(Psi.C_right); 
   Ret.Attr = Psi.Attr;
   Ret.QShift = Psi.QShift;
   Ret.Psi = LinearWavefunction(Psi.Psi.GetSymmetryList());
   for (LinearWavefunction::const_iterator I = Psi.Psi.begin();
        I != Psi.Psi.end(); ++I)
   {
      Ret.Psi.push_back(conj(*I));
   }

   return Ret;
}

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      bool NoTempFile = false;
      bool ShowRealPart = false, ShowImagPart = false, ShowMagnitude = false;
      bool ShowCartesian = false, ShowPolar = false, ShowArgument = false;
      bool ShowRadians = false, ShowCorrLength = false;
      int Rotate = 0;
      int UnitCellSize = 0;
      std::string LhsStr, RhsStr;
      std::vector<std::string> Sector;
      double Tol = 1E-10;
      int Iter = 30;
      bool Sort = false;
      bool Quiet = false;
      bool Reflect = false;
      bool Conj = false;
      bool Print = false;
      std::string String;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("cart,c", prog_opt::bool_switch(&ShowCartesian),
	  "show the result in cartesian coordinates [equivalent to --real --imag]")
	 ("polar,p", prog_opt::bool_switch(&ShowPolar),
	  "show the result in polar coodinates [equivalent to --mag --arg]")
	 ("real,r", prog_opt::bool_switch(&ShowRealPart),
	  "display the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImagPart),
	  "display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("arg,a", prog_opt::bool_switch(&ShowArgument),
          "display the argument (angle) of the result")
	 ("radians,r", prog_opt::bool_switch(&ShowRadians),
	  "display the argument in radians instead of degrees")
	 ("corr,x", prog_opt::bool_switch(&ShowCorrLength),
	  "display the equivalent correlation length")
	 ("unitcell,u", prog_opt::value(&UnitCellSize),
	  "scale the results to use this unit cell size [default wavefunction unit cell]")
         ("notempfile", prog_opt::bool_switch(&NoTempFile),
          "don't use a temporary data file, keep everything in RAM "
          "(faster, but needs enough RAM)")
         ("rotate", prog_opt::value(&Rotate),
          "rotate the unit cell of psi2 this many sites to the left before calculating the overlap [default 0]")
         ("reflect", prog_opt::bool_switch(&Reflect),
          "reflect psi2 (gives parity eigenvalue)")
         ("string", prog_opt::value(&String),
          "use this unit cell operator as a string operator for the overlap")
         ("conj", prog_opt::bool_switch(&Conj),
          "complex conjugate psi2")
         ("q,quantumnumber", prog_opt::value(&Sector),
          "calculate the overlap only in this quantum number sector, "
          "can be used multiple times [default is to calculate all sectors]")
         ("sort,s", prog_opt::bool_switch(&Sort),
          "sort the eigenvalues by magnitude")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
         ("iter", prog_opt::value(&Iter),
          FormatDefault("Maximum subspace size in the Arnoldi basis", Iter).c_str())
         ("quiet", prog_opt::bool_switch(&Quiet),
          "don't show the column headings")
         ("print", prog_opt::bool_switch(&Print), "with --string, Print the MPO to standard output")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose),
          "extra debug output [can be used multiple times]")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lhs", prog_opt::value<std::string>(&LhsStr), "psi1")
         ("rhs", prog_opt::value<std::string>(&RhsStr), "psi2")
         ;

      prog_opt::positional_options_description p;
      p.add("lhs", 1);
      p.add("rhs", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("rhs") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: mp-ioverlap [options] <psi1> <psi2>\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // If no output switches are used, default to showing everything
      if (!ShowRealPart && !ShowImagPart && !ShowMagnitude
	  && !ShowCartesian && !ShowPolar && !ShowArgument
	  && !ShowRadians && !ShowCorrLength)
      {
	 ShowCartesian = true;
	 ShowPolar = true;
	 ShowCorrLength = true;
      }
      
      if (ShowCartesian)
      {
	 ShowRealPart = true;
	 ShowImagPart = true;
      }
      if (ShowPolar)
      {
	 ShowMagnitude = true;
	 ShowArgument = true;
      }
      if (ShowRadians)
	 ShowArgument = true;      

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));

      if (Verbose)
         std::cout << "Loading LHS wavefunction...\n";

      pvalue_ptr<InfiniteWavefunction> Psi1;
      if (NoTempFile)
         Psi1 = pheap::OpenPersistent(LhsStr, mp_pheap::CacheSize(), true);
      else
      {
         mp_pheap::InitializeTempPHeap(Verbose);
         Psi1 = pheap::ImportHeap(LhsStr);
      }

      if (Verbose)
         std::cout << "Loading RHS wavefunction...\n";

      pvalue_ptr<InfiniteWavefunction> Psi2 = RhsStr.empty() ? Psi1 :
         pheap::ImportHeap(RhsStr);

      if (Verbose)
         std::cout << "Calculating overlap...\n";

      // Rotate as necessary.  Do this BEFORE determining the quantum number sectors!
      if (Rotate > 0)
         *Psi2.mutate() = rotate_left(*Psi2, Rotate);

      UnitCell Cell;
      LatticeSite Site;

      if (Reflect)
      {
         if (Verbose)
            std::cout << "Reflecting psi2..." << std::endl;
	 *Psi2.mutate() = reflect(*Psi2);
      }

      if (Conj)
      {
         if (Verbose)
            std::cout << "Conjugating psi2..." << std::endl;
         *Psi2.mutate() = conj(*Psi2);
      }

      // If the wavefunctions don't match in size, then increase them to the lowest common multiple
      if (Psi1->Psi.size() != Psi2->Psi.size())
      {
	 int Size = statistics::lcm(Psi1->Psi.size(), Psi2->Psi.size());
	 if (Verbose)
	 {
	    std::cerr << "Wavefunction unit cells differ, extending wavefunctions to size " << Size << '\n';
	 }
	 *Psi1.mutate() = repeat(*Psi1, Size/Psi1->size());
	 *Psi2.mutate() = repeat(*Psi2, Size/Psi2->size());
	 DEBUG_CHECK_EQUAL(Psi1->size(), Psi2->size());
      }

      // The default UnitCellSize for output is the wavefunction size
      if (UnitCellSize == 0)
      {
	 UnitCellSize = Psi1->Psi.size();
      }
      double ScaleFactor = double(UnitCellSize) / double(Psi1->Psi.size());
	 
      ProductMPO StringOp;
      if (vm.count("string"))
      {
	 InfiniteLattice Lattice;
	 boost::tie(StringOp, Lattice) = ParseProductOperatorAndLattice(String);
	 if (Print)
	 {
	    std::cout << "String MPO is:\n" << StringOp << '\n';
	 }
	 CHECK(Psi1->size() % StringOp.size() == 0)
	    ("Wavefunction size must be a multiple of the string operator size")
	    (Psi1->size())(StringOp.size());
	 StringOp = repeat(StringOp, Psi1->size() / StringOp.size());
      }
      else
      {
         StringOp = ProductMPO::make_identity(ExtractLocalBasis(Psi2->Psi));
      }

      // get the list of quantum number sectors
      std::set<QuantumNumber> Sectors;

      if (!Sector.empty())
      {
         for (std::vector<std::string>::const_iterator I = Sector.begin(); I != Sector.end(); ++I)
         {
            Sectors.insert(QuantumNumber(Psi1->Psi.GetSymmetryList(),*I));
         }
      }
      else
      {
         BasisList B1 = Psi1->Psi.Basis1().Basis();
         BasisList B2 = Psi2->Psi.Basis1().Basis();
         for (unsigned i = 0; i < B1.size(); ++i)
         {
            for (unsigned j = 0; j < B2.size(); ++j)
            {
               inverse_transform_targets(B2[j], B1[i], std::inserter(Sectors, Sectors.begin()));
            }
         }
      }

      if (!Quiet)
      {
	 std::cout << "#" << argv[0];
	 for (int i = 1; i < argc; ++i)
	    std::cout << ' ' << argv[i];
	 std::cout << "\n#quantities are calculated per unit cell size of " << UnitCellSize 
		   << (UnitCellSize == 1 ? " site\n" : " sites\n");
         std::cout << "#sector     ";
         if (ShowRealPart)
            std::cout << "#real                   ";
         if (ShowImagPart)
            std::cout << "#imag                   ";
         if (ShowCorrLength)
            std::cout << "#corr_length            ";
         if (ShowMagnitude)
            std::cout << "#magnitude              ";
         if (ShowArgument)
            std::cout << "#argument" << (ShowRadians ? "(rad)" : "(deg)") << "          ";
         std::cout << '\n';
         std::cout << std::left;
      }

      // Calculate the actual overlaps
      std::vector<TransEigenInfo> EigenList;
      for (std::set<QuantumNumber>::const_iterator I = Sectors.begin(); I != Sectors.end(); ++I)
      {
	 //FiniteMPO StringOp = FiniteMPO::make_identity(ExtractLocalBasis(Psi2->Psi));
         TransEigenInfo Info(*I, overlap(*Psi1, StringOp, *Psi2, *I, Iter, Tol, Verbose));
         if (Sort)
            EigenList.push_back(Info);
         else
            PrintFormat(Info, ShowRealPart, ShowImagPart, ShowCorrLength, ShowMagnitude,
			ShowArgument, ShowRadians, ScaleFactor);
      }

      if (Sort)
      {
         std::sort(EigenList.begin(), EigenList.end());
         for (unsigned i = 0; i < EigenList.size(); ++i)
         {
            PrintFormat(EigenList[i], ShowRealPart, ShowImagPart, ShowCorrLength, ShowMagnitude, ShowArgument,
			ShowRadians, ScaleFactor);
         }
      }


     pheap::Shutdown();

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
