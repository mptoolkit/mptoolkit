
#include "mps/infinitewavefunction.h"
#include "mp/copyright.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "common/terminal.h"
#include "common/environment.h"
#include "common/prog_options.h"
#include "models/hubbard-u1su2.h"
#include "models/hubbard-u1u1.h"
#include "models/hubbard-u1u1-old.h"

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

void PrintFormat(TransEigenInfo const& x, bool ShowRealPart, bool ShowImagPart, bool ShowMagnitude)
{
   std::string SectorStr = boost::lexical_cast<std::string>(x.q);
   std::cout << std::setw(11) << SectorStr << ' ';
   if (ShowRealPart || ShowImagPart || ShowMagnitude)
   {
      if (ShowRealPart)
      {
         std::cout << std::setw(20) << x.x.real();
         if (ShowImagPart || ShowMagnitude)
            std::cout << "    ";
      }
      if (ShowImagPart)
      {
         std::cout << std::setw(20) << x.x.imag();
         if (ShowMagnitude)
            std::cout << "    ";
      }
      if (ShowMagnitude)
      {
         std::cout << std::setw(20) << LinearAlgebra::norm_frob(x.x);
      }
   }
   else // default to C++ complex output
      std::cout << std::setw(40) << x.x;
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
      bool Verbose = false;
      bool NoTempFile = false;
      bool ShowRealPart = false, ShowImagPart = false, ShowMagnitude = false;
      int Rotate = 0;
      std::string LhsStr, RhsStr;
      std::vector<std::string> Sector;
      double Tol = 1E-10;
      int Iter = 30;
      bool Sort = false;
      bool Quiet = false;
      bool Reflect = false;
      bool Conj = false;
      std::string Model;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("real,r", prog_opt::bool_switch(&ShowRealPart),
	  "display the real part of the result")
	 ("imag,i", prog_opt::bool_switch(&ShowImagPart),
	  "display the imaginary part of the result")
         ("mag,m", prog_opt::bool_switch(&ShowMagnitude),
          "display the magnitude of the result")
         ("notempfile", prog_opt::bool_switch(&NoTempFile),
          "don't use a temporary data file, keep everything in RAM "
          "(faster, but needs enough RAM)")
         ("rotate", prog_opt::value(&Rotate),
          "rotate the unit cell of psi2 this many sites to the left before calculating the overlap [default 0]")
         ("reflect", prog_opt::bool_switch(&Reflect),
          "reflect psi2 (gives parity eigenvalue)")
	 ("model", prog_opt::value(&Model),
	  "with the --reflect option, also apply the local reflection operator for this model")
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
         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "extra debug output")
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
         std::cerr << desc << '\n' 
                   << "If none of --real, --imag, --mag are specified, then the default output is a C++ formatted complex number.\n";
         return 1;
      }

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
#if 0
      while (Rotate > 0)
      {
         *Psi2.mutate() = rotate_left(*Psi2, 1);
         --Rotate;
      }
#else
      if (Rotate > 0)
         *Psi2.mutate() = rotate_left(*Psi2, Rotate);
#endif

      if (Reflect)
      {
         if (Verbose)
            std::cout << "Reflecting psi2..." << std::endl;
	 if (!Model.empty())
	 {
	    if (Verbose)
	       std::cout << "Reflecting using model " << Model << std::endl;
	    std::vector<SimpleOperator> ReflectionOp;
	    if (Model == "hubbard-u1su2")
	    {
	       LatticeSite Site = CreateU1SU2HubbardSite();
	       ReflectionOp = std::vector<SimpleOperator>(Psi2->size(), Site["R"]);
	    }
	    else if (Model == "hubbard-u1u1")
	    {
	       LatticeSite Site = CreateU1U1HubbardSite();
	       ReflectionOp = std::vector<SimpleOperator>(Psi2->size(), Site["R"]);
	    }
	    else if (Model == "hubbard-u1u1-old")
	    {
	       LatticeSite Site = CreateU1U1HubbardOldOrderingSite();
	       ReflectionOp = std::vector<SimpleOperator>(Psi2->size(), Site["R"]);
	    }
	    else
	    {
	       PANIC("Unknown model")(Model);
	    }
	    *Psi2.mutate() = reflect(*Psi2, ReflectionOp);
	 }
	 else
	 {
	    *Psi2.mutate() = reflect(*Psi2);
	 }
      }

      if (Conj)
      {
         if (Verbose)
            std::cout << "Conjugating psi2..." << std::endl;
         *Psi2.mutate() = conj(*Psi2);
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
         std::cout << "#sector     ";
         if (!ShowRealPart && !ShowImagPart && !ShowMagnitude)
            std::cout << "#value";
         if (ShowRealPart)
            std::cout << "#real                   ";
         if (ShowImagPart)
            std::cout << "#imag                   ";
         if (ShowMagnitude)
            std::cout << "#magnitude              ";
         std::cout << '\n';
         std::cout << std::left;
      }

      // Calculate the actual overlaps
      std::vector<TransEigenInfo> EigenList;
      for (std::set<QuantumNumber>::const_iterator I = Sectors.begin(); I != Sectors.end(); ++I)
      {
         TransEigenInfo Info(*I, overlap(*Psi1, *Psi2, *I, Iter, Tol, Verbose));
         if (Sort)
            EigenList.push_back(Info);
         else
            PrintFormat(Info, ShowRealPart, ShowImagPart, ShowMagnitude);
      }

      if (Sort)
      {
         std::sort(EigenList.begin(), EigenList.end());
         for (unsigned i = 0; i < EigenList.size(); ++i)
         {
            PrintFormat(EigenList[i], ShowRealPart, ShowImagPart, ShowMagnitude);
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
