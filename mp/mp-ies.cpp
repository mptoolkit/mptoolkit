// -*- C++ -*- $Id$

#include "wavefunction/mpwavefunction.h"
#include "mps/packunpack.h"
#include "lattice/latticesite.h"
#include "wavefunction/operator_actions.h"
#include "mp/copyright.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"
#include "lattice/unitcell.h"
#include "lattice/infinite-parser.h"
#include "linearalgebra/arpack_wrapper.h"
#include <boost/algorithm/string/predicate.hpp>
#include <fstream>
#include "common/formatting.h"
#include <tuple>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

namespace prog_opt = boost::program_options;

// returns true if Name exists and is a regular file
bool FileExists(std::string const& Name)
{
   struct stat buf;
   return stat(Name.c_str(), &buf) != -1 && S_ISREG(buf.st_mode);
}

template <typename Func>
struct PackApplyFunc
{
   PackApplyFunc(PackStateComponent const& Pack_, Func f_) : Pack(Pack_), f(f_) {}

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      StateComponent x = Pack.unpack(In);
      x = f(x);
      Pack.pack(x, Out);
   } 
   PackStateComponent const& Pack;
   Func f;
};

template <typename Func>
PackApplyFunc<Func>
MakePackApplyFunc(PackStateComponent const& Pack_, Func f_)
{
   return PackApplyFunc<Func>(Pack_, f_);
}

std::tuple<std::complex<double>, int, StateComponent>
get_right_eigenvector(LinearWavefunction const& Psi1, QuantumNumber const& QShift1,
		      LinearWavefunction const& Psi2, QuantumNumber const& QShift2, 
		      ProductMPO const& StringOp,
		      double tol = 1E-14, int Verbose = 0)
{
   int ncv = 0;
   int Length = statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size());
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > LeftEigen = 
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
							    RightMultiplyOperator(Psi1, QShift1,
										 StringOp, 
										 Psi2, QShift2, Length, Verbose-1)),
					  n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));
      
   return std::make_tuple(LeftEigen[0], Length, LeftVector);
}


std::tuple<std::complex<double>, int, StateComponent>
get_left_eigenvector(LinearWavefunction const& Psi1, QuantumNumber const& QShift1,
		     LinearWavefunction const& Psi2, QuantumNumber const& QShift2, 
                     ProductMPO const& StringOp,
                     double tol = 1E-14, int Verbose = 0)
{
   int ncv = 0;
   int Length = statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size());
   PackStateComponent Pack(StringOp.Basis1(), Psi1.Basis2(), Psi2.Basis2());
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > LeftEigen = 
         LinearAlgebra::DiagonalizeARPACK(MakePackApplyFunc(Pack,
							    LeftMultiplyOperator(Psi1, QShift1,
										 StringOp, 
										 Psi2, QShift2, Length, Verbose-1)),
					  n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));
      
   return std::make_tuple(LeftEigen[0], Length, LeftVector);
}

// The eigenvalues of the operator, quantum number sector, magnitude and
// phase angle
struct EValueRecord
{
   QuantumNumber q;
   std::complex<double> Eigenvalue;
};

bool CompareMagnitude(EValueRecord const& x, EValueRecord const& y)
{
   return std::abs(x.Eigenvalue) > std::abs(y.Eigenvalue);
}

bool CompareQMagnitude(EValueRecord const& x, EValueRecord const& y)
{
   return x.q < y.q || (x.q == y.q && std::abs(x.Eigenvalue) > std::abs(y.Eigenvalue));
}

void PrintHeading(std::ostream& out, bool ShowRadians)
{
   out << "#Number #Sector     #Degen  #Angle" << (ShowRadians ? "(rad)" : "(deg)")
       <<  "          #Eigenvalue             #Energy\n";
}

void PrintFormat(std::ostream& out, int n, EValueRecord const& x, std::complex<double> Normalizer, bool ShowRadians)
{
   std::string SectorStr = boost::lexical_cast<std::string>(x.q);
   out << std::left << std::setw(7) << n << ' ';
   out << std::left << std::setw(11) << SectorStr << ' ';
   out << std::left << std::setw(6) << degree(x.q) << ' ';
   double Arg = std::arg(x.Eigenvalue/Normalizer);
   if (!ShowRadians)
      Arg *= 180.0 / math_const::pi;
   out << std::setw(20) << std::right << std::fixed << Arg << "  ";
   out << std::setw(20) << std::scientific << std::abs(x.Eigenvalue/Normalizer) << ' ';
   out << std::setw(20) << std::fixed << (-log(std::abs(x.Eigenvalue/Normalizer)));
   out << std::setfill(' ') << '\n';
}

int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      std::string OpStr;
      int Verbose = 0;
      double Tol = 1E-14;
      bool ShowRadians = false;
      bool Quiet = false;
      bool NoGaugeFix = false;
      double GaugeAngle = 0.0;
      bool SplitBySectors = false;
      std::string SplitOutputPrefix;
      

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("quiet,q", prog_opt::bool_switch(&Quiet),
	  "don't show column headings")
	 ("radians", prog_opt::bool_switch(&ShowRadians),
	  "display the angle in radians instead of degrees (also for the --gauge angle, if specified)")
	 ("nogauge", prog_opt::bool_switch(&NoGaugeFix),
	  "don't gauge fix the phase angles, leave it as random output from the eigensolver")
	 ("gauge", prog_opt::value(&GaugeAngle), "gauge fix the largest eigenvalue "
	  "to this phase angle")
	 ("sectors", prog_opt::bool_switch(&SplitBySectors), "sort output by quantum number sector")
	 ("sectors-prefix", prog_opt::value(&SplitOutputPrefix),
	  "Write the output to a separate file for each quantum number sector prefix.sector.dat")
	 ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("wavefunction", prog_opt::value<std::string>(&PsiStr), "psi")
         ("operator", prog_opt::value<std::string>(&OpStr), "operator")
         ;

      prog_opt::positional_options_description p;
      p.add("wavefunction", 1);
      p.add("operator", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1 || vm.count("operator") < 1)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <operator>\n";
         std::cerr << desc << '\n';
	 std::cerr << "This tool constructs the entanglement spectrum of a wavefunction,\n"
		   << "labelling states by the phase angle of a unitary symmetry operator.\n"
		   << "The operator should be a global symmetry of the wavefunction.\n"
		   << "The phase angles are only determined up to a global phase (which\n"
		   << "arises from the global phase of the eigenmatrix of the symmetry\n"
		   << "operator).  By default, the phase is 'gauge-fixed' by setting the\n"
		   << "phase of the largest eigenvalue of the density matrix to be 0.\n"
		   << "The absolute phase will generally need to be corrected in\n"
		   << "post-processing, but the absolute phase of the largest eigenvalue\n"
		   << "can also be set (if known in advance) using the\n"
		   << "--gauge=<angle> parameter.\n"
		   << "\nBy default the entanglement spectrum is printed to standard output,\n"
		   << "ordered by eigenvalue.  Alternatively, the spectrum can be displayed\n"
		   << "for each quantum number sector separately, with the --sectors option.\n"
		   << "The --sectors-prefix=<fileprefix> option implies --sectors, and writes\n"
		   << "each sector to a file of the form fileprefix-sector.dat\n";
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      // if we have --gauge and --nogauge simultaneously then bug out before the dogfight starts
      if (vm.count("gauge") && NoGaugeFix)
      {
	 std::cerr << basename(argv[0]) << ": fatal: --gauge option conflicts with --nogauge\n";
	 exit(1);
      }

      // convert the gauge angle to radians if necessary
      if (!ShowRadians)
	 GaugeAngle *= math_const::pi / 180.0;

      if (vm.count("sectors-prefix"))
	 SplitBySectors = true;

      if (SplitBySectors)
      {
	 std::cerr << "not yet implemented!\n";
	 exit(1);
      }

      // Load the wavefunction
      pvalue_ptr<MPWavefunction> Psi
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      InfiniteWavefunctionLeft InfPsi = Psi->get<InfiniteWavefunctionLeft>();

      // Load the lattice, if it was specified
      //      pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);

      // orthogonalize the wavefunction
      LinearWavefunction Psi1;
      RealDiagonalOperator D;
      std::tie(Psi1, D) = get_left_canonical(InfPsi);
      MatrixOperator Rho = D;
      Rho = scalar_prod(Rho, herm(Rho));
      MatrixOperator Identity = MatrixOperator::make_identity(Psi1.Basis1());
      double Dim = Psi1.Basis1().total_degree();

#if 0
      UnitCell Cell = Lattice->GetUnitCell();
      int UnitCellSize = Cell.size();
      int const NumUnitCells = Psi1.size() / UnitCellSize;
#endif


      //      for Operators  // only one operator supported at this time
      {
	 ProductMPO StringOperator = ParseProductOperatorAndLattice(OpStr).first;

	 if (Psi1.size() % StringOperator.size() != 0)
	 {
	    std::cerr << "mp-ies: error: string operator size is incompatible with the wavefunction size for operator "
		      << OpStr << ", ignoring this operator.\n";
	    exit(1);
	 }

	 StringOperator = repeat(StringOperator, Psi1.size() / StringOperator.size());

	 // Get the matrix
         std::complex<double> e;
	 int n;
         StateComponent v;
         std::tie(e, n, v) = get_right_eigenvector(Psi1, InfPsi.qshift(), Psi1, InfPsi.qshift(), StringOperator, Tol, Verbose);

	 // Now we want the eigenvalues of v
	 // These will be of the form a * exp(i*theta) where a is the density matrix eigenvalue and
	 // theta is the phase (up to an overall phase ambiguity, since it is a constant)

	 std::vector<EValueRecord> EValues;

	 for (std::size_t i = 0; i < v.Basis1().size(); ++i)
	 {
	    if (iterate_at(v[0].data(), i, i))
	    {
	       LinearAlgebra::Vector<std::complex<double>> Eig = 
		  LinearAlgebra::EigenvaluesComplex(v[0](i,i));

	       for (unsigned j = 0; j < size(Eig); ++j)
	       {
		  EValues.push_back({v.Basis1()[i], Eig[j]});
	       }
	    }
	 }

	 // get the sum for normalization, and also the largest eigenvalue
	 double Sum = 0;
	 std::complex<double> LargestEValue = EValues[0].Eigenvalue;
	 for (auto const& e : EValues)
	 {
	    Sum += std::abs(e.Eigenvalue) * degree(e.q);
	    if (std::abs(e.Eigenvalue) > std::abs(LargestEValue))
	       LargestEValue = e.Eigenvalue;
	 }

	 // sort
	 if (SplitBySectors)
	    std::sort(EValues.begin(), EValues.end(), &CompareQMagnitude);
	 else
	    std::sort(EValues.begin(), EValues.end(), &CompareMagnitude);

	 // normalizer, which is a combination of the magnitude and phase correction
	 std::complex<double> Normalizer = Sum;
	 if (!NoGaugeFix)
	    Normalizer = Normalizer * (LargestEValue / 
				       (std::exp(std::complex<double>(0.0,GaugeAngle)) * std::abs(LargestEValue)));

	 // output
	 if (!Quiet)
	 {
	    print_preamble(std::cout, argc, argv);
	    std::cout << "#Eigenvalue of operator is " << format_complex(e) << '\n';
	    PrintHeading(std::cout, ShowRadians);
	 }
	 int k = 0; // eigenvalue number
	 for (auto const& e : EValues)
	 {
	    PrintFormat(std::cout, k++, e, Normalizer, ShowRadians);
	 }
      }

      pheap::Shutdown();
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << '\n';
      return 1;
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
