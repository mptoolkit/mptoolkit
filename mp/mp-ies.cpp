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
   double Magnitude;
   double Radians;
};

bool CompareMagnitude(EValueRecord const& x, EValueRecord const& y)
{
   return x.Magnitude > y.Magnitude;
}

bool CompareQMagnitude(EValueRecord const& x, EValueRecord const& y)
{
   return x.q < y.q || (x.q == y.q && x.Magnitude > y.Magnitude);
}

int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      std::string OpStr;
      int Verbose = 0;
      std::string QSector;
      double Tol = 1E-14;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
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
	 std::cerr << "This tool constructs the entanglement spectrum, labelling states by the phase angle of"
		   << "a unitary symmetry operator.\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));


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

	 if (!QSector.empty())
	 {
	    QuantumNumber q(Psi1.GetSymmetryList(), QSector);
	    StringOperator = StringOperator * ProductMPO::make_identity(StringOperator.LocalBasis2List(), q);
	 }

	 // Get the matrix
         std::complex<double> e;
	 int n;
         StateComponent v;
         std::tie(e, n, v) = get_right_eigenvector(Psi1, InfPsi.qshift(), Psi1, InfPsi.qshift(), StringOperator, Tol, Verbose);

	 TRACE(e);

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
		  EValues.push_back({v.Basis1()[i], std::abs(Eig[j]), std::arg(Eig[j])});
	       }
	    }
	 }

	 // get the sum for normalization
	 double Sum = 0;
	 for (auto const& e : EValues)
	 {
	    Sum += e.Magnitude * degree(e.q);
	 }

	 // sort
	 std::sort(EValues.begin(), EValues.end(), &CompareMagnitude);

	 // output
	 for (auto const& e : EValues)
	 {
	    std::cout << (e.Magnitude/Sum) << ' ' << e.Radians << ' ' << e.q << '\n';
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
