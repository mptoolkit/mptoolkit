
#include "linearalgebra/matrix_utility.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/scalarmatrix.h"
#include "mps/state_component.h"
#include "mps/infinitewavefunction.h"
#include "quantumnumbers/null-quantumnumber.h"
#include "pheap/pheap.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;


using namespace LinearAlgebra;

typedef Matrix<std::complex<double> > ComplexMatrix;

SymmetryList SL("N:Null");
QuantumNumber Q(SL);

inline
Matrix<std::complex<double> > random_unitary(size_type Size)
{
   Matrix<std::complex<double> > M = nrandom_matrix<std::complex<double> >(Size, Size);
   Matrix<std::complex<double> > R = M;
   Matrix<std::complex<double> > Q = QR_Factorize(R);
   Matrix<std::complex<double> > Lambda(Size, Size, 0.0);
   for (unsigned i = 0; i < Size; ++i)
   {
      Lambda(i,i) = R(i,i) / norm_frob(R(i,i));
   }
   return Q*Lambda;
}

inline
Matrix<std::complex<double> > random_unitary(size_type Size1, size_type Size2)
{
   int sz = std::max(Size1, Size2);
   Matrix<std::complex<double> > Mat = random_unitary(sz);
   return Mat(range(0,Size1), range(0,Size2));
}

// Generate a sequentially-generated MPS with a d-dimensional local basis,
// of matrix size m
StateComponent
RandomState(int d, int m)
{
   ComplexMatrix U = random_unitary(d*m, m);

   BasisList LocalBasis(SL);
   for (int i = 0; i < d; ++i)
      LocalBasis.push_back(Q);

   VectorBasis B(SL);
   B.push_back(Q, m);

   StateComponent Result(LocalBasis, B, B);
   for (int i = 0; i < d; ++i)
   {
      Result[i](0,0) = U(range(m*i, m*i+m), all);
   }

   return Result;
}

int main(int argc, char** argv)
{
   try
   {

      int d = 2;
      int m = 3;

      std::string FName;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help,h", "show this help message")
         ("out,o", prog_opt::value(&FName), "Output file (required)")
         ("seed,s", prog_opt::value<unsigned int>(), 
          ("Random seed [range 0.."+boost::lexical_cast<std::string>(RAND_MAX)+"]").c_str())
	 ("max-states,m", prog_opt::value(&m), FormatDefault("Number of states", m).c_str())
	 ("localdimension,d", prog_opt::value(&d), FormatDefault("Local dimension", d).c_str())
         ;

      prog_opt::options_description opt;
      opt.add(desc);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("out") == 0)
      {
         print_copyright(std::cerr);
         std::cerr << "usage: seqgen2 [options]\n";
         std::cerr << desc << '\n';
         return 1;
      }

      unsigned int RandSeed = vm.count("seed") ? (vm["seed"].as<unsigned long>() % RAND_MAX) : (ext::get_unique() % RAND_MAX);
      srand(RandSeed);

      pheap::Initialize(FName, 1, mp_pheap::PageSize(), mp_pheap::CacheSize());
   

      StateComponent A = RandomState(d,m);

      LinearWavefunction Lin(SL);
      Lin.push_back(A);

      InfiniteWavefunction Psi;
      Psi.QShift = Q;
      Psi.C_right = MatrixOperator::make_identity(A.Basis1());
      Psi.C_old = Psi.C_right;
      Psi.Psi = Lin;

      pvalue_ptr<InfiniteWavefunction> PsiPtr = new InfiniteWavefunction(Psi);
      pheap::ShutdownPersistent(PsiPtr);
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
