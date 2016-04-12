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

// Write a matrix as sparse format.
// Force option allows forcing the last element to appear, evn if it is zero, 
// which is needed in the Matlab sparse format because it doesn't list separately the dimensions of the matrix.

void
WriteMatrixAsSparse(std::ostream& out, LinearAlgebra::Matrix<std::complex<double> > const& Op, 
		    double Epsilon = 0.0,
		    int iOffset = 0, int jOffset = 0, bool ForceLast = false,
		    bool Polar = false, bool Radians = false)
{
   for (unsigned i = 0; i < size1(Op); ++i)
   {
      for (unsigned j = 0; j < size2(Op); ++j)
      {
	 if (LinearAlgebra::norm_frob(Op(i,j)) > Epsilon || (ForceLast && i == size1(Op)-1 && j == size2(Op)-1))
	 {
	    out << (i+iOffset) << ' ' << (j+jOffset) << ' ';
	    if (Polar)
	    {
	       double Arg = std::arg(Op(i,j));
	       if (!Radians)
		  Arg *= 180.0 / math_const::pi;
	       out << std::abs(Op(i,j)) << ' ' << Arg << '\n';
	    }
	    else
	    {
	       out << Op(i,j).real() << ' ' << Op(i,j).imag() << '\n';
	    }
	 }
      }
   }
}

void
WriteMatrixOperatorAsSparse(std::ostream& out, MatrixOperator const& Op, double Epsilon = 0.0,
			    bool Polar = false, bool Radians = false)
{
   // We map the basis into a linear index.  So we need to get the offset of each subspace
   std::vector<int> Basis1Offset, Basis2Offset;
   Basis1Offset.push_back(0);
   for (unsigned i = 0; i < Op.Basis1().size(); ++i)
   {
      Basis1Offset.push_back(Basis1Offset.back() + Op.Basis1().dim(i));
   }
   CHECK_EQUAL(Basis1Offset.back(), Op.Basis1().total_dimension());
   Basis2Offset.push_back(0);
   for (unsigned i = 0; i < Op.Basis2().size(); ++i)
   {
      Basis2Offset.push_back(Basis2Offset.back() + Op.Basis2().dim(i));
   }
   CHECK_EQUAL(Basis2Offset.back(), Op.Basis2().total_dimension());

   for (LinearAlgebra::const_iterator<MatrixOperator>::type I = iterate(Op); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<MatrixOperator>::type J = iterate(I); J; ++J)
      {
	 // Offset + 1 to use 1-based addressing
	 WriteMatrixAsSparse(out, *J, Epsilon, Basis1Offset[J.index1()]+1, Basis2Offset[J.index2()]+1,
			     J.index1() == Op.Basis1().size()-1 && J.index2() == Op.Basis2().size()-1,
			     Polar, Radians);
      }
   }

   // make sure that the final element is included, if the bottom-right component of the matrix doesn't exist
   if (!iterate_at(Op.data(), Op.Basis1().size()-1, Op.Basis2().size()-1))
   {
      out << Op.Basis1().total_dimension() << ' ' << Op.Basis2().total_dimension() << " 0.0 0.0\n";
   }
}

void
WriteMatrixOperatorAsSparseStates(std::ostream& out, MatrixOperator const& Op, 
				  std::vector<std::pair<int,int>> const& WhichStates,
				  double Epsilon = 0.0,
				  bool Polar = false, bool Radians = false)
{
   for (unsigned i = 0; i < WhichStates.size(); ++i)
   {
      for (unsigned j = 0; j < WhichStates.size(); ++j)
      {
	 std::complex<double> x(0.0, 0.0);
	 LinearAlgebra::const_inner_iterator<MatrixOperator>::type J 
	    = iterate_at(Op.data(), WhichStates[i].first, WhichStates[j].first);
	 if (J)
	 {
	    x = (*J)(WhichStates[i].second, WhichStates[j].second);
	 }
	 if (LinearAlgebra::norm_frob(x) > Epsilon || (i == WhichStates.size()-1 && j == WhichStates.size()-1))
	 {
	    out << (i+1) << ' ' << (j+1) << ' ';
	    if (Polar)
	    {
	       double Arg = std::arg(x);
	       if (!Radians)
		  Arg *= 180.0 / math_const::pi;
	       out << std::abs(x) << ' ' << Arg << '\n';
	    }
	    else
	    {
	       out << x.real() << ' ' << x.imag() << '\n';
	    }
	 }
      }
   }
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
										 Psi2, QShift2, Length, Verbose-2)),
					  n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   StateComponent LeftVector = Pack.unpack(&(OutVec[0]));
      
   return std::make_tuple(LeftEigen[0], Length, LeftVector);
}



int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      std::string LatticeFile;
      int Verbose = 0;
      double Tol = 1e-14;
      std::vector<std::string> ProductOperators;
      std::vector<std::string> TriangularOperators;
      std::vector<std::string> FiniteOperators;
      std::string Partition;
      std::string RhoFile;
      bool Force = false;
      double Epsilon = 1E-10;
      std::string QSector;
      bool Polar = false;
      bool Radians = false;
      std::string WhichEigenvalues;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("wavefunction,w", prog_opt::value(&PsiStr), "Wavefunction [required]")
	 ("lattice,l", prog_opt::value(&LatticeFile), "Lattice file [required]")
	 ("product,p", prog_opt::value(&ProductOperators), 
	  "Construct matrix elements of a product operator -p file=operator")
	 ("triangular,t", prog_opt::value(&TriangularOperators),
	  "Construct atrix elements of a triangular operator -t file=operator (not yet implemented)")
	 ("finite,i", prog_opt::value(&FiniteOperators),
	  "Construct matrix elements of a finite operator -f file=operator (not yet implemented)")
	 ("sector,s", prog_opt::value(&QSector), "select a different quantum number sector [don't use this unless you know what you are doing]")
	 ("rho", prog_opt::value(&RhoFile), "Construct the density matrix --rho <filename>")
	 ("partition", prog_opt::value(&Partition), "Partition of the wavefunction cell,site (not yet implemented)")
	 ("restrict", prog_opt::value(&WhichEigenvalues), "Use only this set of eigenvalues of the density matrix [list of indices]")
	 ("force,f", prog_opt::bool_switch(&Force), "Overwrite files if they already exist")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
	 ("polar", prog_opt::bool_switch(&Polar), "Write the matrix elements in <magnitude> <argument> format")

	 ("radians", prog_opt::bool_switch(&Radians),
	  "print the argument in radians instead of degrees [implies --polar]")
         ("tol", prog_opt::value(&Tol),
          FormatDefault("Tolerance of the Arnoldi eigensolver", Tol).c_str())
	 ("epsilon,e", prog_opt::value(&Epsilon), FormatDefault("ignore matrix elements smaller than this epsilon", Epsilon).c_str())
	 ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description hidden("Hidden options");

      prog_opt::positional_options_description p;
      p.add("operator", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1 || vm.count("lattice") < 1)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] -w <psi> -l <lattice> -p <file>=<operator> [...] \n";
	 std::cerr << "This tool constructs the auxiliary space represenation of lattice operators,\n"
		   << "which are written to a file in a sparse-matrix format suitable for use in eg MATLAB.\n";
         std::cerr << desc << '\n';
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (Radians)
	 Polar = true;

      // Split the file/operator combinations into separate lists
      std::vector<std::string> ProductOpFile;
      std::vector<std::string> ProductOpStr;
      for (unsigned i = 0; i < ProductOperators.size(); ++i)
      {
	 std::string::iterator I = std::find(ProductOperators[i].begin(), ProductOperators[i].end(), '=');
	 if (I == ProductOperators[i].end())
	 {
	    std::cerr << "mp-aux-matrix: error: argument to --product is not of the form file=operator\n";
	    exit(1);
	 }
	 ProductOpFile.push_back(boost::trim_copy(std::string(ProductOperators[i].begin(), I)));
	 ProductOpStr.push_back(boost::trim_copy(std::string(I+1, ProductOperators[i].end())));
      }
      std::vector<std::string> TriangularOpFile;
      std::vector<std::string> TriangularOpStr;
      for (unsigned i = 0; i < TriangularOperators.size(); ++i)
      {
	 std::string::iterator I = std::find(TriangularOperators[i].begin(), TriangularOperators[i].end(), '=');
	 if (I == TriangularOperators[i].end())
	 {
	    std::cerr << "mp-aux-matrix: error: argument to --triangular is not of the form file=operator\n";
	    exit(1);
	 }
	 TriangularOpFile.push_back(boost::trim_copy(std::string(TriangularOperators[i].begin(), I)));
	 TriangularOpStr.push_back(boost::trim_copy(std::string(I+1, TriangularOperators[i].end())));
      }
      std::vector<std::string> FiniteOpFile;
      std::vector<std::string> FiniteOpStr;
      for (unsigned i = 0; i < FiniteOperators.size(); ++i)
      {
	 std::string::iterator I = std::find(FiniteOperators[i].begin(), FiniteOperators[i].end(), '=');
	 if (I == FiniteOperators[i].end())
	 {
	    std::cerr << "mp-aux-matrix: error: argument to --finite is not of the form file=operator\n";
	    exit(1);
	 }
	 FiniteOpFile.push_back(boost::trim_copy(std::string(FiniteOperators[i].begin(), I)));
	 FiniteOpStr.push_back(boost::trim_copy(std::string(I+1, FiniteOperators[i].end())));
      }

      // If the -f option hasn't been supplied, make sure that the output files don't already exist
      // If any files do exist, then exit immediately, so that we don't end up with a mixture of
      // new files with possibly old files
      if (!Force)
      {
	 for (unsigned i = 0; i < ProductOpFile.size(); ++i)
	 {
	    if (FileExists(ProductOpFile[i]))
	    {
	       std::cerr << "mp-aux-matrix: fatal: output file " << ProductOpFile[i]
			 << " already exists.  Use -f option to overwrite.\n";
	       exit(1);
	    }
	 }
	 for (unsigned i = 0; i < TriangularOpFile.size(); ++i)
	 {
	    if (FileExists(TriangularOpFile[i]))
	    {
	       std::cerr << "mp-aux-matrix: fatal: output file " << TriangularOpFile[i]
			 << " already exists.  Use -f option to overwrite.\n";
	       exit(1);
	    }
	 }
	 for (unsigned i = 0; i < FiniteOpFile.size(); ++i)
	 {
	    if (FileExists(FiniteOpFile[i]))
	    {
	       std::cerr << "mp-aux-matrix: fatal: output file " << FiniteOpFile[i]
			 << " already exists.  Use -f option to overwrite.\n";
	       exit(1);
	    }
	 }
	 if (!RhoFile.empty() && FileExists(RhoFile))
	 {
	    std::cerr << "mp-aux-matrix: fatal: output file " << RhoFile
		      << " already exists.  Use -f option to overwrite.\n";
	       exit(1);
	 }
      }

      // Load the wavefunction
      pvalue_ptr<MPWavefunction> Psi
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      InfiniteWavefunctionLeft InfPsi = Psi->get<InfiniteWavefunctionLeft>();

      // Load the lattice, if it was specified
      pvalue_ptr<InfiniteLattice> Lattice = pheap::ImportHeap(LatticeFile);

      // orthogonalize the wavefunction
      LinearWavefunction Psi1;
      RealDiagonalOperator D;
      std::tie(Psi1, D) = get_left_canonical(InfPsi);
      MatrixOperator Rho = D;
      Rho = scalar_prod(Rho, herm(Rho));
      MatrixOperator Identity = MatrixOperator::make_identity(Psi1.Basis1());
      double Dim = Psi1.Basis1().total_degree();

      UnitCell Cell = Lattice->GetUnitCell();
      int UnitCellSize = Cell.size();
      int const NumUnitCells = Psi1.size() / UnitCellSize;

      // Get the list of states in order.
      // If we didn't specify a list with --restrict then take all of them
      
      std::vector<std::pair<int, int>> WhichStates;
 
      // make a new scope since we have some temporary values
      {
	 // Sort the eigenvalues with an index sort
	 std::vector<std::tuple<int,int,double>> EValues;
	 for (int i = 0; i < int(Rho.Basis1().size()); ++i)
	 {
	    for (int j = 0; j < Rho.Basis1().dim(i); ++j)
	    {
	       EValues.push_back(std::make_tuple(i,j,Rho.Basis1()[i].degree()*LinearAlgebra::norm_frob(Rho(i,i)(j,j))));
	    }
	 }
	 std::sort(EValues.begin(), EValues.end(), [](std::tuple<int,int,double> const& x, std::tuple<int,int,double> const& y)
		   { return std::get<2>(x) < std::get<2>(y); });
	 
	 if (WhichEigenvalues.empty())
	 {
	    // take all of them
	    for (unsigned i = 0; i < EValues.size(); ++i)
	    {
	       WhichStates.push_back({std::get<0>(EValues[i]), std::get<1>(EValues[i])});
	    }
	 }
	 else
	 {
	    std::vector<std::string> Which;
	    boost::split(Which, WhichEigenvalues, boost::is_any_of(", \n\t"), boost::token_compress_on);
	    for (std::string s : Which)
	    {
	       int n = boost::lexical_cast<int>(s);
	       if (n < 1 || n > int(EValues.size()))
	       {
		  std::cerr << "mp-aux-matrix: density matrix index must be in range 1.." 
			    << EValues.size() << ", ignoring element " << n << '\n';
		  continue;
	       }
	       WhichStates.push_back({std::get<0>(EValues[n-1]), std::get<1>(EValues[n-1])});
	    }
	 }
      }

      // Now that we have Rho, save it if necessary
      if (!RhoFile.empty())
      {
	 std::ofstream Out(RhoFile.c_str(), std::ios::out | std::ios::trunc);
	 if (!Out.good())
	 {
	    std::cerr << "mp-aux-matrix: failed to open file " << RhoFile << " for density matrix.\n";
	 }
	 else
	 {
	    Out.precision(getenv_or_default("MP_PRECISION", 14));
	    if (WhichStates.empty())
	       WriteMatrixOperatorAsSparse(Out, Rho, Epsilon, Polar, Radians);
	    else
	       WriteMatrixOperatorAsSparseStates(Out, Rho, WhichStates, Epsilon, Polar, Radians);
	 }
      }

      // reflected and conjugated versions of the wavefunction - we leave them as null
      // until they are actually needed
      LinearWavefunction PsiR, PsiC, PsiRC;

      for (unsigned i = 0; i < ProductOpStr.size(); ++i)
      {
	 std::string OpStr = ProductOpStr[i];
	 std::string FileName = ProductOpFile[i];

	 if (Verbose > 0)
	 {
	    std::cout << "Constructing operator " << OpStr << '\n';
	    std::cout << "Writing to file " << FileName << std::endl;
	 }

	 bool Reflect = false;
	 bool Conjugate = false;
	 LinearWavefunction* Psi2 = &Psi1;
	 // Do we have time reversal or reflection?
	 if (boost::starts_with(OpStr, "r&"))
	 {
	    Reflect = true;
	    OpStr = std::string(OpStr.begin()+2, OpStr.end());
	    if (PsiR.empty())
	    {
	       InfiniteWavefunctionLeft PR = InfPsi;
	       inplace_reflect(PR);
	       // The wavefunction must be in a reflection-symmetric basis, or this isn't valid
	       if (PR.Basis1() != InfPsi.Basis1())
	       {
		  std::cerr << "mp-aux-matrix: cannot reflect operator because the basis is not reflection symmetric.\n"
			    << "mp-aux-matrix: ignoring operator " << OpStr << '\n';
		  continue;
	       }
	       PsiR = get_left_canonical(PR).first;
	    }
	    Psi2 = &PsiR;
	 }
	 else if (boost::starts_with(OpStr,"c&"))
	 {
	    Conjugate = true;
	    OpStr = std::string(OpStr.begin()+2, OpStr.end());
	    if (PsiC.empty())
	    {
	       PsiC = conj(Psi1);
	    }
	    Psi2 = &PsiC;
	 }
	 else if (boost::starts_with(OpStr,"rc&") || boost::starts_with(OpStr,"cr&"))
	 {
	    Reflect = true;
	    Conjugate = true;
	    OpStr = std::string(OpStr.begin()+3, OpStr.end());
	    if (PsiR.empty())
	    {
	       InfiniteWavefunctionLeft PR = InfPsi;
	       inplace_reflect(PR);
	       
	       // The wavefunction must be in a reflection-symmetric basis, or this isn't valid
	       if (PR.Basis1() != InfPsi.Basis1())
	       {
		  std::cerr << "mp-aux-matrix: cannot reflect operator because the basis is not reflection symmetric.\n"
			    << "mp-aux-matrix: ignoring operator " << OpStr << '\n';
		  continue;
	       }
	       PsiR = get_left_canonical(PR).first;
	    }
	    if (PsiRC.empty())
	    {
	       PsiRC = conj(PsiR);
	    }
	    Psi2 = &PsiRC;
	 }

	 ProductMPO StringOperator = ParseProductOperator(*Lattice, OpStr);

	 if (!QSector.empty())
	 {
	    QuantumNumber q(Psi1.GetSymmetryList(), QSector);
	    StringOperator = StringOperator * ProductMPO::make_identity(StringOperator.LocalBasis2List(), q);
	 }

	 // Get the matrix
         std::complex<double> e;
	 int n;
         StateComponent v;
         std::tie(e, n, v) = get_left_eigenvector(Psi1, InfPsi.qshift(), *Psi2, InfPsi.qshift(), StringOperator,
						  Tol, Verbose);

	 // Normalization
	 // it might not be unitary, eg anti-unitary.  So we need to take the 4th power
	 std::complex<double> x = inner_prod(scalar_prod(herm(v), operator_prod(herm(v), v, v)), Rho);
	 v *= std::sqrt(std::sqrt(1.0 / x));

	 MatrixOperator M = v[0];

	 // write to file
	 std::ofstream Out(FileName.c_str(), std::ios::out | std::ios::trunc);
	 if (!Out.good())
	 {
	    std::cerr << "mp-aux-matrix: failed to open file " << FileName << " for operator " << OpStr << '\n';
	 }
	 else 
	 {
	    Out.precision(getenv_or_default("MP_PRECISION", 14));
	    if (WhichStates.empty())
	       WriteMatrixOperatorAsSparse(Out, M, Epsilon, Polar, Radians);
	    else
	       WriteMatrixOperatorAsSparseStates(Out, M, WhichStates, Epsilon, Polar, Radians);
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
