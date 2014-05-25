
#include "mps/infinitewavefunction.h"
#include "mps/packunpack.h"
#include "lattice/latticesite.h"
#include "mps/operator_actions.h"
#include "common/environment.h"
#include "common/terminal.h"
#include "common/prog_options.h"
#include "common/prog_opt_accum.h"
#include "common/environment.h"
#include "interface/inittemp.h"
#include "tensor/tensor_eigen.h"
#include "mp-algorithms/arnoldi.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-u1u1.h"
#include "models/spin-z2.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/kondo-u1u1.h"
#include "models/kondo-so4.h"
#include "models/boson-u1.h"
#include "models/boson-2component-u1z2.h"
#include "models/hubbard-u1u1-old.h"
#include "models/hubbard-u1u1.h"
#include "models/hubbard-u1su2.h"
#include "models/hubbard-so4.h"
#include "models/hubbard.h"

#include "linearalgebra/arpack_wrapper.h"

namespace prog_opt = boost::program_options;

struct LeftMultiplyString
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiplyString(LinearWavefunction const& L1_, LinearWavefunction const& L2_, 
                      QuantumNumber const& QShift_,
                      std::vector<SimpleOperator> const& StringOp_) 
      : L1(L1_), L2(L2_), QShift(QShift_), StringOp(StringOp_) 
   {
      CHECK_EQUAL(L1.size(), L2.size())
         ("The two wavefunctions must be the same length");
      CHECK_EQUAL(L1.size(), StringOp.size())
         ("The string operator must be the same length as the unit cell");
   }

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(r, QShift);;
      std::vector<SimpleOperator>::const_iterator OpI = StringOp.begin();
      LinearWavefunction::const_iterator I1 = L1.begin();
      for (LinearWavefunction::const_iterator I2 = L2.begin(); I2 != L2.end(); ++I1, ++I2, ++OpI)
      {
	 r = operator_prod(herm(*OpI), herm(*I1), r, *I2);
      }
      return r;
   }

   LinearWavefunction const& L1;
   LinearWavefunction const& L2;
   QuantumNumber QShift;
   std::vector<SimpleOperator> StringOp;
};
   

struct MultFuncString
{
   MultFuncString(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2,
                  QuantumNumber const& QShift, 
                  QuantumNumbers::QuantumNumber const& q,
                  std::vector<SimpleOperator> const& StringOp)
      : Mult(Psi1, Psi2, QShift, StringOp), Pack(Psi1.Basis2(), Psi2.Basis2(), q) { }

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   LeftMultiplyString Mult;
   PackMatrixOperator Pack;
};


std::pair<std::complex<double>, MatrixOperator>
get_left_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, 
                     QuantumNumber const& QShift, 
                     std::vector<SimpleOperator> StringOp,
                     double tol = 1E-14, int Verbose = 0)
{
   int ncv = 0;
   QuantumNumber Ident(Psi1.GetSymmetryList());
   PackMatrixOperator Pack(Psi1.Basis1(), Psi2.Basis1(), Ident);
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > LeftEigen = 
         LinearAlgebra::DiagonalizeARPACK(MultFuncString(Psi1, Psi2, QShift, Ident, StringOp),
                                          n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   MatrixOperator LeftVector = Pack.unpack(&(OutVec[0]));
      
   return std::make_pair(LeftEigen[0], LeftVector);
}

struct RightMultiplyString
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   RightMultiplyString(LinearWavefunction const& L1_, LinearWavefunction const& L2_, 
                       QuantumNumber const& QShift_,
                       std::vector<SimpleOperator> const& StringOp_) 
      : L1(L1_), L2(L2_), QShift(QShift_) , StringOp(StringOp_)
   {
      CHECK_EQUAL(L1.size(), StringOp.size())("The string operator must be the same length as the unit cell");
   }

   result_type operator()(argument_type const& x) const
   {
      result_type r = delta_shift(x, adjoint(QShift)); 
      LinearWavefunction::const_iterator I1 = L1.end();
      LinearWavefunction::const_iterator I2 = L2.end();
      std::vector<SimpleOperator>::const_iterator OpI = StringOp.end();
      while (I1 != L1.begin())
      {
         --I1;
         --I2;
         --OpI;
	 r = operator_prod(*OpI, *I1, r, herm(*I2));
      }
      return r;
   }

   LinearWavefunction const& L1;
   LinearWavefunction const& L2;
   QuantumNumber QShift;
   std::vector<SimpleOperator> StringOp;
};

struct MultFuncStringTrans
{
   MultFuncStringTrans(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, 
                       QuantumNumber const& QShift, 
                       QuantumNumbers::QuantumNumber const& q,
                       std::vector<SimpleOperator> const& StringOp)
      : Mult(Psi1, Psi2, QShift, StringOp), Pack(Psi1.Basis1(), Psi2.Basis1(), q) { }

   void operator()(std::complex<double> const* In, std::complex<double>* Out) const
   {
      MatrixOperator x = Pack.unpack(In);
      x = Mult(x);
      Pack.pack(x, Out);
   }

   RightMultiplyString Mult;
   PackMatrixOperator Pack;
};



std::pair<std::complex<double>, MatrixOperator>
get_right_eigenvector(LinearWavefunction const& Psi1, LinearWavefunction const& Psi2, 
                      QuantumNumber const& QShift, 
                      std::vector<SimpleOperator> StringOp,
                      double tol = 1E-14, int Verbose = 0)
{
   int ncv = 0;
   QuantumNumber Ident(Psi1.GetSymmetryList());
   PackMatrixOperator Pack(Psi1.Basis1(), Psi2.Basis1(), Ident);
   int n = Pack.size();
   //   double tolsave = tol;
   //   int ncvsave = ncv;
   int NumEigen = 1;

   std::vector<std::complex<double> > OutVec;
      LinearAlgebra::Vector<std::complex<double> > RightEigen = 
         LinearAlgebra::DiagonalizeARPACK(MultFuncStringTrans(Psi1, Psi2, QShift, Ident, StringOp),
                                          n, NumEigen, tol, &OutVec, ncv, false, Verbose);

   MatrixOperator RightVector = Pack.unpack(&(OutVec[0]));
      
   return std::make_pair(RightEigen[0], RightVector);
}

int main(int argc, char** argv)
{
   try
   {
      std::string PsiStr;
      std::string OpL, OpR;
      //      double Tol = 1e-10;
      int Verbose = 0;
      std::string Model;
      double Spin = 0.5;
      int NMax = 3;
      std::vector<std::string> OperatorStr;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("wavefunction,w", prog_opt::value(&PsiStr), "Wavefunction [required]")
	 ("model", prog_opt::value(&Model), "use this model for the operators [required]")
	 ("spin", prog_opt::value(&Spin), "for spin models, the value of the spin [default 0.5]")
         ("nmax", prog_opt::value(&NMax), "for Bose-Hubbard model, the max number of bosons per site")
	 ("verbose,v", prog_opt_ext::accum_value(&Verbose), "increase verbosity")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("operator", prog_opt::value(&OperatorStr), "operator")
         ;

      prog_opt::positional_options_description p;
      p.add("operator", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("wavefunction") < 1)
      {
         std::cerr << "usage: mp-isymmetry [options] Operator1 [Operator2] ...\n";
         std::cerr << "-w (--wavefunction) and --model are required options.\n";
         std::cerr << desc << '\n';
         return 1;
      }

      // if the number of eigenvalues is specified but
      // the cutoff is not, then set the cutoff to zero
      //      if (vm.count("num-eigenvalues") == 1 && vm.count("eigen-cutoff") == 0)
      //      {
      //         EigenCutoff = 0;
      //      }

      SimpleOperator MyOpL, MyOpR, MyOpL2;
      LatticeSite Site;
      if (Model == "spin")
      {
	 Site = CreateSpinSite(Spin);
      }
      else if (Model == "spin-su2")
      {
	 Site = CreateSU2SpinSite(Spin);
      }
      else if (Model == "spin-u1")
      {
	 Site = CreateU1SpinSite(Spin);
      }
      else if (Model == "uls")
      {
	 Site = CreateU1U1SpinSite();
      }
      else if (Model == "spin-z2")
      {
	 Site = CreateZ2SpinSite(Spin);
      }
      else if (Model == "tj-u1")
      {
	 Site = CreateU1tJSite();
      }
      else if (Model == "sf-u1")
      {
	 Site = CreateU1SpinlessFermion();
      }
      else if (Model == "bh-u1")
      {
	 Site = BosonU1(NMax);
      }
      else if (Model == "bh2-u1z2")
      {
	 Site = CreateBoseHubbard2BosonsU1Z2Site(NMax);
      }
      else if (Model == "hubbard")
      {
         Site = CreateHubbardSite();
      }
      else if (Model == "hubbard-u1u1")
      {
         Site = CreateU1U1HubbardSite();
      }
      else if (Model == "hubbard-u1u1-old")
      {
         Site = CreateU1U1HubbardOldOrderingSite();
      }
      else if (Model == "hubbard-u1su2")
      {
         Site = CreateU1SU2HubbardSite();
      }
      else if (Model == "hubbard-so4")
      {
         Site = CreateSO4HubbardSiteA();
      }
      else if (Model == "klm-u1")
      {
         Site = CreateU1KondoSite();
      }
      else if (Model == "klm-u1u1")
      {
         Site = CreateU1U1KondoSite();
      }
      else if (Model == "klm-so4")
      {
         Site = CreateSO4KondoSiteA();
      }
      else if (Model != "")
      {
	 PANIC("Unknown model");
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      pvalue_ptr<InfiniteWavefunction> Psi1 
         = pheap::OpenPersistent(PsiStr, mp_pheap::CacheSize(), true);

      // The wavefunction
      LinearWavefunction Psi = get_orthogonal_wavefunction(*Psi1);
      double Dim = Psi.Basis1().total_degree();

      // The list of U matrices, for each operator
      std::vector<MatrixOperator> U;

      for (unsigned i = 0; i < OperatorStr.size(); ++i)
      {
         std::vector<SimpleOperator> StringOperator(Psi.size(), Site[OperatorStr[i]]);

         std::complex<double> e;
         MatrixOperator v;
         boost::tie(e, v) = get_left_eigenvector(Psi, Psi, Psi1->shift(), StringOperator);

         std::cout << "Eigenvalue " << OperatorStr[i] << " is " << e << std::endl;

         v *= std::sqrt(Dim); // make it properly unitary

         std::cout << "Trace " << (trace(v) / Dim) << std::endl;
         U.push_back(v);
      }

      // Now go through each operator pair
      std::cout << "#Op1       #Op2       #Commutator-Real  #Commutator-Imag  #Trace-Real    #Trace-Imag\n";
      for (unsigned i = 0; i < U.size(); ++i)
      {
         for (unsigned j = i+1; j < U.size(); ++j)
         {
            std::complex<double> x = inner_prod(U[i]*U[j], U[j]*U[i]) / Dim;
            std::complex<double> tr = trace(U[i]*U[j]) / Dim;
            std::cout << std::setw(10) << std::left << OperatorStr[i] << " "
                      << std::setw(10) << std::left << OperatorStr[j] << " "
                      << std::setw(17) << std::right << std::fixed << x.real() << " "
                      << std::setw(17) << std::right << std::fixed << x.imag() << " "
                      << std::setw(17) << std::right << std::fixed << tr.real() << " "
                      << std::setw(17) << std::right << std::fixed << tr.imag() << std::endl;
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
