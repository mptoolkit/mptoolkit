// -*- C++ -*- $Id$

#include "mpo/triangularoperator.h"
#include "mps/infinitewavefunction.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include "common/environment.h"
#include "interface/inittemp.h"
#include "mp-algorithms/gmres.h"
#include "mp-algorithms/arnoldi.h"

#include "tensor/tensor_eigen.h"

#include "models/spin.h"
#include "models/spin-u1.h"
#include "models/spin-su2.h"
#include "models/tj-u1su2.h"
#include "models/tj-u1.h"
#include "models/spinlessfermion-u1.h"
#include "models/kondo-u1su2.h"
#include "models/kondo-u1.h"
#include "models/hubbard-so4.h"

double const SolverTol = 1e-13;

template <typename CoefficientField>
class Polynomial
{
   public:
      typedef CoefficientField coefficient_type;
      typedef coefficient_type data_type;
      typedef int key_type;

      typedef std::map<int, coefficient_type> container_type;

      typedef typename container_type::iterator       iterator;
      typedef typename container_type::const_iterator const_iterator;
      typedef typename container_type::value_type value_type;

      Polynomial() {}

      // initialize as a constant polynomial
      explicit Polynomial(coefficient_type const& c);
      
      coefficient_type& operator[](int n) { return data_[n]; }
 
      // returns the coefficient of the term x^n, or a default constructed coefficient
      // if the term is zero
      coefficient_type coefficient(int n) const;

      bool empty() const { return data_.empty(); }
      int degree() const;

      // returns the number of non-zero coefficients
      int size() const { return int(data_.size()); }

      iterator begin() { return data_.begin(); }
      iterator end() { return data_.end(); }

      const_iterator begin() const { return data_.begin(); }
      const_iterator end() const { return data_.end(); }

      std::map<int, coefficient_type> data_;
};

template <typename CF>
Polynomial<CF>::Polynomial(coefficient_type const& c)
{
   data_[0] = c;
}

template <typename CF>
int Polynomial<CF>::degree() const
{
   if (data_.empty())
      return 0;
   // else
   const_iterator i = data_.end();
   --i;
   return i->first;
}

template <typename CF>
CF Polynomial<CF>::coefficient(int n) const
{
   const_iterator I = data_.find(n);
   if (I == data_.end())
      return coefficient_type();
   return I->second;
}

template <typename CF>
std::ostream&
operator<<(std::ostream& out, Polynomial<CF> const& x)
{
   if (x.empty())
   {
      return out << "zero polynomial\n";
   }

   for (typename Polynomial<CF>::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      out << "Coefficient of x^" << I->first << " = " << I->second << '\n';
   }
   return out;
}

long Binomial(int n, int k)
{
   if (k > n/2)
      k = n-k;     // take advantage of symmetry
   double r = 1.0;
   for (int i = 1; i <= k; ++i)
   {
      r *= double(n-k+i) / double(i);
   }
   return long(r+0.5); // round to nearest
}

struct OneMinusTransfer
{
   OneMinusTransfer(double ScaleFactor, StateComponent const& Psi, MatrixOperator const& Rho,
                    MatrixOperator const& Identity)
      : Scale_(ScaleFactor), Psi_(Psi), Rho_(Rho), Identity_(Identity) {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-Scale_*operator_prod(herm(Psi_), x, Psi_);
      //      MatrixOperator r = x-Scale_*0.5*(operator_prod(herm(Psi_), x, Psi_) 
      // + operator_prod(Psi_, x, herm(Psi_)));
      r -= inner_prod(r, Rho_) * Identity_; // orthogonalize to the identity
      return r;
   }

   double Scale_;
   StateComponent const& Psi_;
   MatrixOperator const& Rho_;
   MatrixOperator const& Identity_;
};

template <typename Func>
MatrixOperator
LinearSolve(Func F, MatrixOperator Rhs)
{
   MatrixOperator Guess = Rhs;
   int m = 30;
   int max_iter = 10000;
   double tol = SolverTol;
   GmRes(Guess, F, Rhs, m, max_iter, tol, LinearAlgebra::Identity<MatrixOperator>());
   return Guess;
}


Polynomial<MatrixOperator>
SolveMPO_Left(StateComponent const& Psi,
              MpOpTriangular const& Op, MatrixOperator const& Rho,
              MatrixOperator const& Identity, bool Verbose = false)
{
   typedef Polynomial<MatrixOperator> PolyType;
   typedef Polynomial<std::complex<double> > ComplexPolyType;
   int Dim = Op.Basis1().size();  // dimension of the MPO
   std::vector<PolyType> EMat(Dim);  // The E-matrices

   // Initialize the first E matrix
   EMat[Dim-1] = PolyType(Identity);

   // solve recursively
   int Col = Dim-2;
   while (Col >= 0)
   {
      if (Verbose)
      {
         std::cerr << "Solving column " << Col << '\n';
      }
      // Generate the next C matrices, C(n) = sum_{j>Col} Op(j,Col) E_j(n)
      PolyType C;
      for (int j = Col+1; j < Dim; ++j)
      {
         SimpleRedOperator const M = Op(j,Col);
         if (!M.is_null())
         {
            for (PolyType::const_iterator I = EMat[j].begin(); I != EMat[j].end(); ++I)
            {
               for (SimpleRedOperator::const_iterator MIter = M.begin(); 
                    MIter != M.end(); ++MIter)
               {
                  // TODO: add a variant for a SimpleRedOperator to state_component.h
		  if (!I->second.is_null())
		     C[I->first] += operator_prod(herm(*MIter), herm(Psi), I->second, Psi,
						  Op.Basis1()[Col]);
               }
            }
         }
      }

      // Now do the classification, based on the properties of the diagonal operator
      SimpleRedOperator Diag = Op(Col, Col);
      if (norm_frob(scalar(Diag)) < 1E-10)
      {
         // the operator is zero.  In this case the solution is simply
         // E(n) = C(n-1)
         DEBUG_TRACE("Zero diagonal element")(Col)(Diag);
         PolyType E;
         int MaxDegree = C.degree();
         for (int i = MaxDegree; i >= 0; --i)
         {
            E[i] = C[i];
            for (int j = i+1; j <= MaxDegree; ++j)
            {
               E[i] -= double(Binomial(j,i)) * E[j];
            }
         }
         EMat[Col] = E;
      }
      else
      {
         // operator is not zero, assume it is proportional to the identity
         DEBUG_TRACE("Non-zero diagonal element")(Col)(Diag);

         ComplexPolyType EParallel;  // components parallel to the identity, may be zero
         double Fac = norm_frob_sq(scalar(Diag));
         Fac = std::sqrt(Fac / Diag.Basis1().total_degree());
         // is Fac == 1 ?
         if (fabs(Fac - 1.0) < 1E-10)
         {
            // diagonal element is the identity
            DEBUG_TRACE("Unit diagonal element")(Col)(Diag);
         
            // decompose C into components parallel and perpendicular to the identity
            ComplexPolyType CParallel;
            for (PolyType::iterator I = C.begin(); I != C.end(); ++I)
            {
               std::complex<double> Overlap = inner_prod(I->second, Rho);
               I->second -= Overlap*Identity;
               if (norm_frob(Overlap) > 1E-16)
                  CParallel[I->first] = Overlap;
            }

            // components parallel to the identity satisfy equation (23) of the notes
            for (int m = CParallel.degree(); m >= 0; --m)
            {
               EParallel[m+1] = CParallel[m];
               for (int k = m+2; k <= CParallel.degree()+1; ++k)
               {
                  EParallel[m+1] -= double(Binomial(k,m)) * EParallel[k];
               }
               EParallel[m+1] *= 1.0 / (1.0 + m);
            }
         }

         // if we get here, then Fac <= 1
         // do the components perpendicular to the identity
      
         // Components perpendicular to the identity satisfy equation (24)
         PolyType E;
         for (int m = C.degree(); m >= 0; --m)
         {
            MatrixOperator Rhs = C[m];
            for (int k = m+1; k <= C.degree(); ++k)
            {
               Rhs -= double(Binomial(k,m)) * E[k];
            }

            // orthogonalize Rhs against the identity again, which is a kind of
            // DGKS correction
            Rhs -= inner_prod(Rho, Rhs) * Identity;

            double RhsNorm2 = norm_frob_sq(Rhs);
            RhsNorm2 = RhsNorm2 / (Rhs.Basis1().total_degree()*Rhs.Basis2().total_degree());
            //TRACE(RhsNorm2);
            if (RhsNorm2 > 1E-22)
            {
               E[m] = LinearSolve(OneMinusTransfer(Fac, Psi, Rho, Identity), Rhs);
               // do another orthogonalization
               E[m] -= inner_prod(Rho, E[m]) * Identity;
            }
         }

         // Reinsert the components parallel to the identity
         for (ComplexPolyType::const_iterator I = EParallel.begin(); I != EParallel.end(); ++I)
         {
            E[I->first] += I->second * Identity;
         }

         // finished this column
         EMat[Col] = E;
      }

      --Col;
   }

   return EMat[0];
}

Polynomial<std::complex<double> >
ExtractOverlap(Polynomial<MatrixOperator> const& E, MatrixOperator const& Rho)
{
   Polynomial<std::complex<double> > Overlap;
   for (Polynomial<MatrixOperator>::const_iterator I = E.begin(); I != E.end(); ++I)
   {
      Overlap[I->first] = inner_prod(I->second, Rho);
   }
   return Overlap;
}

int main(int argc, char** argv)
{
   if (argc < 2 || argc > 2)
   {
      std::cout << "usage: mp-iexpectation <psi>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<InfiniteWavefunction> PsiPtr = pheap::OpenPersistent(argv[1], CacheSize, true);
   InfiniteWavefunction Psi = *PsiPtr;

#if 0
   // ITF
   SiteBlock Site = CreateSpinSite(0.5);
   double J = -1.0;
   double Lambda = 1.0;
   MpOpTriangular Ham = J * 4.0 * TriangularTwoSite(Site["Sz"], Site["Sz"])
      + Lambda * 2.0 * TriangularOneSite(Site["Sx"]);
   MpOpTriangular A = TriangularOneSite(Site["Sz"]);
   MpOpTriangular ADag = TriangularOneSite(Site["Sz"]);
#endif

#if 0
   // Heisenberg
   SiteBlock Site = CreateSpinSite(0.5);
   double J = -1.0;
   double Jz = 1;
   MpOpTriangular Ham = J * 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"])
                                   + TriangularTwoSite(Site["Sm"], Site["Sp"]))
      + Jz * TriangularTwoSite(Site["Sz"], Site["Sz"]);
   MpOpTriangular A = TriangularOneSite(Site["Sp"]);
   MpOpTriangular ADag = TriangularOneSite(Site["Sm"]);
#endif

#if 1
   // heisenberg with doubled unit cell, at momentum pi
   SiteBlock Site = CreateSpinSite(0.5);
   double J = -1.0;
   double Jz = 1;
   MpOpTriangular Ham = J * 0.5 * (TriangularTwoSite(Site["Sp"], Site["Sm"])
                                   + TriangularTwoSite(Site["Sm"], Site["Sp"]))
      + Jz * TriangularTwoSite(Site["Sz"], Site["Sz"]);
   MpOpTriangular A = TriangularOneSite(Site["Sp"]);
   MpOpTriangular ADag = TriangularOneSite(Site["Sm"]);

   MpOpTriangular Ident = 0.0 * TriangularOneSite(Site["I"]);

   Ham = local_tensor_prod(Ham, Ham);

   A = local_tensor_prod(A,Ident) - local_tensor_prod(Ident, A);
   ADag = local_tensor_prod(ADag,Ident) - local_tensor_prod(Ident, ADag);

#endif


   MatrixOperator Identity = MatrixOperator::make_identity(Psi.Psi.Basis1());
   MatrixOperator Rho = scalar_prod(Psi.C_right, herm(Psi.C_right));
   StateComponent Phi = Psi.Psi.get_front(); // no need to bugger around with C_old,C_right
 
   MatrixOperator LambdaSqrt = SqrtDiagonal(Psi.C_old);
   MatrixOperator LambdaInvSqrt = InvertDiagonal(LambdaSqrt, InverseTol);

   bool Verbose = false;

#if 1
   // double the size of the unit cell, by coarse graining
   Phi = local_tensor_prod(Phi, Phi);
#endif


   Phi = prod(prod(LambdaInvSqrt, Phi), LambdaSqrt);
   Rho = Psi.C_old;
   Identity = Rho;

   TRACE(norm_frob(operator_prod(herm(Phi), Identity, Phi) - Identity));
   TRACE(norm_frob(operator_prod(Phi, Rho, herm(Phi)) - Rho));

   std::cout.precision(14);

   std::vector<double> Delta(1,0.0);  // Delta[0] = 0
   std::vector<MpOpTriangular> f(1,A); // f[0] = A
   std::vector<MpOpTriangular> fDag(1,ADag); // f[0] = ADag
   std::vector<double> Inner;

   MpOpTriangular Op = fDag[0]*f[0] + f[0]*fDag[0];
   Polynomial<MatrixOperator> E = SolveMPO_Left(Phi, Op, Rho, Identity, Verbose);
   Polynomial<std::complex<double> > a = ExtractOverlap(E, Rho);
   TRACE(a);
   Inner.push_back(a[1].real());   // Inner[0] = (f[0], f[0])

   // k=1
   f.push_back(std::complex<double>(0.0,1.0)*(Ham*f[0] - f[0]*Ham));
   fDag.push_back(std::complex<double>(0.0,1.0)*(Ham*fDag[0] - fDag[0]*Ham));
   Op = fDag[1]*f[1] + f[1]*fDag[1];
   E = SolveMPO_Left(Phi, Op, Rho, Identity, Verbose);
   a = ExtractOverlap(E, Rho);
   TRACE(a);
   Inner.push_back(a[1].real());   // Inner[1] = (f[1], f[1])
   Delta.push_back(Inner[1] / Inner[0]);  // Delta[1] = (f[1],f[1])/(f[0],f[0])
   TRACE("k=1")(Delta.back());

   int k = 1;
   while (true)
   {
      // f[k+1] = i L f[k] + Delta[k]*f[k-1]
      f.push_back(std::complex<double>(0.0,1.0)*(Ham*f[k] - f[k]*Ham) + Delta[k]*f[k-1]);
      fDag.push_back(std::complex<double>(0.0,1.0)*(Ham*fDag[k] - fDag[k]*Ham) + Delta[k]*fDag[k-1]);
      ++k;
      Op = fDag[k]*f[k] + f[k]*fDag[k];
      E = SolveMPO_Left(Phi, Op, Rho, Identity, Verbose);
      a = ExtractOverlap(E, Rho);
      TRACE(a);
      Inner.push_back(a[1].real());  
      Delta.push_back(Inner[k] / Inner[k-1]); 
      TRACE(k)(Delta[k]);
   }

   pheap::Shutdown();
}
