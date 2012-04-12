// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include <boost/tuple/tuple.hpp>
#include "linearalgebra/matrixtransform.h"

typedef std::complex<double> complex;


struct GetPositivePart
{
   typedef complex argument_type;
   typedef complex result_type;
   result_type operator()(const complex& x) const
   {
      if (x.imag() != 0)
      {
         WARNING("GetPositivePart: ignoring non-zero imaginary part.")(x);
      }
      return x.real() > 0.0 ? x.real() : 0.0;
   }
};

struct GetNegativePart
{
   typedef complex argument_type;
   typedef complex result_type;
   result_type operator()(const complex& x) const
   {
      if (x.imag() != 0)
      {
         WARNING("GetNegativePart: ignoring non-zero imaginary part.")(x);
      }
      return x.real() < 0.0 ? x.real() : 0.0;
   }
};

boost::tuple<MatrixOperator, MatrixOperator>
SplitPositiveNegative(MatrixOperator const& M)
{
   MatrixOperator Pos, Neg;
   typedef MatrixOperator::value_type InnerMat;
   Pos = transform(M, LinearAlgebra::TransformMatrix<InnerMat, GetPositivePart>());
   Neg = transform(M, LinearAlgebra::TransformMatrix<InnerMat, GetNegativePart>());
   return boost::tie(Pos, Neg);
}

boost::tuple<MPStateComponent, MPStateComponent>
SplitPositiveNegative(MPStateComponent const& M)
{
   MPStateComponent Pos=M, Neg=M;
   typedef MatrixOperator::value_type InnerMat;
   for (unsigned i = 0; i < M.size(); ++i)
   {
      boost::tie(Pos[i], Neg[i]) = SplitPositiveNegative(M[i]);
   }
   return boost::tie(Pos, Neg);
}

void Split(LinearWavefunction& Psi, MatrixOperator& CPos, MatrixOperator& CNeg)
{
   MatrixOperator C = MatrixOperator::make_identity(Psi.Basis2());
   MatrixOperator One = MatrixOperator::make_identity(Psi.Basis1());
   MatrixOperator Zero = One; Zero *= 0.0;
   LinearWavefunction::iterator I = Psi.end();
   --I; // pointing to the last valid element
   SumBasis<VectorBasis> P;
   {
      MPStateComponent A = *I;
      SumBasis<VectorBasis> PNew(A.Basis1(), A.Basis1());
      MPStateComponent APos, ANeg;
      boost::tie(APos, ANeg) = SplitPositiveNegative(A);
      ANeg *= -1.0;
      *I = tensor_col_sum(APos, ANeg, PNew);
      P = PNew;
   }
   while (I != Psi.begin())
   {
      --I;
      MPStateComponent A = *I;
      TRACE(norm_frob_sq(A));
      MPStateComponent APos, ANeg;
      boost::tie(APos, ANeg) = SplitPositiveNegative(A);
      ANeg *= -1.0;
      TRACE(inner_prod(APos, ANeg));
      TRACE(norm_frob_sq(APos))(norm_frob_sq(ANeg));
      MPStateComponent ATopRow = tensor_row_sum(APos, ANeg, P);
      MPStateComponent ABotRow = tensor_row_sum(ANeg, APos, P);
      SumBasis<VectorBasis> PNew(A.Basis1(), A.Basis1());
      *I = tensor_col_sum(ATopRow, ABotRow, PNew);
      TRACE(norm_frob_sq(*I));
      //TRACE(scalar_prod(*I, herm(*I)));
      P = PNew;
   }

   CPos = tensor_row_sum(One, Zero, P);
   CNeg = tensor_row_sum(Zero, One, P);

   Psi = inject_right_old_interface(Psi, C);
   CPos = CPos * C;
   CNeg = CNeg * C;

   LinearWavefunction PsiPos = inject_left_old_interface(CPos, Psi);
   TRACE(norm_frob_sq(CPos))(CPos); 
   PsiPos = inject_right_old_interface(PsiPos, CPos);
   PsiPos.set_front(prod(CPos, PsiPos.get_front()));

   LinearWavefunction PsiNeg = inject_left_old_interface(CNeg, Psi);
   PsiNeg = inject_right_old_interface(PsiNeg, CNeg);
   PsiNeg.set_front(prod(CNeg, PsiNeg.get_front()));

   TRACE(norm_2_sq(PsiPos))(norm_2_sq(PsiNeg));
   TRACE(overlap(PsiPos, PsiNeg));
}

//   // now C should be a 1x1 matrix, which gives the normalization
//   TRACE(C);
//   Psi.set_front(prod(C, Psi.get_front()));
//}

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-make-positive <psi>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(argv[1], CacheSize);

   LinearWavefunction Psi = *PsiPtr;

   MatrixOperator CPos, CNeg;
   Split(Psi, CPos, CNeg);

   *PsiPtr.mutate() = Psi;
   pheap::ShutdownPersistent(PsiPtr);
}
