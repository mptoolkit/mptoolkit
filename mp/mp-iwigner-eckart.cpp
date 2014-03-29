// -*- C++ -*- $Id: mp-wigner-eckart.cpp 1149 2012-04-18 03:12:37Z ianmcc $

//
// Project a wavefunction using the Wigner-Eckart theorem
//
// The Regularize option is bugged - somehow the C_right and C_old
// matrices end up with incompatible bases, perhaps due to sorting of quantum numbers?

#include "mps/infinitewavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "tensor/wigner_eckart.h"
#include "tensor/regularize.h"

using QuantumNumbers::QuantumNumber;
using QuantumNumbers::Projection;

InfiniteWavefunction WignerProjectWavefunction(InfiniteWavefunction const& Psi, 
                                               SymmetryList const& FinalSL, bool RegularizeBasis = false)
{
   InfiniteWavefunction Result;
   Result.Psi = LinearWavefunction(FinalSL);
   VectorBasis b1 = Psi.Psi.Basis1();
   WignerEckartBasis<VectorBasis> W1(b1, FinalSL);
   MatrixOperator Reg1;
   if (RegularizeBasis)
      Reg1 = Regularize(W1.AbelianBasis());

   // The identity projection
   QuantumNumbers::ProjectionList PL = enumerate_projections(Psi.C_old.TransformsAs());
   CHECK_EQUAL(PL.size(), 1U);
   Projection IdentP = PL[0];

   QuantumNumbers::ProjectionList QPL = enumerate_projections(Psi.shift());
   Result.QShift = change(QuantumNumber(FinalSL), QPL[0]);
   
   Result.C_old = wigner_eckart(Psi.C_old, IdentP, W1, W1);
   if (RegularizeBasis)
      Result.C_old = triple_prod(Reg1, Result.C_old, herm(Reg1));

   LinearWavefunction::const_iterator I = Psi.Psi.begin();

   while (I != Psi.Psi.end())
   {
      std::cerr << "Working...\n";
      BasisList LocalBasis = adjoint(I->LocalBasis());
      std::set<QuantumNumber> LocalQN(LocalBasis.begin(), LocalBasis.end());
      WignerEckartBasis<VectorBasis> W2(I->Basis2(), FinalSL);

      BasisList AbelianLocalBasis(FinalSL);
      for (unsigned i = 0; i < I->LocalBasis().size(); ++i)
      {
         QuantumNumbers::ProjectionList pl = enumerate_projections(I->LocalBasis()[i]);
         for (unsigned pi = 0; pi < pl.size(); ++pi)
            {
               AbelianLocalBasis.push_back(map_projection_to_quantum(pl[pi], FinalSL));
            }
      }

      StateComponent Next(AbelianLocalBasis, W1.AbelianBasis(), W2.AbelianBasis());
      int k = 0;
      for (unsigned i = 0; i < I->LocalBasis().size(); ++i)
      {
         QuantumNumbers::ProjectionList pl = enumerate_projections(I->LocalBasis()[i]);
         for (unsigned pi = 0; pi < pl.size(); ++pi)
         {
            Next[k++] = wigner_eckart((*I)[i], pl[pi], W1, W2);
         }
      }

      if (RegularizeBasis)
      {
	 MatrixOperator Reg2 = Regularize(W2.AbelianBasis());
	 Next = triple_prod(Reg1,  Next, herm(Reg2));
	 Reg1 = Reg2;
      }

      Result.Psi.push_back(Next);
      W1 = W2;

      ++I;
   }

   Result.C_right = wigner_eckart(Psi.C_right, IdentP, W1, W1);
   if (RegularizeBasis)
      Result.C_right = triple_prod(Reg1, Result.C_right, herm(Reg1));

   // normalize
   //   Result *= double(degree(Psi.Psi.TransformsAs()));

   return Result;
}

int main(int argc, char** argv)
{
   if (argc != 4)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-iwigner-eckart <input-psi> <output-psi> <symmetry-list>\n";
      return 1;
   }

   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pheap::Initialize(argv[2], 1, PageSize, CacheSize);
   pvalue_ptr<InfiniteWavefunction> PsiIn = pheap::ImportHeap(argv[1]);
   std::string FinalSLStr = argv[3];

   SymmetryList FinalSL = SymmetryList(FinalSLStr);

   pvalue_ptr<InfiniteWavefunction> PsiNew = new InfiniteWavefunction(WignerProjectWavefunction(*PsiIn,
                                                                                                FinalSL, false));


   pheap::ShutdownPersistent(PsiNew);
}
