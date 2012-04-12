// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include "tensor/wigner_eckart.h"
#include "tensor/regularize.h"

using QuantumNumbers::QuantumNumber;
using QuantumNumbers::Projection;

MPWavefunction WignerProjectWavefunction(MPWavefunction const& Psi, 
                                         Projection const& p,
                                         SymmetryList const& FinalSL, bool RegularizeBasis = false)
{
   MPWavefunction Result(FinalSL);
   VectorBasis b1 = Psi.Basis1();
   QPSetType AllowedProj;
   CHECK_EQUAL(b1.size(), 1)("this currently supports only one-dimensional targets");
   AllowedProj.insert(std::make_pair(b1[0], p));
   WignerEckartBasis<VectorBasis> W1(b1, AllowedProj, FinalSL);
   MatrixOperator Reg1;
   if (RegularizeBasis)
      Reg1 = Regularize(W1.AbelianBasis());

   MPWavefunction::const_iterator I = Psi.begin();

   while (I != Psi.end())
   {
      std::cerr << "Working...\n";
      BasisList LocalBasis = adjoint(I->SiteBasis());
      std::set<QuantumNumber> LocalQN(LocalBasis.begin(), LocalBasis.end());
      AllowedProj = UpdateAllowedProjections(AllowedProj, LocalQN);
      WignerEckartBasis<VectorBasis> W2(I->Basis2(), AllowedProj, FinalSL);

      BasisList AbelianLocalBasis(FinalSL);
      for (unsigned i = 0; i < I->SiteBasis().size(); ++i)
      {
         ProjectionList pl = enumerate_projections(I->SiteBasis()[i]);
         for (unsigned pi = 0; pi < pl.size(); ++pi)
            {
               AbelianLocalBasis.push_back(map_projection_to_quantum(pl[pi], FinalSL));
            }
      }

      MPStateComponent Next(AbelianLocalBasis, W1.AbelianBasis(), W2.AbelianBasis());
      int k = 0;
      for (unsigned i = 0; i < I->SiteBasis().size(); ++i)
      {
         ProjectionList pl = enumerate_projections(I->SiteBasis()[i]);
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

      Result.push_back(Next);
      W1 = W2;

      ++I;
   }

   // normalize
   Result *= double(degree(Psi.TransformsAs()));

   return Result;
}

int main(int argc, char** argv)
{
   if (argc != 5)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-wigner-eckart <input-psi> <output-psi> <symmetry-list> <projection>\n";
      return 1;
   }

   int PageSize = getenv_or_default("MP_PAGESIZE", 65536);
   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pheap::Initialize(argv[2], 1, PageSize, CacheSize);
   pvalue_ptr<MPWavefunction> PsiIn = pheap::ImportHeap(argv[1]);
   std::string FinalSLStr = argv[3];
   std::string Proj = argv[4];

   SymmetryList FinalSL = SymmetryList(FinalSLStr);
   Projection p(PsiIn->GetSymmetryList(), Proj);

   pvalue_ptr<MPWavefunction> OutPsi = new MPWavefunction(WignerProjectWavefunction(*PsiIn,
                                                                                    p,
                                                                                    FinalSL));
   pheap::ShutdownPersistent(OutPsi);
}
