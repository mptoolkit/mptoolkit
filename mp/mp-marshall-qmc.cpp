// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/random_wavefunc.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include "mp/copyright.h"
#include "common/environment.h"
#include <tuple>

typedef std::complex<double> complex;

class MpsQmc
{
   public:
      MpsQmc(LinearWavefunction const& Psi_, WavefunctionDesc const& Configuration_,
             int NumSites);

      void Flip();

      bool ShiftRight();   // returns false if we cannot shift right
      bool ShiftLeft();    // returns false if we cannot shift left

      std::complex<double> Value() const { return CurrentValue; }

   private:
      LinearWavefunction Psi;
      WavefunctionDesc Configuration;

      std::vector<SiteBasis> LocalBasis;

      std::stack<MatrixOperator> LeftStack;
      std::stack<MatrixOperator> RightStack;
      int iLoc;
      LinearWavefunction::const_iterator pLoc;  // the iterator corresponding to iLoc
      std::complex<double> CurrentValue;
      int NumFlipSites;

      struct FlipInfoType
      {
         std::vector<int> State;
         MatrixOperator Weight;
         QuantumNumber Q;
      };

      std::vector<FlipInfoType> FlipInfoAtSiteType;
      std::vector<FlipInfoAtSiteType> FlipInfo;
};

MpsQmc::MpsQmc(LinearWavefunction const& Psi_, WavefunctionDesc const& Configuration_)
   : Psi(Psi_), Configuration(Configuration_), NumFlipSites(NumSites), FlipInfo(Psi.size())
{
   CHECK_EQUAL(Configuration.State.size(), Psi.size());
   RightStack.push(MatrixOperator::make_identity(Psi.Basis2()));
   pLoc = Psi.end();
   iLoc = Configuration.State.size();
   while (pLoc != Psi.begin())
   {
      --pLoc; --iLoc;
      RightStack.push(prod((*pLoc)[Configuration.State[iLoc]], 
                           RightStack.top(), 
                           Configuration.Height[iLoc]));
      LocalBasis.push_back(pLoc->LocalBasis());
   }
   LeftStack.push(MatrixOperator::make_identity(Psi.Basis1()));
   std::reverse(LocalBasis.begin(), LocalBasis.end());

   CurrentValue = prod(LeftStack.top(), RightStack.top(), Configuration.Height[0]);

   LinearWavefunction::const_iterator rLoc = pLoc;
   std::advance(rLoc, NumFlipSites);
   while (rLoc != Psi.end())
   {
      std::vector<int> State(NumFlipSites, 0);
      bool done = false;
      while (!done)
      {
         LinearWavefunction::const_iterator I = pLoc;
         FlipInfoType FInfo;
         FInfo.Q = QuantumNumber(Psi.GetSymmetryList());
         FInfo.Weight = MatrixOperator::make_identity(I->Basis1());
         FInfo.State = State;
         for (int i = 0; i < NumFlipSites; ++i, ++I)
         {
            Q = transform_targets(Q, I->LocalBasis()[State[i]])[0];
            Weight = Weight * (*I);
         }
         FlipInfo[iLoc].push_back(FInfo);
         int i = NumFlipSites-1;
         while (i >= 0)
         {
            if (++State[i] == LocalBasis[i+iLoc].size())
            {
               State[i] == 0;
               --i;
            }
         }
         done = (i < 0);
      }
      ++rLoc;
      ++pLoc;
      ++iLoc;
   }

   pLoc = Psi.begin();
   iLoc = 0;
}

void MpsQmc::Flip()
{
   // flip the current state without changing the quantum number
   // firstly, figure out what the current quantum number actually is
   QuantumNumber Q(Psi.GetSymmetryList());
   for (int i = 0; i < NumFlipSites; ++i)
   {
      Q = transform_targets(Q, LocalBasis[iLoc+i])[0];
   }
   // get a list of the possible flips with that quantum number
   std::vector<int> PossibleFlips;
   for (int 

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-marshall-qmc <psi>\n";
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
