// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp/copyright.h"
#include "common/gmpint.h"
#include "common/environment.h"

using namespace Tensor;

typedef std::map<QuantumNumber, gmp::bigint> SymmetrySizes;

void AddSite(SymmetrySizes const& S, BasisList const& Basis, SymmetrySizes& Result)
{
   using QuantumNumbers::QuantumNumberList;
   Result.clear();
   int const bsize = Basis.size();
   SymmetrySizes::const_iterator const Iend = S.end();
   for (int i = 0; i < bsize; ++i)
   {
      for (SymmetrySizes::const_iterator I = S.begin(); I != Iend; ++I)
      {
         QuantumNumberList QList = transform_targets(I->first, Basis[i]);
         QuantumNumberList::const_iterator const Qend = QList.end();
	 for (QuantumNumberList::const_iterator s = QList.begin(); s != Qend; ++s)
         {
	    Result[*s] += I->second;
	 }
      }
   }
}

int main(int argc, char** argv)
{
   if (argc != 2)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "usage: mp-dimensions <lattice>\n";
      return 1;
   }

   long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
   pvalue_ptr<OperatorList> OpList = pheap::OpenPersistent(argv[1], CacheSize, true);

   Lattice L = OpList->GetLattice();
   SymmetrySizes Size;
   Size[QuantumNumber(L.GetSymmetryList())] = 1L;  // 'seed' it with a vacuum state
   for (int i = L.size(); i > 0; --i)
   {
      SymmetrySizes Size2;
      AddSite(Size, L[i].Basis1().Basis(), Size2);
      Size.swap(Size2);
      std::cerr << "Working (" << i << ") n_subspaces = " << Size.size() << '\n';
   }

   gmp::bigint Total = 0;
   std::cout << "Symmetry list is " << L.GetSymmetryList() << '\n';
   std::cout << "Symmetry Sector Dimension\n";
   for (SymmetrySizes::const_iterator I = Size.begin(); I != Size.end(); ++I)
   {
      std::cout << std::setw(15) << boost::lexical_cast<std::string>(I->first) << " "
                << I->second << '\n';
      Total += I->second;
   }
   std::cout << "\nTotal: " << Total << '\n';

   pheap::Shutdown();
}

