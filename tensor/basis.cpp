// -*- C++ -*- $Id$

#include "basis.h"

namespace Tensor
{

//
// BasisList
//

bool
BasisList::is_identity() const
{
   return this->size() == 1 && is_scalar(this->front());
}

bool
BasisList::is_regular() const
{
   std::set<QuantumNumber> Used;
   for (const_iterator I = this->begin(); I != this->end(); ++I)
   {
      if (Used.find(*I) != Used.end())
         return false;

      Used.insert(*I);
   }
   return true;
}

int BasisList::total_degree() const
{
   int Deg = 0;
   for (std::size_t i = 0; i < this->size(); ++i)
   {
      Deg += degree(Q_[i]);
   }
   return Deg;
}

int
BasisList::find_first(QuantumNumber const& q) const
{
   unsigned n = 0;
   while (n < Q_.size() && Q_[n] != q)
      ++n;

   return (n < Q_.size()) ? int(n) : -1;
}

int
BasisList::find_next(QuantumNumber const& q, int n) const
{
   ++n;
   while (n < int(Q_.size()) && Q_[n] != q)
      ++n;

   return (n < int(Q_.size())) ? n : -1;
}

std::ostream& operator<<(std::ostream& Out, BasisList const& B)
{
   Out << "Basis has symmetry " << B.GetSymmetryList()
       << ", size = " << B.size() << '\n';
   Out << "   N  QuantumNumber\n";
   int N = 0;
   for (BasisList::const_iterator I = B.begin(); I != B.end(); ++I)
   {
      Out << std::setw(4) << N << "  "
          << std::setw(13) << *I << '\n';
      ++N;
   }
   return Out;
}

std::string show_projections(BasisList const& B)
{
   using QuantumNumbers::ProjectionList;

   std::ostringstream Out;
   Out << "Basis has symmetry " << B.GetSymmetryList()
      << ", size = " << B.size() << '\n';
   Out << "   N  QuantumNumber  Degree  Projection\n";
   int N = 0;
   for (BasisList::const_iterator I = B.begin(); I != B.end(); ++I)
   {
      bool First = true;
      ProjectionList Pr = enumerate_projections(*I);
      for (ProjectionList::const_iterator PrI = Pr.begin(); PrI != Pr.end(); ++PrI)
      {
         if (First) 
         {
            Out << std::setw(4) << N << "  "
                << std::setw(13) << *I << "  "
                << std::setw(6) << degree(*I) << "  ";
            First = false;
         }
         else
            Out << "                             ";

         Out << std::setw(10) << *PrI << '\n';
      }
      ++N;
   }
   return Out.str();
}

BasisList delta_shift(BasisList const& Orig, QuantumNumbers::Projection const& p)
{
   BasisList Result(Orig.GetSymmetryList());
   for (BasisList::const_iterator I = Orig.begin(); I != Orig.end(); ++I)
      Result.push_back(change(*I, p));

   return Result;
}

#if 0
BasisList RenameSymmetry(BasisList const& BL, SymmetryList const& NewSL)
{
   BasisList Result(NewSL);
   for (unsigned i = 0; i < BL.size(); ++i)
      Result.push_back(RenameSymmetry(BL[i], NewSL));
   return Result;
}
#endif

std::map<int, int>
LinearizeQuantumNumberSubspace(BasisList const& b, QuantumNumbers::QuantumNumber const& q)
{
   std::map<int, int> Result;
   int n = 0;
   for (unsigned i = 0; i < b.size(); ++i)
   {
      if (b[i] == q)
	 Result[i] = n++;
   }
   return Result;
}

//
// VectorBasis
//

int VectorBasis::total_degree() const
{
   int Deg = 0;
   for (std::size_t i = 0; i < this->size(); ++i)
   {
      Deg += degree(Basis_[i]) * Dimension_[i];
   }
   return Deg;
}

std::ostream& operator<<(std::ostream& Out, VectorBasis const& B)
{
   Out << "Basis has symmetry " << B.GetSymmetryList()
       << ", subspace size = " << B.size() 
       << ", dimension = " << B.total_dimension()
       << ", degree = " << B.total_degree() << '\n';
   Out << "   N  QuantumNumber  Dimension\n";
   for (std::size_t N = 0; N < B.size(); ++N)
   {
      Out << std::setw(4) << N << "  "
          << std::setw(13) << B[N] << "  " 
          << std::setw(9) << B.dim(N) << '\n';
   }
   return Out;
}

std::string show_projections(VectorBasis const& B)
{
   using QuantumNumbers::ProjectionList;

   std::ostringstream Out;
   Out << "Basis has symmetry " << B.GetSymmetryList()
       << ", subspace size = " << B.size() 
       << ", dimension = " << B.total_dimension()
       << ", degree = " << B.total_degree() << '\n';
   Out << "   N  QuantumNumber  Dimension  Degree  Projection\n";
   for (std::size_t N = 0; N < B.size(); ++N)
   {
      bool First = true;
      ProjectionList Pr = enumerate_projections(B[N]);
      for (ProjectionList::const_iterator PrI = Pr.begin(); PrI != Pr.end(); ++PrI)
      {
         if (First) 
         {
            Out << std::setw(4) << N << "  "
                << std::setw(13) << B[N] << "  " 
                << std::setw(9) << B.dim(N) << "  "
                << std::setw(6) << degree(B[N]) << "  ";
            First = false;
         }
         else
            Out << "                                        ";

         Out << std::setw(10) << *PrI << '\n';
      }
   }
   return Out.str();
}

VectorBasis delta_shift(VectorBasis const& Orig, QuantumNumbers::Projection const& p)
{
   return VectorBasis(delta_shift(Orig.Basis(), p), 
		      Orig.Dimension_.begin(),
		      Orig.Dimension_.end());
}

VectorBasis
CoerceSymmetryList(VectorBasis& b, SymmetryList const& sl)
{
   VectorBasis Result;
   Result.Basis_ = CoerceSymmetryList(b.Basis_, sl);
   Result.Dimension_ = b.Dimension_;
   return Result;
}

#if 0
VectorBasis RenameSymmetry(VectorBasis const& BL, SymmetryList const& NewSL)
{
   return VectorBasis(RenameSymmetry(BL.Basis_, NewSL), BL.Dimension_);
}
#endif

} // namespace Tensor
