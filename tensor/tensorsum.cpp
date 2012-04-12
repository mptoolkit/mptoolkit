// -*- C++ -*- $Id$

#include "tensorsum.h"
#include <boost/type_traits.hpp>

namespace Tensor
{

void SumBasis<BasisList>::AddBasis(BasisList const& B)
{
   PRECONDITION(B.GetSymmetryList() == this->GetSymmetryList());
   BList_.push_back(B);
   std::vector<int> Map(B.size());
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Map[i] = this->size();
      this->push_back(B[i]);
   }
   Mapping_.push_back(Map);
}

void SumBasis<VectorBasis>::AddBasis(VectorBasis const& B)
{
   PRECONDITION(B.GetSymmetryList() == This_.GetSymmetryList());
   BList_.push_back(B);
   std::vector<int> Map(B.size());
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Map[i] = This_.size();
      This_.push_back(B[i], B.dim(i));
   }
   Mapping_.push_back(Map);
}

} // namespace Tensor
