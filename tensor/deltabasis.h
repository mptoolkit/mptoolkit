/* -*- C++ -*- $Id$

  A version of a BasisList that also contains a delta (projection).
*/

#if !defined(DELTABASIS_H_FHUY4837598YFEH87Y4587OHY)
#define DELTABASIS_H_FHUY4837598YFEH87Y4587OHY

#include "tensor/basis.h"

namespace Tensor
{

class DeltaBasis
{
   public:
      typedef std::pair<QuantumNumber, Projection> value_type;

   private:
      typedef std::vector<value_type> DataType;
      typedef DataType::iterator iterator;

   public:
      typedef DataType::const_iterator const_iterator;

      DeltaBasis() {}

      explicit DeltaBasis(QuantumNumbers::SymmetryList const& S) : S_(S) {}

      const_iterator begin() const { return Q_.begin(); }
      const_iterator end() const { return Q_.end(); }

      value_type const& operator[](int x) const { return Q_[x]; }

      std::size_t size() const { return Q_.size(); }

      bool is_null() const { return S_.is_null(); }

      // returns true if this basis is empty
      bool is_empty() const { return this->size() == 0; }

      void push_back(QuantumNumber const& q, Projection const& p) 
         { DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), S_); 
	 Q_.push_back(std::make_pair(q,p)); }

      // Since we include the delta, the degree and the dimension are the same
      int total_degree() const
      {
	 return this->size(); 
      }

      SymmetryList const& GetSymmetryList() const { return S_; }

   private:
      QuantumNumbers::SymmetryList S_;
      DataType Q_;

   friend bool operator==(DeltaBasis const& b1, DeltaBasis const& b2)
      { return b1.Q_ == b2.Q_; }

   friend bool operator!=(DeltaBasis const& b1, DeltaBasis const& b2)
      { return b1.Q_ != b2.Q_; }

   friend PStream::opstream& operator<<(PStream::opstream& out, DeltaBasis const& B);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, DeltaBasis& B);
   friend void CoerceSymmetryList(DeltaBasis& b, SymmetryList const& sl);
};

} // namespace Tensor

#endif
