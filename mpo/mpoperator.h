// -*- C++ -*- $Id$

// MPOperator is the basic class for a matrix product operator.
// It provides direct access, unlike the compressed LinearOperator.

#if !defined(MPOPERATOR_H_JDCHJKEHY589758YUER89H489)
#define MPOPERATOR_H_JDCHJKEHY589758YUER89H489

#include "operator_component.h"
#include <vector>

class MPOperator
{
   private:
      typedef std::vector<OperatorComponent> DataType;

   public:
      typedef OperatorComponent value_type;
      typedef DataType::iterator iterator;
      typedef DataType::const_iterator const_iterator;
      typedef value_type::basis1_type basis1_type;
      typedef value_type::basis2_type basis2_type;

      MPOperator() {}

      explicit MPOperator(int Size) : Data_(Size) {}

      explicit MPOperator(OperatorComponent const& x) : Data_(1, x) {}

      // Size repeated copies of x
      MPOperator(int Size, OperatorComponent const& x) : Data_(Size, x) {}

      // from an iterator
      template <typename InIter>
      MPOperator(InIter first, InIter last) : Data_(first, last) {}

      bool empty() const { return Data_.empty(); }
      std::size_t size() const { return Data_.size(); }

      bool is_null() const;

      iterator begin() { return Data_.begin(); }
      iterator end() { return Data_.end(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }

      value_type& operator[](int x) { return Data_[x]; }
      value_type const& operator[](int x) const { return Data_[x]; }
    
      value_type const& front() const { return Data_.front(); }
      value_type& front() { return Data_.front(); }

      value_type const& back() const { return Data_.back(); }
      value_type& back() { return Data_.back(); }

      basis1_type Basis1() const { DEBUG_CHECK(!this->empty()); return Data_.front().Basis1(); }
      basis2_type Basis2() const { DEBUG_CHECK(!this->empty()); return Data_.back().Basis2(); }

      QuantumNumbers::SymmetryList GetSymmetryList() const { return Data_[0].GetSymmetryList(); }

      // returns the list of local hilbert spaces for this operator
      std::vector<BasisList> LocalBasis1List() const;
      std::vector<BasisList> LocalBasis2List() const;

      std::vector<OperatorComponent> const& data() const { return Data_; }

   private:
      DataType Data_;
};

std::ostream&
operator<<(std::ostream& out, MPOperator const& op);

// remove unused matrix elements
void cull_unused_elements(MPOperator& Op);

// As an alternative to cull_unused_elements(), we can use a mask vector which indicates which components are unused
void initialize_mask(MPOperator const& Op, std::vector<std::vector<int> >& Mask);

void mask_unused_elements(MPOperator const& Op, std::vector<std::vector<int> >& Mask);

// Does a 2-1 coarse graining of an operator.  The length must be a multiple of 2
MPOperator coarse_grain(MPOperator const& Op);

struct MPOperatorClassification
{
   // indicates that the operator is zero
   bool is_null() const;

   // indicates that the operator a product state, ie a product of 1x1 MPO's
   bool is_product() const;

   // indicates that the operator is a unitary product state, ie a string operator
   bool is_string() const;

   // indicates that the operator is proportional to a unitary product state,
   // up to some complex factor
   bool is_prop_string() const;

   // indicates that the operator is proportional to the identity operator
   bool is_prop_identity() const;

   // indicates that the operator is equal to the identity
   bool is_identity() const;

   // returns true only if the operator fits into no other classification
   bool is_unclassified() const;

   // for operators that are proportional to the identity, returns the factor
   std::complex<double> factor() const;

   // private use only
   std::complex<double> Factor_;
   bool Product_;
   bool String_;
   bool Identity_;
   bool PropString_;
   bool PropIdentity_;
   bool Null_;

   MPOperatorClassification();
};

std::ostream& operator<<(std::ostream& out, MPOperatorClassification const& Class);

MPOperatorClassification classify(MPOperator const& Op);

// plus various functions for acting on states etc

// TODO: find a better location for this function
// Construct an operator that projects onto a given subset of a basis.
SimpleOperator make_projector_onto(BasisList const& Basis, std::set<int> const& Onto);

std::vector<BasisList>
ExtractLocalBasis1(MPOperator const& Op);

#endif
