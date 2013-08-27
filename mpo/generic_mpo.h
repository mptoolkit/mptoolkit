// -*- C++ -*- $Id$

// GenericMPO represents an MPO that is in no particular form

#if !defined(GENERIC_MPO_H_JDCHJKEHY589758YUER89H489)
#define GENERIC_MPO_H_JDCHJKEHY589758YUER89H489

#include "operator_component.h"
#include <vector>

class GenericMPO
{
   private:
      typedef std::vector<OperatorComponent> DataType;

   public:
      typedef OperatorComponent value_type;
      typedef DataType::iterator iterator;
      typedef DataType::const_iterator const_iterator;
      typedef value_type::basis1_type basis1_type;
      typedef value_type::basis2_type basis2_type;

      GenericMPO() {}

      explicit GenericMPO(int Size) : Data_(Size) {}

      explicit GenericMPO(OperatorComponent const& x) : Data_(1, x) {}

      // Size repeated copies of x
      GenericMPO(int Size, OperatorComponent const& x) : Data_(Size, x) {}

      // from an iterator
      template <typename InIter>
      GenericMPO(InIter first, InIter last) : Data_(first, last) {}

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

   friend PStream::opstream& operator<<(PStream::opstream& out, GenericMPO const& op);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, GenericMPO& op);
};

std::ostream&
operator<<(std::ostream& out, GenericMPO const& op);

PStream::opstream&
operator<<(PStream::opstream& out, GenericMPO const& op);

PStream::ipstream&
operator>>(PStream::ipstream& in, GenericMPO& op);

GenericMPO& operator*=(GenericMPO& x, double a);
GenericMPO& operator*=(GenericMPO& x, std::complex<double> a);

GenericMPO operator*(double a, GenericMPO const& x);
GenericMPO operator*(GenericMPO const& x, double a);
GenericMPO operator*(std::complex<double> a, GenericMPO const& x);
GenericMPO operator*(GenericMPO const& x, std::complex<double> a);

// remove unused matrix elements
void cull_unused_elements(GenericMPO& Op);

// As an alternative to cull_unused_elements(), we can use a mask vector which indicates which components are unused
void initialize_mask(GenericMPO const& Op, std::vector<std::vector<int> >& Mask);

void mask_unused_elements(GenericMPO const& Op, std::vector<std::vector<int> >& Mask);

// Does a 2-1 coarse graining of an operator.  The length must be a multiple of 2
GenericMPO coarse_grain(GenericMPO const& Op);

struct GenericMPOClassification
{
   // indicates that the operator is zero
   bool is_null() const;

   // indicates that the operator a product state, ie a product of 1x1 MPO's
   bool is_product() const;

   // indicates that the operator is a unitary product state, ie a string operator
   bool is_unitary() const;

   // indicates that the operator is proportional to a unitary product state,
   // up to some factor
   bool is_prop_unitary() const;

   // indicates that the operator is proportional to the identity operator
   bool is_prop_identity() const;

   // returns true if the operator is the identity multiplied by a complex phase factor of magnitude 1
   bool is_complex_identity() const
   {
      return this->is_prop_identity() && norm_frob(norm_frob(this->factor())-1.0) < 1E-12;
   }

   // indicates that the operator is equal to the identity
   bool is_identity() const;

   // returns true only if the operator fits into no other classification
   bool is_unclassified() const;

   // for operators that are proportional to the identity, returns the factor
   std::complex<double> factor() const;

   // private use only
   std::complex<double> Factor_;
   bool Product_;
   bool Unitary_;
   bool Identity_;
   bool PropUnitary_;
   bool PropIdentity_;
   bool Null_;

   GenericMPOClassification();
};

std::ostream& operator<<(std::ostream& out, GenericMPOClassification const& Class);

GenericMPOClassification classify(GenericMPO const& Op);

// plus various functions for acting on states etc

// TODO: find a better location for this function
// Construct an operator that projects onto a given subset of a basis.
SimpleOperator make_projector_onto(BasisList const& Basis, std::set<int> const& Onto);

std::vector<BasisList>
ExtractLocalBasis1(GenericMPO const& Op);

#endif