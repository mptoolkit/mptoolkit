// -*- C++ -*- $Id$
//
// Matrix product operator for a translationally invariant (or k-dependent)
// product, such as a unitary or non-unitary evolution operator.
//
// Addition isn't defined for ProductMPO.
// These operators necessarily transform as scalars.  ?maybe?

#if !defined(MPTOOLKIT_MPO_PRODUCT_MPO_H)
#define MPTOOLKIT_MPO_PRODUCT_MPO_H

#include "product_mpo.h"

class ProductMPO
{
   private:
      typedef GenericMPO data_type;

   public:
      typedef OperatorComponent         value_type; 
      typedef data_type::const_iterator const_iterator;
      typedef data_type::iterator       iterator;
      typedef OperatorComponent::basis1_type basis1_type;
      typedef OperatorComponent::basis2_type basis2_type;

      ProductMPO() {}

      ProductMPO(ProductMPO const& Other) : Data(Other.Data) {}

      // Removed this constructor because it doesn't make much sense to define a ProductMPO
      // without specifying the LatticeCommute
      //      explicit ProductMPO(int Size) : Data(Size) {}

      ProductMPO(int Size, LatticeCommute Com) : Data(Size, Com) {}

      // Construction from a generic MPO.  The generic MPO must already be in finite form.
      explicit ProductMPO(GenericMPO const& Other);

      ProductMPO& operator=(ProductMPO const& Other) { Data = Other.Data; return *this; }

      // returns the total number of sites this operator contains
      int size() const { return Data.size(); }

      // returns true if this is a zero operator
      bool empty() const { return Data.empty() || Data.front().Basis1().size() == 0; }
      bool is_null() const { return Data.empty() || Data.front().Basis1().size() == 0; }

      // returns true if the operator transforms irreducibly.  This is true
      // iff the left basis contains only a single state, and the right basis contains the vacuum.
      // Note that finite MPO's are not irreducible tensor operators as such (they are more like
      // Wigner operators)
      bool is_irreducible() const;

      // returns true if the operator transforms as a rotational invariant, ie
      // it is irreducible in the scalar symmetry sector
      bool is_scalar() const;

      // returns true if this MPO is the identity operator, that is, a 1x1 MPO that
      // is a product of identity operators.
      bool is_identity() const;

      // precondition: is_irreducible
      QuantumNumbers::QuantumNumber TransformsAs() const;

      // returns the quantum number in the left basis.  If the right basis is the vacuum
      // then this is also the TransformsAs()
      // precondition: Basis1().size() == 1
      QuantumNumbers::QuantumNumber qn1() const;

      // returns the quantum number in the right basis.
      // This doesn't have to be the vacuum state.
      QuantumNumbers::QuantumNumber qn2() const;

      // returns the left-most basis.  This is guaranteed to contain each
      // quantum number at most once.
      basis1_type const& Basis1() const { return Data.front().Basis1(); }

      basis2_type const& Basis2() const { return Data.back().Basis2(); }

      value_type& operator[](int n) { return Data[n]; }
      value_type const& operator[](int n) const { return Data[n]; }

      // iterate over the MPOpComponents at each site
      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      iterator begin() { return Data.begin(); }
      iterator end() { return Data.end(); }

      value_type& front() { return Data.front(); }
      value_type const& front() const { return Data.front(); }

      value_type& back() { return Data.back(); }
      value_type const& back() const { return Data.back(); }

      // small problem here: if this->is_null(), this will not work.
      SymmetryList GetSymmetryList() const { return Data.front().GetSymmetryList(); }

      // implicit conversion to a const GenericMPO
      operator GenericMPO const&() const { return Data; }

      // Return the local basis at the n'th site
      BasisList const& LocalBasis1(int n) const
      { return Data.LocalBasis1(n); }
      BasisList const& LocalBasis2(int n) const
      { return Data.LocalBasis2(n); }

      // returns the list of local hilbert spaces for this operator
      std::vector<BasisList> LocalBasis1List() const
      { return Data.LocalBasis1List(); }
      std::vector<BasisList> LocalBasis2List() const
      { return Data.LocalBasis2List(); }

      LatticeCommute Commute() const { return Data.Commute(); }
      void SetCommute(LatticeCommute x) { Data.SetCommute(x); }

      // direct access to the GenericMPO
      data_type& data() { return Data; }
      data_type const& data() const { return Data; }

      void check_structure() const { Data.check_structure(); }
      void debug_check_structure() const { Data.debug_check_structure(); }

      // Make an identity MPO over the given unit cell basis
      static ProductMPO make_identity(std::vector<BasisList> const& Basis);

   private:
      data_type Data;
};

PStream::opstream& operator<<(PStream::opstream& out, ProductMPO const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, ProductMPO& op);

// Returns the MPO that is Op1 \otimes Op2.
// PRECONDITION: Op1.Basis2() == Op2.Basis1()
ProductMPO join(ProductMPO const& Op1, ProductMPO const& Op2);

// Repeats Op Count number of times, Op^{\oprod Count}.
// PRECONDITION: Op.Basis2() == Op.Basis1()
ProductMPO repeat(ProductMPO const& Op, int Count);

ProductMPO& operator*=(ProductMPO& x, double a);
ProductMPO& operator*=(ProductMPO& x, std::complex<double> a);

ProductMPO operator*(double a, ProductMPO const& x);
ProductMPO operator*(ProductMPO const& x, double a);
ProductMPO operator*(std::complex<double> a, ProductMPO const& x);
ProductMPO operator*(ProductMPO const& x, std::complex<double> a);

ProductMPO prod(ProductMPO const& x, ProductMPO const& y, QuantumNumbers::QuantumNumber const& q);
ProductMPO prod(ProductMPO const& x, ProductMPO const& y);
ProductMPO operator*(ProductMPO const& x, ProductMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
ProductMPO dot(ProductMPO const& x, ProductMPO const& y);

// cross product (if it exists)
ProductMPO cross(ProductMPO const& x, ProductMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
ProductMPO outer(ProductMPO const& x, ProductMPO const& y);

// project a (reducible) operator onto an irreducible component
ProductMPO project(ProductMPO const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.
ProductMPO pow(ProductMPO const& x, int n);

// Exponential operator.
ProductMPO exp(ProductMPO const& x);

// Conjugate
ProductMPO conj(ProductMPO const& x);

// Adjoint
ProductMPO adjoint(ProductMPO const& x);

// optimize the representation
void optimize(ProductMPO& Op);

// completely coarse-grain the MPO into a simple operator.
// The dimensions of this operator are exponentially big in the number of sites
// in x, so be careful!
// For non-abelian symmetries, this coarse-grain occurs from left to right.
SimpleRedOperator coarse_grain(ProductMPO const& x);

// The opposite of coarse_grain - decompose an operator acting on the entire Hilbert space
// into a ProductMPO
ProductMPO fine_grain(SimpleOperator const& x, LatticeCommute Com,
		     std::vector<BasisList> const& LocalBasis1,
		     std::vector<BasisList> const& LocalBasis2);

// Make an identity operator that acts on the same local Hilbert space as x
ProductMPO
MakeIdentityFrom(ProductMPO const& x);

// Make an identity operator that acts on the same local Hilbert space as x,
// with the given quantum number in the auxiliary basis
ProductMPO
MakeIdentityFrom(ProductMPO const& x, QuantumNumber const& q);

// output to a stream
std::ostream& operator<<(std::ostream& out, ProductMPO const& x);

#include "product_mpo.cc"

#endif
