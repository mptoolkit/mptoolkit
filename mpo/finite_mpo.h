// -*- C++ -*- $Id$
//
// Matrix product operator defined on finite support.
// The boundary states are normally one dimensional
// (we also allow the operator to be reducible, 
// representing a sum of quantum number components, in which case the Basis1() will have dimension > 1).
// We used to require that the Basis2() was a scalar, but no longer,
// this isn't possible for extracted components of triangular or
// generic operators.  But we can always do a delta_shift to give a scalar (can we??!?)

#if !defined(FINITE_MPO_H_JDCHJKEHY589758YUER89H489)
#define FINITE_MPO_H_JDCHJKEHY589758YUER89H489

#include "generic_mpo.h"

class FiniteMPO
{
   private:
      typedef GenericMPO data_type;

   public:
      typedef OperatorComponent         value_type; 
      typedef data_type::const_iterator const_iterator;
      typedef data_type::iterator       iterator;
      typedef OperatorComponent::basis1_type basis1_type;
      typedef OperatorComponent::basis2_type basis2_type;

      FiniteMPO() {}

      FiniteMPO(FiniteMPO const& Other) : Data(Other.Data) {}

      explicit FiniteMPO(int Size) : Data(Size) {}

      FiniteMPO(int Size, LatticeCommute Com) : Data(Size, Com) {}

      // Construction from a generic MPO.  The generic MPO must already be in finite form.
      explicit FiniteMPO(GenericMPO const& Other);

      FiniteMPO& operator=(FiniteMPO const& Other) { Data = Other.Data; return *this; }

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
      // WARNING: Use qn1() instead of TransformsAs() where appropriate.
      // For debugging purposes, to audit usage of TransformsAs, we have an include guard
#if !defined(DISABLE_FINITE_MPO_TRANSFORMS_AS)
      QuantumNumbers::QuantumNumber TransformsAs() const;
#endif

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

   private:
      data_type Data;
};

PStream::opstream& operator<<(PStream::opstream& out, FiniteMPO const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, FiniteMPO& op);

// Returns the MPO that is Op1 \otimes Op2.
// PRECONDITION: Op1.Basis2() == Op2.Basis1()
FiniteMPO join(FiniteMPO const& Op1, FiniteMPO const& Op2);

// Repeats Op Count number of times, Op^{\oprod Count}.
// PRECONDITION: Op.Basis2() == Op.Basis1()
FiniteMPO repeat(FiniteMPO const& Op, int Count);

FiniteMPO& operator*=(FiniteMPO& x, double a);
FiniteMPO& operator*=(FiniteMPO& x, std::complex<double> a);

FiniteMPO& operator+=(FiniteMPO& x, FiniteMPO const& y);
FiniteMPO& operator-=(FiniteMPO& x, FiniteMPO const& y);

FiniteMPO operator+(FiniteMPO const& x, FiniteMPO const& y);
FiniteMPO operator-(FiniteMPO const& x, FiniteMPO const& y);

FiniteMPO operator-(FiniteMPO const& x);

FiniteMPO operator*(double a, FiniteMPO const& x);
FiniteMPO operator*(FiniteMPO const& x, double a);
FiniteMPO operator*(std::complex<double> a, FiniteMPO const& x);
FiniteMPO operator*(FiniteMPO const& x, std::complex<double> a);

FiniteMPO prod(FiniteMPO const& x, FiniteMPO const& y, QuantumNumbers::QuantumNumber const& q);
FiniteMPO prod(FiniteMPO const& x, FiniteMPO const& y);
FiniteMPO operator*(FiniteMPO const& x, FiniteMPO const& y);

// dot product - takes into account the multiplicity to rescale the result
FiniteMPO dot(FiniteMPO const& x, FiniteMPO const& y);

// cross product (if it exists)
FiniteMPO cross(FiniteMPO const& x, FiniteMPO const& y);

// outer product of tensors.  This is defined as the product to the maximum
// degree quantum number q.  There is also a scaling factor sqrt(degree(q))
FiniteMPO outer(FiniteMPO const& x, FiniteMPO const& y);

// project a (reducible) quantum number onto an irreducible component
FiniteMPO project(FiniteMPO const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.
FiniteMPO pow(FiniteMPO const& x, int n);

// Exponential operator.
FiniteMPO exp(FiniteMPO const& x);

// Conjugate
FiniteMPO conj(FiniteMPO const& x);

// Adjoint
FiniteMPO adjoint(FiniteMPO const& x);

// optimize the representation
void optimize(FiniteMPO& Op);

// completely coarse-grain the MPO into a simple operator.
// The dimensions of this operator are exponentially big in the number of sites
// in x, so be careful!
// For non-abelian symmetries, this coarse-grain occurs from left to right.
SimpleRedOperator coarse_grain(FiniteMPO const& x);

// The opposite of coarse_grain - decompose an operator acting on the entire Hilbert space
// into a FiniteMPO
FiniteMPO fine_grain(SimpleOperator const& x,
		     std::vector<BasisList> const& LocalBasis1,
		     std::vector<BasisList> const& LocalBasis2);

// Make an identity operator that acts on the same local Hilbert space as x
FiniteMPO
MakeIdentityFrom(FiniteMPO const& x);

// output to a stream
std::ostream& operator<<(std::ostream& out, FiniteMPO const& x);

#include "finite_mpo.cc"

#endif
