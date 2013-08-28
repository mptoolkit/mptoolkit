// -*- C++ -*- $Id$
//
// Matrix product operator defined on finite support.
// The boundary states are normally one dimensional
// (we allow the operator to be reducible, representing a sum of quantum number components)

#if !defined(FINITE_MPO_H_JDCHJKEHY589758YUER89H489)
#define FINITE_MPO_H_JDCHJKEHY589758YUER89H489

#include "generic_mpo.h"

class FiniteMPO
{
   private:
   //      typedef std::deque<OperatorComponent> data_type;
      typedef GenericMPO data_type;

   public:
      typedef OperatorComponent         value_type; 
      typedef data_type::const_iterator const_iterator;
      typedef OperatorComponent::basis1_type basis1_type;
      typedef OperatorComponent::basis2_type basis2_type;

      FiniteMPO() {}

      explicit FiniteMPO(int Size) : Data(Size) {}

      // Construction from a generic MPO.  The generic MPO must already be in finite form.
      explicit FiniteMPO(GenericMPO const& Other);

      // returns the total number of sites this operator contains
      int size() const { return Data.size(); }

      // returns true if this is a zero operator
      bool empty() const { return Data.empty() || Data.front().Basis1().size() == 0; }
      bool is_null() const { return Data.empty() || Data.front().Basis1().size() == 0; }

      // returns true if the operator transforms irreducibly.  This is true
      // iff the left basis contains only a single state.
      bool is_irreducible() const;

      // precondition: is_irreducible
      QuantumNumbers::QuantumNumber TransformsAs() const;

      // returns the left-most basis.  This is guaranteed to contain each
      // quantum number at most once.
      basis1_type const& Basis1() const { return Data.front().Basis1(); }

      basis2_type const& Basis2() const { return Data.back().Basis2(); }

      value_type& operator[](int n) { return Data[n]; }
      value_type const& operator[](int n) const { return Data[n]; }

      // iterate over the MPOpComponents at each site
      const_iterator begin() const { return Data.begin(); }
      const_iterator end() const { return Data.end(); }

      value_type& front() { return Data.front(); }
      value_type const& front() const { return Data.front(); }

      value_type& back() { return Data.back(); }
      value_type const& back() const { return Data.back(); }

      // small problem here: if this->is_null(), this will not work.
      SymmetryList GetSymmetryList() const { return Data.front().GetSymmetryList(); }

      // implicit conversion to a const GenericMPO
      operator GenericMPO const&() const { return Data; }

      // direct access to the MPOpCompressed
      data_type& data() { return Data; }
      data_type const& data() const { return Data; }

   private:
      data_type Data;
};

PStream::opstream& operator<<(PStream::opstream& out, FiniteMPO const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, FiniteMPO& op);

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

// project a (reducible) quantum number onto an irreducible component
FiniteMPO project(FiniteMPO const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.  Only useful for n small!
FiniteMPO pow(FiniteMPO const& x, int n);

// Conjugate
FiniteMPO conj(FiniteMPO const& x);

// Adjoint
FiniteMPO adjoint(FiniteMPO const& x);

// optimize the representation
void optimize(FiniteMPO& Op);

// output to a stream
std::ostream& operator<<(std::ostream& out, FiniteMPO const& x);

#include "finite_mpo.cc"

#endif
