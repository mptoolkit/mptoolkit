// -*- C++ -*- $Id$
//

#if !defined(LINEAROPERATOR_H_JDCHJKEHY589758YUER89H489)
#define LINEAROPERATOR_H_JDCHJKEHY589758YUER89H489

#include "operator_component.h"

class LinearOperator
{
   private:
   typedef std::deque<OperatorComponent> data_type;

   public:
      typedef OperatorComponent         value_type; 
      typedef data_type::const_iterator const_iterator;
      typedef OperatorComponent::basis1_type basis1_type;
      typedef OperatorComponent::basis2_type basis2_type;

      LinearOperator() {}

      explicit LinearOperator(int Size) : Data(Size) {}

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

      // small problem here: if this->is_null(), this will not work.
      SymmetryList GetSymmetryList() const { return Data.front().GetSymmetryList(); }

      // direct access to the MPOpCompressed
      data_type& data() { return Data; }
      data_type const& data() const { return Data; }

   private:
      data_type Data;
};

PStream::opstream& operator<<(PStream::opstream& out, LinearOperator const& op);
PStream::ipstream& operator>>(PStream::ipstream& in, LinearOperator& op);

LinearOperator& operator*=(LinearOperator& x, double a);
LinearOperator& operator*=(LinearOperator& x, std::complex<double> a);

LinearOperator& operator+=(LinearOperator& x, LinearOperator const& y);
LinearOperator& operator-=(LinearOperator& x, LinearOperator const& y);

LinearOperator operator+(LinearOperator const& x, LinearOperator const& y);
LinearOperator operator-(LinearOperator const& x, LinearOperator const& y);

LinearOperator operator-(LinearOperator const& x);

LinearOperator operator*(double a, LinearOperator const& x);
LinearOperator operator*(LinearOperator const& x, double a);
LinearOperator operator*(std::complex<double> a, LinearOperator const& x);
LinearOperator operator*(LinearOperator const& x, std::complex<double> a);

LinearOperator prod(LinearOperator const& x, LinearOperator const& y, QuantumNumbers::QuantumNumber const& q);
LinearOperator prod(LinearOperator const& x, LinearOperator const& y);
LinearOperator operator*(LinearOperator const& x, LinearOperator const& y);

// dot product - takes into account the multiplicity to rescale the result
LinearOperator dot(LinearOperator const& x, LinearOperator const& y);

// project a (reducible) quantum number onto an irreducible component
LinearOperator project(LinearOperator const& x, QuantumNumbers::QuantumNumber const& q);

// power of an operator.  Requires n > 1.  Only useful for n small!
LinearOperator pow(LinearOperator const& x, int n);

// Conjugate
LinearOperator conj(LinearOperator const& x);

// Adjoint
LinearOperator adjoint(LinearOperator const& x);

// output to a stream
inline
std::ostream& operator<<(std::ostream& out, LinearOperator const& x)
{
   for (LinearOperator::const_iterator I = x.begin(); I != x.end(); ++I)
   {
      out << *I << '\n';
   }
   return out;
}

#endif
