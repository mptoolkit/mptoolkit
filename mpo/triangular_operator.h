// -*- C++ -*- $Id$
//
// TriangularOperator: a representation for lattice operators that are in triangular form
// This is the new version, not yet completed.  When it is completed, it will obsolete
// triangularoperator.h

// Discarded design idea: keep track of the diagonal components separately.  Doesn't really work, too hard to combine
// operators.

// This is a refinement of an MPOperator that is lower-triangular.
// It is up to the user to ensure that the TriangularOperator stays in lower-triagular form.
// All functions defined in this header are OK though, the only way to generate a non-lower-triangular
// operator is to modify the components by hand.

#if !defined(TRIANGULAR_OPERATOR_H_SDJFHU47894789H9O43)
#define TRIANGULAR_OPERATOR_H_SDJFHU47894789H9O43

#include "mpoperator.h"
#include <ostream>

// represents a lower triangular operator
class TriangularOperator
{
   private:
      typedef MPOperator DataType;

   public:
      typedef OperatorComponent value_type;
      typedef DataType::iterator iterator;
      typedef DataType::const_iterator const_iterator;
      typedef value_type::basis1_type basis_type;

      TriangularOperator() {}

      explicit TriangularOperator(int Size) : Data_(Size) {}

      // construction as a single-site operator
      explicit TriangularOperator(value_type const& x) : Data_(1, x) {}

      explicit TriangularOperator(std::vector<value_type> const& x) : Data_(x.begin(), x.end()) {}

      std::size_t size() const { return Data_.size(); }

      bool empty() const { return Data_.empty(); }

      iterator begin() { return Data_.begin(); }
      iterator end() { return Data_.end(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }
   
      // The Basis1() and Basis2() always coincide for a TriangularOperator
      basis_type Basis() const { return Data_.front().Basis1(); }
      basis_type Basis1() const { return Data_.front().Basis1(); }
      basis_type Basis2() const { return Data_.back().Basis2(); }

      value_type const& operator[](int n) const { return Data_[n]; }
      value_type& operator[](int n) { return Data_[n]; }

      value_type const& front() const { return Data_.front(); }
      value_type& front() { return Data_.front(); }

      value_type const& back() const { return Data_.back(); }
      value_type& back() { return Data_.back(); }

      QuantumNumber const& TransformsAs() const { return this->Basis().front(); }

      // returns the component at entry (i,j).  Result is a 1x1 MPOperator
      MPOperator operator()(int i, int j) const;

      // extracts a diagonal operator.  This is either null, or a product of SimpleRedOperator's.
      std::vector<SimpleRedOperator> diagonal(int i) const;

     MPOperator const& data() const { return Data_; }

     MPOperator& data() { return Data_; }  // dangerous, probably shouldn't exist

      QuantumNumbers::SymmetryList GetSymmetryList() const { return Data_.GetSymmetryList(); }

   private:
      DataType Data_;
};

std::ostream&
operator<<(std::ostream& out, TriangularOperator const& op);

// extracts a single column from a triangular operator.  Result is an Nx1 row-vector operator
MPOperator extract_column(TriangularOperator const& Op, int Col);

// extracts a single column from a triangular operator, excluding the diagonal.  Result is an Nx1 row-vector operator
MPOperator extract_lower_column(TriangularOperator const& Op, int Col);

void mask_lower_column(TriangularOperator const& Op, int Col, std::vector<std::vector<int> >& Mask);

TriangularOperator TriangularOneSite(SimpleOperator const& x);

// A one-site operator with the given momentum, in angular units
TriangularOperator TriangularOneSite(SimpleOperator const& x, double Momentum);

TriangularOperator TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y,
                                     QuantumNumbers::QuantumNumber const& Trans);

TriangularOperator TriangularTwoSite(SimpleOperator const& x, SimpleOperator const& y);

TriangularOperator TriangularThreeSite(SimpleOperator const& x, 
                                       SimpleOperator const& y, 
                                       SimpleOperator const& z);

// a two point interaction between sites n1 and n2 of a lattice.
// This probably should use SiteOperator's rather than SimpleOperator's, so that
// we can handle fermions, especially for the case n2<n1.
// this is assumed to be a bosonic operator; if n2 < n1 we swap the sites
// (which is incorrect for fermions)
TriangularOperator TwoPointOperator(std::vector<BasisList> const& Sites, 
                                    int n1, SimpleOperator const& x1,
                                    int n2, SimpleOperator const& x2);

// a two-point string operator where the String term is inserted
// at sites n1+1, n1+2, ..., n2-1.
// Because this function can be used to implement fermionic operators,
// we demand normal ordering of sites; it is an error to call this function
// with n2<n1.
TriangularOperator TwoPointStringOperator(std::vector<BasisList> const& Sites, 
					  int n1, SimpleOperator const& x1,
					  SimpleOperator const& String,
					  int n2, SimpleOperator const& x2);

// A one-site operator on a lattice
TriangularOperator OnePointOperator(std::vector<BasisList> const& Sites, 
                                    int n, SimpleOperator const& x);

// replicate an operator on a unit cell this many times
TriangularOperator repeat(TriangularOperator const& x, int Count);

// Two triangular operators are *compatible* if they have the same operator on the
// top-left and bottom-right entries.
bool is_compatible(TriangularOperator const& x, TriangularOperator const& y);

// Multiply by scalar
TriangularOperator operator*(TriangularOperator const& Op, double x);
TriangularOperator operator*(TriangularOperator const& Op, std::complex<double> x);
TriangularOperator operator*(double x, TriangularOperator const& Op);
TriangularOperator operator*(std::complex<double> x, TriangularOperator const& Op);

TriangularOperator& operator*=(TriangularOperator& Op, double x);
TriangularOperator& operator*=(TriangularOperator& Op, std::complex<double> x);

// Addition of triangular operators.  This is only possible if the operators
// are compatible.
TriangularOperator& operator+=(TriangularOperator& Op, TriangularOperator const& x);
TriangularOperator& operator-=(TriangularOperator& Op, TriangularOperator const& x);

TriangularOperator operator+(TriangularOperator const& x, TriangularOperator const& y);
TriangularOperator operator-(TriangularOperator const& x, TriangularOperator const& y);

// extends an operator by repeating it count times
TriangularOperator extend(TriangularOperator const& x, int count);

// does a 2-1 coarse-graining of the operator, which must have an even size
TriangularOperator coarse_grain(TriangularOperator const& x);

// Multiplication of triangular MPO's.  This doesn't depend on the
// compatibility of the operators.
TriangularOperator operator*(TriangularOperator const& x, TriangularOperator const& y);

TriangularOperator& operator*=(TriangularOperator& x, TriangularOperator const& y);

// Get the initial (1x1) E and F matrices.  These only make sense
// if the operator is a scalar.
StateComponent Initial_E(TriangularOperator const& m);
StateComponent Initial_F(TriangularOperator const& m);

// initial matrices for a given vector basis
StateComponent Initial_E(TriangularOperator const& m, VectorBasis const& B);
StateComponent Initial_F(TriangularOperator const& m, VectorBasis const& B);

#endif
