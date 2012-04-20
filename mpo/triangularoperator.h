// -*- C++ -*- $Id$
//
// ***OBSOLETE***
//
// infinite matrix product operator.  This is the old version, that will be replaced
// by the triangular_operator.h header later.
//

#if !defined(TRIANGULAROPERATOR_H_JDCHJKEHY589758YUER89H489)
#define TRIANGULAROPERATOR_H_JDCHJKEHY589758YUER89H489

#include "pstream/pstream.h"
#include "quantumnumbers/quantumnumber.h"
#include "operator_component.h"
#include "siteoperator/siteoperator.h"

class MpOpTriangular
{
   public:
      typedef OperatorComponent::basis1_type basis1_type;
      typedef OperatorComponent::basis2_type basis2_type;
      typedef OperatorComponent::value_type value_type;
  
      MpOpTriangular() {}

      MpOpTriangular(OperatorComponent const& Data);

      SymmetryList GetSymmetryList() const { return Data_.GetSymmetryList(); }
      QuantumNumber TransformsAs() const { return Data_.Basis1()[Basis1().size()-1]; }

      OperatorComponent& data() { return Data_; }
      OperatorComponent const& data() const { return Data_; }

      BasisList const& LocalBasis1() const { return Data_.LocalBasis1(); }
      BasisList const& LocalBasis2() const { return Data_.LocalBasis2(); }

      BasisList const& LocalBasis() const 
      { DEBUG_CHECK_EQUAL(Data_.LocalBasis1(), Data_.LocalBasis2()); return Data_.LocalBasis1(); }

      // This is a one-site operator, so we require that Basis1() == Basis2()
      basis1_type Basis1() const { return Data_.Basis1(); }
      basis2_type Basis2() const { return Data_.Basis2(); }

      value_type operator()(int i, int j) const
      { return Data_(i,j); }

      value_type& operator()(int i, int j)
      { return Data_(i,j); }

   private:
      OperatorComponent Data_;
};

std::ostream&
operator<<(std::ostream& out, MpOpTriangular const& op);

MpOpTriangular TriangularOneSite(SiteOperator const& x);

MpOpTriangular TriangularTwoSite(SiteOperator const& x, SiteOperator const& y,
				 QuantumNumbers::QuantumNumber const& Trans);

MpOpTriangular TriangularTwoSite(SiteOperator const& x, SiteOperator const& y);

MpOpTriangular TriangularThreeSite(SiteOperator const& x, 
                                   SiteOperator const& y, 
                                   SiteOperator const& z);

MpOpTriangular TriangularFourSite(SiteOperator const& w, 
                                  SiteOperator const& x, 
                                  SiteOperator const& y, 
                                  SiteOperator const& z);

// a two-site longer-range interaction, with NumNeighbors >= 0 inserted in between
MpOpTriangular TriangularStretchedTwoSite(SiteOperator const& x, int NumNeighbors,
					  SiteOperator const& y);

// a hack to get a zig-zag chain
MpOpTriangular ZigZagChain(SiteOperator const& S, SiteOperator const& T, double J1, double J2);

MpOpTriangular TriangularTwoSitePBC(SiteOperator const& x, SiteOperator const& y,
                                    QuantumNumbers::QuantumNumber const& Trans);

MpOpTriangular TriangularTwoSitePBC(SiteOperator const& x, SiteOperator const& y);

MpOpTriangular TriangularTwoSitePBC_Boundary(SiteOperator const& x, SiteOperator const& y,
                                    QuantumNumbers::QuantumNumber const& Trans);

MpOpTriangular TriangularTwoSitePBC_Boundary(SiteOperator const& x, SiteOperator const& y);

// arithmetic

// multiply by constant
MpOpTriangular& operator*=(MpOpTriangular& Op, double x);
MpOpTriangular& operator*=(MpOpTriangular& Op, std::complex<double> x);
MpOpTriangular operator*(MpOpTriangular const& Op, double x);
MpOpTriangular operator*(MpOpTriangular const& Op, std::complex<double> x);
MpOpTriangular operator*(double x, MpOpTriangular const& Op);
MpOpTriangular operator*(std::complex<double> x, MpOpTriangular const& Op);

// multiplication
MpOpTriangular operator*(MpOpTriangular const& x, MpOpTriangular const& y);

// addition
MpOpTriangular operator+(MpOpTriangular const& x, MpOpTriangular const& y);
MpOpTriangular& operator+=(MpOpTriangular& x, MpOpTriangular const& y);

// subtraction
MpOpTriangular operator-(MpOpTriangular const& x, MpOpTriangular const& y);

inline
MpOpTriangular local_tensor_prod(MpOpTriangular const& x, MpOpTriangular const& y)
{
   return MpOpTriangular(local_tensor_prod(x.data(), y.data()));
}

// Get the initial (1x1) E and F matrices.  These only make sense
// if the operator is a scalar.
StateComponent Initial_E(MpOpTriangular const& m);
StateComponent Initial_F(MpOpTriangular const& m);

// initial matrices for a given vector basis
StateComponent Initial_E(MpOpTriangular const& m, VectorBasis const& B);
StateComponent Initial_F(MpOpTriangular const& m, VectorBasis const& B);

StateComponent Initial_E(OperatorComponent const& m, VectorBasis const& B);
StateComponent Initial_F(OperatorComponent const& m, VectorBasis const& B);

// returns a MpOpTriangular consisting of the sub-matrix (Offset,Offset+n), (Offset,Offset+n)
MpOpTriangular Section(MpOpTriangular const& x, int Offset, int n);

// returns a MpOpTriangular consisting of the top-left n x n section of x
MpOpTriangular SectionTopLeft(MpOpTriangular const& x, int n);

// returns a MpOpTriangular consisting of the bottom-left n x n section of x
MpOpTriangular SectionBottomRight(MpOpTriangular const& x, int n);

typedef std::deque<OperatorComponent> SimpleMPOperator;

// extracts the local basis from the unit cell of an MPOperator
std::vector<BasisList>
ExtractLocalBasis(SimpleMPOperator const& Op);

#endif
