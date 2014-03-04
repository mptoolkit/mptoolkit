// -*- C++ -*- $Id$
//
// A LatticeOperator is any MPO that
// exists on some lattice.  It knows how it commutes (ie, the Jordan-Wigner string)
// so we can adapt the operator to any multiple of the unit cell.

#if !defined(LATTICE_OPERATOR_H_JHD283YHL)
#define LATTICE_OPERATOR_H_JHD283YHL

#include "siteoperator/unitcell.h"
#include "finite_mpo.h"
#include "triangular_mpo.h"
#include "product_operator.h"
#include "pheap/pvalueptr.h"
//#include "periodic_mpo.h"


enum OperatorType { TypeFiniteMPO, TypeTriangularMPO };

class LatticeOperator
{
   public:

      // returns the basic type of the operator (finite, triangular, periodic)
      OperatorType Type() const;

      // returns the number of sites that this operator is defined over.
      // This is always a multiple of the unit cell size.
      int size() const;

      // returns the Jordan-Wigner string associated with this operator
      std::string JordanWignerString() const;

      // Returns a finite version of the operator.  If it is
      // a triangular operator, then close off the boundaries
      // to make it finite on the unit cell.
      // PRECONDITION: operator is finite or triangular
      FiniteMPO AsFiniteMPO() const;

      // returns the UnitCell of this operator
      UnitCell const& GetUnitCell() const { return *pUnitCell; }

      // shorthand for GetUnitCell().size()
      int UnitCellSize() const { return pUnitCell->size(); }

      // Extend an operator to an integral multiple of the unit cell size,
      // offset from the left by some number of sites.
      // PRECONDITION: Size % this->UnitCellSize() == 0 
      // && Offset % this->UnitCellSize() == 0
      // && Size+Offset >= this->size()
      FiniteMPO AsFiniteMPO(int Size, int Offset) const;

      // Returns a triangular version of the operator.  If it is
      // a finite operator, then convert it into a sum of unit cells
      TriangularMPO AsTriangularMPO() const;

      // Returns a triangular version of the operator.  If it is
      // a finite operator, then convert it into a sum across the unit cell,
      // with momentum k per size.  If it is already triangular, then
      // increase the momentum by k
      TriangularMPO AsTriangularMPOMomentum(double k) const;

      // Returns a triangular version of the operator.  If it is
      // a finite operator, then convert it into a sum across the unit cell,
      // with exponential decrease lambda per size.  
      // If it is already triangular, then multiply by an exponential decay lambda
      TriangularMPO AsTriangularMPODecay(double lambda) const;

      // Returns a triangular version of the operator.  If it is
      // a finite operator, then convert it into a sum across the unit cell,
      // with multiplicative factor f per size.
      TriangularMPO AsTriangularMPO(std::complex<double> f) const;

   private:
      pvalue_ptr<UnitCell> pUnitCell;
      std::string JWString;
      boost::variant<FiniteMPO, TriangularMPO> Operator;
};

LatticeOperator MakeFinite(LatticeOperator const& Op);
LatticeOperator MakeFinite(LatticeOperator const& Op, int Size, int Offset);

LatticeOperator SumUnit(LatticeOperator const& Op);
LatticeOperator SumUnitFactor(LatticeOperator const& Op, std::complex<double> f);

LatticeOperator LocalProduct(LatticeOperator const& A, LatticeOperator const& B);

LatticeOperator operator+(LatticeOperator const& A, LatticeOperator const& B);

LatticeOperator operator*(LatticeOperator const& A, LatticeOperator const& B);

LatticeOperator operator*(double x, LatticeOperator const& Op);

LatticeOperator operator*(std::complex<double> x, LatticeOperator const& Op);


#endif

