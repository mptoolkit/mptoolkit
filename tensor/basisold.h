/*

   basis.h, version 3

   Modified from version 2, 2000-09-05 Ian McCulloch

   TODO: the naming of some functions is a bit confusing.. this is hopefully fixed

   Modified 2003-10-12, Ian McCulloch: removed the 'extra' label.  changed the QuantumNumberType
   to a template parameter.  Hopefully this covers and possible use of the 'extra' label.

*/

#if !defined(BASIS_H_567FG65K342JH7FG687GH768)
#define BASIS_H_567FG65K342JH7FG687GH768

#include "common/trace.h"
//#include "common/math_const.h"
#include "pstream/pstream.h"
#include "pheap/pvalueptr.h"
#include "pheap/prefptr.h"
#include "quantumnumbers/quantumnumber.h"
#include "quantumnumbers/symmetrylist.h"
#include <vector>
#include <iostream>

using QuantumNumbers::SymmetryList;
using QuantumNumbers::QuantumNumber;

class SimpleBasis
{
   public:
      SimpleBasis() {}
      explicit SimpleBasis(SymmetryList const& SList_) : SList(SList_) {}

      int Append(QuantumNumber const& Q_);

      // returns the number of states in the basis
      int size() const { return BasisList.size(); }

      // returns true if this basis is empty
      bool is_null() const { return this->size() == 0; }

      // for a SimpleBasis, the number of subspaces is the same as the size
      int NumSubspaces() const { return this->size(); }

      // for SimpleBasis, the Dimension is always 1
      int Dimension(int s) const { DEBUG_RANGE_CHECK_OPEN(s, 0, this->NumSubspaces()); return 1; }

      QuantumNumber const& qn(int Subspace) const { return BasisList[Subspace]; }

      SymmetryList const& GetSymmetryList() const { return SList; }

   private:
      SymmetryList SList;
      QuantumNumbers::QuantumNumberList BasisList;

   friend PStream::opstream& operator<<(PStream::opstream& out, SimpleBasis const& Basis);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, SimpleBasis& Basis);

   friend bool operator==(SimpleBasis const& b1, SimpleBasis const& b2);
   friend bool operator!=(SimpleBasis const& b1, SimpleBasis const& b2);
};

// generates a 'vacuum basis' from the given symmetry list - that is,
// a basis with a single state that transforms as a scalar.
inline
SimpleBasis MakeVacuumBasis(SymmetryList const& SList)
{
   SimpleBasis Result(SList);
   Result.Append(QuantumNumber(SList));
   return Result;
}

SimpleBasis adjoint(SimpleBasis const& S);

class VectorBasisImpl
{
   public:
      VectorBasisImpl() {}  // needed for streaming

      explicit VectorBasisImpl(SymmetryList const& SList_) : SList(SList_), LinDimension(0) {}

      // appends a list of states with the given quantum number and dimension.  Returns the subspace number.
      int Append(QuantumNumber const& Q_, int Dimension_);

      // returns the number of subspaces in the basis
      int NumSubspaces() const { return QList.size(); }

      int size() const { return QList.size(); }

      // returns the quantum number of the given subspace
      QuantumNumber const& qn(int Subspace) const
      {return QList[Subspace]; }

      // returns the dimension of the given subspace
      int Dimension(int Subspace) const
      { DEBUG_RANGE_CHECK(Subspace, 0, NumSubspaces()-1); return DimensionList[Subspace]; }

      // returns the full dimension of the basis
      int LinearDimension() const { return LinDimension; }

      SymmetryList const& GetSymmetryList() const { return SList; }

      bool operator==(VectorBasisImpl const& Other) const
      { return LinDimension == Other.LinDimension && SList == Other.SList && 
           DimensionList == Other.DimensionList &&
	   QList == Other.QList; }

   private:
      SymmetryList SList;

      // Each subspace is a pair (QuantumNumber, Dimension)
      // but stored separately, so we can use a QuantumNumberList
      QuantumNumbers::QuantumNumberList QList;
      std::vector<int> DimensionList;
      int LinDimension;   // linear dimension, calculated by dead reckening from Append()

   friend PStream::opstream& operator<<(PStream::opstream& out, VectorBasisImpl const& Basis);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, VectorBasisImpl& Basis);
};

std::ostream& operator<<(std::ostream& out, SimpleBasis const& B);

//
// VectorBasis
//

class VectorBasis
{
   public:
      typedef pvalue_ptr<VectorBasisImpl>         pImplType;
      typedef pImplType::lock_type                LockType;

      VectorBasis() : pImpl(new VectorBasisImpl()) {}

      // Constructs an empty VectorBasis object
      explicit VectorBasis(SymmetryList const& QList) : pImpl(new VectorBasisImpl(QList)) { }

      // conversion from SimpleBasis
      VectorBasis(SimpleBasis const& b);

      // appends a group of states.
      int Append(QuantumNumber const& Q, int Dimension) { return pImpl.mutate()->Append(Q, Dimension); }

      // returns the number of subspaces in the basis
      int NumSubspaces() const { return pImpl->NumSubspaces(); }

      // returns the number of subspaces in the basis
      int size() const { return pImpl->size(); }

      // returns true if this basis is empty
      bool is_null() const { return this->size() == 0; }

      // returns the quantum number of the given subspace
      QuantumNumber const& qn(int Subspace) const { return pImpl->qn(Subspace); }

      // returns the dimension of the given subspace
      int Dimension(int Subspace) const { return pImpl->Dimension(Subspace); }

      // returns the full dimension of the basis
      int LinearDimension() const { return pImpl->LinearDimension(); }

      SymmetryList const& GetSymmetryList() const { return pImpl->GetSymmetryList(); }

      bool operator==(VectorBasis const& Other) const { return pImpl == Other.pImpl || *pImpl == *Other.pImpl; }
      bool operator!=(VectorBasis const& Other) const { return pImpl != Other.pImpl && !(*pImpl == *Other.pImpl); }

      // returns a lock on the basis, so that we can do a sequence of mods without calling mutate() every time
      LockType Lock() { return pImpl.lock(); }

 private:
      pImplType pImpl;

   friend PStream::opstream& operator<<(PStream::opstream& out, VectorBasis const& Bi);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, VectorBasis& Bi);
};

std::ostream& operator<<(std::ostream& out, VectorBasis const& B);

//
// inlines
//

#include "basis.cc"

#endif
