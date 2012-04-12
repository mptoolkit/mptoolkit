
#include <iomanip>

//
// SimpleBasis
//

inline
int
SimpleBasis::Append(QuantumNumber const& q)
{
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), this->GetSymmetryList());
   BasisList.push_back(q);
   return BasisList.size() - 1;
}

inline
bool operator==(SimpleBasis const& b1, SimpleBasis const& b2)
{
   return b1.SList == b2.SList && b1.BasisList == b2.BasisList;
}

inline
bool operator!=(SimpleBasis const& b1, SimpleBasis const& b2)
{
   return b1.SList != b2.SList || b1.BasisList != b2.BasisList;
}

//
// VectorBasisImpl
//

inline
int
VectorBasisImpl::Append(QuantumNumber const& q, int Dimension)
{
   DEBUG_PRECONDITION_EQUAL(q.GetSymmetryList(), this->GetSymmetryList());
   QList.push_back(q);
   DimensionList.push_back(Dimension);
   DEBUG_CHECK_EQUAL(QList.size(), DimensionList.size());
   LinDimension += Dimension;
   return QList.size()-1;
}

//
// VectorBasis
//

inline
PStream::opstream& operator<<(PStream::opstream& out, VectorBasis const& B)
{
  return out << B.pImpl;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, VectorBasis& B)
{
  return in >> B.pImpl;
}
