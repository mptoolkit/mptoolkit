// -*- C++ -*- $Id$

namespace Tensor
{

template <typename BLType, typename BRType>
PStream::opstream& operator<<(PStream::opstream& out, ProductBasis<BLType, BRType> const& B)
{
   return out << B.Basis() << B.B1 << B.B2 << B.TransformData << B.ReverseMapping;
}

template <typename BLType, typename BRType>
PStream::ipstream& operator>>(PStream::ipstream& in, ProductBasis<BLType, BRType>& B)
{
   return in >> B.Basis() >> B.B1 >> B.B2 >> B.TransformData >> B.ReverseMapping;
   return in;
}

template <typename BLType, typename BRType>
ProductBasis<BLType, BRType>::ProductBasis(LeftBasisType Basis1_, RightBasisType Basis2_)
  : VectorBasis(Basis1_.GetSymmetryList()), B1(Basis1_), B2(Basis2_), TransformData(B1.size(), B2.size())
{
   DEBUG_PRECONDITION(B1.GetSymmetryList() == B2.GetSymmetryList())
                     (B1.GetSymmetryList())
                     (B2.GetSymmetryList());

   // get a lock on the basis type so we can modify it fast
   VectorBasis::LockType BasisLock(this->Lock());

   // construct the reverse mapping into a std::vector, so we can use push_back.
   // copy it later into a LinearAlgebra::Vector so we have reference counting.
   std::vector<source_type> RMap;

   for (int s1 = 0; s1 < B1.NumSubspaces(); ++s1)
   {
      for (int s2 = 0; s2 < B2.NumSubspaces(); ++s2)
      {
	 std::list<QuantumNumber> QNList;
	 transform_targets(B1.qn(s1), B2.qn(s2), std::back_inserter(QNList));

	 for (std::list<QuantumNumber>::iterator I = QNList.begin(); I != QNList.end(); ++I)
	 {
	    int s = BasisLock->Append(*I, B1.Dimension(s1) * B2.Dimension(s2));
	    DEBUG_CHECK(RMap.size() == s)(RMap.size())(s);
	    RMap.push_back(source_type(s1, s2));
	    TransformData(s1, s2).push_back(s);
	 }
      }
   }
   ReverseMapping = SourceListType(RMap.begin(), RMap.end());
}

template <typename BLType, typename BRType>
ProductBasis<BLType, BRType>::ProductBasis(LeftBasisType Basis1_, RightBasisType Basis2_, QuantumNumber const& Target)
  : VectorBasis(Basis1_.GetSymmetryList()), B1(Basis1_), B2(Basis2_), TransformData(B1.size(), B2.size())
{
   DEBUG_PRECONDITION(B1.GetSymmetryList() == B2.GetSymmetryList())
                     (B1.GetSymmetryList())
                     (B2.GetSymmetryList());
   DEBUG_PRECONDITION(B1.GetSymmetryList() == Target.GetSymmetryList())
                     (B1.GetSymmetryList())
                     (Target.GetSymmetryList());

   // get a lock on the basis type so we can modify it fast
   VectorBasis::LockType BasisLock(this->Lock());

   // construct the reverse mapping into a std::vector, so we can use push_back.
   // copy it later into a LinearAlgebra::Vector so we have reference counting.
   std::vector<source_type> RMap;

   for (int s1 = 0; s1 < B1.NumSubspaces(); ++s1)
   {
      for (int s2 = 0; s2 < B2.NumSubspaces(); ++s2)
      {
	 if (is_transform_target(B1.qn(s1), B2.qn(s2), Target))
	 {
	    int s = BasisLock->Append(Target, B1.Dimension(s1) * B2.Dimension(s2));
	    DEBUG_CHECK(RMap.size() == s)(RMap.size())(s);
	    RMap.push_back(source_type(s1, s2));
	    TransformData(s1, s2).push_back(s);
	 }
      }
   }
   ReverseMapping = SourceListType(RMap.begin(), RMap.end());
}

} // namespace Tensor
