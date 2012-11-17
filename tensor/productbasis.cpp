// -*- C++ -*- $Id$

#include "productbasis.h"

//
// ProductBasis<SimpleBasis, SimpleBasis>
//

ProductBasis<SimpleBasis, SimpleBasis>::ProductBasis(SimpleBasis Basis1_, SimpleBasis Basis2_)
  : SimpleBasis(Basis1_.GetSymmetryList()), B1(Basis1_), B2(Basis2_), TransformData(B1.size(), B2.size())
{
   DEBUG_PRECONDITION(B1.GetSymmetryList() == B2.GetSymmetryList())
                     (B1.GetSymmetryList())
                     (B2.GetSymmetryList());

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
	    int s = this->Append(*I);
	    DEBUG_CHECK(RMap.size() == s)(RMap.size())(s);
	    RMap.push_back(source_type(s1, s2));
	    TransformData(s1, s2).push_back(s);
	 }
      }
   }
   ReverseMapping = SourceListType(RMap.begin(), RMap.end());
}

ProductBasis<SimpleBasis, SimpleBasis>::ProductBasis(SimpleBasis Basis1_, SimpleBasis Basis2_, 
						     QuantumNumber const& Target)
  : SimpleBasis(Basis1_.GetSymmetryList()), B1(Basis1_), B2(Basis2_), TransformData(B1.size(), B2.size())
{
   DEBUG_PRECONDITION(B1.GetSymmetryList() == B2.GetSymmetryList())
                     (B1.GetSymmetryList())
                     (B2.GetSymmetryList());
   DEBUG_PRECONDITION(B1.GetSymmetryList() == Target.GetSymmetryList())
                     (B1.GetSymmetryList())
                     (Target.GetSymmetryList());

   // construct the reverse mapping into a std::vector, so we can use push_back.
   // copy it later into a LinearAlgebra::Vector so we have reference counting.
   std::vector<source_type> RMap;

   for (int s1 = 0; s1 < B1.NumSubspaces(); ++s1)
   {
      for (int s2 = 0; s2 < B2.NumSubspaces(); ++s2)
      {
	 if (is_transform_target(B1.qn(s1), B2.qn(s2), Target))
	 {
	    int s = this->Append(Target);
	    DEBUG_CHECK(RMap.size() == s)(RMap.size())(s);
	    RMap.push_back(source_type(s1, s2));
	    TransformData(s1, s2).push_back(s);
	 }
      }
   }
   ReverseMapping = SourceListType(RMap.begin(), RMap.end());
}

ProductBasis<SimpleBasis, SimpleBasis>
ProductBasis<SimpleBasis, SimpleBasis>::MakeTriangularProjected(LeftBasisType Basis1_, 
								RightBasisType Basis2_, 
								QuantumNumber const& q)
{
   ProductBasis<SimpleBasis, SimpleBasis> Result = ProductBasis(Basis1_.GetSymmetryList());

   Result.B1 = Basis1_;
   Result.B2 = Basis2_;

   // construct the reverse mapping into a std::vector, so we can use push_back.
   // copy it later into a LinearAlgebra::Vector so we have reference counting.
   std::vector<source_type> RMap;

   // Add all targets except for the final quantum number in B1
   for (int s1 = 0; s1 < B1.NumSubspaces()-1; ++s1)
   {
      for (int s2 = 0; s2 < B2.NumSubspaces(); ++s2)
      {
	 std::list<QuantumNumber> QNList;
	 transform_targets(B1.qn(s1), B2.qn(s2), std::back_inserter(QNList));

	 for (std::list<QuantumNumber>::iterator I = QNList.begin(); I != QNList.end(); ++I)
	 {
	    int s = this->Append(*I);
	    DEBUG_CHECK(RMap.size() == s)(RMap.size())(s);
	    RMap.push_back(source_type(s1, s2));
	    TransformData(s1, s2).push_back(s);
	 }
      }
   }
   // Add the remaining targets for the last B1, up to (but not including) the last B2
   int s1 = B1.NumSubspaces()-1;
   for (int s2 = 0; s2 < B2.NumSubspaces()-1; ++s2)
   {
      std::list<QuantumNumber> QNList;
      transform_targets(B1.qn(s1), B2.qn(s2), std::back_inserter(QNList));
      
      for (std::list<QuantumNumber>::iterator I = QNList.begin(); I != QNList.end(); ++I)
      {
	 int s = this->Append(*I);
	 DEBUG_CHECK(RMap.size() == s)(RMap.size())(s);
	 RMap.push_back(source_type(s1, s2));
	 TransformData(s1, s2).push_back(s);
      }
   }

   // Finally, for the last B1 and last B2 add ONLY the quantum number q, which must be possible
   int s2 = B2.NumSubspaces()-1;
   CHECK(is_transform_target(B1.qn(s1), B2.qn(s2), q));
   int s = this->Append(Target);
   DEBUG_CHECK(RMap.size() == s)(RMap.size())(s);
   RMap.push_back(source_type(s1, s2));
   TransformData(s1, s2).push_back(s);

   ReverseMapping = SourceListType(RMap.begin(), RMap.end());
}
