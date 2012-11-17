// -*- C++ -*- $Id$

namespace Tensor
{

template <typename BL, typename BR>
PStream::opstream& operator<<(PStream::opstream& out, ProductBasis<BL, BR> const& B)
{
   return out << B.Basis_ << B.Left_ << B.Right_ << B.TransformData_ << B.ReverseMapping_;
}

template <typename BL, typename BR>
PStream::ipstream& operator>>(PStream::ipstream& in, ProductBasis<BL, BR>& B)
{
   return in >> B.Basis_ >> B.Left_ >> B.Right_ >> B.TransformData_ >> B.ReverseMapping_;
}

template <typename BL, typename BR, typename Base>
ProductBasisBase<BL, BR, Base>::ProductBasisBase(left_basis_type const& Left, 
                                                 right_basis_type const& Right)
   : Basis_(Left.GetSymmetryList()), Left_(Left), Right_(Right), 
     TransformData_(Left_.size(), Right_.size())
{
   DEBUG_PRECONDITION_EQUAL(Left.GetSymmetryList(), Right.GetSymmetryList());

   // construct the reverse mapping into a std::vector, so we can use push_back.
   // copy it later into a LinearAlgebra::Vector so we have reference counting.
   std::vector<source_type> RMap;

   for (std::size_t s1 = 0; s1 < Left_.size(); ++s1)
   {
      for (std::size_t s2 = 0; s2 < Right_.size(); ++s2)
      {
	 std::list<QuantumNumber> QNList;
	 transform_targets(Left_[s1], Right_[s2], std::back_inserter(QNList));

	 for (std::list<QuantumNumber>::iterator I = QNList.begin(); I != QNList.end(); ++I)
	 {
            int s = Basis_.size();
	    Basis_.push_back(*I);
	    RMap.push_back(source_type(s1, s2));
	    TransformData_(s1, s2).push_back(s);
	 }
      }
   }
   ReverseMapping_ = SourceListType(RMap.begin(), RMap.end());
}

template <typename BL, typename BR, typename Base>
ProductBasisBase<BL, BR, Base>::ProductBasisBase(left_basis_type const& Left, 
                                                 right_basis_type const& Right,
                                      QuantumNumber const& Target)
   : Basis_(Left.GetSymmetryList()), Left_(Left), Right_(Right), 
     TransformData_(Left_.size(), Right_.size())
{
   DEBUG_PRECONDITION_EQUAL(Left.GetSymmetryList(), Right.GetSymmetryList());

   // construct the reverse mapping into a std::vector, so we can use push_back.
   // copy it later into a LinearAlgebra::Vector so we have reference counting.
   std::vector<source_type> RMap;

   for (std::size_t s1 = 0; s1 < Left_.size(); ++s1)
   {
      for (std::size_t s2 = 0; s2 < Right_.size(); ++s2)
      {
	 if (is_transform_target(Left_[s1], Right_[s2], Target))
	 {
            int s = Basis_.size();
	    Basis_.push_back(Target);
	    RMap.push_back(source_type(s1, s2));
	    TransformData_(s1, s2).push_back(s);
	 }
      }
   }
   ReverseMapping_ = SourceListType(RMap.begin(), RMap.end());
}

template <typename BL, typename BR, typename Base>
ProductBasisBase<BL, BR, Base>
ProductBasisBase<BL, BR, Base>::MakeTriangularProjected(left_basis_type const& Left_, 
							right_basis_type const& Right_,
							QuantumNumbers::QuantumNumber const& Target)
{
   DEBUG_PRECONDITION_EQUAL(Left_.GetSymmetryList(), Right_.GetSymmetryList());

   ProductBasisBase<BL, BR, Base> Result;
   Result.Basis_ = basis_type(Left_.GetSymmetryList());
   Result.Left_ = Left_;
   Result.Right_ = Right_;
   Result.TransformData_ = LinearAlgebra::Matrix<TargetListType>(Left_.size(), Right_.size());

   // construct the reverse mapping into a std::vector, so we can use push_back.
   // copy it later into a LinearAlgebra::Vector so we have reference counting.
   std::vector<source_type> RMap;

   for (std::size_t s1 = 0; s1 < Left_.size()-1; ++s1)
   {
      for (std::size_t s2 = 0; s2 < Right_.size(); ++s2)
      {
	 std::list<QuantumNumber> QNList;
	 transform_targets(Left_[s1], Right_[s2], std::back_inserter(QNList));

	 for (std::list<QuantumNumber>::iterator I = QNList.begin(); I != QNList.end(); ++I)
	 {
            int s = Result.Basis_.size();
	    Result.Basis_.push_back(*I);
	    RMap.push_back(source_type(s1, s2));
	    Result.TransformData_(s1, s2).push_back(s);
	 }
      }
   }

   std::size_t s1 = Left_.size()-1;
   for (std::size_t s2 = 0; s2 < Right_.size()-1; ++s2)
   {
      std::list<QuantumNumber> QNList;
      transform_targets(Left_[s1], Right_[s2], std::back_inserter(QNList));

      for (std::list<QuantumNumber>::iterator I = QNList.begin(); I != QNList.end(); ++I)
      {
	 int s = Result.Basis_.size();
	 Result.Basis_.push_back(*I);
	 RMap.push_back(source_type(s1, s2));
	 Result.TransformData_(s1, s2).push_back(s);
      }
   }
   std::size_t s2 = Right_.size()-1;
   CHECK(is_transform_target(Left_[s1], Right_[s2], Target));
   int s = Result.Basis_.size();
   Result.Basis_.push_back(Target);
   RMap.push_back(source_type(s1, s2));
   Result.TransformData_(s1, s2).push_back(s);

   Result.ReverseMapping_ = SourceListType(RMap.begin(), RMap.end());

   return Result;
}

} // namespace Tensor
