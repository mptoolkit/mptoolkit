// -*- C++ -*-

#include "density.h"
#include "common/proccontrol.h"
#include "linearalgebra/eigen.h"

//typedef LinearAlgebra::Vector<std::pair<int, LinearAlgebra::Range> > BasisMappingType;


LinearBasis<VectorBasis>::LinearBasis(VectorBasis const& B)
   : VectorBasis(B.GetSymmetryList()), Orig(B), Mapping(B.size())
{
   // Calculate the total dimension of each quantum number in B, and
   // fill in the Range component of the Mapping.
   std::map<QuantumNumber, int> QNDimension;
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      int Dim = B.dim(i);
      QuantumNumber q = B[i];

      int& MappedDim(QNDimension[q]);
      Mapping[i].second = LinearAlgebra::Range(MappedDim, MappedDim + Dim);
      MappedDim += Dim;
   }

   // Now we know the dimensions of the quantum numbers, we can add the subspaces the basis
   std::map<QuantumNumber, int> QNSubspace;
   for (std::map<QuantumNumber, int>::const_iterator I = QNDimension.begin(); 
        I != QNDimension.end(); ++I)
   {
      QNSubspace[I->first] = this->size();
      this->push_back(I->first, I->second);
   }

   // finally, fill in the subspace component of Mapping
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Mapping[i].first = QNSubspace[B[i]];
   }
}

LinearBasis<BasisList>::LinearBasis(BasisList const& B)
   : VectorBasis(B.GetSymmetryList()), Orig(B), Mapping(B.size())
{
   // Calculate the total dimension of each quantum number in B, and
   // fill in the Range component of the Mapping.
   std::map<QuantumNumber, int> QNDimension;
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Mapping[i].second = QNDimension[B[i]]++;

   }

   // Now we know the dimensions of the quantum numbers, we can add the subspaces the basis
   std::map<QuantumNumber, int> QNSubspace;
   for (std::map<QuantumNumber, int>::const_iterator I = QNDimension.begin(); 
        I != QNDimension.end(); ++I)
   {
      QNSubspace[I->first] = this->size();
      this->push_back(I->first, I->second);
   }

   // finally, fill in the subspace component of Mapping
   for (std::size_t i = 0; i < B.size(); ++i)
   {
      Mapping[i].first = QNSubspace[B[i]];
   }
}

int
LinearBasis<BasisList>::ReverseLookup(int s, int index) const
{
   for (std::size_t i = 0; i < Mapping.size(); ++i)
   {
      if (Mapping[i].first == s && Mapping[i].second == index) return i;
   }
   return -1;
}

std::ostream& 
DensityMatrixBase::DensityMatrixReport(std::ostream& outstream, int MaxEigenvalues, bool Base2, bool ShowDegen)
{
   std::ostringstream out;
   out.precision(12);
   out << std::scientific;
   if (MaxEigenvalues < 0) MaxEigenvalues = EigenInfoList.size();
   out << "#Eigenvalue sum = " << this->EigenSum() << '\n';
   out << "#von Neumann Entropy " << (Base2 ? "(base 2)" : "(base e)") << " = " << this->Entropy() << '\n';
   if (MaxEigenvalues > 0)
      out << "#Number    #Eigenvalue         #Degen    #Weight               #Energy             #QuantumNumber\n";
   int n = 0;
   int TotalDegree = 0;
   double EShift = 0;
   double EScale = 1.0;
#if 0
   if (EigenInfoList.size() >= 2)
   {
      EShift = -log(EigenInfoList[0].Eigenvalue / this->EigenSum());
      double E2 = -log(EigenInfoList[1].Eigenvalue / this->EigenSum());
      EScale = E2 - EShift;
   }
#endif
   for (const_iterator Iter = this->begin(); Iter != this->end() && n < MaxEigenvalues; ++Iter)
   {
      double EVal = Iter->Eigenvalue / this->EigenSum();
      double Energy = EVal > 0 ? ((-log(EVal) - EShift) / EScale) : 0.0;

      int OuterDegen = ShowDegen ? Iter->Degree : 1;
      int DisplayDegen = ShowDegen ? 1 : Iter->Degree;

      for (int i = 0; i < OuterDegen; ++i)
      {
         double Weight = EVal * Iter->Degree / OuterDegen;
         TotalDegree += Iter->Degree;
         ++n;
         out << std::right << std::setw(7) << n << "  " 
             << std::right << std::setw(20) << Iter->Eigenvalue 
             << "  " << std::setw(6) << DisplayDegen
             << "  " << std::setw(20) << Weight
             << "  " << std::setw(20) << Energy
             << "  " << std::left 
             << this->Lookup(Iter->Subspace) << '\n';
      }
   }
   out << '#' << n << " out of " << (ShowDegen ? TotalDegree : EigenInfoList.size()) << " eigenvalues shown.  ";
   out << "Total degree = " << TotalDegree << '\n';
   outstream << out.str();
   return outstream;
}

double
DensityMatrixBase::Entropy(bool Base2) const
{
   double x = 0;
   for (const_iterator Iter = begin(); Iter != end(); ++Iter)
   {
      double EVal = Iter->Eigenvalue / this->EigenSum();
      if (EVal > 0)
	 x -= EVal * (Base2 ? log2(EVal) : log(EVal)) * Iter->Degree;
   }
   return x;
}

double
DensityMatrixBase::EvaluateCasimir(int n) const
{
   double x = 0;
   double ESum = this->EigenSum();
   for (const_iterator Iter = begin(); Iter != end(); ++Iter)
   {
      x += (Iter->Eigenvalue / ESum) * Iter->Degree * casimir(this->Lookup(Iter->Subspace), n);
   }
   return x;
}

void DensityMatrixBase::DiagonalizeDMHelper(bool Sort)
{
   //   double Point1 = ProcControl::GetCPUTime();  // Point1 is between construction & diagonalization

   ESum = 0;  // running sum of the eigenvalues
   // diagonalize the DM
   LinearAlgebra::Vector<double> Eigenvalues(MaxLinearDimension);
   for (std::size_t q1 = 0; q1 < RawDMList.size(); ++q1)
   {
      int CurrentDegree = degree(this->Lookup(q1));
      //      std::cout << "Raw DM is\n" << RawDMList[q1] << std::endl;
      //      TRACE(RawDMList[q1])(RawDMList[q1].size1())(RawDMList[q1].size2());
      Eigenvalues = DiagonalizeHermitian(RawDMList[q1]);

      // add the eigenvalues and eigenvector pointers to EigenInfoList
      for (std::size_t i = 0; i < RawDMList[q1].size1(); ++i)
      {
 	 EigenInfoList.push_back(EigenInfo(Eigenvalues[i], CurrentDegree, q1, i));
	 ESum += Eigenvalues[i] * CurrentDegree;
      }
   }

   // sort the eigenvalues
   if (Sort)
      std::sort(EigenInfoList.begin(), EigenInfoList.end(), EigenCompare);
}


//
// DensityMatrix<MatrixOperator>
//

DensityMatrix<MatrixOperator>::DensityMatrix(MatrixOperator const& Op)
   : B(Op.Basis1())
{
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), Op.Basis2());
   // initialize the RawDMList.  At the same time, get the maximum linear size of the subspaces
   RawDMList.resize(B.size());
   MaxLinearDimension = 0;
   for (std::size_t q = 0; q < B.size(); ++q)
   {
      int Dim = B.dim(q);
      //      std::cout << "dimension of " << q1 << " is " << Dim << std::endl;
      RawDMList[q] = RawDMType(Dim, Dim);
      LinearAlgebra::fill(RawDMList[q], std::complex<double>(0));
      MaxLinearDimension = std::max<int>(MaxLinearDimension, Dim);
   }

   // Fill the raw density matrices
   for (LinearAlgebra::const_iterator<MatrixOperator>::type I = iterate(Op); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<MatrixOperator>::type J = iterate(I); J; ++J)
      {
         int tp;
         LinearAlgebra::Range rtp;
         boost::tie(tp, rtp) = B.Lookup(J.index1());
	 int t;
	 LinearAlgebra::Range rt;
	 boost::tie(t, rt) = B.Lookup(J.index2());
	 CHECK_EQUAL(tp,t)("The density matrix must be block-diagonal")(B[tp])(B[t])(Op)(B);
	 RawDMList[tp](rtp, rt) = *J;
      }
   }
   // diagonalize them
   this->DiagonalizeDMHelper();
}

DensityMatrix<MatrixOperator>::DensityMatrix(MatrixOperator const& Op, MatrixOperator const& WavefunctionDM)
   : B(Op.Basis1())
{
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), Op.Basis2());
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), WavefunctionDM.Basis1());
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), WavefunctionDM.Basis2());

   // We have a second RawDMType, for the WavefunctionDM
   std::vector<RawDMType> PsiDMList(B.size());

   // initialize the RawDMList.  At the same time, get the maximum linear size of the subspaces
   RawDMList.resize(B.size());
   MaxLinearDimension = 0;
   for (std::size_t q = 0; q < B.size(); ++q)
   {
      int Dim = B.dim(q);
      //      std::cout << "dimension of " << q1 << " is " << Dim << std::endl;
      RawDMList[q] = RawDMType(Dim, Dim);
      LinearAlgebra::fill(RawDMList[q], std::complex<double>(0));
      PsiDMList[q] = RawDMType(Dim, Dim);
      LinearAlgebra::fill(PsiDMList[q], std::complex<double>(0));
      MaxLinearDimension = std::max<int>(MaxLinearDimension, Dim);
   }

   // Fill the raw density matrices
   for (LinearAlgebra::const_iterator<MatrixOperator>::type I = iterate(Op); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<MatrixOperator>::type J = iterate(I); J; ++J)
      {
         int tp;
         LinearAlgebra::Range rtp;
         boost::tie(tp, rtp) = B.Lookup(J.index1());
	 int t;
	 LinearAlgebra::Range rt;
	 boost::tie(t, rt) = B.Lookup(J.index2());
	 CHECK_EQUAL(tp,t)("The density matrix must be block-diagonal")(B[tp])(B[t])(Op)(B);
	 RawDMList[tp](rtp, rt) = *J;
      }
   }

   // Fill the wavefunction DM
   for (LinearAlgebra::const_iterator<MatrixOperator>::type I = iterate(WavefunctionDM); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<MatrixOperator>::type J = iterate(I); J; ++J)
      {
         int tp;
         LinearAlgebra::Range rtp;
         boost::tie(tp, rtp) = B.Lookup(J.index1());
	 int t;
	 LinearAlgebra::Range rt;
	 boost::tie(t, rt) = B.Lookup(J.index2());
	 CHECK_EQUAL(tp,t)("The density matrix must be block-diagonal")(B[tp])(B[t])(WavefunctionDM)(B);
	 PsiDMList[tp](rtp, rt) = *J;
      }
   }

   // diagonalize the density matrix
   this->DiagonalizeDMHelper(false);

   // Get the eigenvalues with respect to WavefunctionDM
   std::size_t i = 0;
   for (std::size_t j = 0; j < PsiDMList.size(); ++j)
   {
      RawDMType M = RawDMList[j] * PsiDMList[j] * herm(RawDMList[j]);
      for (std::size_t k = 0; k < M.size1(); ++k)
      {
         EigenInfoList[i].Eigenvalue = LinearAlgebra::real(M(k,k));
         //         TRACE(EigenInfoList[i].Eigenvalue);
         ++i;
      }
   }
   CHECK_EQUAL(i, EigenInfoList.size());
   std::sort(EigenInfoList.begin(), EigenInfoList.end(), EigenCompare);
}

//
// DensityMatrix<SimpleOperator>
//

DensityMatrix<SimpleOperator>::DensityMatrix(SimpleOperator const& Op)
   : B(Op.Basis1())
{
   DEBUG_PRECONDITION_EQUAL(Op.Basis1(), Op.Basis2());
   // initialize the RawDMList.  At the same time, get the maximum linear size of the subspaces
   RawDMList.resize(B.size());
   MaxLinearDimension = 0;
   for (std::size_t q = 0; q < B.size(); ++q)
   {
      int Dim = B.dim(q);
      //      std::cout << "dimension of " << q << " is " << Dim << std::endl;
      RawDMList[q] = LinearAlgebra::Matrix<double>(Dim, Dim, 0);
      MaxLinearDimension = std::max<int>(MaxLinearDimension, Dim);
   }
   // Fill the raw density matrices
   for (LinearAlgebra::const_iterator<SimpleOperator>::type I = iterate(Op); I; ++I)
   {
      for (LinearAlgebra::const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         int tp, rtp;
         boost::tie(tp, rtp) = B.Lookup(J.index1());
	 int t, rt;
	 boost::tie(t, rt) = B.Lookup(J.index2());
	 CHECK_EQUAL(tp,t)("The density matrix must be block-diagonal");
 	 RawDMList[tp](rtp, rt) = *J;
      }
   }
   // diagonalize them
   this->DiagonalizeDMHelper();
}

void UpdateKeepList(KeepListType& KeepList,
		    std::set<QuantumNumbers::QuantumNumber> const& SiteQN,
		    VectorBasis const& FullBasis,
		    std::list<EigenInfo>& KeepStates, 
		    std::list<EigenInfo>& DiscardStates,
		    TruncationInfo& Info)
{
   typedef std::set<QuantumNumbers::QuantumNumber> SiteQNType;
   typedef std::list<EigenInfo> StatesListType;

   // Make a set of all quantum numbers that are already kept
   typedef std::set<QuantumNumbers::QuantumNumber> qnType;
   qnType KeptQN;
   for (StatesListType::const_iterator I = KeepStates.begin(); I != KeepStates.end(); ++I)
      KeptQN.insert(FullBasis[I->Subspace]);

   // Iterate through the KeepList and get the transformed states, make sure they
   // exist in KeepStates
   for (KeepListType::const_iterator I = KeepList.begin(); I != KeepList.end(); ++I)
   {
      qnType qn;  // the set of updated quantum numbers corresponding to keep state I
      for (SiteQNType::const_iterator Q = SiteQN.begin(); Q != SiteQN.end(); ++Q)
	 transform_targets(*I, *Q, std::inserter(qn, qn.end()));
      // make sure that at least one of qn exists in KeptQN
      SiteQNType Intersect;
      std::set_intersection(KeptQN.begin(), KeptQN.end(), qn.begin(), qn.end(), 
			    std::inserter(Intersect, Intersect.end()));
      if (Intersect.empty())
      {
	 // No kept state, go through the discard list and resurrect a state
	 StatesListType::iterator Piv = DiscardStates.begin();
	 while (Piv != DiscardStates.end() && qn.count(FullBasis[Piv->Subspace]) == 0)
	    ++Piv;
	 if (Piv == DiscardStates.end())
	 {
	    WARNING("Unexpected: cannot force quantum number into basis")(*I);
	 }
	 else
	 {
	    DEBUG_WARNING("Forcing keep state")(*I)(FullBasis[Piv->Subspace])(Piv->Weight());
	    ++Info.KeptStates_;
	    Info.KeptWeight_ += Piv->Weight();
	    Info.KeptEntropy_ += Piv->Entropy(Info.TotalWeight_);
	    KeptQN.insert(FullBasis[Piv->Subspace]);
	    KeepStates.push_back(*Piv);
	    DiscardStates.erase(Piv);
	 }
      }
   }

   // Now construct the new KeepList.
   KeepList.clear();
   for (StatesListType::const_iterator I = KeepStates.begin(); I != KeepStates.end(); ++I)
   {
      KeepList.insert(FullBasis[I->Subspace]);
   }
   //std::copy(KeepList.begin(),KeepList.end(),std::ostream_iterator<QuantumNumber>(std::cout, " "));
   //std::cout << '\n';
}

//
// Singular value decomposition
//

SingularDecompositionBase::SingularDecompositionBase()
{
}

SingularDecompositionBase::~SingularDecompositionBase()
{
}

void SingularDecompositionBase:: Diagonalize(std::vector<RawDMType> const& M)
{
   ESum = 0;  // Running sum of squares of the singular values
   for (std::size_t i = 0; i < M.size(); ++i)
   {
      LinearAlgebra::Matrix<std::complex<double> > U, Vh;
      LinearAlgebra::Vector<double> D;
      SingularValueDecomposition(M[i], U, D, Vh);

      LeftVectors.push_back(U);
      RightVectors.push_back(Vh);
      SingularValues.push_back(D);
      
      int CurrentDegree = degree(this->Lookup(i));
      for (unsigned j = 0; j < size(D); ++j)
      {
         double Weight = D[j]*D[j]; 
         EigenInfoList.push_back(EigenInfo(Weight, CurrentDegree, i, j));
         ESum += Weight * CurrentDegree;
      }
   }

   std::sort(EigenInfoList.begin(), EigenInfoList.end(), EigenCompare);
}


SingularDecomposition<MatrixOperator, MatrixOperator>::SingularDecomposition(MatrixOperator const& M)
   : B1(M.Basis1()), B2(M.Basis2()), IndexOfi(B1.size(), -1)
{
   // Assemble the raw matrices

   // Iterate through Basis1 and Basis2 and assemble the list of used quantum numbers, in order,
   // and how they map onto the components in Basis1 and Basis2.
   // We can make use of the fact that each quantum number only appears at most once in 
   // the LinearBasis.
   for (unsigned i = 0; i < B1.size(); ++i)
   {
      QuantumNumber q = B1[i];
      int j = B2.find_first(q);
      if (j == -1)
	 continue;

      IndexOfi[i] = UsedQuantumNumbers.size();
      q_iLinear.push_back(i);
      q_jLinear.push_back(j);
      UsedQuantumNumbers.push_back(q);
   }

   // Now assemble the list of raw matrices, initialized to zero
   std::vector<RawDMType> Matrices(q_iLinear.size());
   for (unsigned n = 0; n < q_iLinear.size(); ++n)
   {
      Matrices[n] = RawDMType(B1.dim(q_iLinear[n]), B2.dim(q_jLinear[n]), 0.0);
   }

   // fill the raw matrices with the components in M
   for (MatrixOperator::const_iterator I = iterate(M); I; ++I)
   {
      for (MatrixOperator::const_inner_iterator J = iterate(I); J; ++J)
      {
	 // determine where this (i,j) component fits within the matrices
	 std::pair<int, LinearAlgebra::Range> iIndex = B1.Lookup(J.index1());
	 std::pair<int, LinearAlgebra::Range> jIndex = B2.Lookup(J.index2());
	 // map the index into the used subspaces
	 int Subspace = IndexOfi[iIndex.first];
	 CHECK(Subspace != -1);
	 CHECK_EQUAL(q_jLinear[Subspace], jIndex.first);
	 Matrices[Subspace](iIndex.second, jIndex.second) = *J;
      }
   }

   // do the SVD
   this->Diagonalize(Matrices);
}

QuantumNumber 
SingularDecomposition<MatrixOperator, MatrixOperator>::Lookup(int Subspace) const
{
   return UsedQuantumNumbers[Subspace];
}

void
SingularDecomposition<MatrixOperator, MatrixOperator>::
ConstructOrthoMatrices(std::vector<std::set<int> > const& LinearMapping,
		  MatrixOperator& A, RealDiagonalOperator& C, MatrixOperator& B)
{
   // Construct the truncated basis
   int NumQ = UsedQuantumNumbers.size();
   VectorBasis NewBasis(B1.GetSymmetryList());
   std::vector<int> NewSubspace(NumQ, -1); // maps the q to the new label
   for (int q = 0; q < NumQ; ++q)
   {
      if (LinearMapping[q].size() > 0)
      {
	 NewSubspace[q] = NewBasis.size();
         NewBasis.push_back(UsedQuantumNumbers[q], LinearMapping[q].size());
      }
   }

   // Now we construct the actual matrices
   A = MatrixOperator(B1.MappedBasis(), NewBasis);
   B = MatrixOperator(NewBasis, B2.MappedBasis());
   for (int q = 0; q < NumQ; ++q)
   {
      int ss = NewSubspace[q];
      if (ss == -1) // is this subspace used?
	 continue;

      QuantumNumber Q = UsedQuantumNumbers[q];

      // The set of indices of states we need to keep
      std::vector<int> lm(LinearMapping[q].begin(), LinearMapping[q].end());

      // Now assemble this quantum number component for A
      int i = B1.MappedBasis().find_first(Q);
      while (i != -1)
      {
	 A(i,ss) = LeftVectors[q](B1.Lookup(i).second, lm);
	 i = B1.MappedBasis().find_next(Q, i);
      }

      // and B
      int j = B2.MappedBasis().find_first(Q);
      while (j != -1)
      {
	 B(ss,j) = RightVectors[q](lm, B2.Lookup(j).second);
	 j = B2.MappedBasis().find_next(Q, j);
      }
   }

   // Finally the center matrix
   C = RealDiagonalOperator(NewBasis, NewBasis);
   for (int i = 0; i < NumQ; ++i)
   {
      if (!LinearMapping[i].empty())
      {
	 int b = NewSubspace[i];
	 C(b,b) = LinearAlgebra::DiagonalMatrix<double>(SingularValues[i][std::vector<int>(LinearMapping[i].begin(), LinearMapping[i].end())]);
      }
   }

   A.debug_check_structure();
   B.debug_check_structure();
}

QuantumNumber 
SingularDecomposition<StateComponent, StateComponent>::Lookup(int Subspace) const
{
   return UsedQuantumNumbers[Subspace];
}

SingularDecomposition<StateComponent, StateComponent>::
SingularDecomposition(StateComponent const& A, ProductBasis<BasisList, BasisList> const& Factors_)
   : SingularDecompositionBase(), B1(A.Basis1()), B2(A.Basis2()), Factors(Factors_)
{
   using QuantumNumbers::QuantumNumberList;

   CHECK_EQUAL(A.LocalBasis(), Factors.Basis());

   // denote the operator by <j' || A[k] || j>
   // where A[k] decomposes as S[k_1] T[k_2]
   // The matrix is block diagonal with respect to (j' - k_1) and (j + k_2).
   // So the quantum number for the singular values runs over these numbers.

   // map each distinct quantum number onto successive integers
   int NumQuantum = 0;
   std::map<QuantumNumber, int> QuantumNumberLinearMapping;// maps to the index in UsedQuantumNumbers

   // linear dimensions of the subspaces for quantum number pair (j' - k_1) = (j + k_2) = r
   std::vector<std::pair<int, int> > LinearDimensions;

   // for Basis1, construct the combined basis
   for (unsigned k1 = 0; k1 < Factors.Left().size(); ++k1)
   {
      QuantumNumber k1q_bar = adjoint(Factors.Left()[k1]);

      for (unsigned i = 0; i < A.Basis1().size(); ++i)
      {
	 QuantumNumberList ql = transform_targets(A.Basis1()[i], k1q_bar);
	 for (QuantumNumberList::const_iterator qi = ql.begin(); qi != ql.end(); ++qi)
	 {
	    int LinearIndex;
	    std::map<QuantumNumber, int>::iterator I = QuantumNumberLinearMapping.find(*qi);
	    if (I == QuantumNumberLinearMapping.end())
	    {
	       // not found, add it
	       LinearIndex = UsedQuantumNumbers.size();
	       QuantumNumberLinearMapping[*qi] = LinearIndex;
	       UsedQuantumNumbers.push_back(*qi);
	       LinearDimensions.push_back(std::pair<int, int>(0,0));
	    }
	    else
	       LinearIndex = I->second;

	    int dim = A.Basis1().dim(i);
	    int current = LinearDimensions[LinearIndex].first;
	    LeftSubspaceInfo.push_back(ProductSubspaceInfo(k1, i, LinearIndex,
							   LinearAlgebra::Range(current, current+dim)));
	    LinearDimensions[LinearIndex].first += dim;
	 }
      }
   }

   // same for Basis2

   for (unsigned k2 = 0; k2 < Factors.Right().size(); ++k2)
   {
      QuantumNumber k2q = Factors.Right()[k2];

      for (unsigned i = 0; i < A.Basis2().size(); ++i)
      {
	 QuantumNumberList ql = transform_targets(A.Basis2()[i], k2q);

	 for (QuantumNumberList::const_iterator qi = ql.begin(); qi != ql.end(); ++qi)
	 {
	    int LinearIndex;
	    std::map<QuantumNumber, int>::iterator I = QuantumNumberLinearMapping.find(*qi);
	    if (I == QuantumNumberLinearMapping.end())
	    {
	       // not found, add it.  This is a degenerate case in this context, 
	       // since the quantum number doesn't exist in the left basis so
	       // this component is a zero-dimensional matrix
	       LinearIndex = UsedQuantumNumbers.size();
	       QuantumNumberLinearMapping[*qi] = LinearIndex;
	       UsedQuantumNumbers.push_back(*qi);
	       LinearDimensions.push_back(std::pair<int, int>(0,0));
	    }
	    else
	       LinearIndex = I->second;

	    int dim = A.Basis2().dim(i);
	    int current = LinearDimensions[LinearIndex].second;
	    RightSubspaceInfo.push_back(ProductSubspaceInfo(k2, i, LinearIndex,
							   LinearAlgebra::Range(current, current+dim)));
	    LinearDimensions[LinearIndex].second += dim;
	 }
      }
   }

   NumQuantum = UsedQuantumNumbers.size();

   // Now we can construct the actual matrices, initialized to zero
   std::vector<RawDMType> Matrices(NumQuantum);
   for (int i = 0; i < NumQuantum; ++i)
      Matrices[i] = RawDMType(LinearDimensions[i].first, LinearDimensions[i].second, 0.0);

   // and fill them
   for (unsigned lin_1 = 0; lin_1 < LeftSubspaceInfo.size(); ++lin_1)
   {
      for (unsigned lin_2 = 0; lin_2 < RightSubspaceInfo.size(); ++lin_2)
      {
	 // we only care about blocks that are on the diagonal
	 if (LeftSubspaceInfo[lin_1].LinearIndex != RightSubspaceInfo[lin_2].LinearIndex)
	    continue;

	 // and make sure that block also has some non-zero element
	 if (LeftSubspaceInfo[lin_1].LinearRange.size() == 0 || RightSubspaceInfo[lin_2].LinearRange.size() == 0)
	    continue;

	 // sum over k
	 ProductBasis<BasisList, BasisList>::const_iterator klEnd 
	    = Factors.end(LeftSubspaceInfo[lin_1].k,
			  RightSubspaceInfo[lin_2].k);
	 ProductBasis<BasisList, BasisList>::const_iterator klIter
	    = Factors.begin(LeftSubspaceInfo[lin_1].k,
			    RightSubspaceInfo[lin_2].k);
	 for ( ; klIter != klEnd; ++klIter)
	 {
	    const_inner_iterator<MatrixOperator>::type J = iterate_at(A[*klIter].data(), 
								      LeftSubspaceInfo[lin_1].s, 
								      RightSubspaceInfo[lin_2].s);
	    if (J)
	    {
	       double Coeff 
		  = inverse_product_coefficient(Factors.Left()[LeftSubspaceInfo[lin_1].k],
						Factors.Right()[RightSubspaceInfo[lin_2].k],
						Factors.Basis()[*klIter],
						A.Basis1()[LeftSubspaceInfo[lin_1].s],
						A.Basis2()[RightSubspaceInfo[lin_2].s],
						UsedQuantumNumbers[LeftSubspaceInfo[lin_1].LinearIndex]);

	       // The normalization factor is to satisfy the left-normalization constraint
	       // with the quantum-number dependent prefactor.  See section 7.1 of the nonabelianmp,
	       // equation eq:LeftAMatrixNorm.  Given an ordinary unitary matrix, if we insert it
	       // into an A-matrix, the normalization changes due to the conventions of the reduced
	       // matrix elements.  To counteract this, we need to scale the unitary matrix by
	       // a factor as it is inserted into the A-matrix, and apply the inverse of the
	       // scale factor at this point here.
	       double NormFactor = std::sqrt(double(degree(A.Basis1()[LeftSubspaceInfo[lin_1].s]))
				       / degree(UsedQuantumNumbers[LeftSubspaceInfo[lin_1].LinearIndex]));

	       Matrices[LeftSubspaceInfo[lin_1].LinearIndex](LeftSubspaceInfo[lin_1].LinearRange,
							     RightSubspaceInfo[lin_2].LinearRange)
		  += NormFactor * Coeff * (*J);
	    }
	 }
      }
   }

   // do the SVD
   this->Diagonalize(Matrices);
}

void
SingularDecomposition<StateComponent, StateComponent>::
ConstructOrthoMatrices(std::vector<std::set<int> > const& LinearMapping,
		  StateComponent& A, RealDiagonalOperator& C, StateComponent& B)
{
   int NumQ = UsedQuantumNumbers.size();
   VectorBasis NewBasis(B1.GetSymmetryList());

   // Now we construct the truncated basis
   std::vector<int> NewSubspace(NumQ, -1); // maps the q to the new label
   for (int q = 0; q < NumQ; ++q)
   {
      if (LinearMapping[q].size() > 0)
      {
	 NewSubspace[q] = NewBasis.size();
         NewBasis.push_back(this->Lookup(q), LinearMapping[q].size());
      }
   }
   // Now we construct the actual matrices, firstly for the left
   A = StateComponent(Factors.Left(), B1, NewBasis);
   for (unsigned i = 0; i < LeftSubspaceInfo.size(); ++i)
   {
      int ss = LeftSubspaceInfo[i].LinearIndex;
      if (!LinearMapping[ss].empty())
      {
	 std::vector<int> lm(LinearMapping[ss].begin(), LinearMapping[ss].end());

	 // Apply the normalization factor to ensure A satisfies the left-orthonormalization constraint.
	 // This is the reciprocal of the factor we applied previously.
	 double NormFactor = std::sqrt(double(degree(A.Basis2()[NewSubspace[ss]]))
				       / degree(A.Basis1()[LeftSubspaceInfo[i].s]));

	 A[LeftSubspaceInfo[i].k](LeftSubspaceInfo[i].s, NewSubspace[ss])
	    = NormFactor * LeftVectors[ss](LeftSubspaceInfo[i].LinearRange, lm);
      }
   }

   // now the right
   B = StateComponent(Factors.Right(), NewBasis, B2);
   for (unsigned i = 0; i < RightSubspaceInfo.size(); ++i)
   {
      int ss = RightSubspaceInfo[i].LinearIndex;
      std::vector<int> lm(LinearMapping[ss].begin(), LinearMapping[ss].end());
      if (!lm.empty())
      {
	 B[RightSubspaceInfo[i].k](NewSubspace[ss], RightSubspaceInfo[i].s)
	    = RightVectors[ss](lm, RightSubspaceInfo[i].LinearRange);
      }
   }

   // Finally the center matrix
   C = RealDiagonalOperator(NewBasis, NewBasis);
   for (int i = 0; i < NumQ; ++i)
   {
      if (!LinearMapping[i].empty())
      {
	 int b = NewSubspace[i];
	 C(b,b) = LinearAlgebra::DiagonalMatrix<double>(SingularValues[i][std::vector<int>(LinearMapping[i].begin(), LinearMapping[i].end())]);
      }
   }
}
