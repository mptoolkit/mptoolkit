// -*- C++ -*- $Id$

#include "mpexponential.h"
#include "linearalgebra/eigen.h"
#include "linearalgebra/matrix_utility.h"
#include "mpsvd.h"
#include "tensor/tensor_eigen.h"
#include "tensor/tensor_exponential.h"

struct ExpProd
{
   typedef std::complex<double> result_type;
   typedef SimpleOperator const& first_argument_type;
   typedef SimpleOperator const& second_argument_type;

   result_type operator()(first_argument_type x, second_argument_type y) const
   {
      DEBUG_PRECONDITION_EQUAL(x.Basis1().size(), 1);      
      DEBUG_PRECONDITION_EQUAL(y.Basis2().size(), 1);
      return trace(prod(x, y, QuantumNumbers::QuantumNumber(x.GetSymmetryList())));
   }
};

std::pair<MPOpComponent, MPOpComponent>
TwoSiteExponential(SimpleOperator const& A, SimpleOperator const& B, std::complex<double> x)
{
   QuantumNumber Ident(A.GetSymmetryList());
   ProductBasis<BasisList, BasisList> PB(A.Basis1(), B.Basis1());
   SimpleOperator TP = tensor_prod(A, B, PB, PB, Ident);

   //   TRACE(A)(B)(TP);

   // now we apply the exponential to TP
   if (LinearAlgebra::norm_frob(x) == 0.0)
   {
      BasisList BL = make_vacuum_basis(A.GetSymmetryList());
      MPOpComponent MB = MPOpComponent(B.Basis1(), BL, BL);
      MPOpComponent MA = MPOpComponent(A.Basis1(), BL, BL);
      MB.set_operator(0,0, SimpleOperator::make_identity(A.Basis1()));
      MA.set_operator(0,0, SimpleOperator::make_identity(A.Basis1()));
      return std::make_pair(MA, MB);
   }
   // else

   TP = Tensor::Exponentiate(x*TP);

   //TRACE(TP);

   typedef std::map<PartialProdIndex, std::complex<double> > PartialProdType;
   PartialProdType PartialProd = decompose_tensor_prod(TP, PB, PB);

   // rearrange the matrix elements to do a singular value decomposition
   
   BasisList BondBasis(A.GetSymmetryList()); // bond of the MPO
   std::vector<SimpleOperator> MatLeft, MatRight;

   // firstly, arrange by quantum number (which is easy because we have a scalar operator,
   // so qLeft = adjoint(qRight).  So label the quantum numbers by qRight.

   std::set<QuantumNumber> qUsed;
   for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
   {
      qUsed.insert(I->first.qRight);
   }

   // secondly, for each distinct quantum number, put the matrix elements into a matrix
   // indexed by (Left1,Left2), (Right1,Right2)
   for (std::set<QuantumNumber>::const_iterator q = qUsed.begin(); q != qUsed.end(); ++q)
   {

      // First step is to find out which pairs (Left1,Left2), (Right1,Right2) are
      // used and enumerate them
      typedef std::map<std::pair<int, int>, int> BasisEnumerationType;
      std::map<std::pair<int, int>, int> LeftBasis, RightBasis;
      int nLeft = 0, nRight = 0;
      for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
      {
	 if (*q == I->first.qRight) // filter by quantum number
	 {
	    if (LeftBasis.find(std::make_pair(I->first.Left1, I->first.Left2)) == LeftBasis.end())
	       LeftBasis[std::make_pair(I->first.Left1, I->first.Left2)] = nLeft++;

	    if (RightBasis.find(std::make_pair(I->first.Right1, I->first.Right2)) == RightBasis.end())
	       RightBasis[std::make_pair(I->first.Right1, I->first.Right2)] = nRight++;
	 }
      }

      // now make a matrix
      LinearAlgebra::Matrix<std::complex<double> > Mat(nLeft, nRight, 0.0);
      // fill it
      for (PartialProdType::const_iterator I = PartialProd.begin(); I != PartialProd.end(); ++I)
      {
	 if (*q == I->first.qRight) // filter by quantum number
	 {
	    Mat(LeftBasis[std::make_pair(I->first.Left1, I->first.Left2)], 
		RightBasis[std::make_pair(I->first.Right1, I->first.Right2)]) = I->second;
	 }
      }
      // SVD
      LinearAlgebra::Matrix<std::complex<double> > U, Vt;
      LinearAlgebra::Vector<double> D;
      SingularValueDecomposition(Mat, U, D, Vt);

      // the number of non-zero singular values is equal to the number of degenerate
      // states of quantum number *q we need in the MPO bond local basis

      // get the number of non-zero singular values
      int nnz = 0;
      while (nnz < int(size(D)) && fabs(D[nnz]) > 1e-12) ++nnz;

      //      TRACE(size(D))(nnz);
      
      // assemble the components
      for (int i = 0; i < nnz; ++i)
      {
	 BondBasis.push_back(*q);
	 SimpleOperator MR(B.Basis1(), B.Basis2(), *q);
	 SimpleOperator ML(A.Basis1(), A.Basis2(), adjoint(*q));
	 for (BasisEnumerationType::const_iterator I = RightBasis.begin(); I != RightBasis.end(); ++I)
	 {
	    MR(I->first.first, I->first.second) = Vt(i, I->second);
	 }
	 for (BasisEnumerationType::const_iterator I = LeftBasis.begin(); I != LeftBasis.end(); ++I)
	 {
	    ML(I->first.first, I->first.second) = U(I->second, i) * D[i];
	 }
	 //TRACE(ML)(MR)(tensor_prod(ML, MR, Ident));
	 MatRight.push_back(MR);
	 MatLeft.push_back(ML);
      }
   }

   // assemble the MPO.  Is there an additional recoupling coefficient here in going from
   // the sum of operators to an MPO?  Apparantly not...

   MPOpComponent MB = MPOpComponent(B.Basis1(), BondBasis, make_vacuum_basis(A.GetSymmetryList()));
   MPOpComponent MA = MPOpComponent(A.Basis1(), make_vacuum_basis(A.GetSymmetryList()), BondBasis);

   for (unsigned i = 0; i < BondBasis.size(); ++i)
   {
      MA.set_operator(0, i, MatLeft[i]);
      MB.set_operator(i, 0, MatRight[i]);
   }

#if 0 // test
   SimpleOperator TPtest;
   for (MPOpComponent::const_iterator I = MA.begin(); I != MA.end(); ++I)
   {
      for (MPOpComponent::const_iterator J = MB.begin(); J != MB.end(); ++J)
      {
         if (is_transform_target(I->first, J->first, Ident))
	 {
            TPtest += tensor_prod(I->second, J->second, PB, PB, Ident, ExpProd());
	    //TRACE(tensor_prod(I->second, J->second, PB, PB, Ident, ExpProd()));
	 }
      }
   }
   //TRACE(TPtest);
#endif

   return std::make_pair(MA, MB);
}

void TwoSiteExponential(MPOpComponent& A, MPOpComponent& B, std::complex<double> x)
{
   PRECONDITION_EQUAL(A.Basis1().size(), 1);
   PRECONDITION_EQUAL(B.Basis2().size(), 1);

   QuantumNumber Ident(A.GetSymmetryList());
   PRECONDITION_EQUAL(A.Basis1()[0], Ident);
   PRECONDITION_EQUAL(B.Basis2()[0], Ident);

   //TRACE(A)(B);

   SimpleOperator C = ExpandBasis2(A);
   A = prod(A, C);

   C = ExpandBasis1(B);
   B = prod(C, B);

   //TRACE(A)(B);

   ProductBasis<BasisList, BasisList> PB(A.SiteBasis(), B.SiteBasis());

   SimpleOperator TP;
   for (MPOpComponent::const_iterator I = A.begin(); I != A.end(); ++I)
   {
      for (MPOpComponent::const_iterator J = B.begin(); J != B.end(); ++J)
      {
         if (is_transform_target(I->first, J->first, Ident))
            TP += tensor_prod(I->second, J->second, PB, PB, Ident, ExpProd());
      }
   }

   //TRACE(TP);

   // now we apply the exponential to TP
   LinearAlgebra::Matrix<std::complex<double> > M(TP.data());

   LinearAlgebra::Vector<std::complex<double> > EVec = DiagonalizeHermitian(M);
   for (std::size_t i = 0; i < EVec.size(); ++i)
   {
      EVec[i] = exp(EVec[i] * x);
   }

   zero_all(TP.data());
   M = transpose(M) * diagonal_matrix(EVec) * M;
   for (std::size_t i = 0; i < size1(M); ++i)
   {
      for (std::size_t j = 0; j < size2(M); ++j)
      {
         if (LinearAlgebra::norm_2(M(i,j)) > 1E-10) set_element(TP.data(), i, j, M(i,j));
      }
   }
   //TRACE(TP);

   ProductBasis<BasisList, BasisList> alpha(A.SiteBasis(), adjoint(A.SiteBasis()));
   ProductBasis<BasisList, BasisList> beta(B.SiteBasis(), adjoint(B.SiteBasis()));

   SimpleOperator Partial;// = partial_transpose(TP, PB, PB, alpha, beta);

   SimpleOperator U, Vt;
   LinearAlgebra::Vector<double> D;

   SingularValueDecomposition(Partial, U, D, Vt);

   A = MPOpComponent(A.SiteBasis(), A.Basis1(), U.Basis2());
   B = MPOpComponent(B.SiteBasis(), Vt.Basis1(), B.Basis2());

   for (const_iterator<SimpleOperator>::type I = iterate(U); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         int l,m;
         std::tie(l,m) = alpha.rmap(J.index1());
         QuantumNumber q = alpha[J.index1()];
         SimpleOperator M(A.Basis1(), A.Basis2(), q);
         M(0,J.index2()) = *J * D[J.index2()];
         A[q](l,m) += M;
      }
   }

   for (const_iterator<SimpleOperator>::type I = iterate(Vt); I; ++I)
   {
      for (const_inner_iterator<SimpleOperator>::type J = iterate(I); J; ++J)
      {
         int l,m;
         std::tie(l,m) = beta.rmap(J.index2());
         QuantumNumber q = beta[J.index2()];
         SimpleOperator M(B.Basis1(), B.Basis2(), q);
         M(J.index1(),0) = *J / std::sqrt(double(degree(B.SiteBasis()[l])));
         //TRACE(q)(B.SiteBasis()[l])(B.SiteBasis()[m])(B.Basis1()[J.index1()])(l)(m)(J.index1());
         B[q](l,m) += M;
      }
   }
}

#if 0
SplitOperator BondExponential(std::complex<double> x, SplitOperator const& Op)
{
   SplitOperator Result(Op);
   
   std::complex<double> Scale = 1.0;

   QuantumNumber Ident = Op.TransformsAs();
   CHECK(is_scalar(Ident));

   while (IsProportionalIdentity(Result.Left()) && Result.RightSize() > 1)
   {
      Scale *= IdentityScale(Result.Left());
      Result.Left() = ConstructIdentity(Result.Left().SiteBasis());
      Result.RotateRight();
   }

   BasisList VacuumBasis = make_vacuum_basis(Op.GetSymmetryList());

   CHECK_EQUAL(Result.Left().Basis1(), VacuumBasis);
   CHECK_EQUAL(Result.Right().Basis2(), VacuumBasis);

   //TRACE(Result.RightSize());

   for (int n = 0; n < Result.RightSize()-1; ++n)
   {
      CHECK(IsProportionalIdentity(Result.LookupRight(n)))(n)(Result.LookupRight(n));
      Scale *= IdentityScale(Result.LookupRight(n));
      Result.LookupRight(n) = ConstructIdentity(Result.LookupRight(n).SiteBasis());
   }

   Result.Left() = prod(Result.Left(), Result.Center());
   //TRACE(Scale);
   TwoSiteExponential(Result.Left(), Result.Right(), x * Scale);
   Result.Center() = TruncateBasis1(Result.Right());
   //Result.Center() = SimpleOperator::make_identity(Result.Left().Basis2());

   while (Result.LeftSize() > 1) Result.RotateLeft();

   return Result;
}
#endif
