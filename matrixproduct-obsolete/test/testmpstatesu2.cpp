// -*- C++ -*- $Id$

#include "matrixproduct/mpstate.h"
#include "quantumnumbers/su2.h"
#include "quantumnumbers/u1.h"

using namespace QuantumNumbers;
using namespace Tensor;

int main()
{
   // local site basis
   SymmetryList Symmetry("S:SU(2)");
   QuantumNumbers::QNConstructor<QuantumNumbers::SU2> QN(Symmetry);
   BasisList B(Symmetry);
   B.push_back(QN(1));

   // auxiliary basis
   VectorBasis VB(Symmetry);
   VB.push_back(QN(0.5), 1);

   MatrixOperator AuxIdent(VB, VB, QN(0));
   set_element(AuxIdent.data(), 0, 0, LinearAlgebra::identity_matrix<double>(1));

   MPStateComponent x(B, VB, VB);

   MPStateComponent FB1 = MPStateComponent::ConstructFullBasis1(B, VB);

   MatrixOperator FullIdent1(FB1.Basis1(), FB1.Basis1(), QN(0));
   for (std::size_t i = 0; i < FB1.Basis1().size(); ++i)
   {
      set_element(FullIdent1.data(), i, i, 
		  LinearAlgebra::identity_matrix<double>(1));
   }

   CHECK_CLOSE(scalar_prod(FB1, herm(FB1)), FullIdent1);

   MPStateComponent FB2 =  MPStateComponent::ConstructFullBasis2(VB, B);

   MatrixOperator FullIdent2(FB2.Basis2(), FB2.Basis2(), QN(0));
   for (std::size_t i = 0; i < FB2.Basis2().size(); ++i)
   {
      set_element(FullIdent2.data(), i, i, 
		  LinearAlgebra::identity_matrix<double>(1));
   }

   CHECK_CLOSE(scalar_prod(herm(FB2), FB2), FullIdent2);


   SimpleOperator SiteIdent(B, B, QN(0));
   SiteIdent(0,0) = 1;

   CHECK_CLOSE(operator_prod(SiteIdent, FB1, herm(FB1)), FullIdent1);
   CHECK_CLOSE(operator_prod(SiteIdent, herm(FB2), FB2), FullIdent2);

   CHECK_CLOSE(operator_prod(SiteIdent, FB1, AuxIdent, herm(FB1)), FullIdent1);
   CHECK_CLOSE(operator_prod(SiteIdent, herm(FB2), AuxIdent, FB2), FullIdent2);

   MatrixOperator M = 0.5 * operator_prod(2.0*SiteIdent, herm(FB2), FB2);
   CHECK_CLOSE(M, scalar_prod(herm(FB2), FB2));

   M = 0.5 * operator_prod(2.0*SiteIdent, FB1, herm(FB1));
   CHECK_CLOSE(M, scalar_prod(FB1, herm(FB1)));


   VectorBasis i(Symmetry), j(Symmetry), ip(Symmetry), jp(Symmetry);
   BasisList s(Symmetry), sp(Symmetry);
   i.push_back(QN(3), 1);
   ip.push_back(QN(9), 1);
   j.push_back(QN(1), 1);
   jp.push_back(QN(3), 1);
   s.push_back(QN(2));
   sp.push_back(QN(7));

   MPStateComponent TestA(sp, ip, i);
   set_element(TestA[0].data(), 0, 0,  0.5 * LinearAlgebra::identity_matrix<double>(1));

   MPStateComponent TestB(s, jp, j);
   set_element(TestB[0].data(), 0, 0,  3.0 * LinearAlgebra::identity_matrix<double>(1));

   MatrixOperator TestE(i, j, QN(4));
   set_element(TestE.data(), 0, 0, 4.0 * LinearAlgebra::identity_matrix<double>(1));

   SimpleOperator TestM(sp, s, QN(5));
   TestM(0,0) = 0.166666666666666666666;
   
   MatrixOperator TestRes = operator_prod(TestM, TestA, TestE, herm(TestB), QN(6));

   CHECK_EQUAL(TestRes.TransformsAs(), QN(6));
   CHECK_EQUAL(TestRes.Basis1(), ip);
   CHECK_EQUAL(TestRes.Basis2(), jp);
   // This number is sqrt(7) * sqrt(7) * sqrt(13) * sqrt(15) * 9j(1,2,3, 4,5,6, 3,7,9)
   CHECK_CLOSE(TestRes(0,0)(0,0), 0.159719141249985);


   TestE = MatrixOperator(ip, jp, QN(6));
   set_element(TestE.data(), 0, 0, 4.0 * LinearAlgebra::identity_matrix<double>(1));

   TestRes = operator_prod(TestM, herm(TestA), TestE, TestB, QN(4));

   // This number is the previous * 19 / 7
   CHECK_CLOSE(TestRes(0,0)(0,0), 0.4335233833928164285);

#if 0
   TRACE(TestA);
   M = ExpandBasis1(TestA);
   TRACE(M)(TestA);
   TRACE(prod(M, TestA));
#endif
}

