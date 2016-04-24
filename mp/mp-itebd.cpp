// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/mpoperator.h"
#include "matrixproduct/mpoperatorlist.h"
#include "matrixproduct/wavefunc-utils.h"
#include "matrixproduct/triangularoperator.h"
#include "quantumnumbers/all_symmetries.h"
#include "mp-algorithms/lanczos.h"
#include "pheap/pheap.h"
#include "mp/copyright.h"
#include "interface/operator-parser.h"
#include "common/environment.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>
#include <iostream>
#include "common/environment.h"
#include "matrixproduct/mpexponential.h"
#include "interface/inittemp.h"

#include "models/spin-su2.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include "matrixproduct/infinitewavefunction.h"

namespace prog_opt = boost::program_options;

void RandomizeRowSigns(MatrixOperator& TruncL)
{
   for (unsigned j = 0; j < TruncL.Basis1().size(); ++j)
   {
      for (int k = 0; k < TruncL.Basis1().dim(j); ++k)
      {
	 int sign = ((rand() / 1000) % 2) * 2 - 1;
	 for (unsigned l = 0; l < TruncL.Basis2().size(); ++l)
	 {
	    if (iterate_at(TruncL.data(), j,l))
	    {
	       for (int m = 0; m < TruncL.Basis2().dim(l); ++m)
	       {
		  TruncL(j,l)(k,m) *= sign;
	       }
	    }
	 }
      }
   }
}   

TruncationInfo
DoIteration(MPStateComponent& A, MatrixOperator& Center, MPStateComponent& B,
		 MPOpComponent const& H_A, MPOpComponent const& H_B,
	    StatesInfo const& SInfo, bool ShowReport = false)
{
   //   TRACE(norm_frob_sq(Center));

   // expand left and right
   {
      MatrixOperator U = ExpandBasis2(A);
      Center = U * Center;
   }
   {
      MatrixOperator U = ExpandBasis1(B);
      Center = Center * U;
   }

   //   TRACE(norm_frob_sq(Center));

   // construct the E matrices for the evolution operator
   MPStateComponent E(make_vacuum_basis(A.GetSymmetryList()), A.Basis1(), A.Basis1());
   E[0] = MatrixOperator::make_identity(A.Basis1());
   E = operator_prod(herm(H_A), herm(A), E, A);

   MPStateComponent F(make_vacuum_basis(A.GetSymmetryList()), B.Basis2(), B.Basis2());
   F[0] = MatrixOperator::make_identity(B.Basis2());
   F = operator_prod(H_B, B, F, herm(B));

   // if the wavefunction has ~zero weight, things might get screwy
   if (norm_frob(Center) < 1e-10)
   {
      Center = MakeRandomMatrixOperator(Center.Basis1(), Center.Basis2());
      Center *= 1.0 / norm_frob(Center);    // normalize
   }

   // apply the evolution operator
   MatrixOperator Old = Center;
   Center = operator_prod(E, Center, herm(F));

   //   TRACE(inner_prod(Old, Center));

   //TRACE(Center);

   // normalize
   Center *= 1.0 / norm_frob(Center);

   //TRACE(Center);

   // truncate
      MatrixOperator RhoL = scalar_prod(Center, herm(Center));
      DensityMatrix<MatrixOperator> DML(RhoL);

      //DML.DensityMatrixReport(std::cerr);
      TruncationInfo Info;
      MatrixOperator TruncL = DML.ConstructTruncator(DML.begin(), 
						     TruncateFixTruncationErrorAbsolute(DML.begin(),
											DML.end(),
											SInfo,
											Info));

      if (ShowReport)
	 DML.DensityMatrixReport(std::cout);

      //RandomizeRowSigns(TruncL);
      //TRACE(Info.TruncationError());
      //TRACE(TruncL.Basis1())(TruncL.Basis2());

      Center = TruncL * Center;
      A = prod(A, herm(TruncL));

#if 1
      MatrixOperator RhoR = scalar_prod(herm(Center), Center);
      DensityMatrix<MatrixOperator> DMR(RhoR);
      TruncationInfo InfoR;
      MatrixOperator TruncR = DMR.ConstructTruncator(DMR.begin(), 
						     TruncateFixTruncationErrorAbsolute(DMR.begin(),
											DMR.end(),
											SInfo,
											InfoR));
#else
      MatrixOperator TruncR = TruncL;
      //RandomizeRowSigns(TruncR);
      //TRACE(Center*herm(TruncL))(Center*herm(TruncR));
#endif

      Center = Center * herm(TruncR);
      B = prod(TruncR, B);

      //TRACE(Center);

#if 1
      // Make the singular values positive
   MatrixOperator U, D, Vh;
   SingularValueDecomposition(Center, U, D, Vh);
   A = prod(A, U);
   B = prod(Vh, B);
   Center = D;
#endif

   //TRACE(norm_frob_sq(Center));

   // normalize
   Center *= 1.0 / norm_frob(Center);

   return Info;
}

InfiniteWavefunction MakeWavefunction(MPStateComponent const& A, MatrixOperator const& Lambda2,
				      MPStateComponent const& B, MatrixOperator const& Lambda1)
{
   // Convert to an InfiniteWavefunction 
   InfiniteWavefunction Psi;
   Psi.C_old = Lambda1;
   Psi.Psi.push_back(A);
   MPStateComponent BNew = prod(Lambda2, B);
   MatrixOperator Lambda2New = TruncateBasis2(BNew);
   Psi.Psi.push_back(BNew);
   Psi.C_right = Lambda2New;
   Psi.QShift = QuantumNumbers::QuantumNumber(Lambda2.GetSymmetryList());
   return Psi;
}

int main(int argc, char** argv)
{
   mp_pheap::InitializeTempPHeap(false);

   // set up the lattice
   SiteBlock Site = CreateSU2SpinSite(0.5);
   double J = 1;
   double DeltaTau = 0.2;
   double DeltaT = 0.0;
   double Epsilon = 1e-7;
   double Eps = Epsilon;
   //   MpOpTriangular Ham = TriangularTwoSite(-sqrt(3.0)*Site["S"], 
   //   					  Site["S"], 
   //   					  Site["I"].TransformsAs());
   MPOpComponent H_A, H_B;
   std::tie(H_A, H_B) = TwoSiteExponential(-sqrt(3.0) * Site["S"], Site["S"],
					     std::complex<double>(-DeltaTau, -DeltaT));

   MPOpComponent half_H_A, half_H_B;
   std::tie(half_H_A, half_H_B) = TwoSiteExponential(-sqrt(3.0) * Site["S"], Site["S"],
						       std::complex<double>(-DeltaTau/2.0, -DeltaT/2.0));

   //TRACE(H_A)(H_B);

   int MinStates = 1;
   int MaxStates = 200;
   double TruncCutoff = 0;
   double EigenCutoff = 0; // 1e-10;

   StatesInfo SInfo;
   SInfo.MinStates = MinStates;
   SInfo.MaxStates = MaxStates;
   SInfo.TruncationCutoff = TruncCutoff;
   SInfo.EigenvalueCutoff = EigenCutoff;
   std::cout.precision(12);
   std::cout << SInfo << '\n';
   
   // initial state
   VectorBasis B1 = VectorBasis(make_vacuum_basis(Site.GetSymmetryList()));
   VectorBasis B2 = VectorBasis(Site.Basis1().Basis());

   MPStateComponent A(Site.Basis1().Basis(), B1, B2);
   A[0] = MatrixOperator(B1, B2, Site.Basis1().Basis()[0]);
   A[0](0,0) = LinearAlgebra::Matrix<double>(1,1,1.0);

   MPStateComponent B(Site.Basis1().Basis(), B2, B1);
   B[0] = MatrixOperator(B2, B1, adjoint(Site.Basis1().Basis()[0]));
   B[0](0,0) = LinearAlgebra::Matrix<double>(1,1,1.0);
   
   MatrixOperator Lambda1 = MatrixOperator::make_identity(B.Basis2());
   MatrixOperator Lambda2 = MatrixOperator::make_identity(A.Basis2());
   MatrixOperator Lambda1Inv = InvertDiagonal(Lambda1, Eps), Lambda2Inv = InvertDiagonal(Lambda2, Eps);

   double E1, E2;


      Lambda2 = Lambda2Inv;
      DoIteration(A, Lambda2, B, half_H_A, half_H_B, SInfo);
      Lambda2Inv = InvertDiagonal(Lambda2, Eps);

      // energy
      {
	 MatrixOperator E = operator_prod(herm(SimpleOperator(-sqrt(3.0)*Site["S"])), herm(A), A);
	 MatrixOperator F = operator_prod(SimpleOperator(Site["S"]), B, herm(B));
	 double Energy = inner_prod(triple_prod(herm(E), Lambda2, F), Lambda2).real();
	 E1 = Energy;
	 //TRACE(Energy);
      }

      A = prod(A, Lambda2);
      B = prod(Lambda2, B);

      Lambda1 = Lambda1Inv;
      DoIteration(B, Lambda1, A, H_A, H_B, SInfo);
      Lambda1Inv = InvertDiagonal(Lambda1, Eps);

      // energy
      {
	 MatrixOperator E = operator_prod(herm(SimpleOperator(-sqrt(3.0)*Site["S"])), herm(B), B);
	 MatrixOperator F = operator_prod(SimpleOperator(Site["S"]), A, herm(A));
	 double Energy = inner_prod(triple_prod(herm(E), Lambda1, F), Lambda1).real();
	 E2 = Energy;
	 //TRACE(Energy);
      }

      B = prod(B, Lambda1);
      A = prod(Lambda1, A);


      TruncationInfo Info;

   int NumIter = 200;
   for (int iter = 0; iter < NumIter; ++iter)
   {
      //      if (iter == 2000)
      {
	 DeltaTau *= 1.0 - 1e-4;
	 std::tie(H_A, H_B) = TwoSiteExponential(-sqrt(3.0) * Site["S"], Site["S"],
						   std::complex<double>(-DeltaTau, -DeltaT));
      }

      Lambda2 = Lambda2Inv;

      // energy
      {
	 MatrixOperator E = operator_prod(herm(SimpleOperator(-sqrt(3.0)*Site["S"])), herm(A), A);
	 MatrixOperator F = operator_prod(SimpleOperator(Site["S"]), B, herm(B));
	 double Energy = inner_prod(triple_prod(herm(E), Lambda2, F), Lambda2).real();
	 E1 = Energy;
	 //TRACE(Energy);
      }

      if (iter % 1000 == 0)
      {
	 // to a half-step, to simulate 2nd order S-T decomposition
	 MPOpComponent H_Ax, H_Bx;
	 std::tie(H_Ax, H_Bx) = TwoSiteExponential(-sqrt(3.0) * Site["S"], Site["S"],
						   std::complex<double>(-DeltaTau*0.5/(1.0-1e-4), 0));
	 SInfo.MaxStates = 100;
	 MPStateComponent Ax = A, Bx = B;
	 MatrixOperator L2x = Lambda2;
	 Info = DoIteration(Ax, L2x, Bx, H_Ax, H_Bx, SInfo, iter % 1000 == 0);

	 std::string Filename = "tebdx.psi." + boost::lexical_cast<std::string>(iter);
	 pvalue_ptr<InfiniteWavefunction> Ps = new InfiniteWavefunction(MakeWavefunction(Ax, L2x, Bx, Lambda1));
	 pheap::ExportHeap(Filename, Ps);
      }

      SInfo.MaxStates = 100;
      Info = DoIteration(A, Lambda2, B, H_A, H_B, SInfo, iter % 1000 == 0);
      Eps = std::max(Epsilon, Info.LargestDiscardedEigenvalue());
      Lambda2Inv = InvertDiagonal(Lambda2, Eps);

      // energy
      {
	 MatrixOperator E = operator_prod(herm(SimpleOperator(-sqrt(3.0)*Site["S"])), herm(A), A);
	 MatrixOperator F = operator_prod(SimpleOperator(Site["S"]), B, herm(B));
	 double Energy = inner_prod(triple_prod(herm(E), Lambda2, F), Lambda2).real();
#if 0
	 if ()
	 {
	    std::tie(H_A, H_B) = TwoSiteExponential(-sqrt(3.0) * Site["S"], Site["S"],
						      std::complex<double>(0.0, 0.0));
	 }
	 else
	 {
	    //	    DeltaTau *= 0.5;
	    // TRACE("Reducing DeltaTau")(DeltaTau);
	    std::tie(H_A, H_B) = TwoSiteExponential(-sqrt(3.0) * Site["S"], Site["S"],
						      std::complex<double>(-DeltaTau, -DeltaT));
	 }
#endif
	 E1 = Energy;
	 //TRACE(Energy);

	 std::cout << " BE=" << Energy
		   << " Ent=" << Info.TotalEntropy()
		   << " NS=" << Info.KeptStates() 
		   << " TError=" << Info.TruncationError()
		   << " KEigen=" << Info.SmallestKeptEigenvalue()
		   << " DeltaT=" << DeltaTau
		   << std::endl;
      }

      A = prod(A, Lambda2);
      B = prod(Lambda2, B);

      Lambda1 = Lambda1Inv;


      // energy - is this stable???
      {
	 MatrixOperator E = operator_prod(herm(SimpleOperator(-sqrt(3.0)*Site["S"])), herm(B), B);
	 MatrixOperator F = operator_prod(SimpleOperator(Site["S"]), A, herm(A));
	 double Energy = inner_prod(triple_prod(herm(E), Lambda1, F), Lambda1).real();
	 E2 = Energy;
	 //TRACE(Energy);
      }

      SInfo.MaxStates = 100;
      Info = DoIteration(B, Lambda1, A, H_A, H_B, SInfo, iter % 1000 == 0);
      Eps = std::max(Epsilon, Info.LargestDiscardedEigenvalue());
      //TRACE(Eps);
      Lambda1Inv = InvertDiagonal(Lambda1, Eps);
      //      TRACE(norm_frob_sq(Lambda1))(norm_frob_sq(Lambda1Inv));

      // energy
      {
	 MatrixOperator E = operator_prod(herm(SimpleOperator(-sqrt(3.0)*Site["S"])), herm(B), B);
	 MatrixOperator F = operator_prod(SimpleOperator(Site["S"]), A, herm(A));
	 double Energy = inner_prod(triple_prod(herm(E), Lambda1, F), Lambda1).real();
	 E2 = Energy;
	 //TRACE(Energy);

	 std::cout << " BE=" << Energy
		   << " Ent=" << Info.TotalEntropy()
		   << " NS=" << Info.KeptStates() 
		   << " TError=" << Info.TruncationError()
		   << " KEigen=" << Info.SmallestKeptEigenvalue()
		   << " DeltaT=" << DeltaTau
		   << std::endl;
      }

      B = prod(B, Lambda1);
      A = prod(Lambda1, A);
   }
}



