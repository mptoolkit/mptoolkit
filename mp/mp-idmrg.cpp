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

#include "models/spin-su2.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"

namespace prog_opt = boost::program_options;

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   SuperblockMultiply(MPStateComponent const& Left_,
		      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Left, Psi, herm(Right));
   }

   MPStateComponent Left, Right;
};

inline
SuperblockMultiply::SuperblockMultiply(MPStateComponent const& Left_,
				       MPStateComponent const& Right_)
   : Left(Left_), Right(Right_)
{
}

bool ExpandL = true, ExpandR = true;

void DoIteration(MPOpComponent const& Ham, 
		 MPStateComponent& A1, MPStateComponent& B2,
		 MatrixOperator& OldCenter,
		 MatrixOperator& Center,
		 MPStateComponent& E,
		 MPStateComponent& F)
{
   DEBUG_CHECK_EQUAL(A1.Basis2(), Center.Basis1());
   DEBUG_CHECK_EQUAL(E.Basis2(), Center.Basis1());
   DEBUG_CHECK_EQUAL(E.Basis1(), Center.Basis1());

   DEBUG_CHECK_EQUAL(B2.Basis1(), Center.Basis2());
   DEBUG_CHECK_EQUAL(F.Basis1(), Center.Basis2());
   DEBUG_CHECK_EQUAL(F.Basis2(), Center.Basis2());

   // invert the singular values
   MatrixOperator OldCenterInverse = OldCenter;
   OldCenterInverse *= 0.0;
   for (unsigned i = 0; i < OldCenter.Basis1().size(); ++i)
   {
      if (iterate_at(OldCenter.data(), i, i))
      {
	 for (int j = 0; j < std::min(OldCenter.Basis1().dim(i), OldCenter.Basis2().dim(i)); ++j)
	 {
	    double x = OldCenter(i,i)(j,j).real();
	    if (fabs(x) < 1e-8) 
	       x = 0;
	    else 
	       x = 1.0 / x;
	    OldCenterInverse(i,i)(j,j) = x;
	 }
      }
   }

   // rotate to the right
   MPStateComponent A2 = prod(Center, B2);
   MatrixOperator Cr = TruncateBasis2(A2);

   // rotate to the left
   MPStateComponent B1 = prod(A1, Center);
   MatrixOperator Cl = TruncateBasis1(B1);
   
   // insert the sites
   OldCenter = Center;
   // wavefunction transform

#if 1
   Center = Cr * herm(OldCenterInverse) * Cl;  // wavefunction transformation
#else

   // alternative approach of taking the sqrt of the singular values
   MatrixOperator Ur, Dr, Vhr;
   SingularValueDecomposition(Cr, Ur, Dr, Vhr);
   MatrixOperator Ul, Dl, Vhl;
   SingularValueDecomposition(Cl, Ul, Dl, Vhl);

   // take the square root of the singular values
   for (unsigned i = 0; i < Dr.Basis1().size(); ++i)
   {
      if (iterate_at(Dr.data(), i, i))
      {
	 for (int j = 0; j < std::min(Dr.Basis1().dim(i), Dr.Basis2().dim(i)); ++j)
	 {
	    Dr(i,i)(j,j) = std::sqrt(Dr(i,i)(j,j));
	 }
      }
   }

   for (unsigned i = 0; i < Dl.Basis1().size(); ++i)
   {
      if (iterate_at(Dr.data(), i, i))
      {
	 for (int j = 0; j < std::min(Dl.Basis1().dim(i), Dl.Basis2().dim(i)); ++j)
	 {
	    Dl(i,i)(j,j) = std::sqrt(Dl(i,i)(j,j));
	 }
      }
   }

   Center = Ur * Dr * Vhr * Ul * Dl * Vhl;
#endif


   //   TRACE("transformed wavefunction")(norm_frob(Center));

   // if the wavefunction has ~zero weight, things might get screwy
   if (norm_frob(Center) < 1e-10)
   {
      Center = MakeRandomMatrixOperator(Center.Basis1(), Center.Basis2());
   }

   Center *= 1.0 / norm_frob(Center);    // normalize

   if (ExpandL)
   {
      // expand the left basis
      MatrixOperator U = ExpandBasis2(A2);
      Center = U * Center;
   }
   if (ExpandR)
   {
      // expand the right basis
      MatrixOperator U = ExpandBasis1(B1);
      Center = Center * U;
   }

   E = operator_prod(herm(Ham), herm(A2), E, A2);
   F = operator_prod(Ham, B1, F, herm(B1));
   A1 = A2;
   B2 = B1;
}

int main(int argc, char** argv)
{
   int NumIter = 10;

   SiteBlock Boundary = CreateSU2SpinSite(0.5);
   SiteBlock Site = CreateSU2SpinSite(0.5);

   double J1 = 1;
   double J2 = 0.5;

   int MinStates = 1;
   int MaxStates = 200;
   double TruncCutoff = 0;
   double EigenCutoff = 1e-10;

   std::cout.precision(12);

   StatesInfo SInfo;
   SInfo.MinStates = MinStates;
   SInfo.MaxStates = MaxStates;
   SInfo.TruncationCutoff = TruncCutoff;
   SInfo.EigenvalueCutoff = EigenCutoff;
   std::cout << SInfo << '\n';
   
   MpOpTriangular Ham = TriangularTwoSite(-sqrt(3.0)*Site["S"], 
   					  Site["S"], 
   					  Site["I"].TransformsAs());

   MpOpTriangular BoundaryHam = TriangularTwoSite(-sqrt(3.0)*Boundary["S"], 
   						  Boundary["S"], 
   						  Boundary["I"].TransformsAs());

   //MpOpTriangular Ham = ZigZagChain(-sqrt(3.0)*Site["S"], Site["S"], J1, J2);

   //MpOpTriangular BoundaryHam = ZigZagChain(-sqrt(3.0)*Boundary["S"], Boundary["S"], J1, J2);

   MPStateComponent E = Initial_E(BoundaryHam);
   MPStateComponent F = Initial_F(BoundaryHam);

   // initial 2-site block
   MatrixOperator Center = MatrixOperator::make_identity(E.Basis1());
   MPStateComponent A1 = ConstructFromLeftBasis(Boundary.Basis2().Basis(), Center.Basis1());
   MPStateComponent B2 = ConstructFromRightBasis(Boundary.Basis1().Basis(), Center.Basis2());
   E = operator_prod(herm(BoundaryHam.data()), herm(A1), E, A1);
   F = operator_prod(BoundaryHam.data(), B2, F, herm(B2));

   Center = MakeRandomMatrixOperator(A1.Basis2(), B2.Basis1());
   int Iterations = NumIter;
   double Energy = Lanczos(Center, 
			   SuperblockMultiply(E, F),
			   Iterations);
   TRACE(Energy);
   double BondEnergy = inner_prod(Center, triple_prod(E[1], Center, herm(F[1]))).real();

   // Add one unit cell away from the boundary
   A1 = ConstructFromLeftBasis(Site.Basis2().Basis(), A1.Basis2());
   B2 = ConstructFromRightBasis(Site.Basis1().Basis(), B2.Basis1());
   E = operator_prod(herm(Ham.data()), herm(A1), E, A1);
   F = operator_prod(Ham.data(), B2, F, herm(B2));

   MatrixOperator OldCenter = Center;
   Center = MakeRandomMatrixOperator(A1.Basis2(), B2.Basis1());

   Iterations = NumIter;
   Energy = Lanczos(Center, 
		    SuperblockMultiply(E, F),
		    Iterations);
   double WavefunctionDifference = 0;
   TRACE(Energy);

   for (int i = 0; i < 40000; ++i)
   {
#if 0
      if (i > 400)
      {
	 NumIter = 4;
	 SInfo.MaxStates = 100;
      }
      if (i > 600)
      {
	 NumIter = 4;
	 SInfo.MaxStates = 120;
      }
      if (i > 800)
      {
	 NumIter = 4;
	 SInfo.MaxStates = 140;
      }

      if ((i > 420 && i < 600) || (i > 620 && i < 800) || (i > 820 && i < 1000))
      {
	 ExpandL = i%2;
	 ExpandR = !ExpandL;
      }
      else
      {
	 ExpandL = ExpandR = true;
      }
#endif

      //if (i > 1800)
      //{
      //	 ExpandL = true;
      //	 ExpandR = false;
      //}

      if (i > 100)
	 NumIter = 4;

      if (i > 10000)
	 NumIter = 2;

      // truncate
      MatrixOperator RhoL = scalar_prod(Center, herm(Center));
      DensityMatrix<MatrixOperator> DML(RhoL);
      if (i % 1000 == 0 || i % 1000 == 1)
	 DML.DensityMatrixReport(std::cout);
      TruncationInfo Info;
      MatrixOperator TruncL = DML.ConstructTruncator(DML.begin(), 
						     TruncateFixTruncationErrorAbsolute(DML.begin(),
											DML.end(),
											SInfo,
											Info));

      std::cout << "E=" << Energy 
		<< " BE=" << BondEnergy
		<< " Ent=" << Info.TotalEntropy()
		<< " NS=" << Info.KeptStates() 
		<< " TError=" << Info.TruncationError()
		<< " KEigen=" << Info.SmallestKeptEigenvalue()
		<< " WDiff=" << WavefunctionDifference
		<< std::endl;

      Center = TruncL * Center;
      A1 = prod(A1, herm(TruncL));
      E = triple_prod(TruncL, E, herm(TruncL));
      
      MatrixOperator RhoR = scalar_prod(herm(Center), Center);
      DensityMatrix<MatrixOperator> DMR(RhoR);
      TruncationInfo InfoR;
      MatrixOperator TruncR = DMR.ConstructTruncator(DMR.begin(), 
						     TruncateFixTruncationErrorAbsolute(DMR.begin(),
											DMR.end(),
											SInfo,
											InfoR));
      Center = Center * herm(TruncR);
      B2 = prod(TruncR, B2);
      F = triple_prod(TruncR, F, herm(TruncR));

#if 0
      // Make the singular values positive
      MatrixOperator U, D, Vh;
      Tensor::SingularValueDecomposition(Center, U, D, Vh);
      //TRACE(Center)(D);
      Center = D;
      A1 = prod(A1, U);
      B2 = prod(Vh, B2);
      E = triple_prod(herm(U), E, U);
      F = triple_prod(Vh, F, herm(Vh));
#endif

      DoIteration(Ham.data(), A1, B2, OldCenter, Center, E, F);

      // normalize Center, so we get a proper fidelity
      Center *= 1.0 / norm_frob(Center);

      MatrixOperator CenterSave = Center;
      Iterations = NumIter;
      Energy = Lanczos(Center, 
		       SuperblockMultiply(E, F),
		       Iterations);
      WavefunctionDifference = 1.0 - LinearAlgebra::norm_frob(inner_prod(Center, CenterSave));
      //      TRACE(Energy)(inner_prod(Center, CenterSave));
      //      TRACE(norm_frob(Center))(norm_frob(CenterSave));
      //      TRACE(Center)(CenterSave);
      //      TRACE(EigenvaluesHermitian(scalar_prod(Center, herm(Center))))
      //	 (EigenvaluesHermitian(scalar_prod(CenterSave, herm(CenterSave))));
      //      TRACE(Center.Basis1());
      //      TRACE(CenterSave.Basis1());

      // bond energy
      BondEnergy = inner_prod(Center, triple_prod(E[1], Center, herm(F[1]))).real();
      //+ inner_prod(Center, triple_prod(E[2], Center, herm(F[2]))).real();
      //TRACE(BondEnergy);

      //      TRACE(inner_prod(Center, triple_prod(E[0], Center, herm(F[0]))).real());
      //      TRACE(inner_prod(Center, triple_prod(E[1], Center, herm(F[1]))).real());
      //      TRACE(inner_prod(Center, triple_prod(E[2], Center, herm(F[2]))).real());
      //      TRACE(inner_prod(Center, triple_prod(E[3], Center, herm(F[3]))).real());

      // Adjustment for the energy per site
      //      E[0] -= inner_prod(Center, triple_prod(E[0], Center, herm(F[0]))) * MatrixOperator::make_identity(Center.Basis1());
      //      F[3] -= inner_prod(Center, triple_prod(E[3], Center, herm(F[3]))) * MatrixOperator::make_identity(Center.Basis2());
      E[0] -= (0.5 * Energy) * MatrixOperator::make_identity(Center.Basis1());
      F[2] -= (0.5 * Energy) * MatrixOperator::make_identity(Center.Basis2());

      // over-relaxation step (not effective)
      //if (inner_prod(Center, CenterSave).real() < 0)
      //   Center *= -1.0;
      //MatrixOperator Delta = Center - CenterSave;
      //Center = CenterSave + 0.9 * Delta;

   }
}
