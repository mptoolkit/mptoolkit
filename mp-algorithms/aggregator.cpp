// -*- C++ -*- $Id$

#include "aggregator.h"

void Aggregator::ConstructLeft(std::vector<MPWavefunction> const& Psi,
                               std::vector<MPOperator> const& Op, 
                               std::vector<OperatorType>& LeftMap, 
                               std::vector<double> const& Weights,
                               int MaxStates,
                               double MinTrunc)
{
   // Construct the mapping from the vectors to the result space
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      LeftMap[i] = operator_prod(herm(Aggregate_left), LeftMap[i], Psi[i].Left());
   }
                      
   // Construct the density matrix
   OperatorType Rho;
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      OperatorType ThisPart = triple_prod(LeftMap[i], scalar_prod(Psi[i].Center(), herm(Psi[i].Center())),
                                          herm(LeftMap[i]));
      Rho += (Weights[i] / norm_frob(ThisPart)) * ThisPart;
   }

   // Form the density matrix
   DensityMatrix<OperatorType> DM(Rho);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationError(DM.begin(), DM.end(),
                                                                                   0,
                                                                                   MaxStates,
                                                                                   MinTrunc,
                                                                                   Info));
   if (ShowStates)
      std::cerr << "left density matrix at partition (" << Psi[0].LeftSize() << "," << Psi[0].RightSize() 
                << "), states=" << Info.KeptStates() << ", trunc=" << Info.TruncationError() << '\n';

   // Truncate
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      LeftMap[i] = prod(U, LeftMap[i], Ident);
   }
   Aggregate_left = Aggregate_left * herm(U);

   // construct matrix elements of H
   for (unsigned i = 0; i < Op_left.size(); ++i)
   {
      Op_left[i] = operator_prod(herm(Op[i].Left()), herm(Aggregate_left), Op_left[i], Aggregate_left);
   }
}

void Aggregator::ConstructRight(std::vector<MPWavefunction> const& Psi,
                                std::vector<MPOperator> const& Op, 
                                std::vector<OperatorType>& RightMap, 
                                std::vector<double> const& Weights,
                                int MaxStates,
                                double MinTrunc)
{
   // Construct the mapping from the vectors to the result space
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      RightMap[i] = operator_prod(Aggregate_right, RightMap[i], herm(Psi[i].Right()));
   }

   // Construct the density matrix
   OperatorType Rho;
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      OperatorType ThisPart = triple_prod(RightMap[i], scalar_prod(herm(Psi[i].Center()), Psi[i].Center()),
                                          herm(RightMap[i]));
      Rho += (Weights[i] / norm_frob(ThisPart)) * ThisPart;
   }

   // Form the density matrix
   DensityMatrix<OperatorType> DM(Rho);
   TruncationInfo Info;
   MatrixOperator U = DM.ConstructTruncator(DM.begin(), TruncateFixTruncationError(DM.begin(), DM.end(),
                                                                                   0,
                                                                                   MaxStates,
                                                                                   MinTrunc,
                                                                                   Info));
   if (ShowStates)
      std::cerr << "right density matrix at partition (" << Psi[0].LeftSize() << "," << Psi[0].RightSize() 
                << "), states=" << Info.KeptStates() << ", trunc=" << Info.TruncationError() << '\n';

   // Truncate
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      RightMap[i] = prod(U, RightMap[i], Ident);
   }
   Aggregate_right = prod(U, Aggregate_right);

   // construct matrix elements of H
   for (unsigned i = 0; i < Op_right.size(); ++i)
   {
      Op_right[i] = operator_prod(Op[i].Right(), Aggregate_right, Op_right[i], herm(Aggregate_right));
   }
}

void Aggregator::RotateLeft(std::vector<MPWavefunction>& Psi,  std::vector<MPOperator>& Op)
{
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      Psi[i].RotateLeft();
   }
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Op[i].RotateLeft();
   }
   Aggregate_right = Component::ConstructFullBasis1(Psi[0].Right().SiteBasis(),
                                                    Aggregate_right.Basis1());
}

void Aggregator::RotateRight(std::vector<MPWavefunction>& Psi, std::vector<MPOperator>& Op)
{
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      Psi[i].RotateRight();
   }
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Op[i].RotateRight();
   }
   Aggregate_left = Component::ConstructFullBasis2(Aggregate_left.Basis2(),
                                                   Psi[0].Left().SiteBasis());
}

Aggregator::Aggregator(std::vector<MPWavefunction> Psi, 
                       std::vector<MPOperator> Op, 
                       bool ShowStates_,
                       int MaxStates,
                       double MinTrunc,
                       int Location,
                       std::vector<double> Weights)
   : ShowStates(ShowStates_)
{
   // default weights if they were not supplied
   if (Weights.empty())
      Weights = std::vector<double>(Psi.size(), 1.0);

   // Rotate the wavefunctions to left most position
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      while (Psi[i].LeftSize() > 1)
         Psi[i].RotateLeft();
   }
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      while (Op[i].LeftSize() > 1)
         Op[i].RotateLeft();
   }

   // Initialize the result matrices
   OperatorType LeftVac = OperatorType::make_identity(Psi[0].LeftVacuumBasis());
   std::vector<OperatorType> LeftMap(Psi.size(), LeftVac);

   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Op_left.push_back(make_vacuum_state(Psi[0].LookupLeft(0).Basis1()[0]));
   }

   Aggregate_left = Component::ConstructFullBasis2(Psi[0].Left().Basis1(),
                                                   Psi[0].Left().SiteBasis());

   // Construct the left density matrices up to the center position
   this->ConstructLeft(Psi, Op, LeftMap, Weights, MaxStates, MinTrunc);
   while (Psi[0].LeftSize() < Location)
   {
      this->RotateRight(Psi, Op);
      this->ConstructLeft(Psi, Op, LeftMap, Weights, MaxStates, MinTrunc);
   }

   // now shift all the way to the right and construct the right operators
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      while (Psi[i].RightSize() > 1)
         Psi[i].RotateRight();
   }
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      while (Op[i].RightSize() > 1)
         Op[i].RotateRight();
   }

   OperatorType RightVac = OperatorType::make_identity(Psi[0].RightVacuumBasis());
    std::vector<OperatorType> RightMap(Psi.size(), RightVac);

   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Op_right.push_back(make_vacuum_state(Psi[0].GetSymmetryList()));
   }

   Aggregate_right = Component::ConstructFullBasis1(Psi[0].Right().SiteBasis(), 
                                                    Psi[0].Right().Basis2());


   // Construct the right density matrices up to the center position
   this->ConstructRight(Psi, Op, RightMap, Weights, MaxStates, MinTrunc);
   while (Psi[0].LeftSize() > Location)
   {
      this->RotateLeft(Psi, Op);
      this->ConstructRight(Psi, Op, RightMap, Weights, MaxStates, MinTrunc);
   }

   // Now construct the center matrices
   for (unsigned i = 0; i < Psi.size(); ++i)
   {
      Aggregate_center.push_back(triple_prod(LeftMap[i], Psi[i].Center(), herm(RightMap[i])));
   }

   // Operator center matrices
   for (unsigned i = 0; i < Op.size(); ++i)
   {
      Op_center.push_back(Op[i].Center());
   }
}

std::complex<double>
Aggregator::Expectation(MatrixOperator const& Psi, int n) const
{
   return inner_prod(operator_prod(conj(Op_center[n]), Op_left[n], Psi, herm(Op_right[n])),
                     Psi);
}

std::complex<double>
Aggregator::Expectation(MatrixOperator const& Psi1, int n, MatrixOperator const& Psi2) const
{
   return inner_prod(operator_prod(conj(Op_center[n]), Op_left[n], Psi1, herm(Op_right[n])),
                     Psi2);
}

MatrixOperator
Aggregator::Apply(int n, MatrixOperator const& Psi) const
{
   return operator_prod(conj(Op_center[n]), Op_left[n], Psi, herm(Op_right[n]));
}
