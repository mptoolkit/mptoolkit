// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
// Distributed under the GNU General Public License, version 3 or later.
#include "mp-algorithms/expansion.h"
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
int Checks = 0;

void require(bool Condition, std::string const& Message)
{
   ++Checks;
   if (!Condition) throw std::runtime_error(Message);
}

void oversampling_tests()
{
   require(OversamplingInfo(10,1.0)(8) == 18,"additive oversampling must not be capped by a unit factor");
   require(OversamplingInfo(0,2.0)(8) == 16,"multiplicative oversampling must not be capped by zero additive allowance");
   require(OversamplingInfo(3,1.5)(20) == 30,"larger multiplicative allowance");
   require(OversamplingInfo(3,1.5)(3) == 6,"larger additive allowance");
   require(OversamplingInfo(0,1.5)(3) == 5,"multiplicative allowance rounds to the nearest integer");
   require(OversamplingInfo()(8) == 8,"default oversampling leaves the count unchanged");
}

BasisList null_basis(int Size)
{
   SymmetryList SL("Null:Null");
   BasisList B(SL);
   for (int i = 0; i < Size; ++i) B.push_back(QuantumNumber(SL));
   return B;
}

void growth_test(bool Left, PostExpansionAlgorithm Algo, int Requested, int Extra, int Expected)
{
   auto P = null_basis(4), B = null_basis(6);
   VectorBasis V(P.GetSymmetryList());
   V.push_back(QuantumNumber(P.GetSymmetryList()),1);
   StateComponent C(P,V,V), Env(B,V,V);
   C[0] = MatrixOperator::make_identity(V);
   for (int a = 0; a < 6; ++a) Env[a] = MatrixOperator::make_identity(V);
   OperatorComponent H(P,B,B);
   H(0,0) = H(1,1) = H(5,5) = SimpleOperator::make_identity(P);
   // Three independent channels connect the occupied state to its complement.
   for (int a = 1; a < 4; ++a)
   {
      SimpleOperator X(P,P);
      X(a,0) = X(0,a) = 1.;
      H(a+1,a+1) = X;
   }
   auto Original = C;
   StatesInfo States(Requested);
   TruncationInfo Info;
   auto Lambda = Left
      ? TruncateExpandBasis1(C,Env,H,Env,Algo,States,Extra,0,Info,OversamplingInfo(10,1.0))
      : TruncateExpandBasis2(C,Env,H,Env,Algo,States,Extra,0,Info,OversamplingInfo(10,1.0));
   auto Basis = Left ? C.Basis1() : C.Basis2();
   auto Name = std::string(Left ? "left " : "right ")+Algo.Name()
      +" requested="+std::to_string(Requested)+" extra="+std::to_string(Extra);
   require(Basis.total_dimension() == Expected,
           Name+" expanded dimension="+std::to_string(Basis.total_dimension())+" expected="+std::to_string(Expected));
   require(Info.KeptStates()+Info.ExtraStates() == Expected,Name+" state accounting");
   require(norm_frob((Left ? prod(Lambda,C) : prod(C,Lambda))-Original) < 1e-12,Name+" preserves the input tensor");
   auto Gram = Left ? scalar_prod(C,herm(C)) : scalar_prod(herm(C),C);
   require(norm_frob(Gram-MatrixOperator::make_identity(Basis)) < 1e-12,Name+" orthonormal expanded basis");
}

void growth_tests()
{
   for (bool Left : {false,true})
   {
      for (auto Algo : {PostExpansionAlgorithm::RSVD,PostExpansionAlgorithm::RangeFinding,PostExpansionAlgorithm::Mixing})
      {
         growth_test(Left,Algo,3,1,4);
         growth_test(Left,Algo,1,1,2);
         growth_test(Left,Algo,8,1,4);
         growth_test(Left,Algo,3,0,1);
      }
      growth_test(Left,PostExpansionAlgorithm::NoExpansion,3,1,1);
   }
}
}

int main(int argc, char** argv)
{
   try
   {
      if (argc != 2) throw std::runtime_error("Expected oversampling or growth");
      std::string Mode(argv[1]);
      if (Mode == "oversampling") oversampling_tests();
      else if (Mode == "growth") growth_tests();
      else throw std::runtime_error("Unknown expansion test mode");
      std::cout << "PASS " << Checks << " " << Mode << " checks\n";
   }
   catch (std::exception const& Error)
   {
      std::cerr << "FAIL: " << Error.what() << '\n';
      return 1;
   }
}
