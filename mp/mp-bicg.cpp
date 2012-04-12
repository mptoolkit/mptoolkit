// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "mp-algorithms/cg.h"
#include "mp-algorithms/biconjugategradient.h"
#include "mp/copyright.h"
#include <iostream>

class SolverBiconj : public Solver
{
   public:
      SolverBiconj(MPWavefunction const& Psi_, MPOperator const& Op_, 
                   MPWavefunction const& Rhs_)
         : Solver(Psi_, Op_, Rhs_) {}

      virtual double Solve(int Iterations);

      virtual std::complex<double> Overlap() const;
};

struct SuperblockMultiply
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiply(SimpleOperator const& Op_,
		      MPStateComponent const& Left_,
		      MPStateComponent const& Right_);

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op, Left, Psi, herm(Right));
   }

   SimpleOperator const& Op;
   MPStateComponent const& Left;
   MPStateComponent const& Right;
};

inline
SuperblockMultiply::SuperblockMultiply(SimpleOperator const& Op_,
				       MPStateComponent const& Left_,
				       MPStateComponent const& Right_)
   : Op(Op_), Left(Left_), Right(Right_)
{
}

struct SuperblockMultiplyHerm
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   SuperblockMultiplyHerm(SimpleOperator const& Op_,
                          MPStateComponent const& Left_,
                          MPStateComponent const& Right_)
      : Op(Op_), Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op, herm(Left), Psi, Right);
   }

   SimpleOperator const& Op;
   MPStateComponent const& Left;
   MPStateComponent const& Right;
};

struct MapRhsToLhs
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   MapRhsToLhs(MatrixOperator const& Left_, MatrixOperator const& Right_)
      : Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& r) const
   {
      return triple_prod(herm(Left), r, Right);
   }

   MatrixOperator const& Left;
   MatrixOperator const& Right;
};

struct MapLhsToRhs
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   MapLhsToRhs(MatrixOperator const& Left_, MatrixOperator const& Right_)
      : Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& r) const
   {
      return triple_prod(Left, r, herm(Right));
   }

   MatrixOperator const& Left;
   MatrixOperator const& Right;
};

std::complex<double> SolverBiconj::Overlap() const
{
   return inner_prod(operator_prod(conj(A.Center()), 
                                   yprime_A_x.Left(), 
                                   x.Center(), 
                                   herm(yprime_A_x.Right())), 
                     yprime.Center());
}

double SolverBiconj::Solve(int MaxIterations)
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   int Iter = MaxIterations;
   double Tol = 1E-10;

   BiConjugateGradient(x.Center(), 
                       SuperblockMultiply(conj(A.Center()),
                                          yprime_A_x.Left(), 
                                          yprime_A_x.Right()),
                       SuperblockMultiplyHerm(conj(A.Center()),
                                              yprime_A_x.Left(), 
                                              yprime_A_x.Right()),
                       yprime.Center(),
                       Iter, Tol,
                       MapRhsToLhs(yprime_x.Left(), yprime_x.Right()),
                       MapLhsToRhs(yprime_x.Left(), yprime_x.Right()));


   //   TRACE(Iter)(Tol);
   return Tol;
}

int main(int argc, char** argv)
{
   if (argc < 5 || argc > 8)
   {
      print_copyright(std::cerr);
      std::cerr << "usage: mp-bicg <lattice> <operator> <psi> <rhs> [num-iter] [<maxstates>] [corr]\n";
      return 1;
   }

   TRACE(argc);

   pvalue_ptr<MPWavefunction> Psi = pheap::OpenPersistent(argv[3], 655360);
   pvalue_ptr<OperatorList> System = pheap::ImportHeap(argv[1]);
   std::string Operator = argv[2];
   pvalue_ptr<MPWavefunction> Rhs = pheap::ImportHeap(argv[4]);
   int NumIter = 4;
   if (argc >= 6) NumIter = boost::lexical_cast<int>(argv[5]);
   int MaxStates = DefaultMaxStates;
   if (argc >= 7) MaxStates = boost::lexical_cast<int>(argv[6]);
   double Correction = 0;
   if (argc >= 8) Correction = boost::lexical_cast<double>(argv[7]);

   SolverBiconj solver(*Psi, (*System)[Operator], *Rhs);

   std::cout.precision(14);

   bool First = true;
   for (int Sweeps = 0; Sweeps < 4; ++Sweeps)
   {

   {
      solver.ExpandLeft();
      double E = First ? 0.0 : solver.Solve(NumIter);

      TruncationInfo States = solver.TruncateLeft(MaxStates, Correction);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << '\n';
   }

   // sweep right
   while (solver.RightSize() > 1)
   {
      solver.ShiftRightAndExpand();
       solver.ExpandRight();
      double E = First ? 0.0 : solver.Solve(NumIter);

      TruncationInfo States = solver.TruncateLeft(MaxStates, Correction);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << '\n';
   }
   First = false;

   {
      solver.ExpandRight();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MaxStates, Correction);
      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << '\n';
   }

   // sweep left
   while (solver.LeftSize() > 1)
   {
      solver.ShiftLeftAndExpand();
       solver.ExpandLeft();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MaxStates, Correction);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << '\n';
   }

   }

   *Psi.mutate() = solver.Wavefunction();

   pheap::ShutdownPersistent(Psi);
}
