// -*- C++ -*- $Id$

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include "matrixproduct/cg.h"
#include "matrixproduct/gmres.h"
#include "mp/copyright.h"
#include "tensor/tensor_eigen.h"
#include "tensor/regularize.h"
#include <iostream>

class SolverGmres : public Solver
{
   public:
      SolverGmres(MPWavefunction const& Psi_, MPOperator const& Op_, MPWavefunction const& Rhs_)
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
		      MPStateComponent const& Right_)
      : Op(Op_), Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Op, Left, Psi, herm(Right));
   }

   SimpleOperator Op;
   MPStateComponent Left, Right;
};

struct Precondition
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator const& argument_type;

   Precondition(MPStateComponent const& Left_, MPStateComponent const& Right_, 
                MatrixOperator const& Ei, MatrixOperator const& Fi)
      : yprime_A_x_Left(Left_), yprime_A_x_Right(Right_), EInv(Ei), FInv(Fi) {}

   MatrixOperator operator()(MatrixOperator const& r) const
   {
      MatrixOperator Temp = operator_prod(herm(yprime_A_x_Left), r, yprime_A_x_Right);
      return triple_prod(EInv, Temp, herm(FInv));
   }

   MPStateComponent const& yprime_A_x_Left;
   MPStateComponent const& yprime_A_x_Right;
   MatrixOperator const& EInv;
   MatrixOperator const& FInv;
};

struct RhsInner
{
   typedef std::complex<double> result_type;
   typedef MatrixOperator const& first_argument_type;
   typedef MatrixOperator const& second_argument_type;

   RhsInner(MatrixOperator const& Left_, MatrixOperator const& Right_)
      : Left(Left_), Right(Right_) {}

   result_type operator()(MatrixOperator const& l, MatrixOperator const& r) const
   {
      return inner_prod(triple_prod(Left, l, herm(Right)), r);
   }

   MatrixOperator const& Left;
   MatrixOperator const& Right;
};

struct LhsInner
{
   typedef std::complex<double> result_type;
   typedef MatrixOperator const& first_argument_type;
   typedef MatrixOperator const& second_argument_type;

   LhsInner(MatrixOperator const& Left_, MatrixOperator const& Right_)
      : Left(Left_), Right(Right_) {}

   result_type operator()(MatrixOperator const& r, MatrixOperator const& l) const
   {
      return conj(inner_prod(triple_prod(Left, l, herm(Right)), r));
   }

   MatrixOperator const& Left;
   MatrixOperator const& Right;
};

double SolverGmres::Solve(int MaxIterations)
{
   //   TRACE(Psi.Center())(norm_frob(Psi.Center()));
   //   DEBUG_TRACE(Psi.Center().Basis1())(Psi.Center().Basis2());

   int Iter = MaxIterations;
   double Tol = 1E-10;
   int m = 20;

   MatrixOperator E = scalar_prod(herm(yprime_A_x.Left()), yprime_A_x.Left());
   MatrixOperator F = scalar_prod(herm(yprime_A_x.Right()), yprime_A_x.Right());
   E += MatrixOperator::make_identity(E.Basis1());
   F += MatrixOperator::make_identity(F.Basis1());

   MatrixOperator UE = Tensor::Regularize(E.Basis1());
   E = triple_prod(UE, E, herm(UE));
   InvertHPD(E);
   E = triple_prod(herm(UE), E, UE);

   MatrixOperator UF = Tensor::Regularize(F.Basis1());
   F = triple_prod(UF, F, herm(UF));
   InvertHPD(F);
   F = triple_prod(herm(UF), F, UF);

   GmRes(x.Center(), 
         SuperblockMultiply(conj(A.Center()),
                            yprime_A_x.Left(), 
                            yprime_A_x.Right()),
         yprime.Center(),
         m, Iter, Tol,
         Precondition(yprime_A_x.Left(), yprime_A_x.Right(), E, F));

   TRACE(Iter)(Tol);
   return Tol;
}

std::complex<double> SolverGmres::Overlap() const
{
   return inner_prod(operator_prod(conj(A.Center()), 
                                   yprime_A_x.Left(), 
                                   x.Center(), 
                                   herm(yprime_A_x.Right())), 
                     yprime.Center());
}


int main(int argc, char** argv)
{
   if (argc < 5 || argc > 8)
   {
      print_copyright(std::cerr, "tools", basename(argv[0]));
      std::cerr << "GMRES solver A|x> = |y>\n";
      std::cerr << "usage: mp-gmres <lattice> <operator> <psi> <rhs> [num-iter] [<maxstates>] [corr]\n";
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

   SolverGmres solver(*Psi, (*System)[Operator], *Rhs);

   std::cout.precision(14);

   bool First = true;
   for (int Sweeps = 0; Sweeps < 5; ++Sweeps)
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
      //      solver.ExpandRight();
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
      //      solver.ExpandLeft();
      double E = solver.Solve(NumIter);
      TruncationInfo States = solver.TruncateRight(MaxStates, Correction);

      std::cout << '(' << solver.LeftSize() << ',' << solver.RightSize()
		<< ") " << E << ' ' << States.m << '\n';
   }

   }

   *Psi.mutate() = solver.Wavefunction();

   pheap::ShutdownPersistent(Psi);
}
