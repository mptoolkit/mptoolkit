// -*- C++ -*- $Id$

#if !defined(DDMRG_H_HHVERTHYEHYYOIUDLE89P)
#define DDMRG_H_HHVERTHYEHYYOIUDLE89P

#include "matrixproduct/splitoperator.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/operatorstack.h"
#include "common/math_const.h"



struct Solver
{
   typedef OperatorStack<MPStateComponent> SuperblockOperator;
   typedef OperatorStack<MatrixOperator>   TransformOperator;

   Solver(CenterWavefunction const& Psi_, SplitOperator const& A_, CenterWavefunction const& Rhs_, double Freq_, double Broad_);

   virtual ~Solver();

   CenterWavefunction const& Wavefunction() const { return x; }

   int LeftSize() const { return x.LeftSize(); }
   int RightSize() const { return x.RightSize(); }

   // Calculates the wavefunction, using Iterations number
   // of Conjugate-Gradient steps.   Returns the norm of the residual.
   virtual double Solve(int Iterations) = 0;

   // Does a truncation and shifts the Center matrix to the right.
   // The new left basis is automatically expanded.
   void ShiftRightAndExpand();

   // Does a truncation and shifts the Center matrix to the left.
   // The new right basis is automatically expanded.
   void ShiftLeftAndExpand();

   // 'Expands' the left basis to cover all states.
   void ExpandLeft();

   // 'Expands' the right basis to cover all states.
   void ExpandRight();

   // returns the overlap < y | A | x >
   virtual std::complex<double> Overlap() const = 0;

   TruncationInfo TruncateLeft(int MinStates, int MaxStates, double MinTrunc, double CFactor);
   TruncationInfo TruncateRight(int MinStates, int MaxStates, double MinTrunc, double CFactor);

   // debug sanity check that the wavefunction and operator basis agree
   void DebugCheckBasis() const;

   CenterWavefunction x;             // the left hand side that we are solving for
   SplitOperator A;              // the operator (w+E-H)^2+eta^2
   CenterWavefunction y;             // the fixed right hand side

   SuperblockOperator x_A_x;  // matrix elements   | x > A < x |
   TransformOperator x_y;     // matrix elements   | x > < y |

   double Frequency;
   double Broadening;

   QuantumNumber Ident;
};

#endif
