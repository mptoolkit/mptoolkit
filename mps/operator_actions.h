// -*- C++ -*- $Id$
//
// functions for operators acting on wavefunctions

#if !defined(OPERATOR_ACTIONS_H_SDH47589FOIJHO9JEW)
#define OPERATOR_ACTIONS_H_SDH47589FOIJHO9JEW

#include "mps/linearwavefunction.h"
#include "mpo/mpoperator.h"

// Calculates the operator contraction, where the operator acts
// on the right hand side of m.  Op must have 1x1 boundaries
// +-Psi*-
// |  |
// m  Op
// |  |
// +-Psi-
MatrixOperator apply_right(MatrixOperator const& m, 
			   MPOperator const& Op, 
			   LinearWavefunction const& Psi);

MatrixOperator apply_right_qshift(MatrixOperator const& m, 
				  MPOperator const& Op, 
				  LinearWavefunction const& Psi,
				  QuantumNumber const& QShift);

// Functor to ApplyRight with a fixed operator and wavefunction.
// The fixed operator and wavefunction are stored by reference,
// so make sure they stay in scope while the functor is active.
struct ApplyRight
{
   ApplyRight(MPOperator const& Op, LinearWavefunction const& Psi)
      : Op_(Op), Psi_(Psi)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return apply_right(m, Op_, Psi_);
   }
   
   MPOperator const& Op_;
   LinearWavefunction const& Psi_;
};

// ApplyRight with a quantum number shift.
// The convention is that the initial operator is in the left-side basis,
// and the QShift transforms from the right-side basis to the left side basis.
struct ApplyRightQShift
{
   ApplyRightQShift(MPOperator const& Op, LinearWavefunction const& Psi, QuantumNumber const& QShift)
      : Op_(Op), Psi_(Psi), QShift_(QShift)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return delta_shift(apply_right(m, Op_, Psi_), QShift_);
   }
   
   MPOperator const& Op_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
};


// Calculates the operator contraction, Op must have 1x1 boundaries
// -Psi--+
//   |   |
//   Op* m
//   |   |
// -Psi*-+
MatrixOperator apply_left(MatrixOperator const& m, 
			  MPOperator const& Op, 
			  LinearWavefunction const& Psi);

MatrixOperator apply_left_qshift(MatrixOperator const& m, 
				 MPOperator const& Op, 
				 LinearWavefunction const& Psi,
				 QuantumNumber const& QShift);

// Functor to ApplyLeft with a fixed operator and wavefunction.
// The fixed operator and wavefunction are stored by reference,
// so make sure they stay in scope while the functor is active.
struct ApplyLeft
{
   ApplyLeft(MPOperator const& Op, LinearWavefunction const& Psi)
      : Op_(Op), Psi_(Psi)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return apply_left(m, Op_, Psi_);
   }
   
   MPOperator const& Op_;
   LinearWavefunction const& Psi_;
};

// ApplyLeft with a quantum number shift.
// The convention is that the initial operator is in the left-side basis,
// and the QShift transforms from the right-side basis to the left side basis.
struct ApplyLeftQShift
{
   ApplyLeftQShift(MPOperator const& Op, LinearWavefunction const& Psi, QuantumNumber const& QShift)
      : Op_(Op), Psi_(Psi), QShift_(QShift)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return apply_left(delta_shift(m, adjoint(QShift_)), Op_, Psi_);
   }
   
   MPOperator const& Op_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
};


// Calculates the operator contraction
// +-Psi1*-
// |  |
// In-Op
// |  |
// +-Psi2-
StateComponent 
apply_right(StateComponent const& In, 
            LinearWavefunction const& Psi1, 
            MPOperator const& Op,
            LinearWavefunction const& Psi2);

#endif
