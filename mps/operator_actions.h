// -*- C++ -*- $Id$
//
// functions for operators acting on wavefunctions

#if !defined(OPERATOR_ACTIONS_H_SDH47589FOIJHO9JEW)
#define OPERATOR_ACTIONS_H_SDH47589FOIJHO9JEW

#include "mps/linearwavefunction.h"
#include "mpo/generic_mpo.h"

// Calculates the operator contraction, with a matrix
// actong on the left hand side of the wavefunction.
// Op must have 1x1 boundaries
// +-Psi1*-... Psi1*-
// |  |         |
// m  Op*- ...  Op*
// |  |         |
// +-Psi2--... Psi2--
MatrixOperator 
inject_left(MatrixOperator const& m, 
            GenericMPO const& Op, 
            LinearWavefunction const& Psi);

MatrixOperator 
inject_left(MatrixOperator const& m, 
            LinearWavefunction const& Psi,
            GenericMPO const& Op, 
            LinearWavefunction const& Psi2);

MatrixOperator 
inject_left_qshift(MatrixOperator const& m, 
                   GenericMPO const& Op, 
                   LinearWavefunction const& Psi,
                   QuantumNumber const& QShift);

// Functor to inject_left with a fixed operator and wavefunction.
// The fixed operator and wavefunction are stored by reference,
// so make sure they stay in scope while the functor is active.
struct InjectLeft
{
   InjectLeft(GenericMPO const& Op, LinearWavefunction const& Psi)
      : Op_(Op), Psi_(Psi)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return inject_left(m, Op_, Psi_);
   }
   
   GenericMPO const& Op_;
   LinearWavefunction const& Psi_;
};

// InjectLeft with a quantum number shift.
// The convention is that the initial operator is in the left-side basis,
// and the QShift transforms from the right-side basis to the left side basis.
struct InjectLeftQShift
{
   InjectLeftQShift(GenericMPO const& Op, LinearWavefunction const& Psi, QuantumNumber const& QShift)
      : Op_(Op), Psi_(Psi), QShift_(QShift)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return delta_shift(inject_left(m, Op_, Psi_), QShift_);
   }
   
   GenericMPO const& Op_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
};

// Generalization of the operator contraction
//
// +-Psi1*-     Psi1*-
// |  |          |
// E--Op*-- ...  Op*--
// |  |          |
// +-Psi2-- ... Psi2--

StateComponent 
inject_left(StateComponent const& E, 
            LinearWavefunction const& Psi1, 
            GenericMPO const& Op,
            LinearWavefunction const& Psi2);

//
// inject right variants
//

MatrixOperator 
inject_right(MatrixOperator const& m, 
             GenericMPO const& Op, 
             LinearWavefunction const& Psi);

MatrixOperator
inject_right(MatrixOperator const& m, 
             LinearWavefunction const& Psi1,
             GenericMPO const& Op, 
             LinearWavefunction const& Psi2);

MatrixOperator 
inject_right_qshift(MatrixOperator const& m, 
                    GenericMPO const& Op, 
                    LinearWavefunction const& Psi,
                    QuantumNumber const& QShift);

// Functor to inject_right with a fixed operator and wavefunction.
// The fixed operator and wavefunction are stored by reference,
// so make sure they stay in scope while the functor is active.
struct InjectRight
{
   InjectRight(GenericMPO const& Op, LinearWavefunction const& Psi)
      : Op_(Op), Psi_(Psi)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return inject_right(m, Op_, Psi_);
   }
   
   GenericMPO const& Op_;
   LinearWavefunction const& Psi_;
};

// InjectRight with a quantum number shift.
// The convention is that the initial operator is in the left-side basis,
// and the QShift transforms from the right-side basis to the left side basis.
struct InjectRightQShift
{
   InjectRightQShift(GenericMPO const& Op, LinearWavefunction const& Psi, QuantumNumber const& QShift)
      : Op_(Op), Psi_(Psi), QShift_(QShift)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return inject_right(delta_shift(m, adjoint(QShift_)), Op_, Psi_);
   }
   
   GenericMPO const& Op_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
};


/*
currently not used anywhere, but can be added

StateComponent 
inject_right(LinearWavefunction const& Psi1, 
             MPOperator const& Op,
             LinearWavefunction const& Psi2,
             StateComponent const& In);
*/


// miscellaneous functions

// verify that the local basis of an MPO matches the local basis of a wavefunction
// and that the size is correct
bool
is_local_basis_compatible(LinearWavefunction const& Psi, GenericMPO const& M);
 
// aborts with an error message if !is_local_basis_compatible(Psi,M)
void
local_basis_compatible_or_abort(LinearWavefunction const& Psi, GenericMPO const& M);


struct LeftMultiplyString
{
   typedef MatrixOperator argument_type;
   typedef MatrixOperator result_type;

   LeftMultiplyString(LinearWavefunction const& L1_, 
		      GenericMPO const& StringOp_,
		      LinearWavefunction const& L2_, 
		      QuantumNumber const& QShift_) 
      : L1(L1_), StringOp(StringOp_) , L2(L2_), QShift(QShift_)
   {
      CHECK_EQUAL(L1.size(), L2.size())
         ("The two wavefunctions must be the same length");
      CHECK_EQUAL(L1.size(), StringOp.size())
         ("The string operator must be the same length as the unit cell");
   }

   result_type operator()(argument_type const& x) const
   {
      DEBUG_CHECK_EQUAL(x.Basis1(), L1.Basis1());
      DEBUG_CHECK_EQUAL(x.Basis2(), L2.Basis1());
      //      result_type r = delta_shift(x, QShift);
      return delta_shift(inject_left(x, L1, StringOp, L2), QShift);
   }

   LinearWavefunction const& L1;
   GenericMPO StringOp;
   LinearWavefunction const& L2;
   QuantumNumber QShift;
};



#endif
