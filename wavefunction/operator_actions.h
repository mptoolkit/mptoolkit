// -*- C++ -*- $Id$
//
// functions for operators acting on wavefunctions
//
// These are designed so that the left matrices (E matrices)
// are Hermitian conjugates of the right matrices (F matrices).
//
// If we contact over part of a network using E matrices and
// part-way with F matrices, the final expectation value is
// inner_prod(E, F), or equivalently scalar_prod(herm(E), F).
//
// This means however that if we calculate eigenvalues of
// some MPO from the left, then the resulting eigenvalue
// is the conjgate of what we would expect.
  
#if !defined(OPERATOR_ACTIONS_H_SDH47589FOIJHO9JEW)
#define OPERATOR_ACTIONS_H_SDH47589FOIJHO9JEW

#include "wavefunction/linearwavefunction.h"
#include "mpo/generic_mpo.h"
#include "mpo/product_mpo.h"
#include "common/statistics.h"

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

// Calculates the operator contraction, with a matrix
// acting on the left hand side of the wavefunction.
// Op must have 1x1 boundaries
// --Psi2--... Psi2--+
//    |         |    |
//    Op- ...  Op    m
//    |         |    |
// --Psi1*-... Psi1*-+

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

struct LeftMultiplyOperator
{
   typedef StateComponent argument_type;
   typedef StateComponent result_type;

   LeftMultiplyOperator(LinearWavefunction const& L1_, QuantumNumber const& QShift1_,
			ProductMPO const& StringOp_,
			LinearWavefunction const& L2_, QuantumNumber const& QShift2_,
			int Length_ = 0)
      : Psi1(L1_), QShift1(QShift1_), StringOp(StringOp_) , Psi2(L2_), QShift2(QShift2_), Length(Length_)
   {
      if (Length == 0)
      {
	 Length = statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size());
      }
      else
      {
	 DEBUG_CHECK_EQUAL(Length, statistics::lcm(Psi1.size(), Psi2.size(), StringOp.size()));
      }
   }

   result_type operator()(argument_type const& x) const
   {
      DEBUG_CHECK_EQUAL(x.Basis1(), Psi1.Basis1());
      DEBUG_CHECK_EQUAL(x.Basis2(), Psi2.Basis1());
      StateComponent R = x;
      LinearWavefunction::const_iterator I1 = Psi1.begin();
      LinearWavefunction::const_iterator I2 = Psi2.begin();
      ProductMPO::const_iterator OpIter = StringOp.begin();
      int n = 0;
      QuantumNumber q1 = QuantumNumber(QShift1.GetSymmetryList());
      QuantumNumber q2 = QuantumNumber(QShift2.GetSymmetryList());
      while (I1 != Psi1.end() || I2 != Psi2.end() || OpIter != StringOp.end())
      {
	 if (I1 == Psi1.end())
	 {
	    I1 = Psi1.begin();
	    q1 = delta_shift(q1, QShift1);
	 }
	 if (I2 == Psi2.end())
	 {
	    I2 = Psi2.begin();
	    q2 = delta_shift(q2, QShift2);
	 }
	 if (OpIter == StringOp.end())
	 {
	    OpIter = StringOp.begin();
	 }
	 R = contract_from_left(*OpIter, herm(delta_shift(*I1, adjoint(q1))), R, 
				delta_shift(*I2, adjoint(q2)));
	 ++n;
	 ++I1;
	 ++I2;
	 ++OpIter;
      }
      q1 = delta_shift(q1, QShift1);
      q2 = delta_shift(q2, QShift2);
      CHECK_EQUAL(q1, q2);
      CHECK_EQUAL(n, Length);
      return delta_shift(R, q1);
   }

   LinearWavefunction const& Psi1;
   QuantumNumber QShift1;
   ProductMPO StringOp;
   LinearWavefunction const& Psi2;
   QuantumNumber QShift2;
   int Length;
};


#endif
