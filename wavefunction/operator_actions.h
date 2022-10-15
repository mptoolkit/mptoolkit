// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/operator_actions.h
//
// Copyright (C) 2015-2022 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Reseach publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER
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

#if !defined(MPTOOLKIT_WAVEFUNCTION_OPERATOR_ACTIONS_H)
#define MPTOOLKIT_WAVEFUNCTION_OPERATOR_ACTIONS_H

#include "wavefunction/linearwavefunction.h"
#include "mpo/generic_mpo.h"
#include "mpo/product_mpo.h"
#include "common/statistics.h"

// Calculates the operator contraction, with a matrix
// actong on the left hand side of the wavefunction.
// Op must have 1x1 boundaries
// +-Psi1*-...-Psi1*-
// |  |         |
// m  Op---...--Op
// |  |         |
// +-Psi2--...-Psi2--
MatrixOperator
inject_left(MatrixOperator const& m,
            GenericMPO const& Op,
            LinearWavefunction const& Psi);

MatrixOperator
inject_left(MatrixOperator const& m,
            LinearWavefunction const& Psi1,
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
      : Psi1_(Psi), QShift_(QShift), Op_(Op), Psi2_(Psi), Scale_(1.0)
   {
   }

   InjectLeftQShift(LinearWavefunction const& Psi1, QuantumNumber const& QShift, GenericMPO const& Op, LinearWavefunction const& Psi2)
      : Psi1_(Psi1), QShift_(QShift), Op_(Op), Psi2_(Psi2), Scale_(1.0)
   {
   }

   InjectLeftQShift(LinearWavefunction const& Psi1, QuantumNumber const& QShift, GenericMPO const& Op, LinearWavefunction const& Psi2, std::complex<double> Scale)
      : Psi1_(Psi1), QShift_(QShift), Op_(Op), Psi2_(Psi2), Scale_(Scale)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return delta_shift(inject_left(m, Psi1_, Op_, Psi2_), QShift_) * Scale_;
   }

   LinearWavefunction const& Psi1_;
   QuantumNumber const& QShift_;
   GenericMPO const& Op_;
   LinearWavefunction const& Psi2_;
   std::complex<double> Scale_;
};

// Generalization of the operator contraction
//
// +-Psi1*-     Psi1*-
// |  |          |
// E--Op--- ...  Op---
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
// --Psi1--... Psi1--+
//    |         |    |
//    Op*-...  Op*   m
//    |         |    |
// --Psi2*-... Psi2*-+

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
      : Psi1_(Psi), QShift_(QShift), Op_(Op), Psi2_(Psi), Scale_(1.0)
   {
   }

   InjectRightQShift(LinearWavefunction const& Psi1, QuantumNumber const& QShift, GenericMPO const& Op, LinearWavefunction const& Psi2)
      : Psi1_(Psi1), QShift_(QShift), Op_(Op), Psi2_(Psi2), Scale_(1.0)
   {
   }

   InjectRightQShift(LinearWavefunction const& Psi1, QuantumNumber const& QShift, GenericMPO const& Op, LinearWavefunction const& Psi2, std::complex<double> Scale)
      : Psi1_(Psi1), QShift_(QShift), Op_(Op), Psi2_(Psi2), Scale_(Scale)
   {
   }

   MatrixOperator operator()(MatrixOperator const& m) const
   {
      return inject_right(delta_shift(m, adjoint(QShift_)), Psi1_, Op_, Psi2_);
   }

   LinearWavefunction const& Psi1_;
   QuantumNumber const& QShift_;
   GenericMPO const& Op_;
   LinearWavefunction const& Psi2_;
   std::complex<double> Scale_;
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

//
// LeftMultiplyOperator
// Functor to contract from the left hand side two wavefunctions with a string operator,
// allowing the wavefunctions and operator to have different unit cell sizes.
//

struct LeftMultiplyOperator
{
   typedef StateComponent argument_type;
   typedef StateComponent result_type;

   LeftMultiplyOperator(LinearWavefunction const& L1_, QuantumNumber const& QShift1_,
                        ProductMPO const& StringOp_,
                        LinearWavefunction const& L2_, QuantumNumber const& QShift2_,
                        int Length_ = 0, int Verbose_ = 0)
      : Psi1(L1_), QShift1(QShift1_), StringOp(StringOp_) , Psi2(L2_), QShift2(QShift2_), Length(Length_),
        Verbose(Verbose_)
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
      if (Verbose > 0)
      {
         std::cerr << "contracting operator, total length is " << Length << " sites.\n";
      }
      while (I1 != Psi1.end() || I2 != Psi2.end() || OpIter != StringOp.end())
      {
         if (Verbose > 1)
         {
            std::cerr << "Site " << n << ", E-matrix dimension " << R.size()
                      << "x" << R.Basis1().total_dimension()
                      << "x" << R.Basis2().total_dimension()
                      << '\n';
         }
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
   int Verbose;
};


struct RightMultiplyOperator
{
   typedef StateComponent argument_type;
   typedef StateComponent result_type;

   RightMultiplyOperator(LinearWavefunction const& L1_, QuantumNumber const& QShift1_,
                         ProductMPO const& StringOp_,
                         LinearWavefunction const& L2_, QuantumNumber const& QShift2_,
                         int Length_ = 0, int Verbose_ = 0)
      : Psi1(L1_), QShift1(QShift1_), StringOp(StringOp_) , Psi2(L2_), QShift2(QShift2_), Length(Length_),
        Verbose(Verbose_)
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
      LinearWavefunction::const_iterator I1 = Psi1.end();
      LinearWavefunction::const_iterator I2 = Psi2.end();
      ProductMPO::const_iterator OpIter = StringOp.end();
      int n = Length;
      QuantumNumber q1 = QuantumNumber(QShift1.GetSymmetryList());
      QuantumNumber q2 = QuantumNumber(QShift2.GetSymmetryList());
      if (Verbose > 0)
      {
         std::cerr << "contracting operator, total length is " << Length << " sites.\n";
      }
      while (I1 != Psi1.begin() || I2 != Psi2.begin() || OpIter != StringOp.begin())
      {
         --n;

         if (I1 == Psi1.begin())
         {
            I1 = Psi1.end();
            q1 = delta_shift(q1, adjoint(QShift1));
         }
         --I1;

         if (I2 == Psi2.begin())
         {
            I2 = Psi2.end();
            q2 = delta_shift(q2, adjoint(QShift2));
         }
         --I2;

         if (OpIter == StringOp.begin())
            OpIter = StringOp.end();
         --OpIter;

         if (Verbose > 1)
         {
            std::cerr << "Site " << n << ", E-matrix dimension " << R.size()
                      << "x" << R.Basis1().total_dimension()
                      << "x" << R.Basis2().total_dimension()
                      << '\n';
         }
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
         R = contract_from_right(herm(*OpIter), delta_shift(*I1, adjoint(q1)), R,
                                 herm(delta_shift(*I2, adjoint(q2))));
      }
      q1 = delta_shift(q1, adjoint(QShift1));
      q2 = delta_shift(q2, adjoint(QShift2));
      CHECK_EQUAL(q1, q2);
      //      CHECK_EQUAL(n, 0);
      return delta_shift(R, q1);
   }

   LinearWavefunction const& Psi1;
   QuantumNumber QShift1;
   ProductMPO StringOp;
   LinearWavefunction const& Psi2;
   QuantumNumber QShift2;
   int Length;
   int Verbose;
};


// version of inject_left operator contraction, but
// with a mask (conceptually a vector<bool>, but implemented as a vector<int>).
// Only calculates components of the MPO where the mask bit is set.
// The length of the mask is Op.size()+1

StateComponent
inject_left_mask(StateComponent const& In,
                 LinearWavefunction const& Psi1,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int> > const& Mask);

StateComponent
inject_right_mask(StateComponent const& In,
                 LinearWavefunction const& Psi1,
                 GenericMPO const& Op,
                 LinearWavefunction const& Psi2,
                 std::vector<std::vector<int> > const& Mask);


// Functor to evaluate (1-T)(x) where T is the generalized transfer matrix.
// x is defined in the Basis1.
//      -Psi*-
//        |
//  T =  Op
//        |
//      -Psi--
struct OneMinusTransferLeft
{
   OneMinusTransferLeft(BasicFiniteMPO const& Op, LinearWavefunction const& Psi, QuantumNumber const& QShift)
      : Op_(Op), Psi_(Psi), QShift_(QShift) {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      return x-delta_shift(inject_left(x, Op_, Psi_), QShift_);
   }

   BasicFiniteMPO const& Op_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
};

// Functor to evaluate (1-T)(x) where T is the generalized transfer matrix.
// x is defined in the Basis1.  Optionally, we also optionally orthogonalize against a
// Unit matrix, which is an eigenvalue 1 of T
//      -Psi*-
//        |
//  T =  Op
//        |
//      -Psi--
struct OneMinusTransferLeft_Ortho
{
   OneMinusTransferLeft_Ortho(LinearWavefunction const& Psi1, QuantumNumber const& QShift,
                              BasicFiniteMPO const& Op, LinearWavefunction const& Psi2,
                              MatrixOperator const& LeftUnit,
                              MatrixOperator const& RightUnit, bool Orthogonalize)
      : Psi1_(Psi1), QShift_(QShift), Op_(Op), Psi2_(Psi2),
        LeftUnit_(LeftUnit),
        RightUnit_(RightUnit), Scale_(1.0), Orthogonalize_(Orthogonalize) { }

   OneMinusTransferLeft_Ortho(LinearWavefunction const& Psi1, QuantumNumber const& QShift,
                              BasicFiniteMPO const& Op, LinearWavefunction const& Psi2,
                              MatrixOperator const& LeftUnit,
                              MatrixOperator const& RightUnit, std::complex<double> Scale, bool Orthogonalize)
      : Psi1_(Psi1), QShift_(QShift), Op_(Op), Psi2_(Psi2),
        LeftUnit_(LeftUnit),
        RightUnit_(RightUnit), Scale_(Scale), Orthogonalize_(Orthogonalize) { }

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-delta_shift(inject_left(x, Psi1_, Op_, Psi2_), QShift_)*Scale_;
      if (Orthogonalize_ && r.TransformsAs() == RightUnit_.TransformsAs())
         {
            DEBUG_TRACE(inner_prod(r, RightUnit_))("this should be small");
            DEBUG_TRACE(inner_prod(LeftUnit_, r));
            r -= std::conj(inner_prod(r, RightUnit_)) * LeftUnit_; // orthogonalize to the identity
            DEBUG_TRACE(inner_prod(r, RightUnit_))("this should be zero");
         }
      return r;
   }

   LinearWavefunction const& Psi1_;
   QuantumNumber const& QShift_;
   BasicFiniteMPO const& Op_;
   LinearWavefunction const& Psi2_;
   MatrixOperator const& LeftUnit_;
   MatrixOperator const& RightUnit_;
   std::complex<double> Scale_;
   bool Orthogonalize_;
};

// Functor to evaluate (1-T)(x) where T is the generalized transfer matrix.
// x is defined in the Basis2.
//      -Psi-
//        |
//  T =  Op*
//        |
//      -Psi*--
struct OneMinusTransferRight
{
   OneMinusTransferRight(BasicFiniteMPO const& Op, LinearWavefunction const& Psi, QuantumNumber const& QShift)
      : Op_(Op), Psi_(Psi), QShift_(QShift) {}

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      return x-delta_shift(inject_right(x, Op_, Psi_), adjoint(QShift_));
   }

   BasicFiniteMPO const& Op_;
   LinearWavefunction const& Psi_;
   QuantumNumber const& QShift_;
};

struct OneMinusTransferRight_Ortho
{
   OneMinusTransferRight_Ortho(LinearWavefunction const& Psi1, QuantumNumber const& QShift,
                              BasicFiniteMPO const& Op, LinearWavefunction const& Psi2,
                              MatrixOperator const& LeftUnit,
                              MatrixOperator const& RightUnit, bool Orthogonalize)
      : Psi1_(Psi1), QShift_(QShift), Op_(Op), Psi2_(Psi2),
        LeftUnit_(LeftUnit),
        RightUnit_(RightUnit), Scale_(1.0), Orthogonalize_(Orthogonalize) { }

   OneMinusTransferRight_Ortho(LinearWavefunction const& Psi1, QuantumNumber const& QShift,
                              BasicFiniteMPO const& Op, LinearWavefunction const& Psi2,
                              MatrixOperator const& LeftUnit,
                              MatrixOperator const& RightUnit, std::complex<double> Scale, bool Orthogonalize)
      : Psi1_(Psi1), QShift_(QShift), Op_(Op), Psi2_(Psi2),
        LeftUnit_(LeftUnit),
        RightUnit_(RightUnit), Scale_(Scale), Orthogonalize_(Orthogonalize) { }

   MatrixOperator operator()(MatrixOperator const& x) const
   {
      MatrixOperator r = x-delta_shift(inject_right(x, Psi1_, Op_, Psi2_), adjoint(QShift_))*Scale_;
      if (Orthogonalize_ && r.TransformsAs() == LeftUnit_.TransformsAs())
         {
            r -= std::conj(inner_prod(r, LeftUnit_)) * RightUnit_; // orthogonalize to the identity
         }
      return r;
   }

   LinearWavefunction const& Psi1_;
   QuantumNumber const& QShift_;
   BasicFiniteMPO const& Op_;
   LinearWavefunction const& Psi2_;
   MatrixOperator const& LeftUnit_;
   MatrixOperator const& RightUnit_;
   std::complex<double> Scale_;
   bool Orthogonalize_;
};

#endif
