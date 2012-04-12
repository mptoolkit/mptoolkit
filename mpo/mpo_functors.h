// -*- C++ -*- $Id$
//
// Collection of functors for multiplication of MPO and MPS
//

#if !defined(MPO_FUNCTORS_H_JK98390U8UJ895TU8PGJ5P89TJP)
#define MPO_FUNCTORS_H_JK98390U8UJ895TU8PGJ5P89TJP

#include "mps/state_component.h"
#include "mpo/operator_component.h"
#include "mp-algorithms/gmres.h"

struct OperatorProd_AxBH
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   OperatorProd_AxBH(StateComponent const& Left_,
                     StateComponent const& Right_)
      : Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(Left, Psi, herm(Right));
   }
   
   StateComponent Left, Right;
};

struct OperatorProd_MAxBH
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   OperatorProd_MAxBH(SimpleOperator const& M_, 
                      StateComponent const& Left_,
                      StateComponent const& Right_)
      : M(M_), Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(M, Left, Psi, herm(Right));
   }

   SimpleOperator M;
   StateComponent Left, Right;
};

struct OperatorProd_MAHxB
{
   typedef MatrixOperator result_type;
   typedef MatrixOperator value_type;
   typedef MatrixOperator argument_type;

   OperatorProd_MAHxB(SimpleOperator const& M_, 
                      StateComponent const& Left_,
                      StateComponent const& Right_)
      : M(M_), Left(Left_), Right(Right_) {}

   MatrixOperator operator()(MatrixOperator const& Psi) const
   {
      return operator_prod(herm(M), herm(Left), Psi, Right);
   }

   SimpleOperator M;
   StateComponent Left, Right;
};

struct OperatorProdInner_AxBH
{
   typedef StateComponent result_type;
   typedef StateComponent value_type;
   typedef StateComponent argument_type;

   OperatorProdInner_AxBH(OperatorComponent const& Op_,
                          StateComponent const& Left_,
                          StateComponent const& Right_)
      : Op(Op_), Left(Left_), Right(Right_) {}

   StateComponent operator()(StateComponent const& Psi) const
   {
      return operator_prod_inner(Op, Left, Psi, herm(Right));
   }
   
   OperatorComponent Op;
   StateComponent Left, Right;
};

struct OperatorProdInner_AxBH_Normal
{
   typedef StateComponent result_type;
   typedef StateComponent value_type;
   typedef StateComponent argument_type;

   OperatorProdInner_AxBH_Normal(OperatorComponent const& Op_,
                                 StateComponent const& Left_,
                                 StateComponent const& Right_, MatrixOperator const& Normalizer_)
      : Op(Op_), Left(Left_), Right(Right_), Normalizer(Normalizer_) {}

   StateComponent operator()(StateComponent const& Psi) const
   {
      return triple_prod(Normalizer, operator_prod_inner(Op, Left, Psi, herm(Right)), herm(Normalizer));
   }
   
   OperatorComponent Op;
   StateComponent Left, Right;
   MatrixOperator Normalizer;
};

struct TripleProd_State_AHxB
{
   typedef StateComponent result_type;
   typedef StateComponent value_type;
   typedef StateComponent argument_type;

   TripleProd_State_AHxB(MatrixOperator const& Left_, MatrixOperator const& Right_)
      : Left(Left_), Right(Right_) {}

   StateComponent operator()(StateComponent const& Psi) const
   {
      return triple_prod(herm(Left), Psi, Right);
   }

   MatrixOperator Left, Right;
};


struct OperatorProdInner_AxBH_Gen
{
   typedef StateComponent result_type;
   typedef StateComponent value_type;
   typedef StateComponent argument_type;

   OperatorProdInner_AxBH_Gen(OperatorComponent const& Op_,
                              StateComponent const& Left_,
                              StateComponent const& Right_, MatrixOperator const& Normalizer_)
      : Op(Op_), Left(Left_), Right(Right_), Normalizer(Normalizer_) {}

   StateComponent operator()(StateComponent const& Psi) const
   {
      StateComponent R = operator_prod_inner(Op, Left, Psi, herm(Right));
      int m = 30;
      int max_iter = 10000;
      double tol = 1e-7;
      MatrixOperator Inv = InvertDiagonal(Normalizer, 1E-8);
      StateComponent Result = triple_prod(herm(Inv), R, Inv);
      GmRes(Result, TripleProd_State_AHxB(Normalizer, Normalizer), R, m, max_iter, tol,  LinearAlgebra::Identity<StateComponent>());
      return Result;
   }
   
   OperatorComponent Op;
   StateComponent Left, Right;
   MatrixOperator Normalizer;
};

template <typename T>
struct InnerProdNonOrtho
{
   typedef std::complex<double> result_type;
   typedef T first_argument_type;
   typedef T second_argument_type;

   InnerProdNonOrtho(MatrixOperator const& Ident_) : Ident(Ident_) {}

   std::complex<double> 
   operator()(T const& x, T const& y) const
   {
      return inner_prod(x, triple_prod(Ident, y, herm(Ident)));
   }

   MatrixOperator Ident;
};

template <typename T>
struct NormFrobNonOrtho
{
   typedef double result_type;
   typedef T first_argument_type;
   typedef T second_argument_type;

   NormFrobNonOrtho(MatrixOperator const& Ident_) : Ident(Ident_) {}

   double operator()(T const& x) const
   {
      double r = real(inner_prod(x, triple_prod(Ident, x, herm(Ident))));
      if (r < 0) r = 0;
      return std::sqrt(r);
   }

   MatrixOperator Ident;
};

#endif
