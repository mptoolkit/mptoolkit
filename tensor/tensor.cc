// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// tensor/tensor.cc
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "linearalgebra/matrix.h"

namespace Tensor
{

//
// IrredTensor
//

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>::IrredTensor(basis1_type const& Basis, QuantumNumber const& Trans)
   : Basis1_(Basis), Basis2_(Basis), Trans_(Trans), Data_(Basis.size(), Basis.size())
{
   DEBUG_CHECK_EQUAL(Basis.GetSymmetryList(), Trans.GetSymmetryList());
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>::IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2, 
				QuantumNumber const& Trans)
   : Basis1_(Basis1), Basis2_(Basis2), Trans_(Trans), Data_(Basis1.size(), Basis2.size())
{
   DEBUG_CHECK_EQUAL(Basis1.GetSymmetryList(), Basis2.GetSymmetryList());
   DEBUG_CHECK_EQUAL(Basis1.GetSymmetryList(), Trans.GetSymmetryList());
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>::IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2)
   : Basis1_(Basis1), Basis2_(Basis2), Trans_(Basis1.GetSymmetryList()), 
     Data_(Basis1.size(), Basis2.size())
{
   DEBUG_CHECK_EQUAL(Basis1.GetSymmetryList(), Basis2.GetSymmetryList());
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>::IrredTensor(basis1_type const& Basis1, basis2_type const& Basis2, 
			       QuantumNumber const& Trans, MatrixType const& Data)
   : Basis1_(Basis1), Basis2_(Basis2), Trans_(Trans), Data_(Data)
{
   DEBUG_CHECK_EQUAL(Basis1.GetSymmetryList(), Basis2.GetSymmetryList());
   DEBUG_CHECK_EQUAL(Basis1.GetSymmetryList(), Trans.GetSymmetryList());
   using LinearAlgebra::size1;
   using LinearAlgebra::size2;
   DEBUG_PRECONDITION_EQUAL(size1(Data), Basis1.size());
   DEBUG_PRECONDITION_EQUAL(size2(Data), Basis2.size());
}

template <typename T, typename B1, typename B2, typename S>
inline
typename LinearAlgebra::MatrixBracket<typename IrredTensor<T, B1, B2, S>::MatrixType&, 
				      size_type, size_type>::result_type
IrredTensor<T, B1, B2, S>::operator()(size_type i, size_type j) 
{ 
   DEBUG_CHECK(i < this->Basis1().size())(i)(this->Basis1());
   DEBUG_CHECK(j < this->Basis2().size())(j)(this->Basis2());
   DEBUG_PRECONDITION(is_transform_target(this->qn2(j), this->TransformsAs(), this->qn1(i)))
      (i)(j)(this->TransformsAs())(this->qn1(i))(this->qn2(j));
   return Data_(i,j); 
}

template <typename T, typename B1, typename B2, typename S>
inline
typename LinearAlgebra::MatrixBracket<typename IrredTensor<T, B1, B2, S>::MatrixType, 
				      size_type, size_type>::result_type
IrredTensor<T, B1, B2, S>::operator()(size_type i, size_type j) const
{
   DEBUG_CHECK(i < this->Basis1().size())(i)(this->Basis1());
   DEBUG_CHECK(j < this->Basis2().size())(j)(this->Basis2());
   DEBUG_PRECONDITION(is_transform_target(this->qn2(j), this->TransformsAs(), this->qn1(i)))
      (this->qn2(j))(this->TransformsAs())(this->qn1(i));
   return Data_(i,j); 
}



template <typename T, typename B1, typename B2, typename S>
void implement_check_structure(IrredTensor<T, B1, B2, S> const& x)
{
   for (typename LinearAlgebra::const_iterator<IrredTensor<T, B1, B2, S> >::type 
           I = LinearAlgebra::iterate(x); I; ++I)
   {
      for (typename LinearAlgebra::const_inner_iterator<IrredTensor<T, B1, B2, S> >::type 
              J = LinearAlgebra::iterate(I); J; ++J)
      {
         CHECK(is_transform_target(x.Basis2()[J.index2()], x.TransformsAs(),
				   x.Basis1()[J.index1()]))(x.Basis2()[J.index2()])
            (x.Basis1()[J.index1()])(x.TransformsAs())(*J);
      }
   }
}

template <typename T, typename B1, typename B2, typename S>
void implement_check_structure(IrredTensor<LinearAlgebra::Matrix<T>, B1, B2, S> const& x)
{
   for (auto I = LinearAlgebra::iterate(x); I; ++I)
   {
      for (auto J = LinearAlgebra::iterate(I); J; ++J)
      {
         CHECK(is_transform_target(x.Basis2()[J.index2()], x.TransformsAs(),
				   x.Basis1()[J.index1()]))(x.Basis2()[J.index2()])
            (x.Basis1()[J.index1()])(x.TransformsAs())(*J);

	 CHECK_EQUAL(int(size1(*J)), x.Basis1().dim(J.index1()));
	 CHECK_EQUAL(int(size2(*J)), x.Basis2().dim(J.index2()));
      }
   }
}


template <typename T, typename B1, typename B2, typename S>
void IrredTensor<T, B1, B2, S>::check_structure() const
{
   implement_check_structure(*this);
}

template <typename T, typename B1, typename B2, typename S>
void
IrredTensor<T, B1, B2, S>::CoerceSymmetryList(SymmetryList const& sl)
{
   using Tensor::CoerceSymmetryList;
   Basis1_.CoerceSymmetryList(sl);
   Basis2_.CoerceSymmetryList(sl);
   Trans_.CoerceSymmetryList(sl);
   Data_ = CoerceSymmetryList(Data_, sl);
   this->debug_check_structure();
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>
CoerceSymmetryList(IrredTensor<T, B1, B2, S> const& t, SymmetryList const& sl)
{
   IrredTensor<T, B1, B2, S> Result;
   Result.Basis1_ = CoerceSymmetryList(t.Basis1_, sl);
   Result.Basis2_ = CoerceSymmetryList(t.Basis2_, sl);
   Result.Trans_ = CoerceSymmetryList(t.Trans_, sl);
   Result.Data_ = CoerceSymmetryList(t.Data_, sl);
   Result.debug_check_structure();
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
void
CoerceSymmetryListInPlace(IrredTensor<T, B1, B2, S>& t, SymmetryList const& sl)
{
   t.CoerceSymmetryList(sl);
}

template <typename T, typename B1, typename B2, typename S>
inline
void IrredTensor<T, B1, B2, S>::debug_check_structure() const
{
#if !defined(NDEBUG)
   this->check_structure();
#endif
}

template <typename T, typename B1, typename B2, typename S>
PStream::opstream& operator<<(PStream::opstream& out, IrredTensor<T, B1, B2, S> const& Op)
{
   return out << Op.Basis1_ << Op.Basis2_ << Op.Trans_ << Op.Data_;
}

template <typename T, typename B1, typename B2, typename S>
PStream::ipstream& operator>>(PStream::ipstream& in, IrredTensor<T, B1, B2, S>& Op)
{
   return in >> Op.Basis1_ >> Op.Basis2_ >> Op.Trans_ >> Op.Data_;
}

template <typename T, typename B1, typename B2, typename S>
std::ostream& operator<<(std::ostream& out, IrredTensor<T, B1, B2, S> const& Op)
{
   return out << "Operator transforms with symmetry " << Op.GetSymmetryList()
              << " as " << Op.TransformsAs() << ":\n" << Op.data();
}

template <typename T, typename B1, typename B2, typename S>
std::string show_projections(IrredTensor<T, B1, B2, S> const& Op)
{
   std::ostringstream out;
   out << "Operator transforms with symmetry " << Op.GetSymmetryList() 
       << " as " << Op.TransformsAs() << ":\n";

   std::vector<QuantumNumbers::Projection> Projections;
   enumerate_projections(Op.TransformsAs(), std::back_inserter(Projections));
   for (size_type km = 0; km < Projections.size(); ++km)
   {
      out << "Projection " << std::setw(10) << Projections[km] << " :\n";

      typename IrredTensor<T, B1, B2, S>::const_iterator I = iterate(Op);
      while (I)
      {
         typename IrredTensor<T, B1, B2, S>::const_inner_iterator J = iterate(I);
         while (J)
         {
            QuantumNumber qi = Op.qn1(J.index1());
	    QuantumNumber qj = Op.qn2(J.index2());

	    std::vector<QuantumNumbers::Projection> mi, mj;
	    enumerate_projections(qi, std::back_inserter(mi));
	    enumerate_projections(qj, std::back_inserter(mj));

	    for (size_type ii = 0; ii < mi.size(); ++ii)
	    {
	       for (size_type jj = 0; jj < mj.size(); ++jj)
	       {
                  T elem = *J * clebsch_gordan(qj,    Op.TransformsAs(),  qi,
                                               mj[jj], Projections[km], mi[ii]);

		  if (LinearAlgebra::norm_frob(elem) > 1E-10)
		  {
		     out << "   " << elem
			 << " |" << J.index1() << ": " << qi
			 << ", " << mi[ii]
			 << "> <"
			 << J.index2() << ": " << qj
			 << ", " << mj[jj]
			 << "|\n";
		  }
	       }
	    }
            ++J;
	 }
         ++I;
      }
   }
   return out.str();
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>&
IrredTensor<T, B1, B2, S>::operator+=(IrredTensor const& Op)
{
   // quick return if the right hand side is zero
   if (Op.is_null()) return *this;

   if (this->is_null()) 
   {
      *this = Op;
      return *this;
   }

   PRECONDITION(Basis1_ == Op.Basis1_)(Basis1_)(Op.Basis1_);
   PRECONDITION(Basis2_ == Op.Basis2_)(Basis2_)(Op.Basis2_);
   PRECONDITION(Trans_ == Op.Trans_)(Trans_)(Op.Trans_);
   Data_ += Op.Data_;

   this->debug_check_structure();
   return *this;
}

template <typename T, typename B1, typename B2, typename S>
IrredTensor<T, B1, B2, S>&
IrredTensor<T, B1, B2, S>::operator-=(IrredTensor const& Op)
{
   // quick return if the right hand side is zero
   if (Op.is_null()) return *this;

   if (this->is_null()) 
   {
      *this = -Op;
      return *this;
   }

   PRECONDITION(Basis1_ == Op.Basis1_)(Basis1_)(Op.Basis1_);
   PRECONDITION(Basis2_ == Op.Basis2_)(Basis2_)(Op.Basis2_);
   PRECONDITION(Trans_ == Op.Trans_)(Trans_)(Op.Trans_);
   Data_ -= Op.Data_;

   this->debug_check_structure();
   return *this;
}

//
// triple_prod
//

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, 
          typename T3, typename B4, typename S3>
IrredTensor<typename Prod3Type<T1, T2, T3>::type, B2, B4>
triple_prod(HermitianProxy<IrredTensor<T1, B1, B2, S1> > const& x, 
	    IrredTensor<T2, B1, B3, S2> const& E,
	    IrredTensor<T3, B3, B4, S3> const& y,
	    QuantumNumber qxy,
	    QuantumNumber qEp)
{
   if (x.base().is_null() || y.is_null() || E.is_null()) 
      return IrredTensor<typename Prod3Type<T1, T2, T3>::type, B2, B4>();

   DEBUG_PRECONDITION_EQUAL(y.GetSymmetryList(), x.base().GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(y.GetSymmetryList(), qxy.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(y.GetSymmetryList(), E.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(y.GetSymmetryList(), qEp.GetSymmetryList());
   QuantumNumber Transx = x.base().TransformsAs();
   QuantumNumber Transy = y.TransformsAs();
   QuantumNumber TransE = E.TransformsAs();
   PRECONDITION(is_transform_target(Transy, qxy, Transx))(Transy)(qxy)(Transx);
   PRECONDITION(is_transform_target(qEp, qxy, TransE))(qEp)(qxy)(TransE);
   DEBUG_PRECONDITION_EQUAL(x.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), y.Basis1());

   IrredTensor<T1, B1, B2, S1> const& xb = x.base();  // shortcut

   IrredTensor<typename Prod3Type<T1, T2, T3>::type, B2, B4> 
      Result(x.base().Basis2(), y.Basis2(), qEp);

   typedef typename const_iterator< IrredTensor<T1, B1, B2, S1> >::type       x_iter;
   typedef typename const_inner_iterator< IrredTensor<T1, B1, B2, S1> >::type x_inner;

   typedef typename const_iterator< IrredTensor<T2, B1, B3, S2> >::type       E_iter;
   typedef typename const_inner_iterator< IrredTensor<T2, B1, B3, S2> >::type E_inner;

   typedef typename const_iterator< IrredTensor<T3, B3, B4, S3> >::type       y_iter;
   typedef typename const_inner_iterator< IrredTensor<T3, B3, B4, S3> >::type y_inner;

   // iterate over elements in E first
   for (E_iter Ei = iterate(E.data()); Ei; ++Ei)
   {
      x_iter xi = iterate(xb.data());
      xi += Ei.index();
      QuantumNumber qEi = E.Basis1()[Ei.index()];
      double DegEi = degree(qEi);
      for (E_inner Ein = iterate(Ei); Ein; ++Ein)
      {
         QuantumNumber qEj = E.Basis2()[Ein.index2()];

         y_iter yi = iterate(y.data());
         yi += Ein.index2();
 
         for (x_inner xin = iterate(xi); xin; ++xin)
         {
            QuantumNumber qxip = xb.Basis2()[xin.index2()];
            double DegXi = degree(qxip);
            for (y_inner yin = iterate(yi); yin; ++yin)
            {
               QuantumNumber qyj = y.Basis2()[yin.index2()];
               if (!(is_transform_target(qyj, qEp, qxip))) continue;

               double Coefficient = (DegEi / DegXi) *
                  tensor_coefficient(qyj, Transy, qEj,
       				     qEp, qxy, TransE,
                                     qxip, Transx, qEi);

	       if (fabs(Coefficient) > 1E-10)
	       {
		  add_element(Result.data(), xin.index2(), yin.index2(), Coefficient *
                              herm(*xin) * (*Ein) * (*yin));
	       }
	    }
	 }
      }
   }

   Result.debug_check_structure();
   return Result;
}

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, 
          typename T3, typename B4, typename S3>
IrredTensor<typename Prod3Type<T1, T2, T3>::type, B1, B4>
triple_prod(IrredTensor<T1, B1, B2, S1> const& x, 
	    IrredTensor<T2, B2, B3, S2> const& E,
	    HermitianProxy<IrredTensor<T3, B4, B3, S3> > const& y,
	    QuantumNumber qxy,
	    QuantumNumber qEp)
{
   if (x.is_null() || E.is_null() || y.base().is_null())
      return IrredTensor<typename Prod3Type<T1, T2, T3>::type, B1, B4>();

   //DEBUG_PRECONDITION_EQUAL(x.GetSymmetryList(), y.base().GetSymmetryList());
   //DEBUG_PRECONDITION_EQUAL(x.GetSymmetryList(), qxy.GetSymmetryList());
   //DEBUG_PRECONDITION_EQUAL(x.GetSymmetryList(), E.GetSymmetryList());
   //DEBUG_PRECONDITION_EQUAL(x.GetSymmetryList(), qEp.GetSymmetryList());
   PRECONDITION_EQUAL(x.GetSymmetryList(), y.base().GetSymmetryList());
   PRECONDITION_EQUAL(x.GetSymmetryList(), qxy.GetSymmetryList());
   PRECONDITION_EQUAL(x.GetSymmetryList(), E.GetSymmetryList());
   PRECONDITION_EQUAL(x.GetSymmetryList(), qEp.GetSymmetryList());
   QuantumNumber Transx = x.TransformsAs();
   QuantumNumber Transy = y.base().TransformsAs();
   QuantumNumber TransE = E.TransformsAs();
   PRECONDITION(is_transform_target(Transy, qxy, Transx))(Transy)(qxy)(Transx);
   PRECONDITION(is_transform_target(TransE, qxy, qEp))(TransE)(qxy)(qEp);
   //DEBUG_PRECONDITION_EQUAL(x.Basis2(), E.Basis1());
   //DEBUG_PRECONDITION_EQUAL(E.Basis2(), y.base().Basis2());
   PRECONDITION_EQUAL(x.Basis2(), E.Basis1());
   PRECONDITION_EQUAL(E.Basis2(), y.base().Basis2());

   IrredTensor<typename Prod3Type<T1, T2, T3>::type, B1, B4>
      Result(x.Basis1(), y.base().Basis1(), qEp);

   IrredTensor<T3, B4, B3, S3> const& yb = y.base();  // shortcut

   typedef typename const_iterator< IrredTensor<T1, B1, B2, S1> >::type       x_iter;
   typedef typename const_inner_iterator< IrredTensor<T1, B1, B2, S1> >::type x_inner;

   typedef typename const_iterator< IrredTensor<T2, B2, B3, S2> >::type       E_iter;
   typedef typename const_inner_iterator< IrredTensor<T2, B2, B3, S2> >::type E_inner;

   typedef typename const_iterator< IrredTensor<T3, B4, B3, S3> >::type       y_iter;
   typedef typename const_inner_iterator< IrredTensor<T3, B4, B3, S3> >::type y_inner;

   for (x_iter xi = iterate(x); xi; ++xi)
   {
      QuantumNumber qxp = x.Basis1()[xi.index()];
      for (y_iter yi = iterate(yb); yi; ++yi)
      {
         QuantumNumber qyp = yb.Basis1()[yi.index()];
         if (!is_transform_target(qyp, qEp, qxp)) continue;

         for (x_inner xin = iterate(xi); xin; ++xin)
         {
            for (y_inner yin = iterate(yi); yin; ++yin)
            {
               E_inner Ein = iterate_at(E.data(), xin.index2(), yin.index2());
               if (Ein)
               {
                  /*
                  double Coefficient = tensor_coefficient(x.Basis2()[xin.index2()], Transx, qxp,
							  TransE, qxy, qEp, 
							  yb.Basis2()[yin.index2()], Transy, qyp);
                  */
                  double Coefficient = tensor_coefficient(yb.Basis2()[yin.index2()], Transy, qyp,
							  TransE, qxy, qEp, 
                                                          x.Basis2()[xin.index2()], Transx, qxp);

		  if (fabs(Coefficient) > 1E-10)
		  {
		     add_element(Result.data(), xin.index1(), yin.index1(),
                                 Coefficient * (*xin) * E(xin.index2(), yin.index2()) * herm(*yin));
                  }
	       }
	    }
	 }
      }
   }

   Result.debug_check_structure();
   return Result;
}

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, 
          typename T3, typename B4, typename S3>
inline
IrredTensor<typename Prod3Type<T1, T2, T3>::type, B2, B4>
triple_prod(HermitianProxy<IrredTensor<T1, B1, B2, S1> > const& x, 
	    IrredTensor<T2, B1, B3, S2> const& E,
	    IrredTensor<T3, B3, B4, S3> const& y)
{
   return triple_prod(x, E, y, QuantumNumber(y.GetSymmetryList()), E.TransformsAs());
}

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, 
          typename T3, typename B4, typename S3>
inline
IrredTensor<typename Prod3Type<T1, T2, T3>::type, B1, B4>
triple_prod(IrredTensor<T1, B1, B2, S1> const& x, 
	    IrredTensor<T2, B2, B3, S2> const& E,
	    HermitianProxy<IrredTensor<T3, B4, B3, S3> > const& y)
{
   return triple_prod(x, E, y, QuantumNumber(x.GetSymmetryList()), E.TransformsAs());
}

// add_triple_prod

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, 
          typename T3, typename B4, typename S3, typename T4, typename S4>
void
add_triple_prod(Tensor::IrredTensor<T4, B1, B4, S4>& Result,
                std::complex<double> Factor,
                HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> > const& x, 
                Tensor::IrredTensor<T2, B1, B3, S2> const& E,
                Tensor::IrredTensor<T3, B3, B4, S3> const& y,
                QuantumNumber qxy,
                QuantumNumber qEp)
{
   if (x.base().is_null() || y.is_null() || E.is_null()) 
      return;

   DEBUG_PRECONDITION_EQUAL(y.GetSymmetryList(), x.base().GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(y.GetSymmetryList(), qxy.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(y.GetSymmetryList(), E.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(y.GetSymmetryList(), qEp.GetSymmetryList());
   QuantumNumber Transx = x.base().TransformsAs();
   QuantumNumber Transy = y.TransformsAs();
   QuantumNumber TransE = E.TransformsAs();
   PRECONDITION(is_transform_target(Transy, qxy, Transx))(Transy)(qxy)(Transx);
   PRECONDITION(is_transform_target(qEp, qxy, TransE))(qEp)(qxy)(TransE);
   DEBUG_PRECONDITION_EQUAL(x.base().Basis1(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), y.Basis1());
   DEBUG_PRECONDITION_EQUAL(Result.Basis1(), x.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(Result.Basis2(), y.Basis2());
   DEBUG_PRECONDITION_EQUAL(Result.TransformsAs(), qEp);

   IrredTensor<T1, B1, B2, S1> const& xb = x.base();  // shortcut

   typedef typename const_iterator< IrredTensor<T1, B1, B2, S1> >::type       x_iter;
   typedef typename const_inner_iterator< IrredTensor<T1, B1, B2, S1> >::type x_inner;

   typedef typename const_iterator< IrredTensor<T2, B1, B3, S2> >::type       E_iter;
   typedef typename const_inner_iterator< IrredTensor<T2, B1, B3, S2> >::type E_inner;

   typedef typename const_iterator< IrredTensor<T3, B3, B4, S3> >::type       y_iter;
   typedef typename const_inner_iterator< IrredTensor<T3, B3, B4, S3> >::type y_inner;

   // iterate over elements in E first
   for (E_iter Ei = iterate(E.data()); Ei; ++Ei)
   {
      x_iter xi = iterate(xb.data());
      xi += Ei.index();
      QuantumNumber qEi = E.Basis1()[Ei.index()];
      double DegEi = degree(qEi);
      for (E_inner Ein = iterate(Ei); Ein; ++Ein)
      {
         QuantumNumber qEj = E.Basis2()[Ein.index2()];

         y_iter yi = iterate(y.data());
         yi += Ein.index2();
 
         for (x_inner xin = iterate(xi); xin; ++xin)
         {
            QuantumNumber qxip = xb.Basis2()[xin.index2()];
            double DegXi = degree(qxip);
            for (y_inner yin = iterate(yi); yin; ++yin)
            {
               QuantumNumber qyj = y.Basis2()[yin.index2()];
               if (!(is_transform_target(qyj, qEp, qxip))) continue;

               double Coefficient = (DegEi / DegXi) *
                  tensor_coefficient(qyj, Transy, qEj,
       				     qEp, qxy, TransE,
                                     qxip, Transx, qEi);

	       if (fabs(Coefficient) > 1E-10)
	       {
		  add_element(Result.data(), xin.index2(), yin.index2(), (Factor * Coefficient) *
                              herm(*xin) * (*Ein) * (*yin));
	       }
	    }
	 }
      }
   }

   Result.debug_check_structure();
}


template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, 
          typename T3, typename B4, typename S3, typename T4, typename S4>
void
add_triple_prod(Tensor::IrredTensor<T4, B1, B4, S4>& Result,
                std::complex<double> Factor,
                Tensor::IrredTensor<T1, B1, B2, S1> const& x, 
                Tensor::IrredTensor<T2, B2, B3, S2> const& E,
                HermitianProxy<Tensor::IrredTensor<T3, B4, B3, S3> > const& y,
                QuantumNumber qxy,
                QuantumNumber qEp)
{
   if (x.is_null() || E.is_null() || y.base().is_null())
      return;

   DEBUG_PRECONDITION_EQUAL(x.GetSymmetryList(), y.base().GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(x.GetSymmetryList(), qxy.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(x.GetSymmetryList(), E.GetSymmetryList());
   DEBUG_PRECONDITION_EQUAL(x.GetSymmetryList(), qEp.GetSymmetryList());
   QuantumNumber Transx = x.TransformsAs();
   QuantumNumber Transy = y.base().TransformsAs();
   QuantumNumber TransE = E.TransformsAs();
   DEBUG_PRECONDITION(is_transform_target(Transy, qxy, Transx))(Transy)(qxy)(Transx);
   DEBUG_PRECONDITION(is_transform_target(TransE, qxy, qEp))(TransE)(qxy)(qEp);
   DEBUG_PRECONDITION_EQUAL(x.Basis2(), E.Basis1());
   DEBUG_PRECONDITION_EQUAL(E.Basis2(), y.base().Basis2());
   DEBUG_PRECONDITION_EQUAL(Result.Basis1(), x.Basis1());
   DEBUG_PRECONDITION_EQUAL(Result.Basis2(), y.base().Basis1());
   DEBUG_PRECONDITION_EQUAL(Result.TransformsAs(), qEp);

   IrredTensor<T3, B4, B3, S3> const& yb = y.base();  // shortcut

   typedef typename const_iterator< IrredTensor<T1, B1, B2, S1> >::type       x_iter;
   typedef typename const_inner_iterator< IrredTensor<T1, B1, B2, S1> >::type x_inner;

   typedef typename const_iterator< IrredTensor<T2, B2, B3, S2> >::type       E_iter;
   typedef typename const_inner_iterator< IrredTensor<T2, B2, B3, S2> >::type E_inner;

   typedef typename const_iterator< IrredTensor<T3, B4, B3, S3> >::type       y_iter;
   typedef typename const_inner_iterator< IrredTensor<T3, B4, B3, S3> >::type y_inner;

   for (x_iter xi = iterate(x); xi; ++xi)
   {
      QuantumNumber qxp = x.Basis1()[xi.index()];
      for (y_iter yi = iterate(yb); yi; ++yi)
      {
         QuantumNumber qyp = yb.Basis1()[yi.index()];
         if (!is_transform_target(qyp, qEp, qxp)) continue;

         for (x_inner xin = iterate(xi); xin; ++xin)
         {
            for (y_inner yin = iterate(yi); yin; ++yin)
            {
               E_inner Ein = iterate_at(E.data(), xin.index2(), yin.index2());
               if (Ein)
               {
                  /*
                  double Coefficient = tensor_coefficient(x.Basis2()[xin.index2()], Transx, qxp,
							  TransE, qxy, qEp, 
							  yb.Basis2()[yin.index2()], Transy, qyp);
                  */
                  double Coefficient = tensor_coefficient(yb.Basis2()[yin.index2()], Transy, qyp,
							  TransE, qxy, qEp, 
                                                          x.Basis2()[xin.index2()], Transx, qxp);

		  if (fabs(Coefficient) > 1E-10)
		  {
		     add_element(Result.data(), xin.index1(), yin.index1(),
                                 (Factor * Coefficient) * (*xin) * E(xin.index2(), yin.index2()) * herm(*yin));
                  }
	       }
	    }
	 }
      }
   }

   Result.debug_check_structure();
}


template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, 
          typename T3, typename B4, typename S3, typename T4, typename S4>
inline
void
add_triple_prod(Tensor::IrredTensor<T4, B1, B4, S4>& Result,
                std::complex<double> Factor,
                HermitianProxy<Tensor::IrredTensor<T1, B1, B2, S1> > const& x, 
                Tensor::IrredTensor<T2, B1, B3, S2> const& E,
                Tensor::IrredTensor<T3, B3, B4, S3> const& y)
{
   add_triple_prod(Result, Factor, x, E, y, QuantumNumber(y.GetSymmetryList()), E.TransformsAs());
}

template <typename T1, typename B1, typename B2, typename S1, 
          typename T2, typename B3, typename S2, 
          typename T3, typename B4, typename S3, typename T4, typename S4>
inline
void
add_triple_prod(Tensor::IrredTensor<T4, B1, B4, S4>& Result,
                std::complex<double> Factor,
                Tensor::IrredTensor<T1, B1, B2, S1> const& x, 
                Tensor::IrredTensor<T2, B2, B3, S2> const& E,
                HermitianProxy<Tensor::IrredTensor<T3, B4, B3, S3> > const& y)
{
   add_triple_prod(Result, Factor, x, E, y, QuantumNumber(x.GetSymmetryList()), E.TransformsAs());
}

//
// prod for hermitian
//

template <typename T1, typename B1, typename B2, typename S1,
           typename T2, typename B3, typename S2>
typename IrredProd_Herm<IrredTensor<T2, B1, B2, S2>, 
               HermitianProxy<IrredTensor<T1, B3, B2, S1> > >::result_type
IrredProd_Herm<IrredTensor<T2, B1, B2, S2>, 
               HermitianProxy<IrredTensor<T1, B3, B2, S1> > >
::operator()(first_argument_type x, second_argument_type y) const
{
   DEBUG_PRECONDITION_EQUAL(x.Basis2(), y.base().Basis2());
   DEBUG_CHECK(is_scalar(y.base().TransformsAs()));

   result_type Result(x.Basis1(), y.base().Basis1(), x.TransformsAs());
   QuantumNumber Trans = x.TransformsAs();
   QuantumNumber Ident(x.GetSymmetryList());
   IrredTensor<T1, B3, B2, S1> const& yb = y.base();  // shortcut

   typedef typename const_iterator< IrredTensor<T1, B1, B2, S1> >::type       x_iter;
   typedef typename const_inner_iterator< IrredTensor<T1, B1, B2, S1> >::type x_inner;

   typedef typename const_iterator< IrredTensor<T2, B1, B3, S2> >::type       y_iter;
   typedef typename const_inner_iterator< IrredTensor<T2, B1, B3, S2> >::type y_inner;

   for (x_iter xi = iterate(x.data()); xi; ++xi)
   {
      QuantumNumber qp = x.Basis1()[xi.index1()];
      for (x_inner xin = iterate(xi); xin; ++xin)
      {
         QuantumNumber qpp = x.Basis2()[xin.index2()];

         for (y_iter yi = iterate(yb); yi; ++yi)
         {
            for (y_inner yin = iterate(yi); yin; ++yin)
            {
               if (yin.index2() != xin.index2())
                  continue;

               // because y is a scalar, q == qpp always
               // QuantumNumber q = y.Basis2()[yin.index2()];

               // There is no extra coefficient from the transpose of a scalar operator
               double Coefficient = product_coefficient(Trans, Ident, Trans,
                                                     qp, qpp, qpp);
               if (fabs(Coefficient) > 1E-10)
               {
                  add_element(Result.data(), xin.index1(), yin.index1(), Coefficient *
                              (*xin) * herm(*yin));
               }
            }
         }
      }
   }
   Result.debug_check_structure();
   return Result;
}

template <typename T1, typename B1, typename B2, typename S1,
           typename T2, typename B3, typename S2>
typename IrredProd_Herm<HermitianProxy<IrredTensor<T1, B1, B2, S1> >, 
               IrredTensor<T2, B1, B3, S2> >::result_type
IrredProd_Herm<HermitianProxy<IrredTensor<T1, B1, B2, S1> >, 
               IrredTensor<T2, B1, B3, S2> >
::operator()(first_argument_type x, second_argument_type y) const
{
   DEBUG_PRECONDITION_EQUAL(x.base().Basis1(), y.Basis1());
   DEBUG_PRECONDITION(is_scalar(x.base().TransformsAs()));

   result_type Result(x.base().Basis2(), y.Basis2(), y.TransformsAs());
   QuantumNumber Trans = y.TransformsAs();
   QuantumNumber Ident(y.GetSymmetryList());
   IrredTensor<T1, B1, B2, S1> const& xb = x.base();  // shortcut

   typedef typename const_iterator< IrredTensor<T1, B1, B2, S1> >::type       x_iter;
   typedef typename const_inner_iterator< IrredTensor<T1, B1, B2, S1> >::type x_inner;

   typedef typename const_iterator< IrredTensor<T2, B1, B3, S2> >::type       y_iter;
   typedef typename const_inner_iterator< IrredTensor<T2, B1, B3, S2> >::type y_inner;

   y_iter yi = iterate(y.data());  // always in sync with xi
   for (x_iter xi = iterate(xb.data()); xi; ++xi, ++yi)
   {
      QuantumNumber qpp = y.Basis1()[yi.index1()];
      for (x_inner xin = iterate(xi); xin; ++xin)
      {
         // because x is a scalar, qp == qpp always
         // QuantumNumber qp = x.base().Basis2()[xin.index2()];
         for (y_inner yin = iterate(yi); yin; ++yin)
         {
            QuantumNumber q = y.Basis2()[yin.index2()];

            //if (!(is_transform_target(q, Result.TransformsAs(), qpp))) continue;

            // There is no extra coefficient from the transpose of a scalar operator
            double Coefficient = product_coefficient(Ident, Trans, Trans,
                                                     qpp, q, qpp);
            if (fabs(Coefficient) > 1E-10)
	    {
               add_element(Result.data(), xin.index2(), yin.index2(), Coefficient *
                           herm(*xin) * (*yin));
            }
         }
      }
   }
   Result.debug_check_structure();
   return Result;
}

//
// delta_shift
//

template <typename T, typename B1, typename B2, typename S>
void
Tensor::IrredTensor<T, B1, B2, S>::delta_shift(QuantumNumbers::QuantumNumber const& q)
{
   Basis1_.delta_shift(q);
   Basis2_.delta_shift(q);
}

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
delta_shift(Tensor::IrredTensor<T, B1, B2, S> const& x, 
	    QuantumNumbers::QuantumNumber q,
	    QuantumNumbers::Projection p,
	    B1 const& NewBasis1, B2 const& NewBasis2)
{
   if (x.is_null())
   {
      DEBUG_TRACE("delta_shift: null tensor");
      return x;
   }

   DEBUG_CHECK_EQUAL(NewBasis1, delta_shift(x.Basis1(), p));
   DEBUG_CHECK_EQUAL(NewBasis2, delta_shift(x.Basis2(), p));

   typedef typename IrredTensor<T, B1, B2, S>::const_iterator       const_iterator;
   typedef typename IrredTensor<T, B1, B2, S>::const_inner_iterator const_inner_iterator;

   QuantumNumber Ident(x.GetSymmetryList());

   IrredTensor<T, B1, B2, S> Result(NewBasis1, NewBasis2, x.TransformsAs());

   Result.data() = x.data();

#if 0
   // This was a misguided attempt for non-abelian shifts (which dont actually exist)
   for (const_iterator I = iterate(x); I; ++I)
   {
      for (const_inner_iterator J = iterate(I); J; ++J)
      {
	 // The scale factor here depends on the normalization of the coupling
	 // coefficients.

	 //double Scale = std::sqrt(double(degree(x.qn1(J.index1())))
	 // degree(NewBasis1[J.index1()]));
		 //double Scale = 1.0;
	 double Scale = delta_shift_coefficient(x.qn1(J.index1()), x.TransformsAs(), 
                                                x.qn2(J.index2()), q);

	 Result(J.index1(), J.index2()) = Scale * (*J);
      }
   }
#endif
   return Result;
}

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
delta_shift(Tensor::IrredTensor<T, B1, B2, S> const& x, 
	    QuantumNumbers::QuantumNumber q, 
	    QuantumNumbers::Projection p)
{
   return delta_shift(x, q, p, delta_shift(x.Basis1(), p), delta_shift(x.Basis2(), p));
}

template <typename T, typename B1, typename B2, typename S>
Tensor::IrredTensor<T, B1, B2, S>
delta_shift(Tensor::IrredTensor<T, B1, B2, S> const& x, 
	    QuantumNumbers::QuantumNumber q)
{
   QuantumNumbers::ProjectionList PL = enumerate_projections(q);
   DEBUG_PRECONDITION_EQUAL(q.degree(), 1);

   return delta_shift(x, q, PL[0], delta_shift(x.Basis1(), PL[0]), delta_shift(x.Basis2(), PL[0]));
}

} // namespace Tensor
