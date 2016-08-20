// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/mpopcompressed.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "mpopcompressed.h"
#include "common/trace.h"

//
// do_prod
//

// Visitor to implement the prod() function
struct DoProd : public boost::static_visitor<MPOpCompressed>
{
   DoProd(SimpleOperator& C_, ProductBasis<BasisList, BasisList>& B_)
      : C(C_), B(B_)
   {
      DEBUG_CHECK_EQUAL(C.Basis2(), B.Basis());
   }

   MPOpCompressed operator()(MPOpRepeat const& a, MPOpRepeat const& b) const;
   MPOpCompressed operator()(MPOpRepeat const& a, MPOpArray const& b) const;
   MPOpCompressed operator()(MPOpArray const& a, MPOpRepeat const& b) const;
   MPOpCompressed operator()(MPOpArray const& a, MPOpArray const& b) const;
   MPOpCompressed operator()(MPOpComponent const& a, MPOpComponent const& b) const;

   // The remaining cases are only possible if the representation isnt canonical.  But
   // for the time being, we allow the non-canonical representation
#if 0
   template <typename T, typename U>
   MPOpCompressed operator()(T const& a, U const& b) const
   {
      PANIC("DoProd: representation is not canonical")
         (tracer::typeid_name(a))(tracer::typeid_name(b))(a)(b);
      return MPOpCompressed(); // suppress a warning
   }
#endif

   MPOpCompressed operator()(MPOpComponent const& a, MPOpArray const& b) const;
   MPOpCompressed operator()(MPOpComponent const& a, MPOpRepeat const& b) const;
   MPOpCompressed operator()(MPOpArray const& a, MPOpComponent const& b) const;
   MPOpCompressed operator()(MPOpRepeat const& a, MPOpComponent const& b) const;

   SimpleOperator& C;
   ProductBasis<BasisList, BasisList>& B;
};

MPOpCompressed
DoProd::operator()(MPOpComponent const& a, MPOpComponent const& b) const
{
   DEBUG_CHECK_EQUAL(a.Basis1(), B.Left());
   DEBUG_CHECK_EQUAL(b.Basis1(), B.Right());

   ProductBasis<BasisList, BasisList> B2 = make_product_basis(a.Basis2(), b.Basis2());
   MPOpComponent Next = prod(C, mp_prod(a, b, B, B2));
   C = TruncateBasis2(Next);
   B = B2;
   return MPOpCompressed(Next);
}

MPOpCompressed
DoProd::operator()(MPOpComponent const& a, MPOpArray const& b) const
{
   CHECK(b.size() == 1);
   return apply_visitor(*this, MPOpCompressed(a), b.data().front());
}

MPOpCompressed
DoProd::operator()(MPOpComponent const& a, MPOpRepeat const& b) const
{
   CHECK(b.size() == 1);
   return apply_visitor(*this, MPOpCompressed(a), b.nested());
}

MPOpCompressed
DoProd::operator()(MPOpArray const& a, MPOpComponent const& b) const
{
   CHECK(a.size() == 1);
   return apply_visitor(*this, a.data().front(), MPOpCompressed(b));
}

MPOpCompressed
DoProd::operator()(MPOpRepeat const& a, MPOpComponent const& b) const
{
   CHECK(a.size() == 1);
   return apply_visitor(*this, a.nested(), MPOpCompressed(b));
}

MPOpCompressed
DoProd::operator()(MPOpRepeat const& a, MPOpRepeat const& b) const
{
   //DEBUG_TRACE(a.size())(a.nested().size())(b.size())(b.nested().size());
   CHECK_EQUAL(a.size(), b.size());
   CHECK_EQUAL(a.nested().size(), b.nested().size());

   // iterate until we get a fixed point
   MPOpCompressed Result;
   int Sz = 0;
   while (Sz < a.size())
   {
      SimpleOperator Cold = C;
      MPOpCompressed v = do_prod(C, B, a.nested(), b.nested());
      // do we have a cycle?
      if (C.Basis1() == Cold.Basis1() && C.Basis2() == Cold.Basis2() && norm_frob_sq(C-Cold)<1E-12)
      {
         Result.push_back(MPOpRepeat(a.size()-Sz, v));
         return Result;
      }
      // else
      Result.push_back(v);
      ++Sz;
   }
   return Result;
}

MPOpCompressed
DoProd::operator()(MPOpRepeat const& a, MPOpArray const& b) const
{
   //DEBUG_TRACE(a.size())(a.nested().size())(b.size())(a.logical_size())(b.logical_size());
   MPOpArray Result;
   int const a_nested_size = a.nested().size();
   // iterate through b
   MPOpArray::const_iterator I = b.begin();
   while (I != b.end())
   {
      // if I is a multiple of the nested unit
      if (I->size() % a_nested_size == 0)
      {
         Result.push_back(do_prod(C, B, MPOpRepeat(I->size()/a.nested().size(), a.nested()), *I));
         ++I;
      }
      else
      {
         // in this case, the array must be a split unit cell
         MPOpArray::const_iterator First = I;
         int Sz = I->size();
         ++I;
         while (Sz < a_nested_size)
         {
            CHECK(I != b.end());
            Sz += I->size();
            ++I;
         }
         Result.push_back(do_prod(C, B, a.nested(), MPOpArray(First, I)));
      }
   }
   return Result;
}

MPOpCompressed
DoProd::operator()(MPOpArray const& a, MPOpRepeat const& b) const
{
   //DEBUG_TRACE(a.size())(b.size())(b.nested().size());
   MPOpArray Result;
   int const b_nested_size = b.nested().size();
   // iterate through a
   MPOpArray::const_iterator I = a.begin();
   while (I != a.end())
   {
      // if I is a multiple of the nested unit
      if (I->size() % b_nested_size == 0)
      {
         Result.push_back(do_prod(C, B, *I, MPOpRepeat(I->size()/b.nested().size(), b.nested())));
         ++I;
      }
      else
      {
         // in this case, the array must be a split unit cell
         MPOpArray::const_iterator First = I;
         int Sz = I->size();
         ++I;
         while (Sz < b_nested_size)
         {
            CHECK(I != a.end());
            Sz += I->size();
            ++I;
         }
         Result.push_back(do_prod(C, B, MPOpArray(First, I), b.nested()));
      }
   }
   return Result;
}

MPOpCompressed
DoProd::operator()(MPOpArray const& a, MPOpArray const& b) const
{
   //DEBUG_TRACE(a.size())(b.size());
   // this is the messy case
   MPOpCompressed Result;
   MPOpArray::const_iterator Ia = a.begin(), Ib = b.begin();
   DEBUG_CHECK(Ia != a.end() && Ib != b.end());
   MPOpCompressed Va = *Ia;
   MPOpCompressed Vb = *Ib;

   while (true)
   {
      while (Va.size() == Vb.size())
      {
         Result.push_back(do_prod(C, B, Va, Vb));
         ++Ia;
         ++Ib;
         if (Ia == a.end())
         {
            CHECK(Ib == b.end());
            return MPOpCompressed(Result);
         }
         else
         {
            Va = *Ia;
            Vb = *Ib;
         }
      }

      while (Va.size() < Vb.size())
      {
         MPOpCompressed Temp;
         std::tie(Temp, Vb) = split_operator(Vb, Va.size());
         //TRACE(Va.size())(Temp.size())(Vb.size());
         Result.push_back(do_prod(C, B, Va, Temp));
         ++Ia; Va = *Ia;  // Ia is never singular as the operators are the same overall size
      }

      while (Vb.size() < Va.size())
      {
         MPOpCompressed Temp;
         std::tie(Temp, Va) = split_operator(Va, Vb.size());
         //TRACE(Va.size())(Temp.size())(Vb.size());
         Result.push_back(do_prod(C, B, Temp, Vb));
         ++Ib; Vb = *Ib;  // Ib is never singular as the operators are the same overall size
      }
   }
}

MPOpCompressed do_prod(SimpleOperator& C, ProductBasis<BasisList, BasisList>& B,
                       MPOpCompressed const& a, MPOpCompressed const& b)
{
   return apply_visitor(DoProd(C, B), a, b);
}

//
// do_sum
//

// Visitor to implement the prod() function
struct DoSum : public boost::static_visitor<MPOpCompressed>
{
   DoSum(SimpleOperator& C_, SumBasis<BasisList>& B_)
      : C(C_), B(B_)
   {
      DEBUG_CHECK_EQUAL(C.Basis2(), B);
   }

   MPOpCompressed operator()(MPOpRepeat const& a, MPOpRepeat const& b) const;
   MPOpCompressed operator()(MPOpRepeat const& a, MPOpArray const& b) const;
   MPOpCompressed operator()(MPOpArray const& a, MPOpRepeat const& b) const;
   MPOpCompressed operator()(MPOpArray const& a, MPOpArray const& b) const;
   MPOpCompressed operator()(MPOpComponent const& a, MPOpComponent const& b) const;

#if 0
   // The remaining cases are only possible if the representation isnt canonical
   template <typename T, typename U>
   MPOpCompressed operator()(T const& a, U const& b) const
   {
      PANIC("DoSum: representation is not canonical")
         (tracer::typeid_name(a))(tracer::typeid_name(b));
      return MPOpCompressed(); // suppress a warning
   }
#endif

   MPOpCompressed operator()(MPOpComponent const& a, MPOpArray const& b) const;
   MPOpCompressed operator()(MPOpComponent const& a, MPOpRepeat const& b) const;
   MPOpCompressed operator()(MPOpArray const& a, MPOpComponent const& b) const;
   MPOpCompressed operator()(MPOpRepeat const& a, MPOpComponent const& b) const;

   SimpleOperator& C;
   SumBasis<BasisList>& B;
};

MPOpCompressed
DoSum::operator()(MPOpComponent const& a, MPOpComponent const& b) const
{
   DEBUG_CHECK_EQUAL(a.Basis1(), B.Basis(0));
   DEBUG_CHECK_EQUAL(b.Basis1(), B.Basis(1));

   SumBasis<BasisList> B2(a.Basis2(), b.Basis2());
   MPOpComponent Next = prod(C, tensor_sum(a, b, B, B2));
   C = TruncateBasis2(Next);
   B = B2;
   return MPOpCompressed(Next);
}

MPOpCompressed
DoSum::operator()(MPOpComponent const& a, MPOpArray const& b) const
{
   CHECK(b.size() == 1);
   return apply_visitor(*this, MPOpCompressed(a), b.data().front());
}

MPOpCompressed
DoSum::operator()(MPOpComponent const& a, MPOpRepeat const& b) const
{
   CHECK(b.size() == 1);
   return apply_visitor(*this, MPOpCompressed(a), b.nested());
}

MPOpCompressed
DoSum::operator()(MPOpArray const& a, MPOpComponent const& b) const
{
   CHECK(a.size() == 1);
   return apply_visitor(*this, a.data().front(), MPOpCompressed(b));
}

MPOpCompressed
DoSum::operator()(MPOpRepeat const& a, MPOpComponent const& b) const
{
   CHECK(a.size() == 1);
   return apply_visitor(*this, a.nested(), MPOpCompressed(b));
}


MPOpCompressed
DoSum::operator()(MPOpRepeat const& a, MPOpRepeat const& b) const
{
   CHECK_EQUAL(a.size(), b.size());
   CHECK_EQUAL(a.nested().size(), b.nested().size());

   // iterate until we get a fixed point
   MPOpCompressed Result;
   int Sz = 0;
   while (Sz < a.size())
   {
      SimpleOperator Cold = C;
      MPOpCompressed v = do_sum(C, B, a.nested(), b.nested());
      DEBUG_CHECK_EQUAL(C.Basis1(), Cold.Basis1());
      if (C.Basis2() == Cold.Basis2() && norm_frob_sq(C-Cold)<1E-12)
      {
         Result.push_back(MPOpRepeat(a.size()-Sz, v));
         return Result;
      }
      // else
      Result.push_back(v);
      ++Sz;
   }
   return Result;
}

MPOpCompressed
DoSum::operator()(MPOpRepeat const& a, MPOpArray const& b) const
{
   MPOpArray Result;
   int const a_nested_size = a.nested().size();
   // iterate through b
   MPOpArray::const_iterator I = b.begin();
   while (I != b.end())
   {
      //      TRACE(I->size());
      // if I is a multiple of the nested unit
      if (I->size() % a_nested_size == 0)
      {
         Result.push_back(do_sum(C, B, MPOpRepeat(I->size()/a.nested().size(), a.nested()), *I));
         ++I;
      }
      else
      {
         // in this case, the array must be a split unit cell
         MPOpArray::const_iterator First = I;
         int Sz = I->size();
         ++I;
         while (Sz < a_nested_size)
         {
            CHECK(I != b.end());
            Sz += I->size();
            ++I;
         }
         Result.push_back(do_sum(C, B, a.nested(), MPOpArray(First, I)));
      }
   }
   return Result;
}

MPOpCompressed
DoSum::operator()(MPOpArray const& a, MPOpRepeat const& b) const
{
   MPOpArray Result;
   int const b_nested_size = b.nested().size();
   // iterate through a
   MPOpArray::const_iterator I = a.begin();
   while (I != a.end())
   {
      // if I is a multiple of the nested unit
      if (I->size() % b_nested_size == 0)
      {
         Result.push_back(do_sum(C, B, *I, MPOpRepeat(I->size()/b.nested().size(), b.nested())));
         ++I;
      }
      else
      {
         // in this case, the array must be a split unit cell
         MPOpArray::const_iterator First = I;
         int Sz = I->size();
         ++I;
         while (Sz < b_nested_size)
         {
            CHECK(I != a.end());
            Sz += I->size();
            ++I;
         }
         Result.push_back(do_sum(C, B, MPOpArray(First, I), b.nested()));
      }
   }
   return Result;
}

MPOpCompressed
DoSum::operator()(MPOpArray const& a, MPOpArray const& b) const
{
   // this is the messy case.  We need to match up elements where the nested sizes may be
   // different.
   MPOpCompressed Result;
   MPOpArray::const_iterator Ia = a.begin(), Ib = b.begin();
   DEBUG_CHECK(Ia != a.end() && Ib != b.end());
   MPOpCompressed Va = *Ia;  // get the values, as we might need to split them
   MPOpCompressed Vb = *Ib;

   while (true)
   {
      while (Va.size() == Vb.size())
      {
         Result.push_back(do_sum(C, B, Va, Vb));
         ++Ia;
         ++Ib;
         if (Ia == a.end())
         {
            CHECK(Ib == b.end());
            return MPOpCompressed(Result);
         }
         else
         {
            Va = *Ia;
            Vb = *Ib;
         }
      }

      while (Va.size() < Vb.size())
      {
         MPOpCompressed Temp;
         std::tie(Temp, Vb) = split_operator(Vb, Va.size());
         Result.push_back(do_sum(C, B, Va, Temp));
         ++Ia; Va = *Ia;  // Ia is never singular as the operators are the same overall size
      }

      while (Vb.size() < Va.size())
      {
         MPOpCompressed Temp;
         std::tie(Temp, Va) = split_operator(Va, Vb.size());
         Result.push_back(do_sum(C, B, Temp, Vb));
         ++Ib; Vb = *Ib;  // Ib is never singular as the operators are the same overall size
      }
   }
}

MPOpCompressed do_sum(SimpleOperator& C, SumBasis<BasisList>& B,
                      MPOpCompressed const& a, MPOpCompressed const& b)
{
   return apply_visitor(DoSum(C, B), a, b);
}

//
// inject_right
//

struct DoInjectRight : boost::static_visitor<MPOpCompressed>
{
   DoInjectRight(SimpleOperator& C_) : C(C_) {}

   MPOpCompressed operator()(MPOpComponent const& a) const
   {
      MPOpComponent Result = prod(a, C);
      C = TruncateBasis1(Result);
      return MPOpCompressed(Result);
   }

   MPOpCompressed operator()(MPOpArray const& a) const
   {
      MPOpArray Result;
      MPOpArray::const_iterator I = a.end();
      while (I != a.begin())
      {
         --I;
         Result.push_front(I->apply_visitor(*this));
      }
      return MPOpCompressed(Result);
   }

   MPOpCompressed operator()(MPOpRepeat const& a) const
   {
      // iterate until we get a fixed point
      MPOpCompressed Result;
      int Sz = 0;
      while (Sz < a.size())
      {
         SimpleOperator Cold = C;
         MPOpCompressed v = a.nested().apply_visitor(*this);
         if (C.Basis1() == Cold.Basis1() && C.Basis2() == Cold.Basis2() && norm_frob_sq(C-Cold)<1E-12)
         {
            Result.push_front(MPOpRepeat(a.size()-Sz, v));
            return Result;
         }
         // else
         Result.push_front(v);
         ++Sz;
      }
      return Result;
   }

   SimpleOperator& C;
};

MPOpCompressed inject_right(MPOpCompressed const& a, SimpleOperator& C)
{
   return a.apply_visitor(DoInjectRight(C));
}

//
// inject_left
//

struct DoInjectLeft : boost::static_visitor<MPOpCompressed>
{
   DoInjectLeft(SimpleOperator& C_) : C(C_) {}

   MPOpCompressed operator()(MPOpComponent const& a) const
   {
      MPOpComponent Result = prod(C, a);
      C = TruncateBasis2(Result);
      return MPOpCompressed(Result);
   }

   MPOpCompressed operator()(MPOpArray const& a) const
   {
      MPOpArray Result;
      std::transform(a.begin(), a.end(), std::back_inserter(Result), boost::apply_visitor(*this));
      return MPOpCompressed(Result);
   }

   MPOpCompressed operator()(MPOpRepeat const& a) const
   {
      // iterate until we get a fixed point
      MPOpCompressed Result;
      int Sz = 0;
      while (Sz < a.size())
      {
         SimpleOperator Cold = C;
         MPOpCompressed v = a.nested().apply_visitor(*this);
         if (C.Basis1() == Cold.Basis1() && C.Basis2() == Cold.Basis2() && norm_frob_sq(C-Cold)<1E-12)
         {
            Result.push_back(MPOpRepeat(a.size()-Sz, v));
            return Result;
         }
         // else
         Result.push_back(v);
         ++Sz;
      }
      return Result;
   }

   SimpleOperator& C;
};

MPOpCompressed inject_left(SimpleOperator& C, MPOpCompressed const& a)
{
   return a.apply_visitor(DoInjectLeft(C));
}

//
// multiply_left
//

struct DoMultiplyLeft : boost::static_visitor<MPOpCompressed>
{
   DoMultiplyLeft(SimpleOperator const& C_) : C(C_) {}

   MPOpCompressed operator()(MPOpComponent const& x) const
   {
      return MPOpCompressed(prod(C, x));
   }

   MPOpCompressed operator()(MPOpArray const& x) const
   {
      MPOpArray Temp(x);
      Temp.data().front() = Temp.data().front().apply_visitor(*this);
      return Temp;
   }

   MPOpCompressed operator()(MPOpRepeat const& x) const
   {
      CHECK(x.size() > 0)(x.size());
      MPOpCompressed Front = x.nested().apply_visitor(*this);
      if (x.size() > 1)
         Front.push_back(MPOpRepeat(x.size()-1, x.nested()));
      return Front;
   }

   SimpleOperator const& C;
};

MPOpCompressed multiply_left(SimpleOperator const& C, MPOpCompressed const& a)
{
   return a.apply_visitor(DoMultiplyLeft(C));
}

// operator<< for MPOpCompressed

// this is actually more general than MPOperator, it probably should be somewhere else
template <typename T>
struct PrintRunLengthCompressed : public boost::static_visitor<>
{
   PrintRunLengthCompressed(std::ostream& out_, int NestLevel_) : out(out_), NestLevel(NestLevel_) {}

   void DoNest() const
   {
      out << std::string(NestLevel*4, ' ');
   }

   void operator()(T const& x) const
   {
      this->DoNest();
      out << "Value:\n";
      this->DoNest();
      out << x << '\n';
   };

   void operator()(run_length_array<T> const& x) const
   {
      this->DoNest();
      out << "Array: size " << x.size() << '\n';
      int i = 0;
      for (typename run_length_array<T>::const_iterator I = x.begin(); I != x.end(); ++I)
      {
         this->DoNest();
         out << "Array element " << i++ << '\n';
         I->apply_visitor(PrintRunLengthCompressed<T>(out, NestLevel+1));
      }
   }

   void operator()(run_length_repeat<T> const& x) const
   {
      this->DoNest();
      out << "Repeat: size " << x.size() << '\n';
      x.nested().apply_visitor(PrintRunLengthCompressed<T>(out, NestLevel+1));
   }

   std::ostream& out;
   int NestLevel;
};

std::ostream& operator<<(std::ostream& out, MPOpCompressed const& x)
{
   x.apply_visitor(PrintRunLengthCompressed<MPOpComponent>(out, 0));
   return out;
}

SimpleOperator CollapseBasis(BasisList const& b)
{
   std::set<QuantumNumber> QN(b.begin(), b.end());
   BasisList NewB(b.GetSymmetryList(), QN.begin(), QN.end());
   SimpleOperator C(NewB, b);
   for (unsigned j = 0; j < b.size(); ++j)
   {
      unsigned i = std::find(NewB.begin(), NewB.end(), b[j]) - NewB.begin();
      C(i,j) = 1.0;
   }
   return C;
}

SimpleOperator RemoveEmptyRows(SimpleOperator const& c)
{
   BasisList b(c.GetSymmetryList());
   const_iterator<SimpleOperator>::type I = iterate(c);
   while (I)
   {
      const_inner_iterator<SimpleOperator>::type J = iterate(I);
      if (J) // do we have any elements?
         b.push_back(c.Basis1()[J.index1()]);
      ++I;
   }

   SimpleOperator Result(b, c.Basis2());
   I = iterate(c);
   iterator<SimpleOperator>::type I2 = iterate(Result);
   while (I)
   {
      const_inner_iterator<SimpleOperator>::type J = iterate(I);
      if (J) // do we have any elements?
      {
         *I2 = *I;
         ++I2;
      }
      ++I;
   }
   return Result;
}

MPOpCompressed project_component(MPOpCompressed const& Res, int n1, int n2)
{
   BasisList NewBasis(Res.front().GetSymmetryList());
   NewBasis.push_back(QuantumNumbers::QuantumNumber(Res.front().Basis1()[n1]));
   SimpleOperator C(NewBasis, Res.front().Basis1());
   C(0,n1) = 1.0;
   MPOpCompressed Temp = inject_left(C, Res);
   BasisList NewRightBasis(C.GetSymmetryList());
   NewRightBasis.push_back(C.Basis2()[n2]);
   SimpleOperator D(C.Basis2(), NewRightBasis);
   D(n2, 0) = 1.0;
   C = C*D;
   Temp = inject_right(Temp, C);
   Temp = multiply_left(C, Temp);
   return Temp;
}

struct DoConj : public boost::static_visitor<void>
{
   void operator()(MPOpRepeat& a) const
   {
      a.nested().apply_visitor(*this);
   }

   void operator()(MPOpArray& a) const
   {
      for (MPOpArray::iterator I = a.begin(); I != a.end(); ++I)
      {
         I->apply_visitor(*this);
      }
   }

   void operator()(MPOpComponent& a) const
   {
      a = conj(a);
   }
};

MPOpCompressed conj(MPOpCompressed const& x)
{
   TRACE("conj")(&x);
   MPOpCompressed Result = x;
   Result.apply_visitor(DoConj());
   return Result;
}

struct DoAdjoint : public boost::static_visitor<void>
{
   void operator()(MPOpRepeat& a) const
   {
      a.nested().apply_visitor(*this);
   }

   void operator()(MPOpArray& a) const
   {
      for (MPOpArray::iterator I = a.begin(); I != a.end(); ++I)
      {
         I->apply_visitor(*this);
      }
   }

   void operator()(MPOpComponent& a) const
   {
      a = local_adjoint(a);
   }
};

MPOpCompressed adjoint(MPOpCompressed const& x)
{
   MPOpCompressed Result = x;
   Result.apply_visitor(DoAdjoint());
   return Result;
}
