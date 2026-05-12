// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/zn.cpp
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "zn.h"
#include "symmetrybase.h"
#include <algorithm>
#include <map>

namespace QuantumNumbers
{

namespace detail
{

using ZnTemplateRegistryType = std::map<int, ZnSymmetryCreator>;

ZnTemplateRegistryType& ZnTemplateRegistry()
{
   static ZnTemplateRegistryType Registry;
   return Registry;
}

void RegisterZnTemplate(int n, ZnSymmetryCreator Creator)
{
   if (n >= 2)
      ZnTemplateRegistry()[n] = Creator;
}

SymmetryBase* CreateRegisteredZnTemplate(int n)
{
   auto const i = ZnTemplateRegistry().find(n);
   return i == ZnTemplateRegistry().end() ? nullptr : i->second();
}

} // namespace detail

class RuntimeZnFactory;

class RuntimeZn
{
   public:
      typedef RuntimeZn        ProjectionType;
      typedef RuntimeZn        QuantumNumberType;
      typedef RuntimeZnFactory FactoryType;

      RuntimeZn() : RuntimeZn(2, 0) {}
      explicit RuntimeZn(int N_, half_int x_ = 0);
      explicit RuntimeZn(int N_, int const* InIter);
      explicit RuntimeZn(int N_, std::string const& s);

      std::string ToString() const { return ConvertToString(x); }

      int* Convert(int* OutIter) const
      { *OutIter = x.twice(); return OutIter+1; }

      std::string Type() const { return "Z_" + std::to_string(N); }

      static int Size() { return 1; }
      static int num_casimir() { return 1; }
      static std::string casimir_name(std::string const& QName, int) { return QName; }
      static char const* Suffix() { return ""; }

      bool operator==(RuntimeZn const& y) const { return N == y.N && x == y.x; }
      bool operator!=(RuntimeZn const& y) const { return !(*this == y); }

      int N;
      half_int x;
};

class RuntimeZnFactory
{
   public:
      typedef RuntimeZn QuantumNumberType;
      typedef RuntimeZn ProjectionType;

      explicit RuntimeZnFactory(int N_ = 2) : N(N_) {}

      QuantumNumberType MakeQuantumNumber() const { return QuantumNumberType(N); }
      QuantumNumberType MakeQuantumNumber(int const* q) const { return QuantumNumberType(N, q); }
      QuantumNumberType MakeQuantumNumber(std::string const& s) const { return QuantumNumberType(N, s); }

      ProjectionType MakeProjection(int const* p) const { return ProjectionType(N, p); }
      ProjectionType MakeProjection(std::string const& s) const { return ProjectionType(N, s); }

      std::string Type() const { return "Z_" + std::to_string(N); }
      char const* ProjectionSuffix() const { return ProjectionType::Suffix(); }

      int N;
};

half_int NormalizeZnValue(half_int x, int n)
{
   CHECK(n >= 2);
   int const Modulus = 2*n;
   int r = x.twice() % Modulus;
   if (r < 0)
      r += Modulus;
   return from_twice(r);
}

RuntimeZn::RuntimeZn(int N_, half_int x_)
   : N(N_), x(NormalizeZnValue(x_, N_))
{
}

RuntimeZn::RuntimeZn(int N_, int const* InIter)
   : RuntimeZn(N_, from_twice(*InIter))
{
}

RuntimeZn::RuntimeZn(int N_, std::string const& s)
   : RuntimeZn(N_, convert_string<int>(s.begin(), s.end()))
{
}

RuntimeZn operator+(RuntimeZn x, RuntimeZn y)
{
   CHECK_EQUAL(x.N, y.N);
   return RuntimeZn(x.N, x.x + y.x);
}

RuntimeZn operator-(RuntimeZn x, RuntimeZn y)
{
   CHECK_EQUAL(x.N, y.N);
   return RuntimeZn(x.N, x.x - y.x);
}

RuntimeZn operator-(RuntimeZn x)
{
   return RuntimeZn(x.N, -x.x);
}

std::ostream& operator<<(std::ostream& out, RuntimeZn const& s)
{
   return out << s.x;
}

int degree(RuntimeZn const&)
{
   return 1;
}

double trace(RuntimeZn const&)
{
   return 1;
}

int multiplicity(RuntimeZn const& q1, RuntimeZn const& q2, RuntimeZn const& q)
{
   DEBUG_PRECONDITION_EQUAL(q1+q2, q);
   return 1;
}

bool cross_product_exists(RuntimeZn const&, RuntimeZn const&)
{
   return true;
}

RuntimeZn cross_product_transforms_as(RuntimeZn const& q1, RuntimeZn const& q2)
{
   return q1 + q2;
}

std::complex<double> cross_product_factor(RuntimeZn const&, RuntimeZn const&)
{
   return 1.0;
}

double clebsch_gordan(RuntimeZn const& q1, RuntimeZn const& q2, RuntimeZn const& q,
                      RuntimeZn const& m1, RuntimeZn const& m2, RuntimeZn const& m)
{
   DEBUG_CHECK_EQUAL(q1, m1);
   DEBUG_CHECK_EQUAL(q2, m2);
   DEBUG_CHECK_EQUAL(q, m);
   return (q1+q2 == q) ? 1 : 0;
}

double product_coefficient(RuntimeZn const& k1, RuntimeZn const& k2, RuntimeZn const& k,
                           RuntimeZn const& qp, RuntimeZn const& q, RuntimeZn const& qpp)
{
   return ((k1 + k2 == k) && (qp == k1 + qpp) && (qpp == k2 + q)) ? 1 : 0;
}

double inverse_product_coefficient(RuntimeZn const& k1, RuntimeZn const& k2, RuntimeZn const& k,
                                   RuntimeZn const& qp, RuntimeZn const& q, RuntimeZn const& qpp)
{
   return product_coefficient(k1, k2, k, qp, q, qpp);
}

double tensor_coefficient(RuntimeZn const& k1,  RuntimeZn const& k2,  RuntimeZn const& k,
                          RuntimeZn const& q1p, RuntimeZn const& q2p, RuntimeZn const& qp,
                          RuntimeZn const& q1,  RuntimeZn const& q2,  RuntimeZn const& q)
{
   PRECONDITION(k1 + k2 == k)(k1)(k2)(k);
   PRECONDITION(q1p + q2p == qp)(q1p)(q2p)(qp);
   PRECONDITION(q1 + q2 == q)(q1)(q2)(q);
   PRECONDITION(k1 + q1p == q1)(q1p)(k1)(q1);
   PRECONDITION(k2 + q2p == q2)(q2p)(k2)(q2);
   PRECONDITION(k + qp == q)(qp)(k)(q);
   return 1;
}

double inverse_tensor_coefficient(RuntimeZn const& k1,  RuntimeZn const& k2,  RuntimeZn const& k,
                                  RuntimeZn const& q1p, RuntimeZn const& q2p, RuntimeZn const& qp,
                                  RuntimeZn const& q1,  RuntimeZn const& q2,  RuntimeZn const& q)
{
   return tensor_coefficient(k1, k2, k, q1p, q2p, qp, q1, q2, q);
}

double recoupling(RuntimeZn const& q1, RuntimeZn const& q3, RuntimeZn const& q13,
                  RuntimeZn const& q2, RuntimeZn const& q, RuntimeZn const& q23)
{
   DEBUG_PRECONDITION_EQUAL(q13, q1+q3);
   DEBUG_PRECONDITION_EQUAL(q23, q2+q3);
   DEBUG_PRECONDITION_EQUAL(q, q1+q23);
   return 1;
}

double recoupling_12_3__13_2(RuntimeZn const& q1, RuntimeZn const& q2, RuntimeZn const& q12,
                             RuntimeZn const& q3, RuntimeZn const& q, RuntimeZn const& q13)
{
   DEBUG_PRECONDITION_EQUAL(q12, q1+q2);
   DEBUG_PRECONDITION_EQUAL(q13, q1+q3);
   DEBUG_PRECONDITION_EQUAL(q, q12+q3);
   return 1;
}

RuntimeZn adjoint(RuntimeZn const& q)
{
   return RuntimeZn(q.N, -q.x);
}

double adjoint_coefficient(RuntimeZn const& qp, RuntimeZn const& k, RuntimeZn const& q)
{
   PRECONDITION_EQUAL(qp, q-k);
   return 1;
}

double conj_phase(RuntimeZn const& qp, RuntimeZn const& k, RuntimeZn const& q)
{
   PRECONDITION_EQUAL(qp, q-k);
   return 1;
}

bool is_transform_target(RuntimeZn const& q1, RuntimeZn const& q2, RuntimeZn const& q)
{
   return q == q1+q2;
}

int num_transform_targets(RuntimeZn const&, RuntimeZn const&)
{
   return 1;
}

template <typename OutIter>
void transform_targets(RuntimeZn const& q1, RuntimeZn const& q2, OutIter Out)
{
   *Out++ = q1+q2;
}

int num_inverse_transform_targets(RuntimeZn const&, RuntimeZn const&)
{
   return 1;
}

template <typename OutIter>
void inverse_transform_targets(RuntimeZn const& q1, RuntimeZn const& q, OutIter Out)
{
   *Out++ = q1+q;
}

template <typename OutIter>
void enumerate_projections(RuntimeZn const& q, OutIter Out)
{
   *Out++ = q;
}

bool is_delta(RuntimeZn const& q1, RuntimeZn const& Q, RuntimeZn const& P, RuntimeZn const& q2)
{
   DEBUG_PRECONDITION_EQUAL(Q, P);
   return q1 == P+q2;
}

RuntimeZn difference(RuntimeZn const& q1, RuntimeZn const& q2)
{
   return q1-q2;
}

RuntimeZn negate(RuntimeZn const& p)
{
   return p;
}

RuntimeZn sum(RuntimeZn const& q1, RuntimeZn const& q2)
{
   return q1+q2;
}

bool is_possible(RuntimeZn const&, RuntimeZn const&)
{
   return true;
}

bool is_projection(RuntimeZn const& q, RuntimeZn const& p)
{
   return q == p;
}

RuntimeZn change(RuntimeZn const& q, RuntimeZn const& p)
{
   return q+p;
}

RuntimeZn heighest_weight(RuntimeZn const& p)
{
   return p;
}

double weight(RuntimeZn const& p)
{
   return std::min(p.x.twice(), 2*p.N-p.x.twice());
}

double delta_shift_coefficient(RuntimeZn const& qp, RuntimeZn const& k,
                               RuntimeZn const& q, RuntimeZn const&)
{
   PRECONDITION_EQUAL(qp, q+k);
   return 1;
}

double casimir(RuntimeZn const& s, int i)
{
   CHECK_EQUAL(i, 0);
   return s.x.to_double();
}

class ZnSymmetryFactory : public SymmetryFactory
{
   public:
      SymmetryBase* AttemptCreate(std::string const& Type) override;
};

SymmetryBase* ZnSymmetryFactory::AttemptCreate(std::string const& Type)
{
   if (Type.size() < 3)
      return nullptr;
   if (Type[0] != 'Z' || Type[1] != '_')
      return nullptr;

   int N;
   try
   {
      N = convert_string<int>(Type.begin()+2, Type.end());
   }
   catch (invalid_string_conversion const&)
   {
      return nullptr;
   }

   if (N < 2)
      return nullptr;

   if (SymmetryBase* S = detail::CreateRegisteredZnTemplate(N))
      return S;

   return new BasicSymmetry<RuntimeZn>(RuntimeZnFactory(N));
}

void RegisterZnRuntimeFactory()
{
   RegisterSymmetry(new ZnSymmetryFactory());
}

namespace
{
   NiftyCounter::nifty_counter<RegisterZnRuntimeFactory> ZnRuntimeFactoryRegistrator;
}


} // namespace QuantumNumbers
