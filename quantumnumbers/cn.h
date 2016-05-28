// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/cn.h
//
// Copyright (C) 2003-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  Cn.h

  Created 2003-10-19 Ian McCulloch

  Quantum numbers for the cyclic group C_n.
  Allows for all names of the form "C_n" for integer n.

*/

#if !defined(CN_H_HSDJFH3485789UYHFUIY457UHFUU57YUIY)
#define CN_H_HSDJFH3485789UYHFUIY457UHFUU57YUIY

#include "common/niftycounter.h"
#include "quantumnumber.h"
#include <boost/lexical_cast.hpp>

namespace QuantumNumbers
{

class CnProjection;

class CnFactory;

class Cn
{
   public:
      typedef CnProjection  ProjectionType;
      typedef CnFactory     FactoryType;

      explicit Cn(int N_, int x_ = 0) : N(N_), x((x_ + 10 * N_) % N_) {}
      explicit Cn(int N_, int const* InIter) : N(N_), x((*InIter + 10*N) % N) {}
      explicit Cn(int N_, std::string const& s);

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { *OutIter = x; return OutIter+1; }

      std::string Type() { return "C_" + boost::lexical_cast<std::string>(N); }

      static int Size() { return 1; }

      // Registation is automatic via a nifty counter
      static void Register();

      int N;
      int x;
};

class CnProjection
{
   public:
     typedef Cn QuantumNumberType;

      explicit CnProjection(int N_) : N(N_) {}
      explicit CnProjection(int N_, int const* InIter) : N(N_) {}
      explicit CnProjection(int N_, std::string const& s) : N(N_) {}

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { return OutIter; }

      std::string Type() const { return "C_" + boost::lexical_cast<std::string>(N); }
      static char const* Suffix() { return ""; }

      static int Size() { return 0; }

      int N;
};

struct CnFactory
{
   typedef Cn           QuantumNumberType;
   typedef CnProjection ProjectionType;

    explicit CnFactory(int N_) : N(N_) {}

    QuantumNumberType MakeQuantumNumber() const { return QuantumNumberType(N); }
    QuantumNumberType MakeQuantumNumber(int const* q) const { return QuantumNumberType(N, q); }
    QuantumNumberType MakeQuantumNumber(std::string const& s) const { return QuantumNumberType(N, s); }

    ProjectionType MakeProjection(int const* p) const { return ProjectionType(N, p); }
    ProjectionType MakeProjection(std::string const& s) const { return ProjectionType(N, s); }

   std::string Type() const { return "C_" + boost::lexical_cast<std::string>(N); }
   std::string ProjectionSuffix() const { return ProjectionType::Suffix(); }

   int N;
};

//
// inlines
//

inline
int degree(Cn const&)
{
   return 1;
}

inline
double trace(Cn const&)
{
   return 1;
}

inline
int multiplicity(Cn const& q1, Cn const& q2, Cn const& q)
{
   DEBUG_PRECONDITION(q1.N == q.N);
   DEBUG_PRECONDITION(q2.N == q.N);
   DEBUG_PRECONDITION((q1.x + q2.x)%q1.N == q.x);
   return 1;
}

inline
double CG(Cn const& q1, Cn const& q2, Cn const& q,
	  CnProjection const& m1,  CnProjection const& m2,  CnProjection const& m)
{
   DEBUG_PRECONDITION(q1.N == q.N);
   DEBUG_PRECONDITION(q2.N == q.N);
   DEBUG_PRECONDITION(m1.N == q.N);
   DEBUG_PRECONDITION(m2.N == q.N);
   DEBUG_PRECONDITION(m.N == q.N);
   return ((q1.x + q2.x)%q.N == q.x) ? 1 : 0;
}

inline
double product_coefficient(Cn const& k1, Cn const& k2, Cn const& k,
			  Cn const& qp, Cn const& q, Cn const& qpp)
{
   return ((k1.x + k2.x == k.x) && (qp.x == k1.x + qpp.x) && (qpp.x == k2.x + q.x)) ? 1 : 0;
}

inline
double tensor_coefficient(Cn const& k1,  Cn const& k2,  Cn const& k,
			 Cn const& q1p, Cn const& q2p, Cn const& qp,
			 Cn const& q1,  Cn const& q2,  Cn const& q)
{
   int N = k.N;
   DEBUG_PRECONDITION(k1.N == N);
   DEBUG_PRECONDITION(k2.N == N);
   DEBUG_PRECONDITION(q1p.N == N);
   DEBUG_PRECONDITION(q2p.N == N);
   DEBUG_PRECONDITION(qp.N == N);
   DEBUG_PRECONDITION(q1.N == N);
   DEBUG_PRECONDITION(q2.N == N);
   DEBUG_PRECONDITION(q.N == N);
   PRECONDITION((k1.x + k2.x) % N == k.x);
   PRECONDITION((q1p.x + q2p.x) % N == qp.x);
   PRECONDITION((q1.x + q2.x) % N == q.x);
   PRECONDITION((q1p.x == k1.x) % N + q1.x);
   PRECONDITION((q2p.x == k2.x) % N + q2.x);
   PRECONDITION((qp.x == k.x) % N + q.x);
   return 1;
}

inline
double recoupling(Cn const& q1, Cn const& q3, Cn const& q13,
		  Cn const& q2, Cn const& q, Cn const& q23)
{
   int N = q.N;
   DEBUG_PRECONDITION(q1.N == N);
   DEBUG_PRECONDITION(q2.N == N);
   DEBUG_PRECONDITION(q3.N == N);
   DEBUG_PRECONDITION(q13.N == N);
   DEBUG_PRECONDITION(q23.N == N);
   return (q13.x == (q1.x + q3.x)%N) && (q23.x == (q2.x + q3.x)%N) && (q.x == (q1.x+q2.x+q3.x)%N) && 
     ((q1.x + q23.x)%N == q.x) && ((q13.x + q2.x)%N == q.x) ? 1 : 0;
}

inline
Cn adjoint(Cn const& q)
{
  return Cn(q.N, q.N - q.x);
}

inline
double adjointCoefficient(Cn const& qp, Cn const& k, Cn const& q)
{
   int N = k.N;
   DEBUG_PRECONDITION(qp.N == q.N);
   DEBUG_PRECONDITION(k.N == q.N);
   PRECONDITION(qp.x == (k.x + q.x) % N);
   return 1;
}

inline
bool is_transform_target(Cn const& q1, Cn const& q2, Cn const& q)
{
   DEBUG_PRECONDITION(q1.N == q.N);
   DEBUG_PRECONDITION(q2.N == q.N);
   return (q1.x + q2.x) % q1.N == q.x;
}

inline
int num_transform_targets(Cn const& q1, Cn const& q2)
{
   return 1;
}

template <typename OutIter>
inline
void transform_targets(Cn const& q1, Cn const& q2, OutIter Out)
{
   DEBUG_PRECONDITION(q1.N == q2.N);
   *Out++ = Cn(q1.N, (q1.x + q2.x) % q1.N);
}

inline
int NumInversetransform_targets(Cn const& q1, Cn const& q)
{
   DEBUG_PRECONDITION(q1.N == q.N);
   return 1;
}

template <typename OutIter>
inline
void Inversetransform_targets(Cn const& q1, Cn const& q, OutIter Out)
{
   DEBUG_PRECONDITION(q1.N == q2.N);
   *Out++ = (q.x - q1.x + q1.N) % q1.N;
}

template <typename OutIter>
inline
void enumerate_projections(Cn const& q, OutIter Out)
{
   *Out++ = CnProjection(q.N);
}

inline
bool is_delta(Cn const& q1, Cn const& Q, CnProjection const& P, Cn const& q2)
{
   // note: Q is unused here.
   return (q1.x + Q.x) & Q.N == q2.x;
}

//
// CnProjection
//

namespace
{
   NiftyCounter::nifty_counter<Cn::Register> CnRegistrator;
}

} // namespace QuantumNumbers

#endif
