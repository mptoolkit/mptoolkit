/* -*- C++ -*- $Id$

  iteroperations_blas.h

  overloads to utilize BLAS level 1 for some vector iterator operations.

  Created 2005-01-09 Ian McCulloch
*/

#if !defined(ITEROPERATIONS_BLAS_H_JHCUHIUHIUREYT78348OY)
#define ITEROPERATIONS_BLAS_H_JHCUHIUHIUREYT78348OY

#include "common/blas1f.h"
#include "vectortransformiterator.h"

namespace LinearAlgebra
{

// no iter_fill equivalent in BLAS

// iter_assign

#if 0
template <typename S1, typename S2>
inline
void iter_assign(VecPtrIterator<double, S1> i1, VecPtrIterator<double const, S2> i2)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.size());
   BLAS::dcopy(i2.size(), i2.blas_base(), i2.stride(), i1.blas_base(), i1.stride());
}

template <typename S1, typename S2>
inline
void iter_assign(VecPtrIterator<std::complex<double>, S1> i1, 
		 VecPtrIterator<std::complex<double> const, S2> i2)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.size());
   BLAS::zcopy(i2.size(), i2.blas_base(), i2.stride(), i1.blas_base(), i1.stride());
}

#endif

// iter_inner_prod

template <typename S1, typename S2>
inline
double iter_inner_prod(VecPtrIterator<double const, S1> i1, 
		       VecPtrIterator<double const, S2> i2,
		       InnerProd<double, double>)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.size());
   return BLAS::ddot(i1.size(), i1.blas_base(), i1.stride(), i2.blas_base(), i2.stride());
}

template <typename S1, typename S2>
inline
double iter_inner_prod(VecPtrIterator<double const, S1> i1, 
		       VecPtrIterator<double const, S2> i2,
		       Multiplication<double, double>)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.size());
   return BLAS::ddot(i1.size(), i1.blas_base(), i1.stride(), i2.blas_base(), i2.stride());
}

// complex inner prod - two variants here, for conjugated and non-conjugated

#if defined(HAVE_ZDOTC)

// no conjugation, inner product
template <typename S1, typename S2>
inline
std::complex<double> iter_inner_prod(VecPtrIterator<std::complex<double> const, S1> i1, 
		       VecPtrIterator<std::complex<double> const, S2> i2,
		       InnerProd<std::complex<double>, std::complex<double> >)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.size());
   return BLAS::zdotc(i1.size(), i1.blas_base(), i1.stride(), 
		      i2.blas_base(), i2.stride());
}

// first vector conjugated, ordinary multiplication
template <typename S1, typename S2>
inline
std::complex<double> 
iter_inner_prod(VectorTransformIterator<VecPtrIterator<
		std::complex<double> const, S1>, Conj<std::complex<double> > > i1, 
		VecPtrIterator<std::complex<double> const, S2> i2,
		Multiplication<std::complex<double>, std::complex<double> >)
{
   DEBUG_PRECONDITION_EQUAL(i1.base().size(), i2.size());
   return BLAS::zdotc(i1.base().size(), i1.base().blas_base(), i1.base().stride(), 
		      i2.blas_base(), i2.stride());
}

// second vector conjugated, ordinary multiplication
template <typename S1, typename S2>
inline
std::complex<double> 
iter_inner_prod(VecPtrIterator<std::complex<double> const, S1> i1,
		VectorTransformIterator<VecPtrIterator<
		std::complex<double> const, S2>, Conj<std::complex<double> > > i2, 
		Multiplication<std::complex<double>, std::complex<double> >)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.base().size());
   return std::conj(BLAS::zdotc(i1.size(), i1.blas_base(), i1.stride(), 
				i2.base().blas_base(), i2.base().stride()));
}

// both vectors conjugated, inner product
template <typename S1, typename S2>
inline
std::complex<double> 
iter_inner_prod(VectorTransformIterator<VecPtrIterator<
		std::complex<double> const, S1>, Conj<std::complex<double> > > i1,
		VectorTransformIterator<VecPtrIterator<
		std::complex<double> const, S2>, Conj<std::complex<double> > > i2, 
		InnerProd<std::complex<double>, std::complex<double> >)
{
   DEBUG_PRECONDITION_EQUAL(i1.base().size(), i2.base().size());
   return std::conj(BLAS::zdotc(i1.base().size(), i1.base().blas_base(), i1.base().stride(),
				i2.base().blas_base(), i2.base().stride()));
}
#endif // defined(HAVE_ZDOTC)

#if defined(HAVE_ZDOTU)

// no conjugation, ordinary multiplication
template <typename S1, typename S2>
inline
std::complex<double> 
iter_inner_prod(VecPtrIterator<std::complex<double> const, S1> i1, 
		VecPtrIterator<std::complex<double> const, S2> i2,
		Multiplication<std::complex<double>, std::complex<double> >)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.size());
   return BLAS::zdotu(i1.size(), i1.blas_base(), i1.stride(), 
		      i2.blas_base(), i2.stride());
}

// first vector conjugated, inner product
template <typename S1, typename S2>
inline
std::complex<double> 
iter_inner_prod(VectorTransformIterator<VecPtrIterator<
		std::complex<double> const, S1>, Conj<std::complex<double> > > i1, 
		VecPtrIterator<std::complex<double> const, S2> i2,
		InnerProd<std::complex<double>, std::complex<double> >)
{
   DEBUG_PRECONDITION_EQUAL(i1.base().size(), i2.size());
   return BLAS::zdotu(i1.base().size(), i1.base().blas_base(), i1.base().stride(), 
		      i2.blas_base(), i2.stride());
}

// second vector conjugated, inner product
template <typename S1, typename S2>
inline
std::complex<double> 
iter_inner_prod(VecPtrIterator<std::complex<double> const, S1> i1,
		VectorTransformIterator<VecPtrIterator<
		std::complex<double> const, S2>, Conj<std::complex<double> > > i2, 
		InnerProd<std::complex<double>, std::complex<double> >)
{
   DEBUG_PRECONDITION_EQUAL(i1.size(), i2.base().size());
   return std::conj(BLAS::zdotu(i1.size(), i1.blas_base(), i1.stride(), 
				i2.base().blas_base(), i2.base().stride()));
}

// both vectors conjugated, ordinary multiplication
template <typename S1, typename S2>
inline
std::complex<double> 
iter_inner_prod(VectorTransformIterator<VecPtrIterator<
		std::complex<double> const, S1>, Conj<std::complex<double> > > i1,
		VectorTransformIterator<VecPtrIterator<
		std::complex<double> const, S2>, Conj<std::complex<double> > > i2, 
		Multiplication<std::complex<double>, 
		std::complex<double> >)
{
   DEBUG_PRECONDITION_EQUAL(i1.base().size(), i2.base().size());
   return std::conj(BLAS::zdotu(i1.base().size(), i1.base().blas_base(), i1.base().stride(), 
				i2.base().blas_base(), i2.base().stride()));
}
#endif // defined(HAVE_ZDOTU)

// iter_norm_1

template <typename S1>
inline
double iter_norm_1(VecPtrIterator<double const, S1> i1)
{
   return BLAS::dasum(i1.size(), i1.blas_base(), i1.stride());
}

// complex: BLAS::zamax uses a different defintion of absolute value
// that is incompatible with our definition of norm_1

// iter_norm_2

template <typename S1>
inline
double iter_norm_2(VecPtrIterator<double const, S1> i1)
{
   return BLAS::dnrm2(i1.size(), i1.blas_base(), i1.stride());
}

template <typename S1>
inline
double iter_norm_2(VecPtrIterator<std::complex<double> const, S1> i1)
{
   return BLAS::dznrm2(i1.size(), i1.blas_base(), i1.stride());
}

// iter_norm_2_sq

template <typename S1>
inline
double iter_norm_2_sq(VecPtrIterator<double const, S1> i1)
{
   double x = BLAS::dnrm2(i1.size(), i1.blas_base(), i1.stride());
   return x*x;
}

template <typename S1>
inline
double iter_norm_2_sq(VecPtrIterator<std::complex<double> const, S1> i1)
{
   double x = BLAS::dznrm2(i1.size(), i1.blas_base(), i1.stride());
   return x*x;
}

// iter_norm_inf

template <typename S1>
inline
double iter_norm_inf(VecPtrIterator<double const, S1> i1)
{
   return std::fabs(i1[BLAS::idamax(i1.size(), i1.blas_base(), i1.stride())]);
}

// complex: BLAS::izamax uses a different defintion of absolute value
// that is incompatible with our definition of norm_inf

} // namespace LinearAlgebra

#endif
