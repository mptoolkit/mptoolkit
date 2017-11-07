// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/cuda.h
//
// Copyright (C) 2017 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#if !defined(MPTOOLKIT_ASYNC_ASYNC_MATRIX_H)
#define MPTOOLKIT_ASYNC_ASYNC_MATRIX_H

namespace async
{

template <typename T>
class async_matrix : public blas::BlasMatrix<typename T::value_type, async_matrix<T>, async_tag>
{
   public:
      using base_type          = T;
      using value_type         = typename base_type::value_type;
      using storage_type       = typename base_type::storage_type;
      using const_storage_type = typename base_type::const_storage_type;

      async_matrix(int Rows_, int Cols_) : Base(new base_type(Rows_, Cols_)), Queue(new async::queue())
      {
      }

      async_matrix(async_matrix&& Other) = default;

      template <typename U>
      async_matrix(blas::MatrixRef<value_type, U, async_tag> const& E)
         : async_matrix(E.rows(), E.cols())
      {
         assign(*this, E.as_derived());
      }

      ~async_matrix()
      {
         // schedule the deallocation of Base via the queue
         // self-destruct the queue
      }

      void wait_for(async::queue const& x)
      {
         Queue.wait(x.record());
      }

   private:
      base_type* Base;
      async::queue* Queue;
};

// BLAS-like functions

// **NOTE**:
// This sketch is erroneous, because the nested gemv() call here may involve temporary
// objects which might be destroyed prior to the task getting executed.
// The tasks must involve only raw memory pointers, no proxies.  Therefore we need
// separate async_gpu_matrix and async_matrix classes.

template <typename T, typename U, typename V, typename W>
inline
void
gemv(T alpha, blas::BlasMatrix<T, U, async_tag> const& A,
     blas::BlasVector<T, V, async_tag> const& x, T beta,
     async_vector<W>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(y.size(), A.rows());
   y.wait_for(A.queue());
   y.wait_for(x.queue());
   y.queue().add_task([&]{gemv(alpha, A.base(), x.base(), beta, y.base()); return true;});
   A.wait_for(y.queue());
   x.wait_for(y.queue());
}

// version for an async_gpu_tag
template <typename T, typename U, typename V>
inline
void
gemv(T alpha, blas::BlasMatrix<T, U, async_gpu_tag> const& A,
     blas::BlasVector<T, V, async_gpu_tag> const& x, T beta,
     async_gpu_vector<T>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(y.size(), A.rows());
   y.wait_for(A.queue());
   y.wait_for(x.queue());
   y.queue().add_task(std::bind(cublas::gemv, A.trans(), A.rows(), A.cols(),
                                alpha, A.storage(),
                                A.leading_dimension(), x.storage(), x.stride(),
                                beta, y.storage(), y.stride()));
   A.wait_for(y.queue());
   x.wait_for(y.queue());
}

// another attempt at a generic version; assuming cublas calls have a
// compatible signature to standard blas.
// This is a bit tedious because we need to pass lots of parameters to the lambda function
// by value with automatic variables.

// In async mode we are guaranteed that the objects live longer than this function, because we
// put that intelligence into the destructors of the buffers.  So we don't need reference counting,
// but do need careful handling of temporaries.
template <typename T, typename U, typename V, typename W>
inline
void
gemv(T alpha, blas::BlasMatrix<T, U, async_tag> const& A,
     blas::BlasVector<T, V, async_tag> const& x, T beta,
     async_vector<W>& y)
{
   DEBUG_CHECK_EQUAL(A.cols(), x.size());
   DEBUG_CHECK_EQUAL(y.size(), A.rows());
   y.wait_for(A.queue());
   y.wait_for(x.queue());
   auto Atrans = A.trans();
   auto Arows = A.rows();
   auto Acols = A.cols();
   auto Astorage = A.storage();
   auto Aleading_dimension = A.leading_dimension();
   auto xstorage = x.storage();
   auto xstride = x.stride();
   auto ystorage = y.storage();
   auto ystride = y.stride();
   y.queue().add_task([=](){ gemv(Atrans, Arows, Acols, alpha, Astorage, Aleading_dimension,
                                  xstorage, xstride, beta, ystorage, ystride);});
   A.wait_for(y.queue());
   x.wait_for(y.queue());
}

} // namespace async
