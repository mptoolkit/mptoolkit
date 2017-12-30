// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/cublas.h
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

//
// we possibly need to make the handle reference counted and copyable,
// so that we can pass it by value to cublas calls.  Hence every parameter
// can be passed by value and we don't have any problems with temporaries going
// out of scope if we call functions asynchronously.
//

#if !defined(MPTOOLKIT_CUDA_CUBLAS_H)
#define MPTOOLKIT_CUDA_CUBLAS_H

#include "cuda.h"
#include "gpu_buffer.h"
#include <list>
#include <mutex>
#include <cublas_v2.h>
#include "cub.h"
#include "blas/functors.h"

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/iterator/permutation_iterator.h>

#include <iostream>

namespace cublas
{

// returns the cublas version number
int version();

// helper functions to get a name and error description from a cublasStatus_t
char const* GetErrorName(cublasStatus_t error);

char const* GetErrorString(cublasStatus_t error);

// wrapper for a cublasStatus_t, which can also function as an exception object
class error : public std::runtime_error
{
   public:
      error() = delete;
      error(cublasStatus_t Err) : std::runtime_error(std::string("cuBLAS error: ")
                                                     + cublas::GetErrorString(Err)), err_(Err)
      {std::cerr << "cuBLAS Error " << int(Err) << ' ' << cublas::GetErrorString(Err) << '\n';}

      cublasStatus_t code() const { return err_; }
      operator cublasStatus_t() const { return err_; }

      char const* name() const { return cublas::GetErrorName(err_); }
      char const* string() const { return cublas::GetErrorString(err_); }

   private:
      cublasStatus_t err_;
};

inline
void check_error(cublasStatus_t s)
{
   if (s != CUBLAS_STATUS_SUCCESS)
   {
      throw error(s);
   }
}

// cublas handle.  Moveable, but not copyable.
class handle
{
   public:
      handle() : h_(nullptr) {}
      handle(handle&& other) : h_(other.h_) { other.h_ = nullptr; }
      handle(handle const&) = delete;
      handle& operator=(handle&& other) { std::swap(h_, other.h_); return *this; }
      handle& operator=(handle const&) = delete;
      ~handle() { if (h_) cublasDestroy(h_); }

      cublasHandle_t raw_handle() const { return h_; }

      static handle create() { cublasHandle_t h; cublasCreate(&h); return handle(h); }

      void destroy() { cublasDestroy(h_); h_ = nullptr; }

      // set the stream associated with the handle
      void set_stream(cuda::stream const& s)
      {
         check_error(cublasSetStream(h_, s.raw_stream()));
      }

      void set_pointer_mode(cublasPointerMode_t m)
      {
         check_error(cublasSetPointerMode(h_, m));
      }

   private:
      handle(cublasHandle_t h) : h_(h) {}

      cublasHandle_t h_;
};

// returns the thread-local handle
handle& get_handle();

// utility functions

inline
cublasOperation_t cublas_trans(char c)
{
   switch (c)
   {
   case 'N' : return CUBLAS_OP_N;
   case 'T' : return CUBLAS_OP_T;
   case 'C' : return CUBLAS_OP_C;
   default : throw std::runtime_error("Unsupported TRANS in cublas");
   }
}

template <typename T>
inline
void
setMatrixAsync(int N, int M, T const* A, int lda, cuda::gpu_ptr<T> B, int ldb)
{
   check_error(cublasSetMatrixAsync(N, M, sizeof(T), A, lda, B.device_ptr(), ldb,
                                    B.get_stream().raw_stream()));
}

template <typename T>
inline
void
setVectorAsync(int N, T const* A, int stridea, cuda::gpu_ptr<T> B, int strideb)
{
   check_error(cublasSetVectorAsync(N, sizeof(T), A, stridea, B.device_ptr(), strideb,
                                    B.get_stream().raw_stream()));
}

// blocking
template <typename T>
inline
void
setVector(int N, T const* A, int stridea, cuda::gpu_ptr<T> B, int strideb)
{
   check_error(cublasSetVector(N, sizeof(T), A, stridea, B.device_ptr(), strideb));
}

// helper forwarding functions

inline
void dotc(cublasHandle_t handle, int n,
	  float const* x, int incx,
	  float const* y, int incy,
	  float* result)
{
   cublas::check_error(cublasSdot(handle, n, x, incx, y, incy, result));
}

inline
void dotc(cublasHandle_t handle, int n,
	  double const* x, int incx,
	  double const* y, int incy,
	  double* result)
{
   cublas::check_error(cublasDdot(handle, n, x, incx, y, incy, result));
}

inline
void dotc(cublasHandle_t handle, int n,
	  std::complex<float> const* x, int incx,
	  std::complex<float> const* y, int incy,
	  std::complex<float>* result)
{
   cublas::check_error(cublasCdotc(handle, n,
				   reinterpret_cast<cuComplex const*>(x), incx,
				   reinterpret_cast<cuComplex const*>(y), incy,
				   reinterpret_cast<cuComplex*>(result)));
}

inline
void dotc(cublasHandle_t handle, int n,
	  std::complex<double> const* x, int incx,
	  std::complex<double> const* y, int incy,
	  std::complex<double>* result)
{
   cublas::check_error(cublasZdotc(handle, n,
				   reinterpret_cast<cuDoubleComplex const*>(x), incx,
				   reinterpret_cast<cuDoubleComplex const*>(y), incy,
				   reinterpret_cast<cuDoubleComplex*>(result)));
}

} // namespace cublas

// BLAS functions must go in namespace cuda so they are found during ADL

namespace cuda
{

// BLAS level 1

inline
void
vector_copy(int n, cuda::const_gpu_ptr<double> x, int incx, cuda::gpu_ptr<double> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(y.get_stream());
   y.wait_for(x);
   cublas::check_error(cublasDcopy(H.raw_handle(), n, x.device_ptr(), incx, y.device_ptr(), incy));
   x.wait_for(y);
}

inline
void
vector_copy(int n, cuda::const_gpu_ptr<std::complex<double>> x, int incx, cuda::gpu_ptr<std::complex<double>> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(y.get_stream());
   y.wait_for(x);
   cublas::check_error(cublasZcopy(H.raw_handle(), n,
				   reinterpret_cast<cuDoubleComplex const*>(x.device_ptr()), incx,
				   reinterpret_cast<cuDoubleComplex*>(y.device_ptr()), incy));
   x.wait_for(y);
}

inline
void
vector_scale(int n, double alpha, cuda::gpu_ptr<double> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   cublas::check_error(cublasDscal(H.raw_handle(), n, &alpha, y.device_ptr(), incy));
}

inline
void vector_clear(int N, cuda::gpu_ptr<double> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   if (incy == 1)
   {
      cuda::memset_async(y.get_stream(), y.device_ptr(), 0, N*sizeof(double));
   }
   else
   {
      vector_scale(N, 0, y, incy);
   }
}

inline
void
vector_copy_scaled(int n, double alpha, cuda::const_gpu_ptr<double> x, int incx, cuda::gpu_ptr<double> y, int incy)
{
   vector_copy(n, x, incx, y, incy);
   vector_scale(n, alpha, y, incy);
}

inline
void
vector_add_scaled(int n, double alpha, cuda::const_gpu_ptr<double> x, int incx,
                  cuda::gpu_ptr<double> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   y.wait_for(x);
   cublas::check_error(cublasDaxpy(H.raw_handle(), n, &alpha, x.device_ptr(), incx, y.device_ptr(), incy));
   x.wait_for(y);
}

inline
void
vector_add(int n, cuda::const_gpu_ptr<double> x, int incx,
           cuda::gpu_ptr<double> y, int incy)
{
   vector_add_scaled(n, 1.0, x, incx, y, incy);
}

template <typename T>
inline
std::enable_if_t<is_cuda_floating_point_v<T>, void>
vector_inner_prod(int n,
		  cuda::const_gpu_ptr<T> x, int incx,
                  cuda::const_gpu_ptr<T> y, int incy,
		  cuda::gpu_ref<T>& r)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(r.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_DEVICE);
   r.wait(x.sync());
   r.wait(y.sync());
   cublas::dotc(H.raw_handle(), n, x.device_ptr(), incx, y.device_ptr(), incy, r.device_ptr());
   x.wait(r.sync());
   y.wait(r.sync());
}

template <typename T>
inline
std::enable_if_t<is_cuda_floating_point_v<T>, void>
vector_inner_prod_nested(int n,
			 cuda::const_gpu_ptr<T> x, int incx,
			 cuda::const_gpu_ptr<T> y, int incy,
			 cuda::gpu_ref<T>& r, blas::InnerProd)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(r.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_DEVICE);
   r.wait(x.sync());
   r.wait(y.sync());
   cublas::dotc(H.raw_handle(), n, x.device_ptr(), incx, y.device_ptr(), incy, r.device_ptr());
   x.wait(r.sync());
   y.wait(r.sync());
}

template <typename T>
inline
std::enable_if<is_cuda_floating_point_v<T>, void>
vector_inner_prod_nested(int n,
			 cuda::const_gpu_ptr<T> x, int incx,
			 cuda::const_gpu_ptr<T> y, int incy,
			 cuda::gpu_ref<T>& r, blas::InnerProdNested)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(r.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_DEVICE);
   r.wait(x.sync());
   r.wait(y.sync());
   cublas::dotc(H.raw_handle(), n, x.device_ptr(), incx, y.device_ptr(), incy, r.device_ptr());
   x.wait(r.sync());
   y.wait(r.sync());
}

// BLAS level 2

inline
void
gemv(char Atrans, int M, int N, double alpha,
     cuda::const_gpu_ptr<double> A, int lda, cuda::const_gpu_ptr<double> x, int incx,
     double beta, cuda::gpu_ptr<double> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   y.wait_for(A);
   y.wait_for(x);
   cublas::check_error(cublasDgemv(H.raw_handle(), cublas::cublas_trans(Atrans), M, N,
                           &alpha, A.device_ptr(), lda, x.device_ptr(), incx,
                           &beta, y.device_ptr(), incy));
   A.wait_for(y);
   x.wait_for(y);
}

// BLAS level 2 extensions

// in-place two matrix version implements
// C = alpha*A + beta*C
inline
void
geam(char Atrans, int M, int N,
     double alpha, cuda::const_gpu_ptr<double> A, int lda,
     double beta, cuda::gpu_ptr<double> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   cublas::check_error(cublasDgeam(H.raw_handle(), cublas::cublas_trans(Atrans), CUBLAS_OP_N, M, N,
                           &alpha, A.device_ptr(), lda,
                           &beta, C.device_ptr(), ldc,
                           C.device_ptr(), ldc));
   A.wait_for(C);
}

inline
void
geam(char Atrans, int M, int N,
     std::complex<double> alpha, cuda::const_gpu_ptr<std::complex<double>> A, int lda,
     std::complex<double> beta, cuda::gpu_ptr<std::complex<double>> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   cublas::check_error(cublasZgeam(H.raw_handle(), cublas::cublas_trans(Atrans), CUBLAS_OP_N, M, N,
				   reinterpret_cast<cuDoubleComplex const*>(&alpha),
				   reinterpret_cast<cuDoubleComplex const*>(A.device_ptr()), lda,
				   reinterpret_cast<cuDoubleComplex const*>(&beta),
				   reinterpret_cast<cuDoubleComplex const*>(C.device_ptr()), ldc,
				   reinterpret_cast<cuDoubleComplex*>(C.device_ptr()), ldc));
   A.wait_for(C);
}

// out-of-place version implements
// C = alpha*A + beta*B
inline
void
geam(char Atrans, char Btrans, int M, int N,
     double alpha, cuda::const_gpu_ptr<double> A, int lda,
     double beta,  cuda::const_gpu_ptr<double> B, int ldb,
     cuda::gpu_ptr<double> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   cublas::check_error(cublasDgeam(H.raw_handle(), cublas::cublas_trans(Atrans), cublas::cublas_trans(Btrans), M, N,
				   &alpha, A.device_ptr(), lda,
				   &beta, B.device_ptr(), ldb,
				   C.device_ptr(), ldc));
   A.wait_for(C);
   B.wait_for(C);
}

inline
void
geam(char Atrans, char Btrans, int M, int N,
     std::complex<double> alpha, cuda::const_gpu_ptr<std::complex<double>> A, int lda,
     std::complex<double> beta,  cuda::const_gpu_ptr<std::complex<double>> B, int ldb,
     cuda::gpu_ptr<std::complex<double>> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   cublas::check_error(cublasZgeam(H.raw_handle(), cublas::cublas_trans(Atrans), cublas::cublas_trans(Btrans), M, N,
				   reinterpret_cast<cuDoubleComplex const*>(&alpha),
				   reinterpret_cast<cuDoubleComplex const*>(A.device_ptr()), lda,
				   reinterpret_cast<cuDoubleComplex const*>(&beta),
				   reinterpret_cast<cuDoubleComplex const*>(B.device_ptr()), ldb,
				   reinterpret_cast<cuDoubleComplex*>(C.device_ptr()), ldc));
   A.wait_for(C);
   B.wait_for(C);
}

template <typename T>
inline
void
matrix_copy(char Atrans, int M, int N, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam(Atrans, M, N, 1.0, A, lda, 0.0, B, ldb);
}

template <typename T>
inline
void
matrix_copy_scaled(char Atrans, int M, int N, T alpha, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam(Atrans, M, N, alpha, A, lda, 0.0, B, ldb);
}

template <typename T>
inline
void
matrix_add(char Atrans, int M, int N, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam(Atrans, M, N, 1.0, A, lda, 1.0, B, ldb);
}

template <typename T>
inline
void
matrix_add_scaled(char Atrans, int M, int N, T alpha, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam(Atrans, M, N, alpha, A, lda, 1.0, B, ldb);
}

inline
void
matrix_scale(int M, int N, double alpha, cuda::gpu_ptr<double> A, int lda)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(A.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   double beta = 0.0;
   cublas::check_error(cublasDgeam(H.raw_handle(), CUBLAS_OP_N, CUBLAS_OP_N, M, N,
                           &alpha, A.device_ptr(), lda,
                           &beta, nullptr, 1,
                           A.device_ptr(), lda));
}

inline
void
matrix_scale(int M, int N, std::complex<double> alpha, cuda::gpu_ptr<std::complex<double>> A, int lda)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(A.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   double beta = 0.0;
   cublas::check_error(cublasZgeam(H.raw_handle(), CUBLAS_OP_N, CUBLAS_OP_N, M, N,
				   reinterpret_cast<cuDoubleComplex const*>(&alpha),
				   reinterpret_cast<cuDoubleComplex const*>(A.device_ptr()), lda,
				   reinterpret_cast<cuDoubleComplex const*>(&beta), nullptr, 1,
				   reinterpret_cast<cuDoubleComplex*>(A.device_ptr()), lda));
}

inline
void matrix_clear(int M, int N, cuda::gpu_ptr<double> A, int lda)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(A.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   double beta = 0.0;
   cublas::check_error(cublasDgeam(H.raw_handle(), CUBLAS_OP_N, CUBLAS_OP_N, M, N,
                                   &beta, A.device_ptr(), lda,
                                   &beta, A.device_ptr(), lda,
                                   A.device_ptr(), lda));
}

template <typename T>
std::enable_if_t<is_cuda_floating_point_v<T>, void>
matrix_inner_prod(char Atrans, char Btrans, int M, int N,
		  cuda::const_gpu_ptr<T> A, int ldA,
		  cuda::const_gpu_ptr<T> B, int ldB,
		  cuda::gpu_ref<T>& r);

template <typename T, typename Nested>
inline
std::enable_if_t<is_cuda_floating_point_v<T>, void>
matrix_inner_prod_nested(char Atrans, char Btrans, int M, int N,
			 cuda::const_gpu_ptr<T> A, int ldA,
			 cuda::const_gpu_ptr<T> B, int ldB,
			 cuda::gpu_ref<T>& r, Nested&&)
{
   matrix_inner_prod(Atrans, Btrans, M, N, A, ldA, B, ldB, r);
}

template <typename T>
std::enable_if_t<is_cuda_floating_point_v<T>, void>
matrix_add_inner_prod(char Atrans, char Btrans, int M, int N,
		      cuda::const_gpu_ptr<T> A, int ldA,
		      cuda::const_gpu_ptr<T> B, int ldB,
		      cuda::gpu_ref<T>& r);

template <typename T, typename Nested>
inline
std::enable_if_t<is_cuda_floating_point_v<T>, void>
matrix_add_inner_prod_nested(char Atrans, char Btrans, int M, int N,
			     cuda::const_gpu_ptr<T> A, int ldA,
			     cuda::const_gpu_ptr<T> B, int ldB,
			     cuda::gpu_ref<T>& r, Nested&&)
{
   matrix_add_inner_prod(Atrans, Btrans, M, N, A, ldA, B, ldB, r);
}

// TODO: this implementation is not ideal
inline
void
matrix_norm_frob_sq(int M, int N, const_gpu_ptr<double> A, int lda, gpu_ref<double>& Result)
{
   matrix_inner_prod('N', 'N', M, N, A, lda, A, lda, Result);
}

inline
void
matrix_norm_frob_sq(int M, int N, const_gpu_ptr<std::complex<double>> A, int lda, gpu_ref<double>& Result)
{
   cuda::gpu_ref<std::complex<double>> R;
   matrix_inner_prod('N', 'N', M, N, A, lda, A, lda, R);
   Result.set_wait(get_wait(R).real());
}

// BLAS level 3

inline
void
gemm(char Atrans, char Btrans, int M, int N, int K, double alpha,
     cuda::const_gpu_ptr<double> A, int lda, cuda::const_gpu_ptr<double> B, int ldb,
     double beta, cuda::gpu_ptr<double> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   cublas::check_error(cublasDgemm(H.raw_handle(), cublas::cublas_trans(Atrans), cublas::cublas_trans(Btrans), M, N, K,
                                   &alpha, A.device_ptr(), lda, B.device_ptr(), ldb,
                                   &beta, C.device_ptr(), ldc));
   A.wait_for(C);
   B.wait_for(C);
}

inline
void
gemm(char Atrans, char Btrans, int M, int N, int K, std::complex<double> alpha,
     cuda::const_gpu_ptr<std::complex<double>> A, int lda, cuda::const_gpu_ptr<std::complex<double>> B, int ldb,
     std::complex<double> beta, cuda::gpu_ptr<std::complex<double>> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   cublas::check_error(cublasZgemm(H.raw_handle(), cublas::cublas_trans(Atrans), cublas::cublas_trans(Btrans), M, N, K,
                                   reinterpret_cast<cuDoubleComplex const*>(&alpha),
                                   reinterpret_cast<cuDoubleComplex const*>(A.device_ptr()), lda,
                                   reinterpret_cast<cuDoubleComplex const*>(B.device_ptr()), ldb,
                                   reinterpret_cast<cuDoubleComplex const*>(&beta),
                                   reinterpret_cast<cuDoubleComplex*>(C.device_ptr()), ldc));
   A.wait_for(C);
   B.wait_for(C);
}

} // namespace cuda

#include "cublas.icc"

#endif
