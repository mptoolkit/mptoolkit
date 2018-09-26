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
#include <iostream>

#define USE_GEMM3M

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
   TRACE_CUDA("cublasSetMatrixAsync")(N)(M)(sizeof(T))(A)(lda)(B)(ldb)(B.get_stream().raw_stream());
   check_error(cublasSetMatrixAsync(N, M, sizeof(T), A, lda, B.device_ptr(), ldb,
                                    B.get_stream().raw_stream()));
}

template <typename T>
inline
void
setVectorAsync(int N, T const* A, int stridea, cuda::gpu_ptr<T> B, int strideb)
{
   TRACE_CUDA("cublasSetVectorAsync")(N)(sizeof(T))(A)(stridea)(B.device_pr())(strideb)(B.get_stream().raw_stream());
   check_error(cublasSetVectorAsync(N, sizeof(T), A, stridea, B.device_ptr(), strideb,
                                    B.get_stream().raw_stream()));
}

// blocking
template <typename T>
inline
void
setVector(int N, T const* A, int stridea, cuda::gpu_ptr<T> B, int strideb)
{
   TRACE_CUDA("cublasSetVector")(N)(A)(stridea)(B)(strideb);
   check_error(cublasSetVector(N, sizeof(T), A, stridea, B.device_ptr(), strideb));
}

// helper forwarding functions

inline
void dotc(cublasHandle_t handle, int n,
	  float const* x, int incx,
	  float const* y, int incy,
	  float* result)
{
   TRACE_CUDA("cublasSdot")(n)(x)(incx)(y)(incy)(result);
   cublas::check_error(cublasSdot(handle, n, x, incx, y, incy, result));
}

inline
void dotc(cublasHandle_t handle, int n,
	  double const* x, int incx,
	  double const* y, int incy,
	  double* result)
{
   TRACE_CUDA("cublasDdot")(n)(x)(incx)(y)(incy)(result);
   cublas::check_error(cublasDdot(handle, n, x, incx, y, incy, result));
}

inline
void dotc(cublasHandle_t handle, int n,
	  std::complex<float> const* x, int incx,
	  std::complex<float> const* y, int incy,
	  std::complex<float>* result)
{
   TRACE_CUDA("cublasCdotc")(n)(x)(incx)(y)(incy)(result);
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
   TRACE_CUDA("cublasZdotc")(n)(x)(incx)(y)(incy)(result);
   cublas::check_error(cublasZdotc(handle, n,
				   reinterpret_cast<cuDoubleComplex const*>(x), incx,
				   reinterpret_cast<cuDoubleComplex const*>(y), incy,
				   reinterpret_cast<cuDoubleComplex*>(result)));
}

inline
void copy(cublasHandle_t handle, int n,
          double const* x, int incx,
          double* y, int incy)
{
   TRACE_CUDA("cublasDcopy")(n)(x)(incx)(y)(incy);
   cublas::check_error(cublasDcopy(handle, n, x, incx, y, incy));
}

inline
void copy(cublasHandle_t handle, int n,
          std::complex<double> const* x, int incx,
          std::complex<double>* y, int incy)
{
   TRACE_CUDA("cublasZcopy")(n)(x)(incx)(y)(incy);
   cublas::check_error(cublasZcopy(handle, n,
				   reinterpret_cast<cuDoubleComplex const*>(x), incx,
				   reinterpret_cast<cuDoubleComplex*>(y), incy));
}

inline
void copy(cublasHandle_t handle, int n,
          double const* x, int incx,
          std::complex<double>* y, int incy)
{
   double alpha = 0.0;
   cublasPointerMode_t OldMode;
   cublas::check_error(cublasGetPointerMode(handle, &OldMode));
   if (OldMode != CUBLAS_POINTER_MODE_HOST)
      cublas::check_error(cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_HOST));
   TRACE_CUDA("cublasDscal")(n)(&alpha)(reinterpret_cast<double*>(y)+1)(incy*2);
   cublas::check_error(cublasDscal(handle, n, &alpha,
                                   reinterpret_cast<double*>(y)+1, incy*2));
   if (OldMode != CUBLAS_POINTER_MODE_HOST)
      cublas::check_error(cublasSetPointerMode(handle, OldMode));
   TRACE_CUDA("cublasDcopy")(n)(x)(incx)(y)(incy*2);
   cublas::check_error(cublasDcopy(handle, n, x, incx,
				   reinterpret_cast<double*>(y), incy*2));
}

 inline
void scale(cublasHandle_t handle, int n,
           double const* alpha, double* y, int incy)
{
   TRACE_CUDA("cublasDscal")(n)(alpha)(y)(incy);
   cublas::check_error(cublasDscal(handle, n, alpha, y, incy));
}

inline
void scale(cublasHandle_t handle, int n,
           std::complex<double> const* alpha, std::complex<double>* y, int incy)
{
   TRACE_CUDA("cublasZscal")(n)(alpha)(y)(incy);
   cublas::check_error(cublasZscal(handle, n, reinterpret_cast<cuDoubleComplex const*>(alpha),
                                   reinterpret_cast<cuDoubleComplex*>(y), incy));
}

inline
void axpy(cublasHandle_t handle, int n,
          double* alpha, double const* x, int incx,
          double* y, int incy)
{
   TRACE_CUDA("cublasDaxpy")(n)(alpha)(x)(incx)(y)(incy);
   cublas::check_error(cublasDaxpy(handle, n, alpha, x, incx, y, incy));
}

inline
void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
          int m, int n, double const* alpha, double const* A, int lda,
          double const* beta, double const* B, int ldb, double* C, int ldc)
{
   TRACE_CUDA("cublasDgeam")(transa)(transb)(m)(n)(alpha)(A)(lda)(beta)(B)(ldb)(C)(ldc);
   cublas::check_error(cublasDgeam(handle, transa, transb, m, n,
                                   alpha, A, lda,
                                   beta, B, ldb,
                                   C, ldc));
}

inline
void geam(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
          int m, int n, std::complex<double> const* alpha, std::complex<double> const* A, int lda,
          std::complex<double> const* beta, std::complex<double> const* B, int ldb,
          std::complex<double>* C, int ldc)
{
   TRACE_CUDA("cublasZgeam")(transa)(transb)(m)(n)(alpha)(A)(lda)(beta)(B)(ldb)(C)(ldc);
   cublas::check_error(cublasZgeam(handle, transa, transb, m, n,
				   reinterpret_cast<cuDoubleComplex const*>(alpha),
				   reinterpret_cast<cuDoubleComplex const*>(A), lda,
				   reinterpret_cast<cuDoubleComplex const*>(beta),
				   reinterpret_cast<cuDoubleComplex const*>(B), ldb,
				   reinterpret_cast<cuDoubleComplex*>(C), ldc));
}

inline
void gemm(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
          int m, int n, int k, double const* alpha, double const* A, int lda,
          double const* B, int ldb,
          double const* beta, double* C, int ldc)
{
   TRACE_CUDA("cublasDgemm")(transa)(transb)(m)(n)(k)(alpha)(A)(lda)(B)(ldb)(beta)(C)(ldc);
   cublas::check_error(cublasDgemm(handle, transa, transb, m, n, k,
                                   alpha, A, lda, B, ldb,
                                   beta, C, ldc));
}

inline
void gemm(cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
          int m, int n, int k, std::complex<double> const* alpha,
          std::complex<double> const* A, int lda,
          std::complex<double> const* B, int ldb,
          std::complex<double> const* beta, std::complex<double>* C, int ldc)
{
   TRACE_CUDA("cublasZgemm")(transa)(transb)(m)(n)(k)(alpha)(A)(lda)(B)(ldb)(beta)(C)(ldc);
#if defined(USE_GEMM3M)
   if (cuda::get_device_properties().compute_capability_major() >= 5)
   {
      cublas::check_error(cublasZgemm3m(handle, transa, transb, m, n, k,
                                        reinterpret_cast<cuDoubleComplex const*>(alpha),
                                        reinterpret_cast<cuDoubleComplex const*>(A), lda,
                                        reinterpret_cast<cuDoubleComplex const*>(B), ldb,
                                        reinterpret_cast<cuDoubleComplex const*>(beta),
                                        reinterpret_cast<cuDoubleComplex*>(C), ldc));
   }
   else
#endif
   {
      cublas::check_error(cublasZgemm(handle, transa, transb, m, n, k,
                                      reinterpret_cast<cuDoubleComplex const*>(alpha),
                                      reinterpret_cast<cuDoubleComplex const*>(A), lda,
                                      reinterpret_cast<cuDoubleComplex const*>(B), ldb,
                                      reinterpret_cast<cuDoubleComplex const*>(beta),
                                      reinterpret_cast<cuDoubleComplex*>(C), ldc));
   }
}

inline
void dgmm(cublasHandle_t handle, cublasSideMode_t mode,
          int m, int n,
          double const* A, int lda,
          double const* x, int incx,
          double* C, int ldc)
{
   TRACE_CUDA("cublasDgmm")(mode)(m)(n)(A)(lda)(x)(incx)(C)(ldc);
   cublas::check_error(cublasDdgmm(handle, mode, m, n, A, lda, x, incx, C, ldc));
}

inline
void dgmm(cublasHandle_t handle, cublasSideMode_t mode,
          int m, int n,
          std::complex<double> const* A, int lda,
          std::complex<double> const* x, int incx,
          std::complex<double>* C, int ldc)
{
   TRACE_CUDA("cublasZgmm")(mode)(m)(n)(A)(lda)(x)(incx)(C)(ldc);
   cublas::check_error(cublasZdgmm(handle, mode, m, n,
                                   reinterpret_cast<cuDoubleComplex const*>(A), lda,
                                   reinterpret_cast<cuDoubleComplex const*>(x), incx,
                                   reinterpret_cast<cuDoubleComplex*>(C), ldc));
}

inline
void gemv(cublasHandle_t handle, cublasOperation_t transa,
          int m, int n, double const* alpha, double const* A, int lda,
          double const* x, int incx, double const* beta,
          double* y, int incy)
{
   TRACE_CUDA("cublasDgemv")(transa)(m)(n)(alpha)(A)(lda)(x)(incx)(beta)(y)(incy);
   cublas::check_error(cublasDgemv(handle, transa, m, n, alpha, A, lda,
                                   x, incx, beta, y, incy));
}
} // namespace cublas

// BLAS functions must go in namespace cuda so they are found during ADL

namespace cuda
{

// functions acting on scalars (BLAS 'level 0')

template <typename T>
inline
void
clear(cuda::gpu_ref<T>& r)
{
   r.set_wait(T{});
}

template <typename T>
inline
void
clear(cuda::gpu_ref<T>&& r)
{
   r.set_wait(T{});
}

// BLAS level 1

// this allows for mixed real/complex
template <typename T, typename U>
inline
void
vector_copy(int n, cuda::const_gpu_ptr<T> x, int incx, cuda::gpu_ptr<U> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   y.wait_for(x);
   cublas::copy(H.raw_handle(), n, x.device_ptr(), incx, y.device_ptr(), incy);
   x.wait_for(y);
}

template <typename T, typename U>
inline
void
vector_copy_conj(int n, cuda::const_gpu_ptr<T> x, int incx, cuda::gpu_ptr<U> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   y.wait_for(x);
   cublas::copy(H.raw_handle(), n, x.device_ptr(), incx, y.device_ptr(), incy);
   x.wait_for(y);
   double alpha = -1.0;
   cublas::scale(H.raw_handle(), n, &alpha, reinterpret_cast<double*>(y.device_ptr())+1, incy*2);
}

inline
void
vector_scale(int n, double alpha, cuda::gpu_ptr<double> y, int incy)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(y.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   cublas::scale(H.raw_handle(), n, &alpha, y.device_ptr(), incy);
}

inline
void vector_clear(int N, cuda::gpu_ptr<double> y, int incy)
{
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
   cublas::axpy(H.raw_handle(), n, &alpha, x.device_ptr(), incx, y.device_ptr(), incy);
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
std::enable_if_t<cuda::is_cuda_floating_point_v<T>, void>
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
std::enable_if_t<cuda::is_cuda_floating_point_v<T>, void>
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
std::enable_if<cuda::is_cuda_floating_point_v<T>, void>
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

inline
void
vector_norm_frob_sq(int n, cuda::const_gpu_ptr<double>x, int incx, cuda::gpu_ref<double>& r)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(r.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_DEVICE);
   r.wait(x.sync());
   cublas::dotc(H.raw_handle(), n, x.device_ptr(), incx, x.device_ptr(), incx, r.device_ptr());
   x.wait(r.sync());
}

inline
void vector_conj(int n, cuda::gpu_ptr<std::complex<double>> x, int incx)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(x.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   double alpha = -1.0;
   cublas::scale(H.raw_handle(), n, &alpha, reinterpret_cast<double*>(x.device_ptr())+1, incx*2);
}

//  TODO: vector_norm_frob


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
   cublas::gemv(H.raw_handle(), cublas::cublas_trans(Atrans), M, N,
                &alpha, A.device_ptr(), lda, x.device_ptr(), incx,
                &beta, y.device_ptr(), incy);
   A.wait_for(y);
   x.wait_for(y);
}

// BLAS level 2 extensions

// in-place two matrix version implements
// C = alpha*A + beta*C
template <typename T>
inline
void
geam(char Atrans, int M, int N,
     T alpha, cuda::const_gpu_ptr<T> A, int lda,
     T beta, cuda::gpu_ptr<T> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   cublas::geam(H.raw_handle(), cublas::cublas_trans(Atrans), CUBLAS_OP_N, M, N,
                &alpha, A.device_ptr(), lda,
                &beta, C.device_ptr(), ldc,
                C.device_ptr(), ldc);
   A.wait_for(C);
}

// out-of-place version implements
// C = alpha*A + beta*B
template <typename T>
inline
void
geam(char Atrans, char Btrans, int M, int N,
     T alpha, cuda::const_gpu_ptr<T> A, int lda,
     T beta,  cuda::const_gpu_ptr<T> B, int ldb,
     cuda::gpu_ptr<T> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   cublas::geam(H.raw_handle(), cublas::cublas_trans(Atrans), cublas::cublas_trans(Btrans), M, N,
                &alpha, A.device_ptr(), lda,
                &beta, B.device_ptr(), ldb,
                C.device_ptr(), ldc);
   A.wait_for(C);
   B.wait_for(C);
}

template <typename T>
inline
void
matrix_copy(char Atrans, int M, int N, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam(Atrans, M, N, T(1.0), A, lda, T(0.0), B, ldb);
}

template <typename T>
inline
void
matrix_copy(int M, int N, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam('N', M, N, T(1.0), A, lda, T(0.0), B, ldb);
}

template <typename T>
inline
void
matrix_copy_scaled(char Atrans, int M, int N, T alpha, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam(Atrans, M, N, alpha, A, lda, T(0.0), B, ldb);
}

template <typename T>
inline
void
matrix_add(char Atrans, int M, int N, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam(Atrans, M, N, T(1.0), A, lda, T(1.0), B, ldb);
}

template <typename T>
inline
void
matrix_add_scaled(char Atrans, int M, int N, T alpha, const_gpu_ptr<T> A, int lda, gpu_ptr<T> B, int ldb)
{
   geam(Atrans, M, N, alpha, A, lda, T(1.0), B, ldb);
}

template <typename T>
inline
void
matrix_scale(int M, int N, T alpha, cuda::gpu_ptr<T> A, int lda)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(A.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   T beta = 0.0;
   cublas::geam(H.raw_handle(), CUBLAS_OP_N, CUBLAS_OP_N, M, N,
                &alpha, A.device_ptr(), lda,
                &beta, A.device_ptr(), lda,
                A.device_ptr(), lda);
}

template <typename T>
inline
void matrix_clear(int M, int N, cuda::gpu_ptr<T> A, int lda)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(A.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   T beta = 0.0;
   cublas::geam(H.raw_handle(), CUBLAS_OP_N, CUBLAS_OP_N, M, N,
                &beta, A.device_ptr(), lda,
                &beta, A.device_ptr(), lda,
                A.device_ptr(), lda);
}

inline
void matrix_conj(int M, int N, cuda::gpu_ptr<std::complex<double>> A, int lda)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(A.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   double alpha = -1.0;
   if (lda == M)
   {
      cublas::scale(H.raw_handle(), M*N, &alpha, reinterpret_cast<double*>(A.device_ptr())+1, 2);
   }
   else
   {
      for (int r = 0; r < N; ++r)
      {
	 cublas::scale(H.raw_handle(), M, &alpha, reinterpret_cast<double*>(A.device_ptr())+r*lda+1, 2);
      }
   }
}

template <typename T>
std::enable_if_t<cuda::is_cuda_floating_point_v<T>, void>
matrix_inner_prod(char Atrans, char Btrans, int M, int N,
		  cuda::const_gpu_ptr<T> A, int ldA,
		  cuda::const_gpu_ptr<T> B, int ldB,
		  cuda::gpu_ref<T>& r);

template <typename T, typename Nested>
inline
std::enable_if_t<cuda::is_cuda_floating_point_v<T>, void>
matrix_inner_prod_nested(char Atrans, char Btrans, int M, int N,
			 cuda::const_gpu_ptr<T> A, int ldA,
			 cuda::const_gpu_ptr<T> B, int ldB,
			 cuda::gpu_ref<T>& r, Nested&&)
{
   matrix_inner_prod(Atrans, Btrans, M, N, A, ldA, B, ldB, r);
}

template <typename T>
std::enable_if_t<cuda::is_cuda_floating_point_v<T>, void>
matrix_add_inner_prod(char Atrans, char Btrans, int M, int N,
		      cuda::const_gpu_ptr<T> A, int ldA,
		      cuda::const_gpu_ptr<T> B, int ldB,
		      cuda::gpu_ref<T>& r);

template <typename T, typename Nested>
inline
std::enable_if_t<cuda::is_cuda_floating_point_v<T>, void>
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

template <typename T>
inline
void
gemm(char Atrans, char Btrans, int M, int N, int K, T alpha,
     cuda::const_gpu_ptr<T> A, int lda, cuda::const_gpu_ptr<T> B, int ldb,
     T beta, cuda::gpu_ptr<T> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(B);
   cublas::gemm(H.raw_handle(), cublas::cublas_trans(Atrans), cublas::cublas_trans(Btrans), M, N, K,
                &alpha, A.device_ptr(), lda, B.device_ptr(), ldb,
                &beta, C.device_ptr(), ldc);
   A.wait_for(C);
   B.wait_for(C);
}

inline
void
dgmm(int M, int K,
     cuda::const_gpu_ptr<std::complex<double>> x, int incx,
     cuda::const_gpu_ptr<std::complex<double>> B, int ldb,
     cuda::gpu_ptr<std::complex<double>> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(x);
   C.wait_for(B);
   cublas::dgmm(H.raw_handle(), CUBLAS_SIDE_LEFT, M, K,
                B.device_ptr(), ldb,
                x.device_ptr(), incx,
                C.device_ptr(), ldc);
   x.wait_for(C);
   B.wait_for(C);
}

// mixed real/complex

// dgmm for mixed real/complex can't be implemented as a cublas call,
// instead we have our own kernel, in cub.h

inline
void
gdmm(int M, int K,
     cuda::const_gpu_ptr<std::complex<double>> A, int lda,
     cuda::const_gpu_ptr<std::complex<double>> y, int incy,
     cuda::gpu_ptr<std::complex<double>> C, int ldc)
{
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(y);
   cublas::dgmm(H.raw_handle(), CUBLAS_SIDE_RIGHT, M, K,
                A.device_ptr(), lda,
                y.device_ptr(), incy,
                C.device_ptr(), ldc);
   A.wait_for(C);
   y.wait_for(C);
}

// mixed real/complex
inline
void
gdmm(int M, int K,
     cuda::const_gpu_ptr<std::complex<double>> A, int lda,
     cuda::const_gpu_ptr<double> y, int incy,
     cuda::gpu_ptr<std::complex<double>> C, int ldc)
{
   // We can treat an M*K complex matrix as a real matrix of size (2M)*K.
   cublas::handle& H = cublas::get_handle();
   H.set_stream(C.get_stream());
   H.set_pointer_mode(CUBLAS_POINTER_MODE_HOST);
   C.wait_for(A);
   C.wait_for(y);
   cublas::dgmm(H.raw_handle(), CUBLAS_SIDE_RIGHT, 2*M, K,
                reinterpret_cast<double const*>(A.device_ptr()), 2*lda,
                y.device_ptr(), incy,
                reinterpret_cast<double*>(C.device_ptr()), 2*ldc);
   A.wait_for(C);
   y.wait_for(C);
}

} // namespace cuda

#include "cublas.icc"

#endif
