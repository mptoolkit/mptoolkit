// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/cusovler.cpp
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

#include "cusolver.h"

namespace cusolver
{

char const*
GetErrorName(cusolverStatus_t error)
{
    switch (error)
    {
        case CUSOLVER_STATUS_SUCCESS:
            return "CUSOLVER_STATUS_SUCCESS";

        case CUSOLVER_STATUS_NOT_INITIALIZED:
            return "CUSOLVER_STATUS_NOT_INITIALIZED";

        case CUSOLVER_STATUS_ALLOC_FAILED:
            return "CUSOLVER_STATUS_ALLOC_FAILED";

        case CUSOLVER_STATUS_INVALID_VALUE:
            return "CUSOLVER_STATUS_INVALID_VALUE";

        case CUSOLVER_STATUS_ARCH_MISMATCH:
            return "CUSOLVER_STATUS_ARCH_MISMATCH";

        case CUSOLVER_STATUS_EXECUTION_FAILED:
            return "CUSOLVER_STATUS_EXECUTION_FAILED";

        case CUSOLVER_STATUS_INTERNAL_ERROR:
            return "CUSOLVER_STATUS_INTERNAL_ERROR";

        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
           return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
    }

    return "<unknown>";
}

char const*
GetErrorString(cusolverStatus_t error)
{
    switch (error)
    {
        case CUSOLVER_STATUS_SUCCESS:
            return "cuSOLVER: SUCCESS";

        case CUSOLVER_STATUS_NOT_INITIALIZED:
            return "cuSOLVER: library not initialized";

        case CUSOLVER_STATUS_ALLOC_FAILED:
            return "cuSOLVER: memory allocation failed";

        case CUSOLVER_STATUS_INVALID_VALUE:
            return "cuSOLVER: invalid value or parameter";

        case CUSOLVER_STATUS_ARCH_MISMATCH:
            return "cuSOLVER: feature not supported by this architecture (possible double-precision?)";

        case CUSOLVER_STATUS_EXECUTION_FAILED:
            return "cuSOLVER: GPU kernel execution failed";

        case CUSOLVER_STATUS_INTERNAL_ERROR:
            return "cuSOLVER: internal error";

        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
           return "cuSOLVER: matrix type not supported";

    }

    return "<cusolver-unknown>";
}

} // namespace cusolver

namespace blas
{

namespace detail
{

void DiagonalizeSymmetric(int Size, cuda::gpu_ptr<double> A, int ldA, cuda::gpu_ptr<double> Eigen)
{
   cusolver::handle& H = cusolver::get_handle();
   H.set_stream(A.get_stream());
   A.wait_for(Eigen);
   int lWork;
   cusolver::check_error(cusolverDnDsyevd_bufferSize(H.raw_handle(),
                                                     CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                                     Size, A.device_ptr(), ldA, Eigen.device_ptr(), &lWork));
   double* Work = static_cast<double*>(cuda::allocate_gpu_temp_memory(lWork*sizeof(double)));
   int* DevInfo = static_cast<int*>(cuda::allocate_gpu_temp_memory(sizeof(int)));
   cusolver::check_error(cusolverDnDsyevd(H.raw_handle(),
                                          CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                          Size, A.device_ptr(), ldA, Eigen.device_ptr(), Work, lWork, DevInfo));
   int Info;
   cuda::memcpy_device_to_host(DevInfo, &Info, sizeof(int));
   CHECK_EQUAL(Info, 0);
   Eigen.wait_for(A);
   cuda::free_gpu_temp_memory(DevInfo, sizeof(int));
   cuda::free_gpu_temp_memory(Work, lWork*sizeof(double));
}

void DiagonalizeHermitian(int Size, cuda::gpu_ptr<std::complex<double>> A, int ldA, 
			  cuda::gpu_ptr<double> Eigen)
{
   cusolver::handle& H = cusolver::get_handle();
   H.set_stream(A.get_stream());
   A.wait_for(Eigen);
   int lWork;
   cusolver::check_error(cusolverDnZheevd_bufferSize(H.raw_handle(),
                                                     CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                                     Size, 
						     reinterpret_cast<cuDoubleComplex*>(A.device_ptr()),
						     ldA, Eigen.device_ptr(), &lWork));
   cuDoubleComplex* Work = static_cast<cuDoubleComplex*>(cuda::allocate_gpu_temp_memory(lWork*sizeof(cuDoubleComplex)));
   int* DevInfo = static_cast<int*>(cuda::allocate_gpu_temp_memory(sizeof(int)));
   cusolver::check_error(cusolverDnZheevd(H.raw_handle(),
                                          CUSOLVER_EIG_MODE_VECTOR, CUBLAS_FILL_MODE_LOWER,
                                          Size, 
					  reinterpret_cast<cuDoubleComplex*>(A.device_ptr()), 
					  ldA, Eigen.device_ptr(), Work, lWork, DevInfo));
   int Info;
   cuda::memcpy_device_to_host(DevInfo, &Info, sizeof(int));
   CHECK_EQUAL(Info, 0);
   Eigen.wait_for(A);
   cuda::free_gpu_temp_memory(DevInfo, sizeof(int));
   cuda::free_gpu_temp_memory(Work, lWork*sizeof(cuDoubleComplex));
}

void SingularValueDecomposition(char job, int Rows, int Cols, 
				cuda::gpu_ptr<double> Data, int LeadingDim, 
				cuda::gpu_ptr<double> Dvec,
				cuda::gpu_ptr<double> Umat, int ldU, 
				cuda::gpu_ptr<double> Vmat, int ldV)
{
   int min_mn = std::min(Rows, Cols);
   cusolver::handle& H = cusolver::get_handle();
   H.set_stream(Data.get_stream());
   Data.wait_for(Dvec);
   Data.wait_for(Umat);
   Data.wait_for(Vmat);
   int lWork;
   cusolver::check_error(cusolverDnDgesvd_bufferSize(H.raw_handle(), Rows, Cols, &lWork));

   double* Work = static_cast<double*>(cuda::allocate_gpu_temp_memory(lWork*sizeof(double)));
   double* rWork = static_cast<double*>(cuda::allocate_gpu_temp_memory((min_mn-1)*sizeof(double)));

   cuda::gpu_ref<int> DevInfo;
   cusolver::check_error(cusolverDnDgesvd(H.raw_handle(), job, job, Rows, Cols, Data.device_ptr(), LeadingDim,
					  Dvec.device_ptr(), Umat.device_ptr(), ldU,
					  Vmat.device_ptr(), ldV, Work, lWork, rWork, DevInfo.device_ptr()));

   Dvec.wait_for(Data);
   Umat.wait_for(Data);
   Vmat.wait_for(Data);
   cuda::free_gpu_temp_memory(Work, lWork*sizeof(double));
   cuda::free_gpu_temp_memory(rWork, (min_mn-1)*sizeof(double));
}

void SingularValueDecomposition(int Rows, int Cols, 
				cuda::gpu_ptr<double> Data, int LeadingDim, 
				cuda::gpu_ptr<double> Dvec,
				cuda::gpu_ptr<double> Umat, int ldU, 
				cuda::gpu_ptr<double> Vmat, int ldV)
{
   SingularValueDecomposition('S', Rows, Cols, Data, LeadingDim, Dvec, Umat, ldU, Vmat, ldV);
}

void SingularValueDecompositionFull(int Rows, int Cols, 
				    cuda::gpu_ptr<double> Data, int LeadingDim, 
				    cuda::gpu_ptr<double> Dvec,
				    cuda::gpu_ptr<double> Umat, int ldU, 
				    cuda::gpu_ptr<double> Vmat, int ldV)
{
   SingularValueDecomposition('A', Rows, Cols, Data, LeadingDim, Dvec, Umat, ldU, Vmat, ldV);
}


void SingularValueDecomposition(char job, int Rows, int Cols, 
				cuda::gpu_ptr<std::complex<double>> Data, int LeadingDim, 
				cuda::gpu_ptr<double> Dvec,
				cuda::gpu_ptr<std::complex<double>> Umat, int ldU, 
				cuda::gpu_ptr<std::complex<double>> Vmat, int ldV)
{
   int min_mn = std::min(Rows, Cols);
   cusolver::handle& H = cusolver::get_handle();
   H.set_stream(Data.get_stream());
   Data.wait_for(Dvec);
   Data.wait_for(Umat);
   Data.wait_for(Vmat);
   int lWork;
   cusolver::check_error(cusolverDnZgesvd_bufferSize(H.raw_handle(), Rows, Cols, &lWork));
   
   cuDoubleComplex* Work = static_cast<cuDoubleComplex*>(cuda::allocate_gpu_temp_memory(lWork*sizeof(cuDoubleComplex)));
   double* rWork = static_cast<double*>(cuda::allocate_gpu_temp_memory((min_mn-1)*sizeof(double)));

   cuda::gpu_ref<int> DevInfo;
   cusolver::check_error(cusolverDnZgesvd(H.raw_handle(), job, job, Rows, Cols, 
					  reinterpret_cast<cuDoubleComplex*>(Data.device_ptr()), LeadingDim,
					  Dvec.device_ptr(), 
					  reinterpret_cast<cuDoubleComplex*>(Umat.device_ptr()), ldU,
					  reinterpret_cast<cuDoubleComplex*>(Vmat.device_ptr()), ldV, 
					  Work, lWork, rWork, DevInfo.device_ptr()));

   Dvec.wait_for(Data);
   Umat.wait_for(Data);
   Vmat.wait_for(Data);
   cuda::free_gpu_temp_memory(Work, lWork*sizeof(cuDoubleComplex));
   cuda::free_gpu_temp_memory(rWork, (min_mn-1)*sizeof(cuDoubleComplex));
}

void SingularValueDecomposition(int Rows, int Cols, 
				cuda::gpu_ptr<std::complex<double>> Data, int LeadingDim, 
				cuda::gpu_ptr<double> Dvec,
				cuda::gpu_ptr<std::complex<double>> Umat, int ldU, 
				cuda::gpu_ptr<std::complex<double>> Vmat, int ldV)
{
   SingularValueDecomposition('S', Rows, Cols, Data, LeadingDim, Dvec, Umat, ldU, Vmat, ldV);
}

void SingularValueDecompositionFull(int Rows, int Cols, 
				    cuda::gpu_ptr<std::complex<double>> Data, int LeadingDim, 
				    cuda::gpu_ptr<double> Dvec,
				    cuda::gpu_ptr<std::complex<double>> Umat, int ldU, 
				    cuda::gpu_ptr<std::complex<double>> Vmat, int ldV)
{
   SingularValueDecomposition('A', Rows, Cols, Data, LeadingDim, Dvec, Umat, ldU, Vmat, ldV);
}

} // namespace detail

} // namespace cuda
