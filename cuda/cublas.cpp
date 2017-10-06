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

#include "cublas.h"

namespace cublas
{

char const*
cublasGetErrorName(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}

char const*
cublasGetErrorString(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "cuBLAS library not initialized";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "memory allocation failed";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "invalid value or parameter";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "feature not supported by this architecture (possible double-precision?)";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "invalid GPU memory mapping";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "GPU kernel execution failed";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "internal error";
    }

    return "<cublas-unknown>";
}

} // namespace cubas
