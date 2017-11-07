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

#include "cusovler.h"

namespace cusovler
{

char const*
GetErrorName(cublasStatus_t error)
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
GetErrorString(cublasStatus_t error)
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

} // namespace cusovler
