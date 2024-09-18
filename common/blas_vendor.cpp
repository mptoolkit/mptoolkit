// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/blas_vendor.cpp
//
// Copyright (C) 2024 Ian McCulloch <ian@qusim.net>
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

#if !defined(HAVE_CONFIG_H)
#error "need config.h to proceed!"
#endif
#include "config.h"

#include "blas_vendor.h"

#include <string>
#include <cstring>

#if defined(BLAS_VENDOR_MKL)

// MKL has two versions of the get_version_string function. The version as exported in the library as
// mkl_get_version_string is intended to be called from Fortran, and it completely fills the buffer with blank spaces,
// and doesn't null terminate. The acutal C function is MKL_Get_Version_String, and mkl.h #define's mkl_get_version_string
// to be MKL_Get_Version_String.  This version doesn't do padding and null-terminates the buffer.
extern "C" void MKL_Get_Version_String(char* buf, int len);

std::string BLAS_Vendor()
{
   return "MKL";
}

std::string BLAS_Version()
{
   char buf[198];
   MKL_Get_Version_String(buf, 198);
   return std::string(buf);
}

#elif defined(BLAS_VENDOR_OPENBLAS)

extern "C" char* openblas_get_config();

std::string BLAS_Vendor()
{
   return "OpenBLAS";
}

std::string BLAS_Version()
{
   return openblas_get_config();
}

#else
// generic BLAS
//

std::string BLAS_Vendor()
{
   return "generic";
}

std::string BLAS_Version()
{
   return "(no version information available)";
}

#endif
