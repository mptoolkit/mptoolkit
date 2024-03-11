// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/blas_vendor.h
//
// Copyright (C) 2012-2024 Ian McCulloch <ian@qusim.net>
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

#if !defined(COMMON_BLAS_VENDOR_H)
#define COMMON_BLAS_VENDOR_H

#include <string>

// returns BLAS vendor information
std::string BLAS_Vendor();

std::string BLAS_Version();

#endif
