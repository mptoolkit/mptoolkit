// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// cuda/cuda-setup.h
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
// Basic functions for cuda initialization and device enumeration.
// This header can be used even if there is no cuda environment;
// in that case it implements dummy versions of the interfaces
// that indicate that cuda is not enabled, and no devices are available.

#if !defined(MPTOOLKIT_CUDA_CUDA_SETUP_H)
#define MPTOOLKIT_CUDA_CUDA_SETUP_H

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include <vector>
#include <string>

namespace cuda
{

// returns true if the cuda runtime library is available
bool is_cuda_enabled();

// returns the number of available cuda devices (0 if cuda isn't available)
int num_cuda_devices();

// returns the number of available cuda devices (0 and error message if cuda isn't available)
std::string num_cuda_devices_str();

// returns a vector of names of the available cuda devices (empty
// if cuda isn't available)
std::vector<std::string> get_cuda_device_names();

// returns the device number set by the MP_CUDA_DEVICE environment variable, or
// -1 if the variable isn't set
int mp_cuda_device();

// initializes cuda, intended to be called once by the main thread.
// returns the cuda device that will be used (mp_cuda_device(), or 0 if not set)
// Returns -1 if cuda isn't available, or the MP_CUDA_DEVICE is invalid
int setup_cuda();

// sets up cuda in a new thread.  This must be called by every new thread prior to
// calling any cuda functions.
int setup_cuda_thread();

// dummy implementation for the case where we don't have cuda enabled
#if !defined(HAVE_CUDA)
inline
bool is_cuda_enabled()
{
   return false;
}

inline
int num_cuda_devices()
{
   return 0;
}

inline
std::string num_cuda_devices_str()
{
   return "0 (cuda is not enabled)";
}

inline
std::vector<std::string> get_cuda_device_names()
{
   return std::vector<std::string>();
}

inline
int mp_cuda_device()
{
   return -1;
}

inline
int setup_cuda()
{
   return -1;
}

inline
int setup_cuda_thread()
{
   return -1;
}
#endif // !HAVE_CUDA

} // namespace cuda

#endif
