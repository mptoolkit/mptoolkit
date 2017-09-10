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

#include "cuda-setup.h"

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#if !defined(HAVE_CUDA)
#error "CUDA is required to compile this program!"
#endif

#include "cuda.h"
#include "common/environment.h"

namespace cuda
{

bool is_cuda_enabled()
{
   return true;
}

int num_cuda_devices()
{
   return cuda::num_devices();
}

std::vector<std::string> get_cuda_device_names()
{
   int n = cuda::num_devices();
   std::vector<std::string> Result(n);
   for (int i = 0; i < n; ++i)
   {
      cuda::device_properties d = cuda::get_device_properties(i);
      Result[i] = d.name();
   }
   return Result;
}

int const dev = getenv_or_default("MP_CUDA_DEVICE", 0);
//int const num_dev = cuda::num_devices();

int mp_cuda_device()
{
   int d = getenv_or_default("MP_CUDA_DEVICE", 0);
   if (dev == -1 || dev >= cuda::num_devices())
      return -1;
   return d;
}

int setup_cuda()
{
   if (dev == -1 || dev >= cuda::num_devices())
      return -1;
   cuda::set_device(dev);
   return dev;
}

int setup_cuda_thread()
{
   if (dev == -1 || dev >= cuda::num_devices())
      return -1;
   cuda::set_device(dev);
   return dev;
}
} // namespace cuda
