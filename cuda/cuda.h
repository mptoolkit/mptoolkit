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

#if !defined(MPTOOLKIT_CUDA_CUDA_H)
#define MPTOOLKIT_CUDA_CUDA_H

#include <mutex>
#include <cuda_runtime.h>

namespace cuda
{

// wrapper for a cudaError_t, which can also function as an exception object
class error : public std::runtime_error
{
   public:
      error() = delete;
      error(cudaError_t Err) : err_(Err) {}

      cuda_error_t() error() const { return err_; }
      operator cudaError_t() const { return err_; }

      virtual const char* what() noexcept { return cuda_errorName(err_); }

      char const* name() const { return cudaErrorName(err_); }
      char const* string() const { return cudaErrorString(err_); }

   private:
      cudaError_t err_;
};

// helper function to check error returns and throw an exception on error
inline
void check_error(cudaError_t e)
{
   if (e != cudaSuccess)
      throw error(e);
}
   
// wrapper for a cuda stream.  Moveable, but not copyable.
// Streams are allocated by a pool.
class stream
{
   public:
      stream();
      stream(stream const&) = delete;
      stream(stream&& other);
      stream& operator=(stream const&) = delete;
      stream& operator=(stream&& other);
      ~stream();

      cudaStream_t raw_stream() const { return stream_; }

   private:
      static cudaStream_t Allocate();
      static std::mutex FreeListMutex;
      static std::list<cudaStream_t> FreeList;

      cudaStream_t stream_;
}:

class event
{
   public:
      event();

      // record the event as part of stream s
      void record(stream const& s);

      // returns true if work has been sucessfully completed
      bool is_complete() const;

      cudaEvent_t raw_event() const { return event_; }

   private:
      static cudaStream_t Allocate();
      static std::mutex FreeListMutex;
      static std::list<cudaEvent_t> FreeList;

      cudaEvent_t event_;
};

class timer
{
   public:
      timer();
      
      // record the starting point as part of stream s
      void start(stream const& s);

      // record the stop point as part of stream s
      void stop(stream const& s);

      // returns true if the start point has been triggered
      bool is_started() const;

      // returns true if the start and stop points have been triggered
      bool is_complete() const;

      // precondition: is_complete()
      // returns the time duration
      float elapsed_time_ms() const;

      cudaEvent_t raw_start() const { return start_; }
      cudaEvent_t raw_stop() const { return stop_; }

   private:
      static cudaStream_t Allocate();
      static std::mutex FreeListMutex;
      static std::list<cudaEvent_t> FreeList;

      cudaEvent_t start_;
      cudaEvent_t stop_;
};


} // namespace cuda

#include "cuda.ipp"

#endif
