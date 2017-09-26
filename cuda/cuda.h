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

#include <list>
#include <mutex>
#include <cuda_runtime.h>
#include "common/atomicrefcount.h"

#include <iostream>

namespace cuda
{

// wrapper for a cudaError_t, which can also function as an exception object
class error : public std::runtime_error
{
   public:
      error() = delete;
      error(cudaError_t Err) : std::runtime_error(cudaGetErrorString(Err)), err_(Err)
      {std::cerr << "Error " << int(Err) << '\n';}

      cudaError_t code() const { return err_; }
      operator cudaError_t() const { return err_; }

      virtual const char* what() noexcept { return cudaGetErrorName(err_); }

      char const* name() const { return cudaGetErrorName(err_); }
      char const* string() const { return cudaGetErrorString(err_); }

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

// device enumeration

// returns the number of devices with compute capability >= 2
int num_devices();

// sets the device.  This must be called per host thread.
void set_device(int d);

// returns which device is currently in use
int get_device();

class device_properties
{
   public:
      device_properties() {}
      device_properties(cudaDeviceProp const& p) : p_(p) {}

      std::string name() const { return p_.name; }
      std::size_t total_global_memory() const { return p_.totalGlobalMem; }

      // TODO: many other properties we can add here

   private:
      friend device_properties get_device_properties(int d);
      cudaDeviceProp p_;
};

// get the properties of a specific device
device_properties get_device_properties(int d);

// synchronize on the device, and block until all submitted tasks are finished
void device_synchronize();

class event;

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

      // wait for the given event
      void wait(event const& e);

      friend void swap(stream& a, stream& b)
      {
	 std::swap(a.stream_, b.stream_);
      }

   private:
      static cudaStream_t Allocate();
      static std::mutex FreeListMutex;
      static std::list<cudaStream_t> FreeList;

      cudaStream_t stream_;
      AtomicRefCount count_;
};

class event
{
   public:
      event();
      event(event const&) = delete;
      event(event&& other);
      event& operator=(event const&) = delete;
      event& operator=(event&& other);
      ~event();

      // record the event as part of stream s
      void record(stream const& s);

      // returns true if work has been sucessfully completed
      bool is_complete() const;

      cudaEvent_t raw_event() const { return event_; }

      friend void swap(event& a, event& b)
      {
	 std::swap(a.event_, b.event_);
      }

   private:
      static cudaEvent_t Allocate();
      static std::mutex FreeListMutex;
      static std::list<cudaEvent_t> FreeList;

      cudaEvent_t event_;
};

inline
void
stream::wait(event const& e)
{
   check_error(cudaStreamWaitEvent(stream_, e.raw_event(), 0));
}

inline
void
event::record(stream const& s)
{
   check_error(cudaEventRecord(event_, s.raw_stream()));
}

class timer
{
   public:
      timer();
      timer(timer const&) = delete;
      timer(timer&& other);
      timer& operator=(timer const&) = delete;
      timer& operator=(timer&& other);
      ~timer();

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

      friend void swap(timer& a, timer& b)
      {
	 std::swap(a.start_, b.start_);
	 std::swap(a.stop_, b.stop_);
      }

   private:
      static cudaEvent_t Allocate();
      static std::mutex FreeListMutex;
      static std::list<cudaEvent_t> FreeList;

      cudaEvent_t start_;
      cudaEvent_t stop_;
};

// copy GPU memory asyncronously
void memcpy_device_to_device_async(stream const& s, void const* src, void* dest, std::size_t size);

inline
void memcpy_device_to_device_async(stream const& s, void const* src, void* dest, std::size_t size)
{
   cudaMemcpyAsync(dest, src, size, cudaMemcpyDeviceToDevice, s.raw_stream());
}

void memcpy_host_to_device_async(stream const& s, void const* src, void* dest, std::size_t size);

inline
void memcpy_host_to_device_async(stream const& s, void const* src, void* dest, std::size_t size)
{
   cudaMemcpyAsync(dest, src, size, cudaMemcpyHostToDevice, s.raw_stream());
}


} // namespace cuda

#include "cuda.icc"

#endif
