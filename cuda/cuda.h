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

//
// Debug support:
// define CUDA_TRACE_FUNCTIONS to get detailed information for CUDA API calls.
//
// define CUDA_SYNCHRONIZE to force all streams to use stream 0, and synchronize
// after each API call.
//
// define CUDA_TRACE_STREAMS to print debug info whenever a stream is created or destroyed.
//

#if !defined(MPTOOLKIT_CUDA_CUDA_H)
#define MPTOOLKIT_CUDA_CUDA_H

// cuda.h must be valid C++11, since CUDA 8 doesn't have C++14 support.

#include <list>
#include <mutex>
#include <memory>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include "common/atomicrefcount.h"
#include <iostream>

#if defined(CUDA_TRACE_FUNCTIONS)
#define TRACE_CUDA(Msg) TRACE(Msg)
#else
#define TRACE_CUDA(Msg) DUMMY_TRACE(Msg)
#endif

#if defined(CUDA_TRACE_STREAMS)
#define TRACE_STREAM(Msg) TRACE(Msg)
#else
#define TRACE_STREAM(Msg) DUMMY_TRACE(Msg)
#endif


namespace cuda
{

// wrapper for a cudaError_t, which can also function as an exception object
class error : public std::runtime_error
{
   public:
      error() = delete;
      error(cudaError_t Err) : std::runtime_error(cudaGetErrorString(Err)), err_(Err)
      {std::cerr << "CUDA Error " << int(Err) << ' ' << cudaGetErrorString(Err) << '\n';}

      cudaError_t code() const { return err_; }
      operator cudaError_t() const { return err_; }

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
#if defined(CUDA_SYNCHRONIZE)
   e = cudaDeviceSynchronize();
   if (e != cudaSuccess)
      throw error(e);
#endif
}

// device enumeration

// returns the number of devices with compute capability >= 2
int num_devices();

// returns a string representation of the number of cuda devices;
// this gives more information in the case where there are no cuda
// devices (eg, whether there are no devices or whether there is some driver error).
std::string num_devices_str();

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

// wrapper for a cuda stream.
// Streams are allocated by a pool.
class stream
{
   public:
      stream();
      stream(stream const&);
      stream(stream&& other);
      stream& operator=(stream const&);
      stream& operator=(stream&& other);
      ~stream();

      // Normally the stream destructor will emit a warning if we destroy the
      // stream while it has pending operations, becuase this is often a bug (eg,
      // the memory resources associated with the stream will probably get deallocated
      // at the same time as the stream).  But in some cases this is OK, eg if the stream
      // is part of a gpu_ref that is synchronized to a parent stream.  In that case,
      // we can call this function immediately prior to destroying the stream.
      void safe_to_destroy();

      cudaStream_t raw_stream() const { return stream_; }

      // record the event at the current location in the stream,
      // as a single-use event
      event record() const;

      // blocks until the stream has completed all operations
      void synchronize() const;

      // wait for the given event
      void wait(event const& e) const;

      // returns true if there are pending operations on the stream.
      // returns false if all operations on the stream have completed.
      bool is_running() const;

      friend void swap(stream& a, stream& b)
      {
         using std::swap;
	 swap(a.stream_, b.stream_);
         swap(a.count_, b.count_);
      }

   private:
      void sub_reference();

      cudaStream_t stream_;
      shared_counter count_;
};

//
// event
//
// A reference-counted wrapper around the cuda event type.
// This uses a double-indirection (via shared_ptr<cudaEvent_t>)
//
//

class event
{
   public:
      event();
      event(event const& other);
      event(event&& other);
      event& operator=(event const& other);
      event& operator=(event&& other);
      ~event();

      // clears the event - equivalent to *this = event()
      void clear();

      // returns true if this is a null event
      bool is_null() const { return !event_; }

      // returns true if work has been sucessfully completed
      bool is_complete() const;

      // blocks the caller until the event has triggered
      void wait() const;

      cudaEvent_t raw_event() const { return *event_; }

   private:
      friend class stream;

      // the only constructor (aside from move/copy-construction)
      explicit event(cudaStream_t s);

      void sub_reference();

      std::shared_ptr<cudaEvent_t> event_;
      shared_counter count_;

      friend class event_ref;
};

class event_ref
{
   public:
      event_ref() = delete;

      event_ref(event const& e)
         : event_(e.event_) {}

      event_ref(event_ref&& Other) = default;
      event_ref(event_ref const& other) = default;

      event_ref& operator=(event_ref const& Other) = default;
      event_ref& operator=(event_ref&& Other) = default;

      // returns true if work has been sucessfully completed
      bool is_complete() const;

   private:
      std::shared_ptr<cudaEvent_t> event_;
};

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

// cuda floating point types
// is_cuda_floating_point metafunction is true for floating point and complex
// types that are supported on cuda.

template <typename T>
struct is_cuda_floating_point : std::false_type {};

#if defined(__cpp_variable_templates)
template <typename T>
constexpr bool is_cuda_floating_point_v = is_cuda_floating_point<T>::value;
#endif

template <> struct is_cuda_floating_point<float> : std::true_type {};
template <> struct is_cuda_floating_point<double> : std::true_type {};
template <> struct is_cuda_floating_point<std::complex<float>> : std::true_type {};
template <> struct is_cuda_floating_point<std::complex<double>> : std::true_type {};

// standard complex arithmetic doesn't work on cuda - the std::complex functions
// are not specified to run on the GPU.  instead we need to cast to the types
// cuComplex and cuDoubleComplex.

template <typename T>
struct cuda_complex;

template <typename T>
struct cuda_complex<T const> : cuda_complex<T> {};

template <typename T>
struct cuda_complex<T&> : cuda_complex<T> {};

template <typename T>
struct cuda_complex<T&&> : cuda_complex<T> {};

template <>
struct cuda_complex<float>
{
   using type = float;
};

template <>
struct cuda_complex<double>
{
   using type = double;
};

template <>
struct cuda_complex<std::complex<float>>
{
   using type = cuFloatComplex;
};

template <>
struct cuda_complex<std::complex<double>>
{
   using type = cuDoubleComplex;
};

template <typename T>
using cuda_complex_t = typename cuda_complex<T>::type;

template <typename T>
typename cuda_complex<T>::type*
cuda_complex_cast(T* x)
{
   return reinterpret_cast<cuda_complex_t<T>*>(x);
}

template <typename T>
typename cuda_complex<T>::type const*
cuda_complex_cast(T const* x)
{
   return reinterpret_cast<cuda_complex_t<T> const*>(x);
}

// copy GPU memory syncronously
void memcpy_device_to_host(void const* src, void* dest, std::size_t size);

inline
void memcpy_device_to_host(void const* src, void* dest, std::size_t size)
{
   TRACE_CUDA("cudaMemcpy cudaMemcpyDeviceToHost")(dest)(src)(size);
   check_error(cudaMemcpy(dest, src, size, cudaMemcpyDeviceToHost));
}

inline
void memcpy_host_to_device(void const* src, void* dest, std::size_t size)
{
   TRACE_CUDA("cudaMemcpy cudaMemcpyHostToDevice")(dest)(src)(size);
   check_error(cudaMemcpy(dest, src, size, cudaMemcpyHostToDevice));
}

// copy GPU memory asyncronously
void memcpy_device_to_device_async(stream const& s, void const* src, void* dest, std::size_t size);

inline
void memcpy_device_to_device_async(stream const& s, void const* src, void* dest, std::size_t size)
{
   TRACE_CUDA("cudaMemcpyAsync cudaMemcpyDeviceToDevice")(dest)(src)(size)(s.raw_stream());
   check_error(cudaMemcpyAsync(dest, src, size, cudaMemcpyDeviceToDevice, s.raw_stream()));
}

void memcpy_host_to_device_async(stream const& s, void const* src, void* dest, std::size_t size);

inline
void memcpy_host_to_device_async(stream const& s, void const* src, void* dest, std::size_t size)
{
   TRACE_CUDA("cudaMemcpyAsync cudaMemcpyHostToDevice")(dest)(src)(size)(s.raw_stream());
   check_error(cudaMemcpyAsync(dest, src, size, cudaMemcpyHostToDevice, s.raw_stream()));
}

inline
void memset(void* dest, int x, std::size_t size)
{
   TRACE_CUDA("cudaMemset")(dest)(x)(size);
   check_error(cudaMemset(dest, x, size));
}

inline
void memset_async(stream const& s, void* dest, int x, std::size_t size)
{
   TRACE_CUDA("cudaMemsetAsync")(dest)(x)(size)(s.raw_stream());
   check_error(cudaMemsetAsync(dest, x, size, s.raw_stream()));
}

} // namespace cuda

#include "cuda.icc"

#endif
