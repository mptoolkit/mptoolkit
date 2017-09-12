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

#if !defined(MPTOOLKIT_CUDA_CUBLAS_H)
#define MPTOOLKIT_CUDA_CUDBLAS_H

#include "cuda.h"
#include <list>
#include <mutex>
#include <cublas_v2.h>

#include <iostream>

namespace cublas
{

// returns the cublas version number
int version();

// cublas handle.  Moveable, but not copyable.
class handle
{
   public:
      handle() : h_(nullptr) {}
      handle(handle&& other) : h_(other.h_) { other.h_ = nullptr; }
      handle(handle const&) = delete;
      handle& operator=(handle&& other) { std::swap(h_, other.h_); return *this; }
      hasndle& operator=(handle const&) = delete;
      ~handle() { if (h_) cublasDestroy(h_); }

      cublasHandle_t raw_handle() const { return h_; }

      static handle create() { cublasHandle_t h; cublasCreate(&h); return handle(h); }

      void destroy() { cublasDestroy(h_); h_ = nullptr; }

   private:
      handle(cublasHandle_t h) : h_(h) {}

      cublasHandle_t h_;
};

// set the stream associated with the given cublas handle
void set_stream(handle const& h, cuda::stream const& s);

class gpu_ref;

template <typename T>
class gpu_buffer
{
   public:

      gpu_ref<T> operator[](int n);

      // block modifications to this buffer until event e has been triggered
      void wait(event const& e);

      // 'lock' the buffer, ie forbid read/write operations until it is unlocked.
      // Locking can be done recursively.
      void lock();

      void unlock();

   private:
      T* Ptr;
      cuda::stream Stream;
      cuda::event Sync;
      int LockCount;

      // if we want to reference count the references, then this is what it looks like:
      std::map<size_t, gpu_ref> LiveReferencex;
};

template <typename T>
void
gpu_buffer<T>::lock()
{
   ++LockCount;
}

template <typename T>
void
gpu_buffer<T>::unlock()
{
   --LockCount;
}


template <typename T>
class gpu_ref
{
   public:
      gpu_ref() = delete;

      ~gpu_ref()
      {
	 Sync.record(Stream);
	 Buf->wait(Sync);
      }

      gpu_ref(gpu_ref&& other) { using std::swap; swap(Buf, other.Buf); swap(Ptr, other.Ptr); swap(Stream, other.Stream);
	 swap(Sync, other.Sync); }

      gpu_ref(gpu_ref const&) = delete;

      gpu_ref& operator=(gpu_ref<T>&&) = delete;

      gpu_ref& operator=(gpu_ref<T> const& Other)
      {
	 // wait until the other stream has finished writing
	 Stream.wait(Other.Sync);
	 // do the copy
	 memcpy_device_to_device_async(Stream, Other.Ptr, Ptr, sizeof(T));
	 // signal our event that we have finished writing
	 Sync.record(Stream);
	 // make the other stream wait until the copy is finished
	 Other.Stream.wait(Sync);
	 return *this;
      }

      void set_block(T const& x)
      {
	 memcpy_host_to_device(Stream, &x, Other.Ptr, sizeof(T));
	 Stream.record(Sync);
      }

      // sets the element asynchronously.
      void set(T const& x)
      {
	 T const* Ptr = get_global_constant(x);
	 memcpy_host_to_device_async(Stream, Ptr, Other.Ptr, sizeof(T));
	 Stream.record(Sync);
      }

      // blocking get operation
      T get_block()
      {
	 T x;
	 memcpy_device_to_host(Ptr, &x, sizeof(T));
	 return x;
      }

      // non-blocking get.  Returns a cuda::event, when triggered the host memory location is valid.
      cuda::event get(T* x)
      {
	 cuda::event E;
	 memcpy_device_to_host_async(Stream, Ptr, x, sizeof(T));
	 Stream.record(E);
	 return E;
      }

   private:
      friend gpu_ref<T> gpu_buffer<T>::operator[](int n);

      gpu_ref(gpu_buffer* Buf_, T* Ptr_) : Buf(Buf_), Ptr(Ptr_) { Buf->lock(); }

      gup_buffer* Buf;
      T* Ptr;
      cuda::stream Stream;
      cuda::event Sync;
};



} // namespace cublas

#endif
