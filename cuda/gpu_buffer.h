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

#if !defined(MPTOOLKIT_CUDA_GPU_BUFFER_H)
#define MPTOOLKIT_CUDA_GPU_BUFFER_H

#include "cuda.h"
#include <memory>
#include "common/trace.h"
#include "blas/arena.h"

namespace cuda
{

template <typename T>
class gpu_ref;

inline
std::size_t round_up(std::size_t numToRound, std::size_t multiple)
{
    if (multiple == 0)
        return numToRound;

    std::size_t remainder = numToRound % multiple;
    if (remainder == 0)
        return numToRound;

    return numToRound + multiple - remainder;
}

class AllocationBlock
{
   public:
      AllocationBlock() = delete;

      explicit AllocationBlock(std::size_t Size_, bool FreeOnDest = true);

      AllocationBlock(AllocationBlock&& Other) : Align(Other.Align),
						 Size(Other.Size),
						 BasePtr(Other.BasePtr),
						 CurrentOffset(Other.CurrentOffset),
						 NumAllocations(Other.NumAllocations),
                                                 FreeOnDestructor(Other.FreeOnDestructor)
      {
	 Other.BasePtr = nullptr;
      }

      AllocationBlock(AllocationBlock const&) = delete;

      AllocationBlock& operator=(AllocationBlock&&) = delete;
      AllocationBlock& operator=(AllocationBlock const&) = delete;

      ~AllocationBlock() noexcept;

      // attempts to allocate Size bytes, aligned at Align.  Returns nullptr if no allocation
      // was possible.
      void* try_allocate(size_t RequestSize, size_t RequestAlign)
      {
	 CHECK(Align % RequestAlign == 0)(Align)(RequestAlign);

	 if (RequestSize > Size)
	    return nullptr;

	 size_t NextOffset = round_up(CurrentOffset, RequestAlign);
	 if (NextOffset > Size - RequestSize)
	    return nullptr;

         ++NumAllocations;
	 CurrentOffset = NextOffset + RequestSize;
	 return static_cast<void*>(BasePtr + NextOffset);
      }

      // free memory that was previously allocated.
      // Returns true if this caused the allocation block to become empty
      bool free(void* Ptr, size_t AllocSize)
      {
	 // TOOD: add some debugging stuff
	 if (--NumAllocations == 0)
	 {
	    CurrentOffset = 0;
	    return true;
	 }
         return false;
      }

      bool try_free(void* Ptr, size_t AllocSize)
      {
	 if (static_cast<unsigned char*>(Ptr) < BasePtr)
	    return false;

	 if (BasePtr+Size <= static_cast<unsigned char*>(Ptr))
	    return false;

	 this->free(Ptr, Size);
	 return true;
      }

   private:
      size_t Align;   // minimum alignment of this block
      size_t Size;
      unsigned char* BasePtr;  // base pointer of the allocation
      size_t CurrentOffset;    // current one-past-the-end of allocated blocks
      int NumAllocations;
      bool FreeOnDestructor;
};

inline
AllocationBlock::AllocationBlock(std::size_t Size_, bool Free)
   : Align(256), Size(Size_), CurrentOffset(0), NumAllocations(0), FreeOnDestructor(Free)
{
   void* Ptr;
   // memory alignment from cudaMalloc is always at least 256 bytes
   // TODO: handle out-of-memory conditions properly.
   check_error(cudaMalloc(&Ptr, Size_));
   BasePtr = static_cast<unsigned char*>(Ptr);
}

inline
AllocationBlock::~AllocationBlock() noexcept
{
   if (FreeOnDestructor && BasePtr)
   {
      DEBUG_TRACE_IF(NumAllocations != 0)(NumAllocations);
      if (NumAllocations == 0)
         check_error(cudaFree(static_cast<void*>(BasePtr)));
   }
}

std::size_t const DefaultBlockMultiple = 16777216;  // 2^24 = 16MB

class BlockAllocator : public blas::AllocatorBase
{
   public:
      BlockAllocator(std::size_t Mult = DefaultBlockMultiple, bool Free = true) : BlockMultiple(Mult), FreeOnDestructor(Free) {}

      void* allocate(std::size_t Size, std::size_t Align);

      void* allocate(std::size_t Size);

      void free(void* Ptr, std::size_t Size);

      virtual ~BlockAllocator() noexcept {}

   private:
      std::list<AllocationBlock> Allocations;

      std::size_t BlockMultiple;
      bool FreeOnDestructor;
};

inline
blas::arena make_gpu_block_allocator()
{
   return blas::arena(new BlockAllocator());
}

// The destructor of gpu_buffer must block until the stream is finshed.
// When run in async mode we want the destructor of an async gpu matrix to run as a task.
// The way to handle this is to provide a way to move the memory buffer out of the
// gpu_buffer so that it can be destroyed asynchronously.

// Walk the GpuBufferPendingDelete list and deallocate buffers associated with any streams
// that have finished execution.  Non-blocking.
void TryFlushGpuBuffers();

// Force synchronization and memory deallocation of all streams
void SyncFlushGpuBuffers();

//
void AddToPendingDelete(cuda::stream& Stream, blas::arena& Arena, void* Ptr, std::size_t ByteSize);

template <typename T>
class gpu_ptr;

template <typename T>
class const_gpu_ptr;

template <typename T>
class gpu_ref;

template <typename T>
class gpu_buffer
{
   public:
      using storage_type       = gpu_ptr<T>;
      using const_storage_type = const_gpu_ptr<T>;
      using reference          = gpu_ref<T>;
      using const_reference    = gpu_ref<T>;

      gpu_buffer() noexcept
      : Ptr(nullptr), ByteSize(0) {}

      // Async version
      ~gpu_buffer() noexcept;

      gpu_buffer(gpu_buffer&& other) noexcept
      : Ptr(other.Ptr), ByteSize(other.ByteSize), Stream(std::move(other.Stream)),
	Sync(std::move(other.Sync)), Arena(std::move(other.Arena))
      {
         other.Ptr = nullptr;
      }

      gpu_buffer(gpu_buffer const& other) = delete;

      gpu_buffer& operator=(gpu_buffer&&) = delete;
      gpu_buffer& operator=(gpu_buffer const&) = delete;

      gpu_ref<T> operator[](int n);

      gpu_ptr<T> ptr();
      const_gpu_ptr<T> ptr() const;
      const_gpu_ptr<T> cptr() const;

      gpu_ptr<T> ptr(int Offset);
      const_gpu_ptr<T> ptr(int Offset) const;
      const_gpu_ptr<T> cptr(int Offset) const;

      std::size_t size() const { return ByteSize / sizeof(T); }

      std::size_t byte_size() const { return ByteSize; }

      // block modifications to this buffer until event e has been triggered
      void wait(event const& e) const
      {
         Stream.wait(e);
      }

      template <typename U>
      void wait_for(gpu_ptr<U> const& Other) const;

      template <typename U>
      void wait_for(const_gpu_ptr<U> const& Other) const;

      template <typename U>
      void wait_for(gpu_ref<U> const& Other) const;

      template <typename U>
      void wait_for(gpu_buffer<U> const& Other) const
      {
         this->wait(Other.sync());
         // Generally, calling wait_for() is going to be done just before
         // doing some operation that writes to the buffer, so an invalidate()
         // will follow shortly.  Note that we could force a sync() operation
         // before waiting - currently we don't force a sync() so that means
         // that if something else is going to wait_for() this buffer, then
         // they should do that BEFORE calling this->wait_for(Other).
         // Correct:
         //    A.wait_for(B);
         //    B.wait_for(C);
         // Inefficient (causes A to block for longer than it needs to):
         //    B.wait_for(C);
         //    A.wait_for(B);
      }

      event sync() const
      {
         if (Sync.is_null())
            Sync = Stream.record();
         return Sync;
      }

      // resets the CUDA event associated with the stream
      void invalidate()
      {
         Sync.clear();
      }

      T* device_ptr() { this->invalidate(); return Ptr; }
      T const* device_ptr() const { return Ptr; }

      cuda::stream const& get_stream() const { return Stream; }

      static gpu_buffer allocate(std::size_t Size, blas::arena A)
      {
	 return gpu_buffer(Size, A);
      }

   private:
      gpu_buffer(int Count, blas::arena a)
	 : ByteSize(Count*sizeof(T)), Arena(a)
      {
	 Ptr = static_cast<T*>(Arena.allocate(ByteSize, sizeof(T)));
	 if (!Ptr)
	    throw std::runtime_error("No GPU memory");
      }

      friend class gpu_ref<T>;

      T* Ptr;
      std::size_t ByteSize;
      mutable cuda::stream Stream;
      mutable cuda::event Sync;
      blas::arena Arena;
};

// a gpu_ptr is a weak version of a gpu_buffer - it contains a stream and a device pointer.
// Use when the access to the pointer happens via buffer backing stream, ie no effective
// parallelization over different components of the buffer.

template <typename T>
class const_gpu_ptr
{
   public:
      const_gpu_ptr() = delete;

      const_gpu_ptr(gpu_buffer<T> const& Buf_, int Offset_) : Buf(Buf_), Offset(Offset_) {}

      ~const_gpu_ptr() = default;

      void wait(event const& e) const
      {
         Buf.wait(e);
      }

      template <typename U>
      void wait_for(gpu_buffer<U> const& Other) const
      {
         Buf.wait_for(Other);
      }

      template <typename U>
      void wait_for(gpu_ptr<U> const& Other) const
      {
         Buf.wait_for(Other);
      }

      template <typename U>
      void wait_for(const_gpu_ptr<U> const& Other) const
      {
         Buf.wait_for(Other);
      }

      event sync() const
      {
         return Buf.sync();
      }

      T const* device_ptr() const { return Buf.device_ptr() + Offset; }

      cuda::stream const& get_stream() const { return Buf.get_stream(); }

      gpu_buffer<T> const& buffer() const { return Buf; }
      int offset() const { return Offset; }

   private:
      gpu_buffer<T> const& Buf;
      int Offset;
};

template <typename T>
const_gpu_ptr<T> operator+(const_gpu_ptr<T> const& x, int Offset)
{
   return const_gpu_ptr<T>(x.buffer(), x.offset() + Offset);
}

template <typename T>
class gpu_ptr
{
   public:
      gpu_ptr() = delete;

      gpu_ptr(gpu_buffer<T>& Buf_, int Offset_) : Buf(Buf_), Offset(Offset_) {}

      ~gpu_ptr() = default;

      operator const_gpu_ptr<T>() const { return const_gpu_ptr<T>(Buf, Offset); }

      void wait(event const& e)
      {
         Buf.wait(e);
      }

      template <typename U>
      void wait_for(gpu_buffer<U> const& Other)
      {
         Buf.wait_for(Other);
      }

      template <typename U>
      void wait_for(gpu_ptr<U> const& Other)
      {
         Buf.wait_for(Other);
      }

      template <typename U>
      void wait_for(const_gpu_ptr<U> const& Other)
      {
         Buf.wait_for(Other);
      }

      event sync() const
      {
         return Buf.sync();
      }

      void invalidate()
      {
         Buf.invalidate();
      }

      T const* device_ptr() const { return Buf.device_ptr() + Offset; }
      T* device_ptr() { return Buf.device_ptr() + Offset; }

      cuda::stream const& get_stream() const { return Buf.get_stream(); }

      gpu_buffer<T>& buffer() const { return Buf; }
      int offset() const { return Offset; }

   private:
      gpu_buffer<T>& Buf;
      int Offset;
};

template <typename T>
gpu_ptr<T> operator+(gpu_ptr<T> const& x, int Offset)
{
   return gpu_ptr<T>(x.buffer(), x.offset() + Offset);
}

// a reference to a location within a gpu_buffer
template <typename T>
class gpu_ref
{
   public:
      using remove_proxy_t = T;

      gpu_ref(); // = delete;

      explicit gpu_ref(T* Ptr_, stream* ParentStream_ = nullptr)
         : Ptr(Ptr_), ParentStream(ParentStream_)
      {
         if (ParentStream)
            Stream.wait(ParentStream->record());
      }

      ~gpu_ref()
      {
         // If we have a parent stream then we don't need to synchronize our stream on destruction,
         // since the parent will wait for us.  In that case, we want to destroy our stream
         // without invoking the warning that we are destroying a stream with pending operations.
         if (ParentStream)
         {
            if (Sync.is_null())
               Sync = Stream.record();
            ParentStream->wait(Sync);
            Stream.safe_to_destroy();
         }
         else
         {
            Stream.synchronize();
         }
      }

      gpu_ref(gpu_ref&& other) : Ptr(other.Ptr), ParentStream(other.ParentStream),
                                 Stream(std::move(other.Stream)),
                                 Sync(std::move(other.Sync))
      {
         other.ParentStream = nullptr;
      }

      gpu_ref(gpu_ref const& Other) = delete;

      gpu_ref& operator=(gpu_ref<T>&&) = delete;

      gpu_ref& operator=(gpu_ref<T> const& Other)
      {
	 // wait until the other stream has finished writing
	 Stream.wait(Other.sync());
	 // do the copy
	 memcpy_device_to_device_async(Stream, Other.Ptr, Ptr, sizeof(T));
	 // signal our event that we have finished writing
	 Sync = Stream.record();
	 // make the other stream wait until the copy is finished
	 Other.wait(Sync);
	 return *this;
      }

      event sync() const { if (Sync.is_null()) Sync = Stream.record(); return Sync; }

      void wait(event const& e) const { Stream.wait(e); }

      void set_wait(T const& x)
      {
         Stream.synchronize();
         Sync.clear();
	 memcpy_host_to_device(&x, Ptr, sizeof(T));
      }

      // sets the element asynchronously.
      void set(T const& x)
      {
         Sync.clear();
	 memcpy_host_to_device_async(Stream, &x, Ptr, sizeof(T));
      }

      // blocking get operation
      T get_wait() const
      {
	 T x;
         Stream.synchronize();
	 memcpy_device_to_host(Ptr, &x, sizeof(T));
	 return x;
      }

      // non-blocking get.  Returns a cuda::event, when triggered the host memory location is valid.
      cuda::event get(T* x)
      {
	 cuda::event E;
	 memcpy_device_to_host_async(Stream, Ptr, x, sizeof(T));
	 Sync = Stream.record();
	 return Sync;
      }

      T* device_ptr() { return Ptr; }
      cuda::stream const& get_stream() const { return Stream; }
      cuda::stream& get_stream() { return Stream; }

   private:
      T* Ptr;
      cuda::stream* ParentStream;
      mutable cuda::stream Stream;
      mutable cuda::event Sync;
};

template <typename T>
inline
T
get_wait(gpu_ref<T> const& x)
{
   return x.get_wait();
}

// function to return a newly allocated gpu_ref
template <typename T>
gpu_ref<T>
allocate_gpu_ref()
{
   static gpu_buffer<T> Buf = gpu_buffer<T>::allocate(1, make_gpu_block_allocator());
   return Buf[0];
}

template <typename T>
inline
gpu_ref<T>::gpu_ref()
   : gpu_ref(allocate_gpu_ref<T>())
{
}

template <typename T>
inline
gpu_ref<T>
gpu_buffer<T>::operator[](int n)
{
   return gpu_ref<T>(Ptr+n, &Stream);
}

//
// helper functions for temporary allocations
//

void* allocate_gpu_temp_memory(int Size);

void free_gpu_temp_memory(void* Buf, int Size);

template <typename T>
gpu_buffer<T>
allocate_gpu_temporary(int Size);

} // namsepace cuda

#include "gpu_buffer.icc"

#endif
