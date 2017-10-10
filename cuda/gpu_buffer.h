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

      explicit AllocationBlock(std::size_t Size_, std::size_t Align_, bool FreeOnDest = true);

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

      ~AllocationBlock();

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
AllocationBlock::AllocationBlock(std::size_t Size_, std::size_t Align_, bool Free)
   : Align(Align_), Size(Size_), CurrentOffset(0), NumAllocations(0), FreeOnDestructor(Free)
{
   void* Ptr;
   check_error(cudaMalloc(&Ptr, Size_));
   BasePtr = static_cast<unsigned char*>(Ptr);
}

inline
AllocationBlock::~AllocationBlock()
{
   if (FreeOnDestructor && BasePtr)
   {
      TRACE("Freeing");
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

   private:
      std::list<AllocationBlock> Allocations;

      std::size_t BlockMultiple;
      bool FreeOnDestructor;
};

inline
void* BlockAllocator::allocate(std::size_t Size, std::size_t Align)
{
   // walk the existing allocations and try to find somewhere that fits
   for (auto& a : Allocations)
   {
      void* p = a.try_allocate(Size, Align);
      if (p)
	 return p;
   }
   // failed to allocate, need another block
   size_t BlockSize = round_up(Size, BlockMultiple);
   Allocations.push_back(AllocationBlock(BlockSize, Align, FreeOnDestructor));
   void* p = Allocations.back().try_allocate(Size, Align);
   return p;
}

inline
void* BlockAllocator::allocate(std::size_t Size)
{
   return this->allocate(Size, 1);
}

inline
void BlockAllocator::free(void* Ptr, std::size_t Size)
{
   for (auto& a : Allocations)
   {
      if (a.try_free(Ptr, Size))
	 return;
   }
   PANIC("Attempt to free memory that wasn't allocated by this allocator!")(Ptr);
}

inline
blas::arena make_gpu_block_allocator()
{
   return blas::arena(new BlockAllocator());
}

template <typename T>
class gpu_ptr;

template <typename T>
class const_gpu_ptr;

template <typename T>
class gpu_buffer
{
   public:
      gpu_buffer() = delete;

      ~gpu_buffer() { Arena.free(Ptr, ByteSize); }

      gpu_buffer(gpu_buffer&& other) = default;
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

      // block modifications to this buffer until event e has been triggered
      void wait(event const& e)
      {
         Stream.wait(e);
      }

      template <typename U>
      void wait_for(gpu_ptr<U> const& Other);

      template <typename U>
      void wait_for(const_gpu_ptr<U> const& Other);

      template <typename U>
      void wait_for(gpu_buffer<U> const& Other)
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

      event sync() const
      {
         return Buf.sync();
      }

      T const* device_ptr() const { return Buf.device_ptr() + Offset; }

      cuda::stream const& get_stream() const { return Buf.get_stream(); }

   private:
      gpu_buffer<T> const& Buf;
      int Offset;
};

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

   private:
      gpu_buffer<T>& Buf;
      int Offset;
};

// a reference to a location within a gpu_buffer
template <typename T>
class gpu_ref
{
   public:
      gpu_ref() = delete;

      gpu_ref(T* Ptr_, stream* ParentStream_ = nullptr)
         : Ptr(Ptr_), ParentStream(ParentStream_)
      {
         if (ParentStream)
            Stream.wait(ParentStream->record());
      }

      ~gpu_ref()
      {
         if (ParentStream)
         {
            if (Sync.is_null())
               Sync = Stream.record();
            ParentStream->wait(Sync);
         }
      }

      gpu_ref(gpu_ref&& other) : Ptr(other.Ptr), ParentStream(other.ParentStream), Stream(std::move(other.Stream)),
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
         Sync.clear();
	 memcpy_host_to_device(Stream, &x, Ptr, sizeof(T));
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
   static gpu_buffer<T> Buf = gpu_buffer<T>::allocate(100, make_gpu_block_allocator());
   return Buf[0];
}

template <typename T>
inline
gpu_ref<T>
gpu_buffer<T>::operator[](int n)
{
   return gpu_ref<T>(Ptr+n, &Stream);
}

} // namsepace cuda

#include "gpu_buffer.icc"

#endif
