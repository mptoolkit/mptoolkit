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

namespace cuda
{

template <typename T>
class gpu_ref;

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

      explicit AllocationBlock(std::size_t Size_, std::size_t Align_);

      AllocationBlock(AllocationBlock&& Other) : Align(Other.Align),
						 Size(Other.Size),
						 BasePtr(Other.BasePtr),
						 CurrentOffset(Other.CurrentOffset),
						 NumAllocations(Other.NumAllocations)
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

	 CurrentOffset = NextOffset + RequestSize;
	 return static_cast<void*>(BasePtr + NextOffset);
      }

      // free memory that was previously allocated.
      // Returns true if this caused the allocation block to become empty
      bool free(void* Ptr, size_t Size)
      {
	 // TOOD: add some debugging stuff
	 if (--NumAllocations == 0)
	 {
	    CurrentOffset = 0;
	    return true;
	 }
      }

      bool try_free(void* Ptr, size_t Size)
      {
	 if (static_cast<unsigned char*>(Ptr) < BasePtr)
	    return false;

	 if (BasePtr+Size < static_cast<unsigned char*>(Ptr))
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
};

AllocationBlock::AllocationBlock(std::size_t Size_, std::size_t Align_)
   : Align(Align_), Size(Size_), CurrentOffset(0), NumAllocations(0)
{
   void* Ptr;
   check_error(cudaMalloc(&Ptr, Size_));
   BasePtr = static_cast<unsigned char*>(Ptr);
}

AllocationBlock::~AllocationBlock()
{
   if (BasePtr)
      check_error(cudaFree(static_cast<void*>(BasePtr)));
}

std::size_t const DefaultBlockMultiple = 16777216;  // 2^24 = 16MB

class AllocatorBase
{
   public:
      AllocatorBase() {}

      virtual void* allocate(std::size_t Size, std::size_t Align) = 0;

      virtual void* allocate(std::size_t Size) = 0;

      virtual void free(void* Ptr, std::size_t Size) = 0;

      virtual ~AllocatorBase() = 0;
};

inline
AllocatorBase::~AllocatorBase()
{
}

class BlockAllocator : public AllocatorBase
{
   public:
      BlockAllocator(std::size_t Mult = DefaultBlockMultiple) : BlockMultiple(Mult) {}

      void* allocate(std::size_t Size, std::size_t Align);

      void* allocate(std::size_t Size);

      void free(void* Ptr, std::size_t Size);

   private:
      std::list<AllocationBlock> Allocations;

      std::size_t BlockMultiple;
      std::size_t MinAlign;
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
   Allocations.push_back(AllocationBlock(BlockSize, Align));
   void* p = Allocations.back().try_allocate(Size, Align);
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
}

class arena
{
   public:
      arena() {}

      explicit arena(std::shared_ptr<AllocatorBase> Alloc_) : Alloc(Alloc_) {}

      explicit arena(AllocatorBase* Alloc_) : Alloc(Alloc_) {}

      void* allocate(std::size_t Size, std::size_t Align)
      { return Alloc->allocate(Size, Align); }

      void* allocate(std::size_t Size)
      { return Alloc->allocate(Size); }

      void free(void* Ptr, std::size_t Size)
      {
	 Alloc->free(Ptr, Size);
      }

   private:
      std::shared_ptr<AllocatorBase> Alloc;
};

inline
arena get_block_allocator()
{
   return arena(new BlockAllocator());
}

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

      // block modifications to this buffer until event e has been triggered
      void wait(event const& e);

      T* device_ptr() { return Ptr; }
      T const* device_ptr() const { return Ptr; }

      static gpu_buffer allocate(std::size_t Size, arena A)
      {
	 return gpu_buffer(Size, A);
      }

   private:
      gpu_buffer(int Count, arena a)
	 : ByteSize(Count*sizeof(T)), Arena(a)
      {
	 Ptr = static_cast<T*>(Arena.allocate(ByteSize, sizeof(T)));
	 if (!Ptr)
	    throw std::runtime_error("No GPU memory");
      }

      friend class gpu_ref<T>;

      T* Ptr;
      std::size_t ByteSize;
      cuda::stream Stream;
      cuda::event Sync;
      arena Arena;
};

// a reference to a location within a gpu_buffer
template <typename T>
class gpu_ref
{
   public:
      gpu_ref() = delete;

      ~gpu_ref()
      {
	 Sync.record(Stream);
	 Buf->wait(Sync);
	 Buf->remove(this);
      }

      gpu_ref(gpu_ref&& other) { using std::swap; swap(Buf, other.Buf); swap(Ptr, other.Ptr); swap(Stream, other.Stream);
	 swap(Sync, other.Sync); }

      gpu_ref(gpu_ref const& Other) = delete;

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

      void set_wait(T const& x)
      {
	 memcpy_host_to_device(Stream, &x, Ptr, sizeof(T));
	 Sync.record(Stream);
      }

      // sets the element asynchronously.
      void set(T const& x)
      {
	 T const* Ptr = get_global_constant(x);
	 memcpy_host_to_device_async(Stream, &x, Ptr, sizeof(T));
	 Sync.record(Stream);
      }

      // blocking get operation
      T get_wait()
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
	 Sync.record(Stream);
	 return E;
      }

   private:
      friend gpu_ref<T> gpu_buffer<T>::operator[](int n);

      gpu_ref(gpu_buffer<T>* Buf_, T* Ptr_) : Buf(Buf_), Ptr(Ptr_) { }

      gpu_buffer<T>* Buf;
      T* Ptr;
      AtomicRefCount Count;
      cuda::stream Stream;
      cuda::event Sync;
};

} // namsepace cuda

#endif
