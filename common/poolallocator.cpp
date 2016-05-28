// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/poolallocator.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if defined(HACK_MULTITHREAD) && !defined(MULTITHREAD)
#define MULTITHREAD
#endif

#if defined(HAVE_CONFIG_H)
#include "config.h"         // for configuration options
#endif
#include "poolallocator.h"
#include "mutex.h"
#include <iostream>
#include <iomanip>
#include <set>

namespace
{

// the free list.  We allocate +1 more so we can use it as a 1-based array
// BlockRec stores a linked list of allocated blocks.
void* FreeList[PoolAlloc::MaxBlockAllocatorUnits+1] = {NULL};
void* BlockRecHead[PoolAlloc::MaxBlockAllocatorUnits+1] = {NULL};
void* BlockRecTail[PoolAlloc::MaxBlockAllocatorUnits+1] = {NULL};
pthread::mutex BlockMutex[PoolAlloc::MaxBlockAllocatorUnits+1];

// returns the UnitSize of the corresponding actual size
inline size_t GetUnitSize(size_t size)
{
   return (size == 0) ? 1 : (size-1) / PoolAlloc::MinAlign + 1;
}

// allocates a new pool of memory
void* AllocatePool(size_t UnitSize, size_t NAlloc = PoolAlloc::DefaultAllocationUnits)
{
   size_t AllocSize = UnitSize * PoolAlloc::MinAlign;

   // Initialize the data into a free list, we allocate MinAlign*2 extra bytes,
   // for the bock list and we also store NAlloc.
   char* Array = static_cast<char*>(::operator new(NAlloc*AllocSize + PoolAlloc::MinAlign*2));
   for (int i = 0; i < int(NAlloc)-1; ++i)
   {
      *static_cast<void**>(static_cast<void*>(Array + i*AllocSize +  PoolAlloc::MinAlign*2)) 
	= (Array + (i+1)*AllocSize + PoolAlloc::MinAlign*2);
   }
   *static_cast<void**>(static_cast<void*>(Array + (NAlloc-1)*AllocSize +  PoolAlloc::MinAlign*2)) = NULL;

   // set the first pointer as an index into the blocks
   *static_cast<void**>(static_cast<void*>(Array)) = NULL;
   *reinterpret_cast<size_t*>(Array + PoolAlloc::MinAlign) = NAlloc;
   if (BlockRecTail[UnitSize] != NULL) *static_cast<void**>(BlockRecTail[UnitSize]) = Array;
   BlockRecTail[UnitSize] = Array;
   if (BlockRecHead[UnitSize] == NULL) BlockRecHead[UnitSize] = Array;

   return Array +  PoolAlloc::MinAlign*2;
}

// returns true if p is a valid address of a block of data of size UnitSize.
bool IsValidAddress(void* p, size_t UnitSize)
{ 
   void* Ptr = BlockRecHead[UnitSize];
   while (Ptr)
   {
      char* First = static_cast<char*>(Ptr) +  PoolAlloc::MinAlign * 2;
      int NAlloc = *reinterpret_cast<size_t*>(static_cast<char*>(Ptr) +  PoolAlloc::MinAlign);

      if (First <= p && p < First + NAlloc * UnitSize * PoolAlloc::MinAlign)
      {
	 ptrdiff_t Offset = static_cast<char*>(p) - First;
	 return Offset % (UnitSize * PoolAlloc::MinAlign) == 0;
      }
      Ptr = *static_cast<void**>(Ptr);
   }
   return false;
}

void CheckForAllocatedBlocks(size_t UnitSize, std::set<void*> const& Free)
{ 
   void* Ptr = BlockRecHead[UnitSize];
   while (Ptr)
   {
      char* First = static_cast<char*>(Ptr) +  PoolAlloc::MinAlign * 2;
      int NAlloc = *reinterpret_cast<size_t*>(static_cast<char*>(Ptr) +  PoolAlloc::MinAlign);
      for (int i = 0; i < NAlloc; ++i)
      {
	 if (Free.count(First) == 0) 
	 {
	    std::cerr << "\n      Block size " << (UnitSize * PoolAlloc::MinAlign) << " at " << (void*)First << " has leaked!";
	 }
	 First += UnitSize * PoolAlloc::MinAlign;
      }
      Ptr = *static_cast<void**>(Ptr);
   }
}

// verifies that a pointer is valid in preparation for freeing it.  
// If it is already on the free list, it is reported as being freed twice.
void DebugCheck(void* p, size_t size)
{
   size_t UnitSize = GetUnitSize(size);
   void* Ptr = FreeList[UnitSize];
   while (Ptr != NULL && Ptr != p)
   {
      Ptr = *static_cast<void**>(Ptr);
   }
   if (Ptr == p)
   {
     std::cerr << "PoolAlloc: Block of " << size << " bytes at " << p << " was freed twice!" << std::endl;
     // walk the heap then seg fault
     PoolAlloc::walk_heap();
     PANIC("PoolAlloc: double deallocation of a block")(size)(p);
   }

   // also check that p is a valid pointer
   if (!IsValidAddress(p, UnitSize))
   {
     std::cerr << "PoolAlloc: Attempted deallocate of block of " << size << " bytes at " << p 
	       << " failed - block address is not valid!" << std::endl;
     // walk the heap then seg fault
     PoolAlloc::walk_heap();
     PANIC("PoolAlloc: double deallocation of a block")(size)(p);
   }
}

// Walks the heap for the given UnitSize.  Returns false if 
// there are allocated memory blocks, true if there are no allocated memory blocks.
// If Verbose is true, detailed debug info is written to std::cerr
// If AssumeLeaked is true, then any allocated memory blocks are reported as memory leaks. 
bool WalkHeap(size_t UnitSize, bool Verbose = false, bool AssumeLeaked = true)
{
   // calculate how many total blocks exist of size unitSize
   int CountUsed = 0;
   void* Ptr = BlockRecHead[UnitSize];
   if (Ptr && Verbose > 0)
   {
      std::cerr << "Memory regions in use for block size " << UnitSize * PoolAlloc::MinAlign << '\n';
   }
   while (Ptr)
   {
      int ThisSize = *reinterpret_cast<size_t*>(static_cast<char*>(Ptr) + PoolAlloc::MinAlign);
      if (Verbose > 0) std::cerr << "  Block starting " << Ptr << ", " << ThisSize << " records, " 
		<< (ThisSize * UnitSize * PoolAlloc::MinAlign + PoolAlloc::MinAlign*2) << " bytes\n";
       CountUsed += ThisSize;
      Ptr = *static_cast<void**>(Ptr);
   }

   if (CountUsed == 0) return true;

   bool HeapOk = true; // flag is set false if we detect leaked memory blocks

   // calculate how many free blocks there are
   std::set<void*> FreePointers;
   int CountFree = 0;
   Ptr = FreeList[UnitSize];
   while (Ptr)
   {
      ++CountFree;
      FreePointers.insert(Ptr);
      Ptr = *static_cast<void**>(Ptr);
   }

   if (CountFree != CountUsed) HeapOk = false;

   if (Verbose || CountFree != CountUsed) 
   {
     std::cerr << "Block size " << (UnitSize * PoolAlloc::MinAlign) << ", total=" << CountUsed
	     << ", free=" << CountFree << ", allocated=" << (CountUsed - CountFree);
      if (AssumeLeaked && CountFree != CountUsed)
      {
	 std::cerr << " LEAKED!";
	 // find out which pointers are leaked
	 CheckForAllocatedBlocks(UnitSize, FreePointers);
      }
      std::cerr << std::endl;
   }

   return HeapOk;
}

} // namespace

namespace PoolAlloc
{

void* allocate(size_t size)
{
   size_t UnitSize = GetUnitSize(size);

   PRECONDITION(UnitSize < MaxBlockAllocatorUnits);

   pthread::mutex::sentry MyLock(BlockMutex[UnitSize]);

   if (!(FreeList[UnitSize])) FreeList[UnitSize] = AllocatePool(UnitSize, 
								PoolAlloc::DefaultAllocationUnits/UnitSize+1);

   void* Ret = FreeList[UnitSize];
   FreeList[UnitSize] = *static_cast<void**>(FreeList[UnitSize]);

   //   std::cout << "Allocated " << size << " at " << (void*) Ret << '\n';

   return Ret;
}

void deallocate(void* p, size_t size)
{
  //   std::cout << "Deallocating " << size << " at " << (void*) p << '\n';

   if (!p) return;
   size_t UnitSize = GetUnitSize(size);

#if defined(POOLALLOCATOR_DEBUG)
   DebugCheck(p, size);
#endif

   *static_cast<void**>(p) = FreeList[UnitSize];
   FreeList[UnitSize] = p;
}

void walk_heap()
{
   std::cerr << "PoolAllocator heap status:\n";
   for (int i = 1; i < int(MaxBlockAllocatorUnits); ++i)
   {
      WalkHeap(i, true, false);
   }
}

void PoolAllocatorInit()
{
   // nothing to do here
}

void PoolAllocatorExit()
{
  // only bother freeing the memory if NDEBUG is not defined,
   // if it is, we just let it leak.
#if defined(POOL_ALLOCATOR_LEAK_CHECK)
#if defined(POOL_ALLOCATOR_VERBOSE)
   bool const Verbose = true;
#else
   bool const Verbose = false;
#endif

   for (int i = 1; i < MaxBlockAllocatorUnits; ++i)
   {
       // see if the heap checks out OK, and if it does, free the blocks
       if (WalkHeap(i, Verbose))
       {
	  void* Ptr = BlockRecHead[i];
	  while (Ptr)
	  {
	     void* PtrSave = Ptr;
	     Ptr = *static_cast<void**>(Ptr);
	     ::operator delete(static_cast<char*>(PtrSave));
	  }
       }
   }
#endif
}

} // namespace PoolAlloc
