// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/stackallocator.cpp
//
// Copyright (C) 2012-2016 Ian McCulloch <ian@qusim.net>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

#include "stackallocator.h"
#include "trace.h"
#include <utility>
#include <memory>
#include <stack>

namespace StackAlloc
{

size_t const BlockSize = 1024 * 1024 * 16;

typedef std::stack<void*> BlockStackType;
typedef std::stack<std::pair<void*, void*> > UsedBlockStackType;
BlockStackType* FreeBlocks(NULL);
UsedBlockStackType* UsedBlocks(NULL);
void* Block = NULL;
void* Current = NULL;

#if !defined(NDEBUG)
std::stack<std::pair<void*, size_t> >* DebugTrace;
#endif

namespace detail
{

void StackAllocatorInit()
{
   FreeBlocks = new BlockStackType;
   UsedBlocks = new UsedBlockStackType;
   Block = ::operator new(BlockSize);
   Current = Block;

#if !defined(NDEBUG)
   DebugTrace = new std::stack<std::pair<void*, size_t> >();
#endif
}

void StackAllocatorExit()
{
  // lots of other stuff leaks....
  //   delete FreeBlocks;
  //   delete UsedBlocks;
}

} // namespace detail

inline size_t RoundSize(size_t Size)
{
   return (Size == 0) ? 1 : ((Size-1) / MinAlign + 1) * MinAlign;
}

void* allocate(size_t Size)
{
   if (Size > BlockSize) return ::operator new(Size);

   Size = RoundSize(Size);

   if (Size > BlockSize - (static_cast<char*>(Current) - static_cast<char*>(Block)))
   {
      UsedBlocks->push(std::make_pair(Block, Current));
      if (FreeBlocks->empty())
      {
         Current = Block = ::operator new(BlockSize);
      }
      else
      {
         Current = Block = FreeBlocks->top();
         FreeBlocks->pop();
      }
   }

   void* Ret = Current;
   Current = static_cast<char*>(Current) + Size;

#if !defined(NDEBUG)
   (*DebugTrace).push(std::make_pair(__builtin_return_address(0), Size));
#endif

   return Ret;
}

void deallocate(void* ptr, size_t Size)
{
   if (Size > BlockSize)
   {
      ::operator delete(ptr);
      return;
   }

   Size = RoundSize(Size);

   if (Current == Block)
   {
      FreeBlocks->push(Block);
      std::tie(Block, Current) = UsedBlocks->top();
      UsedBlocks->pop();
   }

   Current = static_cast<char*>(Current) - Size;
#if !defined(NDEBUG)
   CHECK(Current == ptr)(Current)(ptr)(Size)(DebugTrace->top().first)(DebugTrace->top().second);
   DebugTrace->pop();
#else
   CHECK(Current == ptr)(Current)(ptr)(Size);
#endif
}

} // namespace
