// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pheapfsblock.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  Finger grained allocator for the page files.

  Created 2002-07-16 Ian McCulloch

  This is analagous to the old subpageallocator.h, although this version is much easier to use.

  To use:

  Obtain an OutputBuffer* by calling BlockFileSystem::allocate_output_buffer().
  This gives you a memory buffer that you can directly access.  Remember to use bump(N)
  to increment the buf_ptr by N bytes.  The size of the buffer is unpsecified, it may
  be smaller than the page size.

  Once data has been written to the buffer, call OutputBuffer::commit().  This returns
  a BlockRecord that can be used to read back the data later.  The OutputBuffer then
  points after the committed data and can be used for another allocation, if required,
  (But the remaining OutputBuffer might have zero size).

  Once you have finished with a (committed) OutputBuffer, return it to the system with
  BlockFileSystem::return_output_buffer().

  An additional function, BlockFileSystem::cycle_output_buffer() is provided for the
  case that you want to return an exhausted OutputBuffer and immediately allocate another.
*/

#if !defined(PHEAPFSBLOCK_H_JDFHUIERTY78TY78HFIUH35789GHWE87O)
#define PHEAPFSBLOCK_H_JDFHUIERTY78TY78HFIUH35789GHWE87O

#include "pheapfsv4.h"
#include "common/trace.h"

namespace PHeapFileSystem
{

class InputBuffer;
class OutputBuffer;
class BlockRecord;

class BlockFileSystem : public FileSystem
{
   public:
      BlockFileSystem();
      ~BlockFileSystem();

      OutputBuffer* allocate_output_buffer();

      // returns an output buffer that has been obtained previously from
      // allocate_output_buffer() to the system.  It is invalid to call this
      // function with a buffer that has uncommited data (a buffer is uncommited
      // when the buf_ptr() is larger than the begin_ptr()).
      void return_output_buffer(OutputBuffer* Buf);

      // A combination of return and allocate.  Returns Buf, and returns another
      // allocated buf.
      OutputBuffer* cycle_output_buffer(OutputBuffer* Buf);

      // a shortcut for Temp = new BlockFileSystem(); Temp->initialize(...); return Temp;
      //      static BlockFileSystem* Create(std::string const& FileName, int NumFiles,
      //                                     size_t PageSize, size_t PageCacheByteSize);

      // a shortcut for Temp = new BlockFileSystem(); Temp->restart(...); return Temp;
//      static BlockFileSystem* Open(std::string const& FileName, size_t PageCacheByteSize, bool ReadOnly = false);

      void WriteBlockRecord(PStream::opstream& out, BlockRecord const& Rec);

      BlockRecord ReadBlockRecord(PStream::ipstream& in);

   private:
      pthread::mutex BufferListMutex;
      std::vector<OutputBuffer*> Buffers;
};

#if 0
inline
BlockFileSystem* BlockFileSystem::Create(std::string const& FileName, int NumFiles,
                                         size_t PageSize, size_t PageCacheByteSize,
                                         bool Unlink)
{
   BlockFileSystem* FS = new BlockFileSystem();
   FS->create(FileName, NumFiles, PageSize, PageCacheByteSize, Unlink);
   return FS;
}

inline
BlockFileSystem* BlockFileSystem::Open(std::string const& FileName,
                                       size_t PageCacheByteSize,
                                       bool ReadOnly)
{
   BlockFileSystem* FS = new BlockFileSystem();
   FS->open(FileName, PageCacheByteSize, ReadOnly);
   return FS;
}
#endif

// BlockRecord controls the page reference count
class BlockRecord
{
   public:
      BlockRecord();
      // Compiler generated copy ctor, copy assignment are OK

      // returns the number of bytes in the block.
      size_t size() const { return Length; }

   private:
      BlockRecord(PageId const& Page_, size_t Offset_, size_t Length_);

      PageId Page;
      size_t Offset;
      size_t Length;

   friend class OutputBuffer;
   friend class InputBuffer;
   friend class BlockFileSystem;

   friend std::ostream& operator<<(std::ostream& out, BlockRecord const& r);
};

inline
std::ostream& operator<<(std::ostream& out, BlockRecord const& r)
{
   return out << "BlockRecord(Page=" << r.Page << ", Offset=" << r.Offset
              << ", Length=" << r.Length << ")";
}

// these control the page lock count
class OutputBuffer
{
   public:
      OutputBuffer();
      // constructs an output buffer from a page buffer.
      explicit OutputBuffer(WriteBuffer const& Buf);

      // returns the start of the buffer
      unsigned char* buf_begin() const { return BeginPtr; }

      // returns the one-past-end of the buffer
      unsigned char* buf_end() const { return EndPtr; }

      // returns the current pointer into the buffer
      unsigned char* buf_ptr() const { return CurrentPtr; }

      // returns the number of bytes remaining in the buffer
      size_t remain() const { return EndPtr - CurrentPtr; }

      // Increments the curent pointer by Count
      void bump(size_t Count);

      // Forces the current pointer to Ptr.  Ptr must be in the range buf_begin() <= Ptr <= buf_end()
      void set_buf_ptr(unsigned char* Ptr)
      { DEBUG_RANGE_CHECK(Ptr, buf_begin(), buf_end()); CurrentPtr = Ptr; }

      // Commits a buffer and returns the associated block record.
      // This buffer then points AFTER the committed buffer, in the same page.
      BlockRecord commit();

   private:
      ~OutputBuffer();

      OutputBuffer(OutputBuffer const& Buf);              // not implemented
      OutputBuffer& operator=(OutputBuffer const& Buf);   // not implemented

      WriteBuffer WB;
      unsigned char* BeginPtr;
      unsigned char* CurrentPtr;
      unsigned char* EndPtr;

   friend class BlockFileSystem;
};

class InputBuffer
{
   public:
      InputBuffer();
      explicit InputBuffer(BlockRecord const& R);

      // compiler-generated copy ctor and copy assignment are OK

      // returns the start of the buffer
      unsigned char const* buf_begin() const { return BeginPtr; }

      // returns the one-past-end of the buffer
      unsigned char const* buf_end() const { return EndPtr; }

      // returns the current pointer into the buffer
      unsigned char const* buf_ptr() const { return CurrentPtr; }

      // returns the number of bytes remaining in the buffer
      size_t remain() const { return EndPtr - CurrentPtr; }

      // Increments the curent pointer by Count
      void bump(size_t Count) { CurrentPtr += Count; }

      // Forces the current pointer to Ptr.  Ptr must be in the range buf_begin() <= Ptr <= buf_end()
      void set_buf_ptr(unsigned char const* Ptr)
      { DEBUG_RANGE_CHECK(Ptr, buf_begin(), buf_end()); CurrentPtr = Ptr; }

   private:
      ReadBuffer RB;
      unsigned char const* BeginPtr;
      unsigned char const* CurrentPtr;
      unsigned char const* EndPtr;
};

//
// inlines
//

inline
void OutputBuffer::bump(size_t Count)
{
   CurrentPtr += Count;
}

} // namespace PHeapFileSystem

#endif
