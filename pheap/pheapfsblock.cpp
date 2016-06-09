// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pheapfsblock.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

#include "pheapfsblock.h"

#include <iostream>
#include <iomanip>

namespace PHeapFileSystem
{

// MinimumUsefulBufferSize is a cutoff for deciding whether to re-use a buffer.
// Buffers smaller than this value are not re-used as it is assumed that 
// the overhead per buffer is larger than the useful storage space.
size_t const MinimumUsefulBufferSize = 64;

//
// BlockRecord
//

BlockRecord::BlockRecord()
{
}

BlockRecord::BlockRecord(PageId const& Page_, size_t Offset_, size_t Length_)
   : Page(Page_), Offset(Offset_), Length(Length_)
{
}

//
// OutputBuffer
//

OutputBuffer::~OutputBuffer()
{
}

OutputBuffer::OutputBuffer(WriteBuffer const& Buf)
   : WB(Buf), BeginPtr(Buf.buffer()),
     CurrentPtr(Buf.buffer()), EndPtr(Buf.buffer() + Buf.get_fs()->get_page_size())
{
} 

BlockRecord OutputBuffer::commit()
{
   BlockRecord Temp(WB.get_page(), BeginPtr - WB.buffer(), CurrentPtr - BeginPtr);
   BeginPtr = CurrentPtr;
   return Temp;
}

//
// InputBuffer
//

InputBuffer::InputBuffer() 
   : BeginPtr(NULL), CurrentPtr(NULL), EndPtr(NULL)
{
}

InputBuffer::InputBuffer(BlockRecord const& R)
   : RB(R.Page)
{
   CurrentPtr = BeginPtr = RB.buffer() + R.Offset;
   EndPtr = BeginPtr + R.Length;
}   

//
// BlockFileSystem
//

BlockFileSystem::BlockFileSystem()
{
}

BlockFileSystem::~BlockFileSystem()
{
}

OutputBuffer* BlockFileSystem::allocate_output_buffer()
{
   pthread::mutex::sentry Lock(BufferListMutex);
   if (Buffers.empty())
   {
      //     std::cout << "Allocating buffer: new buffer." << std::endl;
      return new OutputBuffer(allocate());
   }
   else
   {
      OutputBuffer* Ret = Buffers.back();
      //      std::cout << "Allocating buffer: recycled buffer has remain=" << Ret->remain() << std::endl;
      Buffers.pop_back();
      return Ret;
   }
}

void BlockFileSystem::return_output_buffer(OutputBuffer* Buf)
{
   DEBUG_PRECONDITION(Buf != NULL);
   DEBUG_PRECONDITION(Buf->buf_ptr() == Buf->buf_begin());

   //   std::cout << "returning buffer to system, used = " << (Buf->buf_ptr() - Buf->buf_begin())
      //	     << " remain = " << Buf->remain() << std::endl;

   size_t Remain = Buf->buf_end() - Buf->buf_begin();
   if (Remain < MinimumUsefulBufferSize)
   {
      delete Buf;
   }
   else
   {
      pthread::mutex::sentry Lock(BufferListMutex);
      Buffers.push_back(Buf);
   }
}

OutputBuffer* BlockFileSystem::cycle_output_buffer(OutputBuffer* Buf)
{
   DEBUG_PRECONDITION(Buf != NULL);
   DEBUG_PRECONDITION(Buf->buf_ptr() == Buf->buf_begin());

   //   std::cout << "cycling buffer, used = " << (Buf->buf_ptr() - Buf->buf_begin())
   //	     << " remain = " << Buf->remain() << std::endl;

   size_t Remain = Buf->buf_end() - Buf->buf_begin();
   if (Remain < MinimumUsefulBufferSize)
   {
      delete Buf;
      return allocate_output_buffer();
   }
   return Buf;
}

void BlockFileSystem::WriteBlockRecord(PStream::opstream& out, BlockRecord const& Rec)
{
   this->write_page_id(out, Rec.Page);
   out.write<uint32>(Rec.Offset);
   out.write<uint32>(Rec.Length);
}

BlockRecord BlockFileSystem::ReadBlockRecord(PStream::ipstream& in)
{
   BlockRecord Rec;
   Rec.Page = this->read_page_id(in);
   Rec.Offset = in.read<uint32>();
   Rec.Length = in.read<uint32>();
   return Rec;
}

} // namespace PHeapFileSystem
