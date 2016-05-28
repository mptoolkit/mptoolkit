// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pheapstream.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "common/trace.h"
#include "pheap.h"
#include <string.h>

namespace PHeapFileSystem
{

// opheapstream

opheapstream::opheapstream(BlockFileSystem* FS_, Descriptor* Descriptor_) 
  : opstream(0, NULL, NULL, NULL), MyDescriptor(Descriptor_), MyFS(FS_)
{
   CHECK(MyFS != NULL);
   // Allocate a filesystem output buffer
   OutBuffer = MyFS->allocate_output_buffer();
   // synchronize our stream output buffer
   this->set_buffer(OutBuffer->buf_begin(), OutBuffer->buf_end(), OutBuffer->buf_ptr());
}

opheapstream::~opheapstream()
{
   // flush the buffer
   if (this->buf_ptr() != this->buf_begin())
   {
      OutBuffer->set_buf_ptr(this->buf_ptr());
      MyDescriptor->append_block(OutBuffer->commit());
   }
   MyFS->return_output_buffer(OutBuffer);
}

void opheapstream::overflow()
{
   // synchronize the OutBuffer pointer to the stream pointer
   OutBuffer->set_buf_ptr(this->buf_ptr());
   MyDescriptor->append_block(OutBuffer->commit());
   OutBuffer = MyFS->cycle_output_buffer(OutBuffer);
   // synchronize the stream buffer to the new OutBuffer
   this->set_buffer(OutBuffer->buf_begin(), OutBuffer->buf_end(), OutBuffer->buf_ptr());
}

void opheapstream::put_id(id_type ID)
{
   MyDescriptor->append_id(ID);
}

// ipheapstream

ipheapstream::ipheapstream(BlockFileSystem* FS_, Descriptor const* Descriptor_)
  : ipstream(0, NULL, NULL, NULL), MyDescriptor(Descriptor_), MyFS(FS_)
{
   CurrentBlock = MyDescriptor->block_begin();
   CurrentId = MyDescriptor->id_begin();
   InBuffer = InputBuffer(*CurrentBlock);
   // synchronize the stream buffer to the filesystem input buffer
   this->set_buffer(InBuffer.buf_begin(), InBuffer.buf_end(), InBuffer.buf_ptr());
   this->get_format();
}

ipheapstream::~ipheapstream()
{
}

void ipheapstream::underflow()
{
   InBuffer.set_buf_ptr(this->buf_ptr());   // synchronize the buffer pointer
   ++CurrentBlock;
   if (CurrentBlock == MyDescriptor->block_end()) return;  // end-of-buffer (EOF)

   InBuffer = InputBuffer(*CurrentBlock);
   this->set_buffer(InBuffer.buf_begin(), InBuffer.buf_end(), InBuffer.buf_ptr());
}

id_type ipheapstream::get_id()
{
   if (CurrentId == MyDescriptor->id_end()) return 0;
   return *CurrentId++;
}

} // namespace PHeapFileSystem
