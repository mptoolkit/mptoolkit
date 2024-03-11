// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/pfilestream.cpp
//
// Copyright (C) 2004-2016 Ian McCulloch <ian@qusim.net>
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

#include "pfilestream.h"

namespace PStream
{

//
// opfstream
//

opfilestream::~opfilestream()
{
   if (FD != -1)
   {
      this->flush();
   }
   delete[] this->get_buffer()->buf_begin();
}

  /*
opfilestream::opfilestream()
  : opstream(format::Current, NULL, NULL, NULL), FD(-1)
{
   byte_type* Buf = new byte_type[DefaultBufferSize];
   this->set_buffer(Buf, Buf+DefaultBufferSize, Buf);
}
  */

opfilestream::opfilestream(int Format, size_t BufferSize)
  : opstream(Format, NULL, NULL, NULL), FD(-1)
{
   byte_type* Buf = new byte_type[BufferSize];
   this->set_buffer(Buf, Buf+BufferSize, Buf);
}

void opfilestream::overflow()
{
   int Size = this->buf_ptr() - this->buf_begin();
   if (Size == 0) return;

   ssize_t Written = ::write(FD, this->buf_begin(), Size);
   if (Written != Size)
   {
      throw_runtime_errno();
   }
   this->set_buf_ptr(this->buf_begin());
}

//
// ipfilestream
//

ipfilestream::~ipfilestream()
{
   delete[] this->get_buffer()->buf_begin();
}

/*
ipfilestream::ipfilestream()
  : ipstream(format::Current, NULL, NULL, NULL), FD(-1), BufferSize(DefaultBufferSize)
{
   byte_type* Buf = new byte_type[BufferSize];
   this->set_buffer(Buf, Buf+BufferSize, Buf+BufferSize);
}
*/

ipfilestream::ipfilestream(int Format, size_t BufferSize_)
   : ipstream(Format, NULL, NULL, NULL), FD(-1), BufferSize(BufferSize_), Eof(false)
{
   byte_type* Buf = new byte_type[BufferSize];
   this->set_buffer(Buf, Buf+BufferSize, Buf+BufferSize);
}

void ipfilestream::underflow()
{
   CHECK(this->buf_ptr() == this->buf_end());
   CHECK(!Eof);

   ssize_t BytesRead = ::read(FD, const_cast<byte_type*>(this->buf_begin()), BufferSize);
   if (BytesRead == -1) throw_runtime_errno();
   if (BytesRead == 0) Eof = true;

   this->set_buf_end(this->buf_begin() + BytesRead);
   this->set_buf_ptr(this->buf_begin());
}

} // namespace PStream
