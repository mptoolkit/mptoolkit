// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/pdummystream.h
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

/*
  A 'dummy' output stream class that does nothing; there is a fixed size buffer
  but overflow() simply resets the buf_ptr to the beginning of the buffer.
  This is not useful, except to test in-cache I/O performance without any buffer allocation overheads.

  Created 2004-03-02 Ian McCulloch
*/

#if !defined(PDUMMYSTREAM_H_DSKJHRIUWY365YR87HF4873O)
#define PDUMMYSTREAM_H_DSKJHRIUWY365YR87HF4873O

#include "pstream.h"

namespace PStream
{

class opdummystream : public opstream
{
   public:
      static size_t const DefaultBufferSize = 4096;

      opdummystream(int Format = format::XDR, size_t BufferSize = DefaultBufferSize);

      void set_format(int Format) { opstream::set_format(Format); }

      size_t bytes_written() const;

   protected:
      virtual void overflow();

   private:
      byte_type* Buffer;
      size_t Written;
};

inline
opdummystream::opdummystream(int Format, size_t BufferSize)
  : opstream(Format, NULL, NULL, NULL), Buffer(new byte_type[BufferSize]), Written(0)
{
   this->set_buffer(Buffer, Buffer + BufferSize, Buffer);
}

inline
void opdummystream::overflow()
{
   Written += this->buf_ptr() - Buffer;
   this->set_buf_ptr(Buffer);
}

inline
size_t opdummystream::bytes_written() const
{
   return Written + (this->buf_ptr() - Buffer);
}

} // namespace PStream

#endif
