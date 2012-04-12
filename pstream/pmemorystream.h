/* -*- C++ -*- $Id$

  A simple vector-based memory stream class

  Created 2004-03-02 Ian McCulloch
*/

#include "pstream.h"

namespace PStream
{

class opmemorystream : public opstream
{
   public:
      static size_t const DefaultBufferSize = 4096;

      opmemorystream(int Format = format::XDR, size_t BufferSize = DefaultBufferSize);

      ~opmemorystream();
      
      void set_format(int Format) { opstream::set_format(Format); }

   protected:
      virtual void overflow();

   private:
      std::vector<unsigned char> Buffer;
};

opmemorystream::~opmemorystream()
{
}

opmemorystream::opmemorystream(int Format, size_t BufferSize)
  : opstream(Format, NULL, NULL, NULL), Buffer(BufferSize)
{
   this->set_buffer(&Buffer[0], &Buffer[0] + BufferSize, &Buffer[0]);
}

void opmemorystream::overflow()
{
   int UsedSize = this->buf_ptr() - this->buf_begin();

   size_t NewBufSize = DefaultBufferSize + UsedSize * 2;
   Buffer.resize(NewBufSize);

   this->set_buffer(&Buffer[0], &Buffer[0] + NewBufSize, &Buffer[0] + UsedSize);
}

};
