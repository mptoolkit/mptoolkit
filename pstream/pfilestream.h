// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// pstream/pfilestream.h
//
// Copyright (C) 2002-2016 Ian McCulloch <ian@qusim.net>
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
  Persistent stream class for a POSIX file descriptor

  Created 2002-06-16 Ian McCulloch

  For opfilestream, the destructor flushes the stream but it does not close the
  file descriptor.  The file descriptor must be closed manually.

*/

#if !defined(PFILESTREAM_H_IUY4387RY743Y7G73T78GG)
#define PFILESTREAM_H_IUY4387RY743Y7G73T78GG

#if defined(HAVE_CONFIG_H)
#include "config.h" // for LARGEFILE configuration
#endif
#include "pstream.h"
#include "common/trace.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <exception>
#include <errno.h>

namespace PStream
{

#if !defined(HAVE_STRERROR_R)
inline
int strerror_r(int errnum, char *buf, size_t buflen)
{
   strncpy(buf, strerror(errnum), buflen);
   return 0;
}
#endif

class runtime_errno : public std::exception
{
   public:
      runtime_errno() : m_Errno(errno) { strerror_r(m_Errno, StrBuffer, StrBufferSize); }
      runtime_errno(int Errno_) : m_Errno(Errno_) { strerror_r(m_Errno, StrBuffer, StrBufferSize); }

      virtual char const* what() const throw() { return StrBuffer; }

      int what_errno() const { return m_Errno; }

   private:
      static size_t const StrBufferSize = 240;
      int m_Errno;
      char StrBuffer[StrBufferSize];
};

inline
void throw_runtime_errno()
{
#if !defined(ERRNO_ERRORS_NON_FATAL)
   int Errno = errno;
   PANIC(Errno) << strerror(Errno);
#else
   throw runtime_errno(errno);
#endif
}

inline
void check_errno(int ReturnCode)
{
   if (ReturnCode == -1) throw_runtime_errno();
}

class opfilestream : public opstream
{
   public:
      static size_t const DefaultBufferSize = 4096;

      opfilestream(int Format = format::XDR, size_t BufferSize = DefaultBufferSize);

      ~opfilestream();

      void open(std::string const& pathname, int flags, mode_t mode = 0666);

      void close();

      off_t lseek(off_t offset, int whence);
      off_t tell() const;

      // truncates the file at the current location
      void truncate();

      // truncates the file to the given size
      void truncate(off_t Loc);

      void put_format(int NewFormat) { opstream::put_format(NewFormat); }
      void put_format() { opstream::put_format(); }

      void set_fd(int FD_) { FD = FD_; }
      int get_fd() { return FD; }

      void set_format(int Format) { opstream::set_format(Format); }

   protected:
      virtual void overflow();

      int FD;
};

class ipfilestream : public ipstream
{
   public:
      static size_t const DefaultBufferSize = 4096;

      ipfilestream(int Format = format::XDR, size_t BufferSize = DefaultBufferSize);

      ~ipfilestream();

      void open(std::string const& pathname, int flags, mode_t mode = 0666);

      void close();

      off_t lseek(off_t offset, int whence);
      off_t tell() const;

      int get_format() { return ipstream::get_format(); }

      void set_fd(int FD_) { FD = FD_; }
      int get_fd() { return FD; }

      void set_format(int Format) { ipstream::set_format(Format); }

      bool eof() { if (!Eof && (this->buf_ptr() == this->buf_end())) this->underflow(); return Eof; }

   protected:
      virtual void underflow();

      int FD;
      size_t BufferSize;
      bool Eof;
};

//
// inlines
//

// opfstream

inline
void opfilestream::open(std::string const& pathname, int flags, mode_t mode)
{
   if (FD != -1) this->flush();
   check_errno((FD = ::open(pathname.c_str(), flags, mode)));
}

inline
void opfilestream::close()
{
   this->flush();
   check_errno(::close(FD));
   FD = -1;
}

inline
off_t opfilestream::lseek(off_t offset, int whence)
{
   this->flush();
   off_t Result = ::lseek(FD, offset, whence);
   if (Result == off_t(-1)) throw_runtime_errno();
   return Result;
}

inline
off_t opfilestream::tell() const
{
   off_t Result = ::lseek(FD, 0, SEEK_CUR);
   if (Result == off_t(-1)) throw_runtime_errno();
   return Result + (this->buf_ptr() - this->buf_begin());
}

inline
void opfilestream::truncate()
{
   this->truncate(this->tell());
}

inline
void opfilestream::truncate(off_t loc)
{
   check_errno(::ftruncate(FD, loc));
}

// ipfstream

inline
void ipfilestream::open(std::string const& pathname, int flags, mode_t mode)
{
   check_errno((FD = ::open(pathname.c_str(), flags, mode)));
   this->set_buf_ptr(this->buf_end());              // force underflow on next read
   Eof = false;
}

inline
void ipfilestream::close()
{
   check_errno(::close(FD));
   Eof = false;
   FD = -1;
}

inline
off_t ipfilestream::tell() const
{
   off_t Where = ::lseek(FD, 0, SEEK_CUR);
   if (Where == off_t(-1)) throw_runtime_errno();
   return Where + (this->buf_ptr() - this->buf_end());
}

inline
off_t ipfilestream::lseek(off_t offset, int whence)
{
   // if we are seeking relative to the current position, adjust the offset for the
   // current position in the buffer
   if (whence == SEEK_CUR)
      offset += this->buf_ptr() - this->buf_end();

   off_t Result = ::lseek(FD, offset, whence);
   if (Result == off_t(-1)) throw_runtime_errno();
   this->set_buf_ptr(this->buf_end());              // force underflow on next read
   Eof = false;                                     // clear EOF flag
   return Result;
}

} // namespace PStream

#endif
