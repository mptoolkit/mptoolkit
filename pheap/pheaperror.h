// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pheap/pheaperror.h
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

// Error classes for pheap

#if !defined(MPTOOLKIT_PHEAP_PHEAPERROR_H)
#define MPTOOLKIT_PHEAP_PHEAPERROR_H

#include <exception>
#include <string>

namespace pheap
{

class PHeapError : public std::runtime_error
{
   public:
      PHeapError(std::string const& s) : std::runtime_error(s) {}
};

class PHeapVersionMismatch : public PHeapError
{
   public:
      PHeapVersionMismatch(std::string const& s) : PHeapError(s) {}
};

class PHeapFileError : public PHeapError
{
   public:
      PHeapFileError(std::string Filename_, std::string const& Why_, std::string const& s) 
	 : PHeapError(s), Filename(Filename_), Why(Why_) {}

      PHeapFileError(std::string Filename_, std::string const& Why_) 
	 : PHeapError("Error occurred processing file " + Filename_ + " : " + Why_), 
	   Filename(Filename_), Why(Why_) {}

      std::string Filename;
      std::string Why;
};


class PHeapCannotCreateFile : public PHeapFileError
{
   public:
      PHeapCannotCreateFile(std::string const& Filename_, std::string const& Why_)
	 : PHeapFileError(Filename_, Why_, "Cannot create file " + Filename_ + "\nReason: " + Why_) {}
};

class PHeapCannotOpenFile : public PHeapFileError
{
   public:
      PHeapCannotOpenFile(std::string const& Filename_, std::string const& Why_)
	 : PHeapFileError(Filename_, Why_, "Cannot open file " + Filename_ + "\nReason: " + Why_) {}
};

} // namespace PHeap

#endif
