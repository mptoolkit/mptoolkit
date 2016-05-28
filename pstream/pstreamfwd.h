// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// pstream/pstreamfwd.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if !defined(PSTREAMFWD_H_HDJYGFYUGT743675887)
#define PSTREAMFWD_H_HDJYGFYUGT743675887

#include "config.h"
#include "common/inttype.h"

namespace PStream
{

typedef unsigned char     byte_type;
typedef inttype::uint64_t id_type;

class opstream;
class ipstream;

template <int Which>
class opstreambuf;

template <int Which>
class ipstreambuf;

}

#endif
