// -*- C++ -*- $Id$

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
