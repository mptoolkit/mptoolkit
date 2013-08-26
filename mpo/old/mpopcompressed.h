// -*- C++ -*- $Id$
//
// Run-length compressed matrix product operator,
// and associated low-level functions
//

#if !defined(MPOPCOMPRESSED_H_JDCHJKEHY589758YUER89H489)
#define MPOPCOMPRESSED_H_JDCHJKEHY589758YUER89H489

#include "common/runlengthcompressed.h"
#include "pstream/pstream.h"
#include "quantumnumbers/quantumnumber.h"
#include "mpopcomponent.h"

typedef run_length_compressed<MPOpComponent> MPOpCompressed;
typedef run_length_repeat<MPOpComponent> MPOpRepeat;
typedef run_length_array<MPOpComponent> MPOpArray;

// Splits a compressed operator into two pieces, the first piece has size
// Loc, the second piece has size x.size()-Loc.
inline
std::pair<MPOpCompressed, MPOpCompressed>
split_operator(MPOpCompressed const& x, int Loc)
{
   return split(x, Loc);
}

#endif
