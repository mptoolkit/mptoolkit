// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/sha256.h
//
// Copyright (C) 2017 Ian McCulloch <ian@qusim.net>
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

#if !defined(MPTOOLKIT_COMMON_SHA256_H)
#define MPTOOLKIT_COMMON_SHA256_H

#include <cstdint>
#include <array>
#include <vector>

namespace SHA256
{

// Returns the SHA256 hash of the given vector
std::array<uint32_t, 8>
hash(std::vector<uint32_t> const& Data);

// Constructs a properly added 512-bit block from
// an array (which must be size less than 13)
std::array<uint32_t, 16>
pad_array(std::vector<uint32_t> const& Data);

// returns the SHA256 of the supplied block
std::array<uint32_t, 8>
hash_block(std::array<uint32_t, 16> const& block);


} // namespace

#endif
