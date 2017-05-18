// -*- C++ -*-
//
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
