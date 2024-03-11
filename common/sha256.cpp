// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// common/sha256.cpp
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

#include "sha256.h"
#include <cstdint>

namespace SHA256
{

#define DBL_INT_ADD(a,b,c) if (a > 0xffffffff - (c)) ++b; a += c;
#define ROTLEFT(a,b) (((a) << (b)) | ((a) >> (32-(b))))
#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ ((x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ ((x) >> 10))

uint32_t k[64] = {
	0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
	0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
	0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
	0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
	0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
	0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
	0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
	0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

std::array<uint32_t, 8>
hash(std::vector<uint32_t> const& Data)
{
   return hash_block(pad_array(Data));
}

std::array<uint32_t, 16> pad_array(std::vector<uint32_t> const& Data)
{
   // Equation for padding is Len + 1 + k = 448 mod 512
   // Simplified to: Len + 1 + k = 448
   //		  448 - 1 - Len = k
   //		  447 - Len = k
   // Len = length of message in bits
   // k = how much zero's to add so new message will be a multiple of 512.
   std::array<uint32_t, 16> Block;
   std::copy(Data.begin(), Data.end(), Block.begin());
   int n = Data.size();
   Block[n++] = 0x80000000;
   while (n < 14)
   {
      Block[n++] = 0;
   }
   // length is 64-bits big-endian
   Block[14] = 0;
   Block[15] = Data.size() * 32;
   return Block;
}

std::array<uint32_t, 8>
hash_block(std::array<uint32_t, 16> const& block, std::array<uint32_t, 8> const& State)
{
   std::array<uint32_t, 64> m;

   for (int i = 0; i < 16; ++i)
   {
      m[i] = block[i];
   }
   for (int i = 16; i < 64; ++i)
   {
      m[i] = SIG1(m[i - 2]) + m[i - 7] + SIG0(m[i - 15]) + m[i - 16];
   }

   uint32_t a = State[0];
   uint32_t b = State[1];
   uint32_t c = State[2];
   uint32_t d = State[3];
   uint32_t e = State[4];
   uint32_t f = State[5];
   uint32_t g = State[6];
   uint32_t h = State[7];

   for (int i = 0; i < 64; ++i)
   {
      uint32_t t1 = h + EP1(e) + CH(e, f, g) + k[i] + m[i];
      uint32_t t2 = EP0(a) + MAJ(a, b, c);
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
   }

   return std::array<uint32_t, 8>{a,b,c,d,e,f,g,h};
}

std::array<uint32_t, 8>
hash_block(std::array<uint32_t, 16> const& block)
{
   return hash_block(block, {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                           0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19});
}

} // namespace SHA256
