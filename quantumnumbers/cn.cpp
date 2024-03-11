// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// quantumnumbers/cn.cpp
//
// Copyright (C) 2003-2016 Ian McCulloch <ian@qusim.net>
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
  cn.cpp

  Implementation of cyclic group symmetry

  Created 2003-10-20 Ian McCulloch

*/

#include "cn.h"
#include "symmetrybase.h"

namespace QuantumNumbers
{

class CnSymmetryFactory : public SymmetryFactory
{
   public:
      virtual SymmetryBase* AttemptCreate(std::string const& Type);
};

SymmetryBase* CnSymmetryFactory::AttemptCreate(std::string const& Type)
{
   // see if Type has the form "C_N" for N an integer
   if (Type.size() < 3) return NULL;
   if (Type[0] != 'C' || Type[1] != '_') return NULL;

   int N;
   try
   {
      N = convert_string<int>(Type.begin()+2, Type.end());
   }
   catch (invalid_string_conversion const&)
   {
      return NULL;
   }

   return new BasicSymmetry<Cn>(CnFactory(N));
}

Cn::Cn(int N_, std::string const& s)
  : N(N_), x((boost::lexical_cast<int>(s) + 10*N_) % N_)
{
}

std::string
Cn::ToString() const
{
   return ConvertToString(x);
}

void Cn::Register()
{
   RegisterSymmetry(new CnSymmetryFactory());
}

std::string
CnProjection::ToString() const
{
   return std::string();
}

} // namespace QuantumNumbers
