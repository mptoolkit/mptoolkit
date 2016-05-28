// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// quantumnumbers/symmetrybase.cpp
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

#include "common/trace.h"
#include "symmetrybase.h"

namespace QuantumNumbers
{

//
// SymmetryBase
//

// global data, initialized by nifty counter
std::map<std::string, SymmetryBase*>* SymmetryBase::CreatedInstances = NULL;
std::vector<SymmetryFactory*>* RegisteredObjects = NULL;

void
SymmetryBase::InitializeInstances()
{
   SymmetryBase::CreatedInstances = new std::map<std::string, SymmetryBase*>();
   RegisteredObjects = new std::vector<SymmetryFactory*>();
}

SymmetryBase* SymmetryBase::Create(std::string const& Type)
{
   if (CreatedInstances->count(Type)) return (*CreatedInstances)[Type];

   SymmetryBase* S = NULL;
   std::size_t i = 0;
   while (S == NULL && i != RegisteredObjects->size())
   {
      S = (*RegisteredObjects)[i]->AttemptCreate(Type);
      ++i;
   }
   if (S == NULL) { PANIC(Type)("Construction of quantum number failed, unrecognized symmetry type."); }

   (*CreatedInstances)[Type] = S;
   return S;
}

SymmetryBase::~SymmetryBase()
{
}

SymmetryFactory::~SymmetryFactory()
{
}

void SymmetryFactory::Register(SymmetryFactory* F)
{
   RegisteredObjects->push_back(F);
}

} // namespace QuantumNumbers
