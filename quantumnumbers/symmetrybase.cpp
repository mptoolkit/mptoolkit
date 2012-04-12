// -*- C++ -*- $Id$

#include "common/trace.h"
#include "symmetrybase.h"

namespace QuantumNumbers
{

namespace
{

// this is a local static singleton, to avoid global
// static construction ordering problems.
// This data is not freed before program exit (deliberate & harmless leak).
std::vector<SymmetryFactory*>& 
RegisteredObjects()
{
   static std::vector<SymmetryFactory*> Instance;
   return Instance;
}

} // namespace

//
// SymmetryBase
//

std::map<std::string, SymmetryBase*> SymmetryBase::CreatedInstances;

SymmetryBase* SymmetryBase::Create(std::string const& Type)
{
   if (CreatedInstances.count(Type)) return CreatedInstances[Type];

   SymmetryBase* S = NULL;
   std::size_t i = 0;
   while (S == NULL && i != RegisteredObjects().size())
   {
      S = RegisteredObjects()[i]->AttemptCreate(Type);
      ++i;
   }
   if (S == NULL) { PANIC(Type)("Construction of quantum number failed, unrecognized symmetry type."); }

   CreatedInstances[Type] = S;
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
   RegisteredObjects().push_back(F);
}

} // namespace QuantumNumbers
