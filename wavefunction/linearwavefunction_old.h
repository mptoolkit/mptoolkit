// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// wavefunction/linearwavefunction_old.h
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
//
// LinearWavefunctionOld - the old format for LinearWavefunction
// This is only to support legacy streaming, to facilitate
// conversion to the new MPWavefunction
// 

#if !defined(MPTOOLKIT_WAVEFUNCTION_LINEARWAVEFUNCTION_OLD_H)
#define MPTOOLKIT_WAVEFUNCTION_LINEARWAVEFUNCTION_OLD_H

#include "state_component.h"
#include "pheap/pvalueptr.h"
#include "pheap/pvalueiterator.h"
#include "interface/attributes.h"

class LinearWavefunctionOld
{
   public:
      typedef StateComponent value_type;
      typedef pvalue_handle<StateComponent> handle_type;
      typedef std::list<handle_type> container_type;
      typedef container_type::iterator base_iterator;
      typedef container_type::const_iterator const_base_iterator;
      typedef pvalue_handle_iterator<base_iterator> iterator;
      typedef const_pvalue_handle_iterator<const_base_iterator> const_iterator;

      LinearWavefunctionOld() {}

      SymmetryList GetSymmetryList() const { return this->begin()->GetSymmetryList(); }

      const_iterator begin() const { return const_iterator(Data.begin()); }
      const_iterator end() const { return const_iterator(Data.end()); }

      base_iterator base_begin() { return Data.begin(); }
      base_iterator base_end() { return Data.end(); }

      const_base_iterator base_begin() const { return Data.begin(); }
      const_base_iterator base_end() const { return Data.end(); }

      AttributeList const& Attributes() const { return Attr; }

   private:
      SymmetryList SList;
      container_type Data;
      AttributeList Attr;

   friend PStream::opstream& operator<<(PStream::opstream& out, LinearWavefunctionOld const& psi);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, LinearWavefunctionOld& psi);
};

inline
PStream::opstream& operator<<(PStream::opstream& out, LinearWavefunctionOld const& psi)
{
   return out << psi.SList << psi.Data << psi.Attr;
}

inline
PStream::ipstream& operator>>(PStream::ipstream& in, LinearWavefunctionOld& psi)
{
   return in >> psi.SList >> psi.Data >> psi.Attr;
}

#endif
