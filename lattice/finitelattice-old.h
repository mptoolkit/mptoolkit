// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// lattice/finitelattice-old.h
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

// A FiniteLattice is an array of Unit Cells.  This uses run-length compression to handle
// large arrays.

#if !defined(MPTOOLKIT_LATTICE_FINITELATTICE_H)
#define MPTOOLKIT_LATTICE_FINITELATTICE_H

#include "lattice/latticesite.h"
#include "common/runlengthcompressed.h"
#include "pheap/pvalueptr.h"

typedef run_length_compressed<UnitCell> UnitCellListType;
typedef pvalue_ptr<UnitCellListType> UnitCellPtrType;

class FiniteLattice
{
   public:
      typedef UnitCell                        value_type;
      typedef run_length_compressed<UnitCell> data_type;
      typedef data_type::const_iterator       const_iterator;

      FiniteLattice();
      FiniteLattice(UnitCell const& s);
      FiniteLattice(UnitCell const& s, UnitCell const& t);
      FiniteLattice(UnitCell const& s, UnitCell const& t, UnitCell const& u);
      FiniteLattice(UnitCell const& s, UnitCell const& t, UnitCell const& u, UnitCell const& v);

      FiniteLattice(SymmetryList const& sl, UnitCell const& s);
      FiniteLattice(SymmetryList const& sl, UnitCell const& s, UnitCell const& t);
      FiniteLattice(SymmetryList const& sl, UnitCell const& s, UnitCell const& t,
              UnitCell const& u);
      FiniteLattice(SymmetryList const& sl, UnitCell const& s, UnitCell const& t,
              UnitCell const& u, UnitCell const& v);

      FiniteLattice(int RepeatCount, FiniteLattice const& l);
      FiniteLattice(FiniteLattice const& x1, FiniteLattice const& x2);
      FiniteLattice(FiniteLattice const& x1, FiniteLattice const& x2, FiniteLattice const& x3);
      FiniteLattice(FiniteLattice const& x1, FiniteLattice const& x2, FiniteLattice const& x3, FiniteLattice const& x4);

      FiniteLattice(int Size, UnitCell const& s);

      // Constructs a singleton lattice, setting the coordinate at the same time
      FiniteLattice(LatticeSite const& s, std::string const& Coord);

      template <typename T>
      FiniteLattice(LatticeSite const& s, T const& Coord);

      run_length_compressed<LatticeSite> const& data() const { return *Data_; }

      SymmetryList GetSymmetryList() const { return Data_->front().GetSymmetryList(); }

      // fowards to run_length_compressed
      bool empty() const { return Data_->empty(); }
      int size() const { return Data_->size(); }
      int leaf_count() const { return Data_->leaf_count(); }
      int node_count() const { return Data_->node_count(); }
      value_type const& front() const { return Data_->front(); }
      value_type const& back() const { return Data_->back(); }

      const_iterator begin() const { return Data_->begin(); }
      const_iterator end() const { return Data_->end(); }

      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v) const;

      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v);

      LatticeSite const& operator[](int i) const;

      // returns the local basis at site i of the lattice
      SiteBasis LocalBasis(int i) const
      { return this->operator[](i)["I"].Basis(); }

   private:
      typedef run_length_compressed<LatticeSite> FiniteLatticeType;
      pvalue_ptr<FiniteLatticeType> Data_;

   friend PStream::opstream& operator<<(PStream::opstream& out, FiniteLattice const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, FiniteLattice& L);

   friend bool operator==(FiniteLattice const& u1, FiniteLattice const& u2);
   friend bool operator!=(FiniteLattice const& u1, FiniteLattice const& u2);
};

// Make a new lattice out of RepeatCount copies of x
FiniteLattice repeat(FiniteLattice const& x, int RepeatCount);

// Make a new lattice out of the join of x and y
FiniteLattice join(FiniteLattice const& x, FiniteLattice const& y);
FiniteLattice join(FiniteLattice const& x, FiniteLattice const& y, FiniteLattice const& z);
FiniteLattice join(FiniteLattice const& x, FiniteLattice const& y, FiniteLattice const& z, FiniteLattice const& w);
FiniteLattice join(FiniteLattice const& x, FiniteLattice const& y, FiniteLattice const& z, FiniteLattice const& w,
             FiniteLattice const& v);

bool
operator==(FiniteLattice const& u1, FiniteLattice const& u2);

bool
operator!=(FiniteLattice const& u1, FiniteLattice const& u2);

std::ostream&
operator<<(std::ostream& out, FiniteLattice const& u);

#include "unitcell.cc"

#endif
