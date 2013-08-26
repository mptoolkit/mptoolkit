// -*- C++ -*- $Id$

// A UnitCell is an array of LatticeSite's.  This uses run-length compression to handle
// large arrays.

// This was formerly known as a Lattice, and used to contain a mechanism to set
// the coordinates of each site.  This is removed now and the coordinates of
// each site are fixed to be 1,2,...

#if !defined(UNITCELL_H_GT3QY87GYO87Y437YWLO87YLO0)
#define UNITCELL_H_GT3QY87GYO87Y437YWLO87YLO0

#include "siteoperator/latticesite.h"
#include "common/runlengthcompressed.h"

class UnitCell
{
   public:
      typedef LatticeSite                        value_type;
      typedef run_length_compressed<LatticeSite> data_type;
      typedef data_type::const_iterator          const_iterator;

      UnitCell();
      UnitCell(LatticeSite const& s);
      UnitCell(LatticeSite const& s, LatticeSite const& t);
      UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u);
      UnitCell(LatticeSite const& s, LatticeSite const& t, LatticeSite const& u, LatticeSite const& v);

      UnitCell(SymmetryList const& sl, LatticeSite const& s);
      UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t);
      UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, 
              LatticeSite const& u);
      UnitCell(SymmetryList const& sl, LatticeSite const& s, LatticeSite const& t, 
              LatticeSite const& u, LatticeSite const& v);

      UnitCell(int RepeatCount, UnitCell const& l);
      UnitCell(UnitCell const& x1, UnitCell const& x2);
      UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3);
      UnitCell(UnitCell const& x1, UnitCell const& x2, UnitCell const& x3, UnitCell const& x4);

      UnitCell(int Size, LatticeSite const& s);

      // Constructs a singleton lattice, setting the coordinate at the same time
      UnitCell(LatticeSite const& s, std::string const& Coord);

      template <typename T>
      UnitCell(LatticeSite const& s, T const& Coord);

      run_length_compressed<LatticeSite> const& data() const { return Data_; }

      SymmetryList GetSymmetryList() const { return Data_.front().GetSymmetryList(); }
   
      // fowards to run_length_compressed
      bool empty() const { return Data_.empty(); }
      int size() const { return Data_.size(); }
      int leaf_count() const { return Data_.leaf_count(); }
      int node_count() const { return Data_.node_count(); }
      value_type const& front() const { return Data_.front(); }
      value_type const& back() const { return Data_.back(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }

      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v) const;

      template <typename Visitor>
      typename Visitor::result_type
      apply_visitor(Visitor const& v);

      LatticeSite const& operator[](int i) const;

   private:
      run_length_compressed<LatticeSite> Data_;

   friend PStream::opstream& operator<<(PStream::opstream& out, UnitCell const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, UnitCell& L);
};

// Make a new lattice out of RepeatCount copies of x
UnitCell repeat(UnitCell const& x, int RepeatCount);

// Make a new lattice out of the join of x and y
UnitCell join(UnitCell const& x, UnitCell const& y);
UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z);
UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w);
UnitCell join(UnitCell const& x, UnitCell const& y, UnitCell const& z, UnitCell const& w,
	     UnitCell const& v);

#include "unitcell.cc"

#endif
