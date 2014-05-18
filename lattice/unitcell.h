// -*- C++ -*- $Id$

// A UnitCell is an array of LatticeSite's, that represents a primitive unit for a lattice.
// The sites in the lattice are denoted [0], [1], .... using square brackets, and zero-based counting.
//
// For translationally invariant systems, this is the basic unit that defines the repeated unit.

// the unit cell defines operators with finite support.  Eg, the UnitCell might contain N LatticeSite's,
// each of which has a spin operator S.  We refer to the n'th LatticeSite operator as S[n].  We might also
// define another operator on the UnitCell as S = S[0]+S[1]+...+S[N].
// Once we get to a Lattice, as a collection of UnitCell's, we would refer to the total spin operator at
// the i'th UnitCell as S(i).  An individual spin at the n'th site of the i'th UnitCell would be referred to
// as S(i)[n].  The Lattice itself will probably define its own operator S as the sum of S(i) over all unit cells.
//
// The unit cell of the operators is allowed to be a multiple of the lattice UnitCell. In this case,
// we would refer to a specific operator on the lattice as the left-most unit cell where it has non-trivial support.
//
// The UnitCell always contains the global identity operator "I".  This is added by the UnitCell constructor.
//
// If the UnitCell is a single site, we allow local operators to act as UnitCell operators too, and
// automatically generate them if necessary.
//
// How to handle the J-W string needs more thought.
//

#if !defined(MPTOOLKIT_LATTICE_UNITCELL_H)
#define MPTOOLKIT_LATTICE_UNITCELL_H

#include "latticesite.h"
#include "mpo/finite_mpo.h"
#include <vector>
#include <map>
#include <string>

class UnitCell
{
   public:
      typedef LatticeSite               value_type;
      typedef std::vector<LatticeSite>  data_type;
      typedef data_type::const_iterator const_iterator;

      typedef std::map<std::string, FiniteMPO>  operator_map_type;
      typedef operator_map_type::const_iterator const_operator_iterator;

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

      SymmetryList GetSymmetryList() const { return Data_.front().GetSymmetryList(); }

      // functions acting on the LatticeSite
   
      bool empty() const { return Data_.empty(); }
      int size() const { return Data_.size(); }
      value_type const& front() const { return Data_.front(); }
      value_type const& back() const { return Data_.back(); }

      const_iterator begin() const { return Data_.begin(); }
      const_iterator end() const { return Data_.end(); }

      // returns the i'th LatticeSite
      LatticeSite const& operator[](int i) const;

#if 0
      // Visitor pattern tools for iterating over the LatticeSite's
      template <typename Visitor>
      typename Visitor::result_type
      apply_for_each_site(Visitor const& v) const;

      template <typename Visitor>
      typename Visitor::result_type
      apply_for_each_site(Visitor const& v);
#endif

      // returns the local basis at site i of the lattice
      SiteBasis LocalBasis(int i) const
      { return this->operator[](i).Basis2(); }

      // Functions for the associated operators

      // returns true if the named operator exists on this unit cell
      bool operator_exists(std::string const& s) const;

      // lookup a unit cell operator, which can be a local operator
      // if the unit cell size is one site.
      FiniteMPO Operator(std::string const& Op) const;

      // lookup a unit cell operator (not a local operator), or
      // adds it if it doesn't already exist.
      FiniteMPO& Operator(std::string const& Op);

      // lookup the Jordan Wigner string of the specified operator
      std::string JWString(std::string const& Op) const;

      // returns true if the specified local operator Op[n] exists at site n of the UnitCell.
      bool operator_exists(std::string const& Op, int n) const;

      // lookup a local operator
      FiniteMPO Operator(std::string const& Op, int n) const;

      // Parse an operator of the form O or O[n]
   //      FiniteMPO Parse(std::string const& s);

   //std::string JWString(std::string const& Op, int n) const;

      // returns a begin() iterator into the unit cell operators (not local operators!)
      const_operator_iterator operator_begin() const;
      const_operator_iterator operator_end() const;

   private:
      data_type Data_;
      operator_map_type OperatorMap_;

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

// Utility function to make the identity MPO
FiniteMPO
identity_mpo(UnitCell const& c);

FiniteMPO
identity_mpo(UnitCell const& c, QuantumNumbers::QuantumNumber const& q);

#if 0 // These operators probably don't make much sense
bool
operator==(UnitCell const& u1, UnitCell const& u2);

bool
operator!=(UnitCell const& u1, UnitCell const& u2);
#endif

std::ostream&
operator<<(std::ostream& out, UnitCell const& u);

#include "unitcell.cc"

#endif
