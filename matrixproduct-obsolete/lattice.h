// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/lattice.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#if !defined(LATTICE_H_GT3QY87GYO87Y437YWLO87YLO0)
#define LATTICE_H_GT3QY87GYO87Y437YWLO87YLO0

#include "block.h"
#include "siteoperator/siteoperator.h"
#include "common/runlengthcompressed.h"
#include "mpopcompressed.h"
#include "mpoperator.h"


// TODO: this stuff doesn't belong here

typedef Block<SiteOperator> SiteBlock;

// do the flip conjugation of a SiteBlock, used for making an ancilla site that
// transforms as the conjugate
SiteBlock flip_conj(SiteBlock const& b);

class Lattice
{
   public:
      typedef SiteBlock value_type;
      typedef run_length_compressed<SiteBlock> data_type;
      typedef data_type::const_iterator const_iterator;

      Lattice();
      Lattice(SiteBlock const& s);
      Lattice(SiteBlock const& s, SiteBlock const& t);
      Lattice(SiteBlock const& s, SiteBlock const& t, SiteBlock const& u);
      Lattice(SiteBlock const& s, SiteBlock const& t, SiteBlock const& u, SiteBlock const& v);

      Lattice(SymmetryList const& sl, SiteBlock const& s);
      Lattice(SymmetryList const& sl, SiteBlock const& s, SiteBlock const& t);
      Lattice(SymmetryList const& sl, SiteBlock const& s, SiteBlock const& t, 
              SiteBlock const& u);
      Lattice(SymmetryList const& sl, SiteBlock const& s, SiteBlock const& t, 
              SiteBlock const& u, SiteBlock const& v);

      Lattice(int RepeatCount, Lattice const& l);
      Lattice(Lattice const& x1, Lattice const& x2);
      Lattice(Lattice const& x1, Lattice const& x2, Lattice const& x3);
      Lattice(Lattice const& x1, Lattice const& x2, Lattice const& x3, Lattice const& x4);

      Lattice(int Size, SiteBlock const& s);

      // Constructs a singleton lattice, setting the coordinate at the same time
      Lattice(SiteBlock const& s, std::string const& Coord);

      template <typename T>
      Lattice(SiteBlock const& s, T const& Coord);

      run_length_compressed<SiteBlock> const& data() const { return Data_; }

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

      SiteBlock const& operator[](int i) const;

      // Return the coordinate of the given site
      std::string const& coordinate_at_site(int Site) const;

       // return the site index of the given coordinate
      int site_at_coordinate(std::string const& s) const;

      // helper function to get the site index with a generic coordinate
      template <typename T1>
      int site_at_coordinate(T1 const& x1) const;

      template <typename T1, typename T2>
      int site_at_coordinate(T1 const& x1, T2 const& x2) const;

      template <typename T1, typename T2, typename T3>
      int site_at_coordinate(T1 const& x1, T2 const& x2, T3 const& x3) const;

      // fixes the coordinates to be an increasing sequence of integers 1,2,...
      void fix_coordinates();

      // fixes the coordinates to be a decreasing sequence of integers ...,2,1
      void fix_coordinates_reverse();

      // fixes the coordinates to be an increasing sequence of integers n, n+1, ...
      void fix_coordinates_starting_from(int n);

      // fixes the coordinates to be an increasing sequence of integers 1,2,...
      // as the first coordinate.  ie, if there are other coordinates, then
      // the new coordinate comes first, 1,X 2,X .... rather than the default
      // X,1 X,2 ....
      void fix_coordinates_prepend();

      void fix_coordinates_prepend_starting_from(int n);

      // assuming that all sites already have distinct coordinates, fix them as they are
      void fix_coordinates_unique();

      // Set the names of the coordinates from an iterator sequence; uses lexical_cast
      // to convert to a string representation
      template <typename FwdIter>
      void fix_coordinates_from_sequence(FwdIter start, FwdIter finish);

      // Shortcut functions for setting the coordinates of small lattices
      template <typename T>
      void fix_coordinates(T const& c1);
      template <typename T>
      void fix_coordinates(T const& c1, T const& c2);
      template <typename T>
      void fix_coordinates(T const& c1, T const& c2, T const& c3);
      template <typename T>
      void fix_coordinates(T const& c1, T const& c2, T const& c3, T const& c4);

      void fix_coordinates(char const* c1);
      void fix_coordinates(char const* c1, char const* c2);
      void fix_coordinates(char const* c1, char const* c2, char const* c3);
      void fix_coordinates(char const* c1, char const* c2, char const* c3, char const* c4);

   private:
      // Sets the 'reverse mapping' SiteAtCoord_ and sets CoordinatesFixed_ = true
      void Fixate(); 

      run_length_compressed<SiteBlock> Data_;
      std::vector<std::string> Coordinates_;
      bool CoordinatesFixed_;
      std::map<std::string, int> SiteAtCoord_;

   friend PStream::opstream& operator<<(PStream::opstream& out, Lattice const& L);
   friend PStream::ipstream& operator>>(PStream::ipstream& in, Lattice& L);
};

// Make a new lattice out of RepeatCount copies of x
Lattice repeat(Lattice const& x, int RepeatCount);

// Make a new lattice out of the join of x and y
Lattice join(Lattice const& x, Lattice const& y);
Lattice join(Lattice const& x, Lattice const& y, Lattice const& z);
Lattice join(Lattice const& x, Lattice const& y, Lattice const& z, Lattice const& w);
Lattice join(Lattice const& x, Lattice const& y, Lattice const& z, Lattice const& w,
	     Lattice const& v);

// Constructs a tensor product operator over the lattice L.
// If Operator does not exist at some site in L, the identity operator "I"
// is substituted.  Operator must transform as a scalar.
// The matrix basis for the operator is b at all sites.
MPOpCompressed 
CreateTensorProductOpCompressed(Lattice const& L, 
                                std::string const& Operator, 
                                BasisList const& b);

// Constructs a local operator acting over a lattice L.
// The bosonic/fermionic character of the operator is properly taken into account.
MPOpCompressed 
CreateMPOpCompressed(Lattice const& L, std::string const& Operator, int Site);

MPOperator CreateMPOperator(Lattice const& L, std::string const& Operator, int Site);

// Returns how the given operator transforms.  This assumes that the
// operator transforms the same way on all sites.
QuantumNumbers::QuantumNumber
LatticeSiteOperatorTransformsAs(Lattice const& L, std::string const& Operator);

// Given a symbolic operator that contains a matrix of strings representing
// local operators, and a basis (which must correspond with how the local operators transform),
// returns the corresponding MPOpComponent.
MPOpComponent ConstructFromSymbolic(LinearAlgebra::SparseMatrix<std::string> const& M,
                                    SiteBlock const& Block,
                                    BasisList const& B);

// Constructs an operator that is the sum of Op1 at every site.
// If Op1 does not exist at some sites, zero is substituted instead.
MPOperator
CreateRepeatedOperator(Lattice const& L, std::string const& Op1);

// Constructs an operator that is the sum of dot(Op1(i), Op2(i+1)) at every site.
// This takes into account the scaling factor for non-abelian quantum numbers.
// If Op1 or Op2 does not exist at some sites, zero is substituted instead.
MPOperator
CreateRepeatedOperator(Lattice const& L, std::string const& Op1, std::string const& Op2);

// Constructs an operator that is the sum of dot(Op1(i), Op2 (i+1) * Op3(i+2)) at every site.
// This takes into account the scaling factor for non-abelian quantum numbers.
// If Op1 or Op3 does not exist at some sites, zero is substituted instead.
// Op2 should be a scalar operator; this doesn't handle triple products.
MPOperator
CreateRepeatedOperator(Lattice const& L, std::string const& Op1, std::string const& Op2,
                       std::string const& Op3);

#include "lattice.cc"

#endif
