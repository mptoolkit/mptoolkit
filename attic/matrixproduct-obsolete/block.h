// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// matrixproduct-obsolete/block.h
//
// Copyright (C) 2004-2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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

/*
  Created 2004-09-05 Ian McCulloch

  Class to hold a group of operators.  Typically, these will be site operators.
  In 'classic' DMRG it could also hold the block operators needed to construct
  the Hamiltonian.
*/

#if !defined(BLOCK_H_853YF987RHVHYC85HYT87HGGO2)
#define BLOCK_H_853YF987RHVHYC85HYT87HGGO2

#include "pstream/pstream.h"
#include "pheap/pvalueptr.h"
#include "quantumnumbers/symmetrylist.h"
#include <map>

template <typename OperatorT>
class Block;

using QuantumNumbers::SymmetryList;

template <typename OperatorT>
PStream::opstream& operator<<(PStream::opstream& out, Block<OperatorT> const& B);

template <typename OperatorT>
PStream::ipstream& operator>>(PStream::ipstream& in, Block<OperatorT>& B);

template <typename OperatorT>
class Block
{
   private:
      typedef std::map<std::string, OperatorT> DataType;

   public:
      typedef typename DataType::value_type     value_type;
      typedef typename DataType::iterator       iterator;
      typedef typename DataType::const_iterator const_iterator;
      typedef typename OperatorT::basis1_type   basis1_type;
      typedef typename OperatorT::basis2_type   basis2_type;

      Block() : Data(new DataType()) {}

      // precondition: !empty()
      SymmetryList GetSymmetryList() const;

      // precondition: !empty()
      basis1_type const& Basis1() const;
      basis2_type const& Basis2() const;

      bool empty() const { return Data->empty(); }

      const_iterator begin() const { return Data->begin(); }
      const_iterator end() const { return Data->end(); }

      OperatorT& operator[](std::string const& s) { return (*Data.mutate())[s]; }
      OperatorT const& operator[](std::string const& s) const;

      const_iterator find(std::string const& s) const { return Data->find(s); }

      bool exists(std::string const& s) const { return Data->find(s) != Data->end(); }

      void CoerceSymmetryList(QuantumNumbers::SymmetryList const& sl);

   private:
      typedef pvalue_ptr<DataType> ptr_type;
      ptr_type Data;

      friend PStream::opstream& operator<< <>(PStream::opstream& out, Block const& B);

      friend PStream::ipstream& operator>> <>(PStream::ipstream& in, Block& B);
};

#include "block.cc"

#endif
