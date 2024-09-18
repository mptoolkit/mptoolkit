// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// attic/matrixproduct-obsolete/periodicoperator.h
//
// Copyright (C) 2004-2017 Ian McCulloch <ian@qusim.net>
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
//
// Periodic matrix product operator.
// Perhaps we want to extend this to a non-trivial unit cell at some stage?
//

#if !defined(MPOPPERIODIC_H_JDCHJKEHY589758YUER89H489)
#define MPOPPERIODIC_H_JDCHJKEHY589758YUER89H489

#include "pstream/pstream.h"
#include "quantumnumbers/quantumnumber.h"
#include "mpopcomponent.h"
#include "linearoperator.h"

class MpOpPeriodic
{
   public:
      typedef MPOpComponent::mapped_type mapped_type;
      typedef MPOpComponent::basis1_type basis_type;

      MpOpPeriodic(int Size, MPOpComponent const& Data);
      MpOpPeriodic(int Size, SimpleOperator const& Q, MPOpComponent const& Data);

      QuantumNumber TransformsAs() const { return adjoint(Q_.TransformsAs()); }

      // size == 0 denotes indefinite size (or infinite)
      int size() const { return Size_; }

      void set_size(int Size) { Size_ = Size; }

      SimpleOperator const& Q() const { return Q_; }
      MpOpComponent const& data() const { return Data_; }

      MpOpComponent& data() { return Data_; }

      mapped_type& operator[](QuantumNumber const& q);
      mapped_type operator[](QuantumNumber const& q) const;

      BasisList const& SiteBasis() const { return Data_.SiteBasis(); }
      basis_type Basis() const { return Data_.Basis1(); }

   private:
      int Size_;
      SimpleOperator Q_;
      MPOpComponent Data_;
};

LinearOperator ConvertPeriodicToOpen(MPOpPeriodic const& Op);


#endif
