// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp-algorithms/ea-dmrg.h
//
// Copyright (C) 2023 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#if !defined(MPTOOLKIT_MP_ALGORITHMS_EA_DMRG_H)
#define MPTOOLKIT_MP_ALGORITHMS_EA_DMRG_H

#include "ef-matrix.h"
#include "mps/packunpack.h"
#include "wavefunction/ea.h"

// Struct to hold the settings to initialize the EA_DMRG class.
struct EA_DMRGSettings
{
   double Tol = 1e-10;
   double GMRESTol = 1e-13;
   double UnityEpsilon = DefaultEigenUnityEpsilon;
   bool Quiet = false;
   int Verbose = 0;
};

class EA_DMRG
{
   public:
      EA_DMRG() {}

      EA_DMRG(EAWavefunction const& Psi_, BasicTriangularMPO const& HamMPO, EA_DMRGSettings Settings_);

      // Return the current wavefunction.
      EAWavefunction Wavefunction() const;

      // Solve for the window components with optimal energy for the current site.
      void SolveCurrentSite();

      // Move the orthogonality center one site to the left/right.
      void IterateLeft();
      void IterateRight();

      // Sweep left and then sweep right.
      void SweepLR();

   private:
      EAWavefunction Psi;
      EFMatrix EF;

      std::vector<LinearWavefunction> WindowVec;
      std::vector<LinearWavefunction::iterator> WIVec;
      int Site;
      int LeftStop;
      int RightStop;
      int Sweep = 0;

      double Tol;
      bool Quiet;
      int Verbose;
};

// Functor for the effective Hamiltonian acting at position i in each window.
class HEff
{
   public:
      HEff(EFMatrix* EF_, int Site_);

      std::deque<StateComponent> operator()(std::deque<StateComponent> WDeque);

   private:
      EFMatrix* EF;
      int Site;
      std::complex<double> ExpIK;
};

// A wrapper for the effective Hamiltonian functor for using with ARPACK.
class PackHEff
{
   public:
      PackHEff(HEff* H_, std::deque<StateComponent> WDeque);

      // Apply HEff to input vector.
      void operator()(std::complex<double> const* In_, std::complex<double>* Out_) const;

      // Convert raw array to deque of StateComponents.
      std::deque<StateComponent> unpack(std::complex<double> const* In_) const;

      // Convert deque of StateComponents to an array readable by the ARPACK wrapper.
      void pack(std::deque<StateComponent> WDeque, std::complex<double>* Out_) const;

      // Size of the raw array.
      int size() const { return Size; }

   private:
      std::deque<PackStateComponent> Pack;
      HEff* H;
      int Size;
};

#endif
