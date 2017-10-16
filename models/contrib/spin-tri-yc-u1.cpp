// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// models/contrib/spin-tri-yc-u1.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
// Copyright (C) 2015,2016 Seyed N. Saadatmand <s.saadatmand@uq.edu.au>
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
// YC configuration of a triangular lattice.
// The default unit-cell size is the width value.
//
// Example for a width-6 lattice (site numbers in brackets are periodic repeats in the vertical
// direction (i.e. top-left (5) is the same site as the bottom-left 5).
// Sites 6,7,8,9,10,11 are the second unit cell, e.g. 6 is (1)[0].
//
//      (17)
//  (11) |
//(5) | 12
// |  6< |
// 0< | 13
// |  7< |
// 1< | 14
// |  8< |
// 2< | 15
// |  9< |
// 3< | 16
// | 10< |
// 4< | 17
// | 11< |
// 5< | (12)
// | (6)
//(0)
//

#include "pheap/pheap.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "mp/copyright.h"
#include "models/spin-u1.h"
#include "common/terminal.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

// struct to get the prefactors of the flux terms.  The complete flux term at angle theta is
// uniform_part + cos(theta)*real_part + i*sin(theta)*imag_part.
// The prefactors give the contribution for hopping from site
struct FluxPhase
{
   FluxPhase(int w_) : w(w_) {}

   double uniform(int s1, int s2)
   {
      CHECK(std::abs(s1-s2) <= 2 || (w - abs(s1-s2)) <= 2);
      return (std::abs(s1-s2) <= 2) ? 1 : 0;
   }

   double real(int s1, int s2)
   {
      // if we don't cross the boundary, then there is no flux contribution
      if (std::abs(s1-s2) <= 2) return 0;
      // else
      return 1;
   }

   double imag(int s1, int s2)
   {
      // if we don't cross the boundary, then there is no flux contribution
      if (std::abs(s1-s2) <= 2) return 0;
      // else if we are going in the +ve direction (ie, s2 < s1, becase we wrapped around)
      // then the flux is +1, otherwise -1
      return (s2 < s1) ? 1 : -1;
   }
   int w;
};

int main(int argc, char** argv)
{
   try
   {
      half_int Spin = 0.5;
      int w = 3;
      std::string FileName;
      bool NoReflect = false;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
         ("Spin,S", prog_opt::value(&Spin), "magnitude of the spin [default 0.5]")
         ("width,w", prog_opt::value(&w), "width of the cylinder [default 3]")
         ("out,o", prog_opt::value(&FileName), "output filename [required]")
         ;

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(desc).style(prog_opt::command_line_style::default_style ^
                                          prog_opt::command_line_style::allow_guessing).
                      run(), vm);
      prog_opt::notify(vm);

      OperatorDescriptions OpDescriptions;
      OpDescriptions.add_operators()
         ("H_J1"                   , "nearest neighbor Heisenberg spin exchange")
         ("H_J1z"                  , "nearest neighbor Sz-Sz coupling")
         ("H_J1p"                  , "nearest neighbor XY coupling")
         ("H_J1p_flux_uniform"     , "nearest neighbor XY spin exchange with flux (uniform part)")
         ("H_J1p_flux_cos"         , "nearest neighbor XY spin exchange with flux (real part)")
         ("H_J1p_flux_sin"         , "nearest neighbor XY spin exchange with flux (imaginary part)")
         ("H_J2"                   , "next-nearest neighbor Heisenberg spin exchange")
         ("H_J2z"                  , "next-nearest neighbor Sz-Sz coupling")
         ("H_J2p"                  , "next-nearest neighbor XY coupling")
         ("H_J2p_flux_uniform"     , "next-nearest neighbor XY spin exchange with flux (uniform part)")
         ("H_J2p_flux_cos"         , "next-nearest neighbor XY spin exchange with flux (real part)")
         ("H_J2p_flux_sin"         , "next-nearest neighbor XY spin exchange with flux (imaginary part)")
         ;

      OpDescriptions.add_functions()
         ("H_J1_flux"                            , "nearest neighbor Heisenberg spin exchange with flux phase")
         ("H_J2_flux"                            , "nearest neighbor Heisenberg spin exchange with flux phase")
         ;

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << basename(argv[0]) << " [options]\n";
         std::cerr << desc << '\n';
         std::cerr << "Constructs a triangular lattice in the YC configuration with lattice vector (0,1).\n";
         std::cerr << OpDescriptions << '\n';
            ;
         return 1;
      }

      LatticeSite Site = SpinU1(Spin);
      UnitCell Cell = repeat(Site, w);
      InfiniteLattice Lattice(&Cell);

      FluxPhase phase(w);

      // Add some operators on the unit cell, for convenience
      UnitCellOperator Sp(Cell, "Sp"), Sm(Cell, "Sm"), Sz(Cell, "Sz");
      UnitCellOperator I(Cell, "I"), Trans(Cell, "Trans");

      for (int i = 0; i < w; ++i)
      {
         Sp += Sp[i];     // total S+ on a leg of cylinder
         Sm += Sm[i];     // total S- on a leg of cylinder
         Sz += Sz[i];     // total Sz on a leg of cylinder
      }

      // Translation operator
      Trans = I(0);
      for (int i = 0; i < w-1; ++i)
      {
         Trans = Trans(0) * Cell.swap_gate_no_sign(i, i+1);
      }

      // Add it to the lattice
      Lattice["Ty"] = prod_unit_left_to_right(UnitCellMPO(Trans(0)).MPO(), w);

      // Reflection operator
      UnitCellMPO Ry = I(0);
      if (!NoReflect)
      {
         for (int c = 0; c < w; ++c)
         {
            UnitCellMPO ThisR = I(0);
            // get the 'pivot' site/bond that we reflect about
            int const p1 = c/2;
            int const p2 = (c+1)/2;

            // if we're reflecting about a bond, do that first
            if (p1 != p2)
               ThisR = ThisR * Cell.swap_gate_no_sign(p1,p2);

            int i1 = (p1+w-1)%w;
            int i2 = (p2+1)%w;

            while (i1 != p1 + w/2)
            {
               ThisR = ThisR * Cell.swap_gate_no_sign(i1,i2);
               i1 = (i1+w-1)%w;
               i2 = (i2+1)%w;
            }

            ThisR.translate(c*w);
            Ry = Ry * ThisR;
         }
      }

      // Add it to the lattice
      Lattice["Ry"] = prod_unit_left_to_right(Ry.MPO(), w*w);

      // Hamiltonian operators, constructed as a sum over terms that act on a unit cell
      UnitCellMPO H1, H1z, H1p, H1p_flux_uniform, H1p_flux_cos, H1p_flux_sin;
      UnitCellMPO H2, H2z, H2p, H2p_flux_uniform, H2p_flux_cos, H2p_flux_sin;

      for (int i = 0; i < w; ++i)
      {
         //
         // Nearest-neighbour bonds
         //

         // Z-component
         UnitCellMPO n1_z_vertical = Sz(0)[i]*Sz(0)[(i+1)%w];    // vertical
         UnitCellMPO n1_z_up       = Sz(0)[i]*Sz(1)[i];          // +60 degrees
         UnitCellMPO n1_z_down     = Sz(0)[i]*Sz(1)[(i+1)%w];    // -60 degrees

         // Sp*Sm components
         UnitCellMPO n1_plus_vertical = 0.5 * Sp(0)[i]*Sm(0)[(i+1)%w];
         UnitCellMPO n1_plus_up       = 0.5 * Sp(0)[i]*Sm(1)[i];
         UnitCellMPO n1_plus_down     = 0.5 * Sp(0)[i]*Sm(1)[(i+1)%w];

         // Sm*Sp components
         UnitCellMPO n1_minus_vertical = 0.5 * Sm(0)[i]*Sp(0)[(i+1)%w];
         UnitCellMPO n1_minus_up       = 0.5 * Sm(0)[i]*Sp(1)[i];
         UnitCellMPO n1_minus_down     = 0.5 * Sm(0)[i]*Sp(1)[(i+1)%w];

         // Hamiltonian

         H1z += n1_z_vertical + n1_z_up + n1_z_down;

         H1p += n1_plus_vertical + n1_minus_vertical
            + n1_plus_up + n1_minus_up
            + n1_plus_down + n1_minus_down;

         H1p_flux_uniform += phase.uniform(i,i+1) * (n1_plus_vertical + n1_minus_vertical)
            + (n1_plus_up + n1_minus_up)
            + phase.uniform(i,i+1) * (n1_plus_down + n1_minus_down);

         H1p_flux_cos += phase.real(i,i+1) *
            (n1_plus_vertical + n1_minus_vertical + n1_plus_down + n1_minus_down);

         H1p_flux_sin += phase.imag(i,i+1) *
            (n1_plus_vertical - n1_minus_vertical + n1_plus_down - n1_minus_down);

         //
         // Next-nearest neighbor bonds
         //

         // Z-component
         UnitCellMPO n2_z_horizontal = Sz(0)[i]*Sz(1)[(i+w-1)%w];
         UnitCellMPO n2_z_up         = Sz(0)[i]*Sz(1)[(i+2)%w];
         UnitCellMPO n2_z_down       = Sz(0)[i]*Sz(1)[(i+w-1)%w];

         // Sp*Sm components
         UnitCellMPO n2_plus_horizontal = 0.5 * Sp(0)[i]*Sm(1)[(i+w-1)%w];
         UnitCellMPO n2_plus_up         = 0.5 * Sp(0)[i]*Sm(1)[(i+2)%w];
         UnitCellMPO n2_plus_down       = 0.5 * Sp(0)[i]*Sm(1)[(i+w-1)%w];

         // Sm*Sp components
         UnitCellMPO n2_minus_horizontal = 0.5 * Sm(0)[i]*Sp(1)[(i+w-1)%w];
         UnitCellMPO n2_minus_up         = 0.5 * Sm(0)[i]*Sp(1)[(i+2)%w];
         UnitCellMPO n2_minus_down       = 0.5 * Sm(0)[i]*Sp(1)[(i+w-1)%w];

         // Hamiltonian

         H2z += n2_z_horizontal + n2_z_up + n2_z_down;

         H2p += n2_plus_horizontal + n2_minus_horizontal
            + n2_plus_up + n2_minus_up
            + n2_plus_down + n2_minus_down;

         H2p_flux_uniform += phase.uniform(i,i-1) * (n2_plus_horizontal + n2_minus_horizontal)
            + phase.uniform(i,i+2) * (n2_plus_up + n2_minus_up)
            + phase.uniform(i, i-1) * (n2_plus_down + n2_minus_down);

         H2p_flux_cos += phase.real(i,i-1) * (n2_plus_horizontal + n2_minus_horizontal)
            + phase.real(i,i+2) * (n2_plus_up + n2_minus_up)
            + phase.real(i, i-1) * (n2_plus_down + n2_minus_down);

         H2p_flux_sin += phase.imag(i,i-1) * (n2_plus_horizontal + n2_minus_horizontal)
            + phase.imag(i,i+2) * (n2_plus_up + n2_minus_up)
            + phase.imag(i, i-1) * (n2_plus_down + n2_minus_down);
      }

      // Add the operators to the lattice

      Lattice["H_J1"] = sum_unit(H1z + H1p);
      Lattice["H_J1z"] = sum_unit(H1z);
      Lattice["H_J1p"] = sum_unit(H1p);
      Lattice["H_J1p_flux_uniform"] = sum_unit(H1p_flux_uniform);
      Lattice["H_J1p_flux_cos"] = sum_unit(H1p_flux_cos);
      Lattice["H_J1p_flux_sin"] = sum_unit(H1p_flux_sin);

      Lattice["H_J2"] = sum_unit(H2z + H2p);
      Lattice["H_J2z"] = sum_unit(H2z);
      Lattice["H_J2p"] = sum_unit(H2p);
      Lattice["H_J2p_flux_uniform"] = sum_unit(H2p_flux_uniform);
      Lattice["H_J2p_flux_cos"] = sum_unit(H2p_flux_cos);
      Lattice["H_J2p_flux_sin"] = sum_unit(H2p_flux_sin);

      Lattice.func("H_J1_flux")(arg("theta")) = "H_J1z + H_J1p_flux_uniform + cos(theta)*H_J1p_flux_cos + i*sin(theta)*H_J1p_flux_sin";
      Lattice.func("H_J2_flux")(arg("theta")) = "H_J2z + H_J2p_flux_uniform + cos(theta)*H_J2p_flux_cos + i*sin(theta)*H_J2p_flux_sin";

      // save the lattice
      pheap::ExportObject(FileName, Lattice);
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << '\n';
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!\n";
      return 1;
   }
}
