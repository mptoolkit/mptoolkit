// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// Copyright (C) 2026 Ian McCulloch <ian@qusim.net>
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

#include "common/terminal.h"
#include "lattice/infinitelattice.h"
#include "lattice/unitcelloperator.h"
#include "models/spin-u1.h"
#include "mp/copyright.h"
#include "pheap/pheap.h"
#include <boost/program_options.hpp>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>

namespace
{

/// \brief Matrix unit |row><column| in the physical spin-2 Sz basis.
SiteOperator matrix_unit(LatticeSite const& site, int row, int column)
{
   QuantumNumbers::QNConstructor<QuantumNumbers::U1> qn(site.GetSymmetryList());
   SiteOperator result(site.Basis1(), qn(row-column), LatticeCommute::Bosonic);
   result(std::to_string(row), std::to_string(column)) = 1.0;
   return result;
}

/// \brief SO(5) vector site with conserved physical Sz = 2 Q + Tz.
LatticeSite so5_site()
{
   LatticeSite site = SpinU1(2);
   site["P_plus"] = matrix_unit(site, 2, 2);
   site["P_minus"] = matrix_unit(site, -2, -2);
   site["P_zero"] = site["I"] - site["P_plus"] - site["P_minus"];
   site["Q"] = site["P_plus"] - site["P_minus"];
   site["Tz"] = matrix_unit(site, 1, 1) - matrix_unit(site, -1, -1);
   site["Tp"] = std::sqrt(2.0) * (matrix_unit(site, 1, 0) + matrix_unit(site, 0, -1));
   site["Tm"] = adjoint(site["Tp"]);
   site["T2"] = 2.0 * site["P_zero"];

   // Spherical vector: [Tz,r_a]=a r_a, [Q,r_a]=-r_a.
   // Their physical Sz charges are -1, -2 and -3, respectively.
   site["r_p"] = matrix_unit(site, 1, 2) + matrix_unit(site, -2, -1);
   site["r_0"] = matrix_unit(site, 0, 2) - matrix_unit(site, -2, 0);
   site["r_m"] = matrix_unit(site, -1, 2) + matrix_unit(site, -2, 1);
   site["r_p_dag"] = adjoint(site["r_p"]);
   site["r_0_dag"] = adjoint(site["r_0"]);
   site["r_m_dag"] = adjoint(site["r_m"]);

   // The on-site part of the collective invariant is nonzero.
   site["c_minus"] = -3.0 * matrix_unit(site, -2, 2);
   site["c_plus"] = adjoint(site["c_minus"]);

   std::array const descriptions = {
      std::pair{"P_plus", "Projector onto Sz=2, the SO(5) highest-weight ray"},
      std::pair{"P_minus", "Projector onto Sz=-2"},
      std::pair{"P_zero", "Projector onto the embedded SO(3) triplet, Sz=-1,0,1"},
      std::pair{"Q", "SO(2) charge (Sz^3-Sz)/6"},
      std::pair{"Tz", "Embedded SO(3) Cartan generator (4*Sz-Sz^3)/3"},
      std::pair{"Tp", "Embedded SO(3) raising operator P_zero*Sp*P_zero/sqrt(3)"},
      std::pair{"Tm", "Embedded SO(3) lowering operator, adjoint(Tp)"},
      std::pair{"T2", "On-site embedded SO(3) Casimir, 2*P_zero"},
      std::pair{"r_p", "Spherical broken lowerer with Tz component +1, physical Sz charge -1"},
      std::pair{"r_0", "Spherical broken lowerer with Tz component 0, physical Sz charge -2"},
      std::pair{"r_m", "Spherical broken lowerer with Tz component -1, physical Sz charge -3"},
      std::pair{"r_p_dag", "Adjoint of r_p"},
      std::pair{"r_0_dag", "Adjoint of r_0"},
      std::pair{"r_m_dag", "Adjoint of r_m"},
      std::pair{"c_minus", "On-site r_0^2-2*r_p*r_m = -3*|-2><2|"},
      std::pair{"c_plus", "Adjoint of c_minus"}
   };
   for (auto const& [name, description] : descriptions)
      site[name].set_description(description);
   return site;
}

/// \brief Add a finite-support wrapper f(cell){L} on L consecutive sites.
void finite_wrapper(UnitCell& cell, std::string const& name, std::string const& infinite_name)
{
   cell.func(name)(arg("L")) = "fsup(0,L," + infinite_name + ")";
   cell.func(name).set_description("Finite restriction of " + infinite_name
                                  + "; f(cell){L} acts on [cell,cell+L), L>=1");
}

}

int main(int argc, char** argv)
{
   namespace po = boost::program_options;
   try
   {
      std::string filename;
      po::options_description options("Allowed options", terminal::columns());
      options.add_options()
         ("help", "Show this help message")
         ("out,o", po::value(&filename), "Output lattice filename [required]");

      po::variables_map vm;
      po::store(po::command_line_parser(argc, argv).options(options)
                   .style(po::command_line_style::default_style
                          ^ po::command_line_style::allow_guessing).run(), vm);
      po::notify(vm);

      OperatorDescriptions descriptions;
      descriptions.set_description("SO(5) vector ferromagnet in the spin-2 Sz U(1) basis");
      descriptions.author("Ian McCulloch", "ian@qusim.net");
      descriptions.add_operators()
         ("H_SO5", "SO(5) ferromagnet; open-chain ground energy -(L-1)")
         ("H", "Alias for H_SO5")
         ("Q_tot", "Total SO(2) charge")
         ("Sz_tot", "Total physical spin-2 Sz = 2*Q_tot+Tz_tot")
         ("Tz_tot", "Total embedded SO(3) Cartan generator")
         ("Tp_tot", "Total embedded SO(3) raising operator")
         ("Tm_tot", "Total embedded SO(3) lowering operator")
         ("T2_tot", "Total SO(3) Casimir, including all on-site and pair terms")
         ("R_p", "Collective sum of r_p")
         ("R_0", "Collective sum of r_0")
         ("R_m", "Collective sum of r_m")
         ("C_minus", "Collective R_0^2-2*R_p*R_m; delta Q=-2, delta Sz=-4")
         ("C_plus", "Adjoint of C_minus");
      descriptions.add_functions()
         ("H_select", "H_SO5+lambda*T2_tot, lambda defaults to 1; finite-state selector");

      if (vm.count("help") || !vm.count("out"))
      {
         print_copyright(std::cerr);
         std::cerr << "usage: " << argv[0] << " [options]\n" << options << '\n'
                   << descriptions << '\n'
                   << "Local basis labels: -2,-1,0,1,2. Highest weight: 2.\n"
                   << "Finite example: fsup(0,16,C_minus), or C_minus_finite(0){L=16}.\n"
                   << "See docs/models/so5-spin2-u1.md for conventions and examples.\n";
         return vm.count("help") ? 0 : 1;
      }

      UnitCell cell(so5_site());
      UnitCellOperator I(cell, "I"), Sp(cell, "Sp"), Sm(cell, "Sm"), Sz(cell, "Sz");
      UnitCellOperator Q(cell, "Q"), Tz(cell, "Tz"), Tp(cell, "Tp"), Tm(cell, "Tm"), T2(cell, "T2");
      UnitCellOperator rp(cell, "r_p"), r0(cell, "r_0"), rm(cell, "r_m"), cm(cell, "c_minus");
      UnitCellOperator SS(cell, "SS"), Tdot(cell, "T_dot"), bond(cell, "h_SO5");
      SS = Sz(0)*Sz(1) + 0.5*(Sp(0)*Sm(1) + Sm(0)*Sp(1));
      SS.set_description("Physical spin-2 exchange S(0).S(1)");
      Tdot = Tz(0)*Tz(1) + 0.5*(Tp(0)*Tm(1) + Tm(0)*Tp(1));
      Tdot.set_description("Embedded SO(3) exchange T(0).T(1)");

      UnitCellMPO h = -Q(0)*Q(1) - Tdot(0);
      for (std::string const name : {"r_p", "r_0", "r_m"})
      {
         UnitCellOperator r(cell, name);
         h -= adjoint(r(0))*r(1) + r(0)*adjoint(r(1));
      }
      bond = h;
      bond.set_description("SO(5) bond: 1+5*SS/6-SS^2/9-SS^3/18; spectrum -1,1,4");

      for (std::string const name : {"H_SO5", "Q_tot", "Sz_tot", "Tz_tot", "Tp_tot", "Tm_tot",
                                     "T2_tot", "R_p", "R_0", "R_m", "C_minus", "C_plus"})
         finite_wrapper(cell, name + "_finite", name);
      finite_wrapper(cell, "H_finite", "H_SO5");
      cell.func("H_select_finite")(arg("L"), arg("lambda")=1)
         = "fsup(0,L,H_SO5)+lambda*fsup(0,L,T2_tot)";
      cell.func("H_select_finite").set_description("Finite H_SO5+lambda*T2_tot on L sites");

      InfiniteLattice lattice(&cell);
      lattice["H_SO5"] = sum_unit(h);
      lattice["H"] = lattice["H_SO5"];
      lattice["Q_tot"] = sum_unit(Q(0));
      lattice["Sz_tot"] = sum_unit(Sz(0));
      lattice["Tz_tot"] = sum_unit(Tz(0));
      lattice["Tp_tot"] = sum_unit(Tp(0));
      lattice["Tm_tot"] = sum_unit(Tm(0));
      lattice["R_p"] = sum_unit(rp(0));
      lattice["R_0"] = sum_unit(r0(0));
      lattice["R_m"] = sum_unit(rm(0));

      lattice["T2_tot"] = sum_unit(T2(0))
         + 2.0*sum_string(Tz(0), I(0), Tz(0))
         + sum_string(Tp(0), I(0), Tm(0))
         + sum_string(Tm(0), I(0), Tp(0));

      // sum_string defaults to a scalar total charge. C_minus instead has
      // physical Sz charge -4, so its string terms need the explicit target.
      QuantumNumbers::QNConstructor<QuantumNumbers::U1> qn(cell.GetSymmetryList());
      BasicTriangularMPO cminus = sum_unit(cm(0))
         + 2.0*sum_string(r0(0), I(0), r0(0), 1, qn(-4))
         - 2.0*sum_string(rp(0), I(0), rm(0), 1, qn(-4))
         - 2.0*sum_string(rm(0), I(0), rp(0), 1, qn(-4));
      lattice["C_minus"] = cminus;
      lattice["C_plus"] = adjoint(cminus);
      lattice.func("H_select")(arg("lambda")=1) = "H_SO5+lambda*T2_tot";

      lattice.set_command_line(argc, argv);
      lattice.set_operator_descriptions(descriptions);
      pheap::ExportObject(filename, lattice);
   }
   catch (std::exception const& error)
   {
      std::cerr << "Exception: " << error.what() << '\n';
      return 1;
   }
}
