// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// mp/mp-entropy.cpp
//
// Copyright (C) 2025 Jesse Osborne <j.osborne@uqconnect.edu.au>
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

#include "common/environment.h"
#include "common/formatting.h"
#include "common/prog_opt_accum.h"
#include "common/prog_options.h"
#include "common/terminal.h"
#include "interface/inittemp.h"
#include "mp/copyright.h"
#include "wavefunction/mpwavefunction.h"
#include "mps/state_component.h"

#include <algorithm>
#include <cctype>
#include <iterator>
#include <set>
#include <type_traits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace
{
   std::string trim(std::string s)
   {
      auto const not_space = [](unsigned char c) { return !std::isspace(c); };
      s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
      s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
      return s;
   }

   std::vector<std::string> split_commas(std::string const& token)
   {
      std::vector<std::string> parts;
      std::size_t start = 0;
      while (start <= token.size())
      {
         std::size_t comma = token.find(',', start);
         if (comma == std::string::npos)
         {
            parts.emplace_back(token.substr(start));
            break;
         }
         parts.emplace_back(token.substr(start, comma - start));
         start = comma + 1;
      }
      if (parts.empty())
         parts.push_back(token);
      return parts;
   }

   // TODO: In C++23, replace int_range with std::views::iota.

   template <typename Int>
   struct int_range
   {
      static_assert(std::is_integral<Int>::value, "int_range requires integral type");

      Int first_;
      Int last_;

      struct iterator
      {
         using iterator_category = std::forward_iterator_tag;
         using value_type        = Int;
         using difference_type   = Int;
         using pointer           = void;
         using reference         = Int;

         Int current_;

         Int operator*() const noexcept { return current_; }
         iterator& operator++() noexcept { ++current_; return *this; }
         bool operator!=(iterator const& other) const noexcept { return current_ != other.current_; }
      };

      iterator begin() const noexcept { return iterator{first_}; }
      iterator end()   const noexcept { return iterator{static_cast<Int>(last_ + 1)}; } // inclusive range
   };

   std::vector<int> parse_site_token(std::string token)
   {
      token = trim(std::move(token));
      if (token.empty())
         return {};

      std::string compact;
      compact.reserve(token.size());
      for (unsigned char c : token)
      {
         if (!std::isspace(c))
            compact.push_back(static_cast<char>(c));
      }

      if (compact.empty())
         return {};

      auto const parse_integer = [](std::string const& value) -> int
      {
         std::size_t processed = 0;
         int const result = std::stoi(value, &processed, 10);
         if (processed != value.size())
            throw std::runtime_error("Invalid site number: " + value);
         return result;
      };

      std::size_t range_pos = std::string::npos;
      for (std::size_t i = 1; i < compact.size(); ++i)
      {
         if (compact[i] == '-')
         {
            range_pos = i;
            break;
         }
      }

      if (range_pos == std::string::npos)
      {
         return {parse_integer(compact)};
      }

      std::string const start_text = compact.substr(0, range_pos);
      std::string const end_text = compact.substr(range_pos + 1);

      if (start_text.empty() || end_text.empty())
         throw std::runtime_error("Invalid site specification: " + token);

      int const start = parse_integer(start_text);
      int const end = parse_integer(end_text);

      if (end < start)
      {
         throw std::runtime_error(
            "Invalid site range: " + std::to_string(start) + "-" + std::to_string(end) +
            " (end is less than start)");
      }

      std::vector<int> result;
      result.reserve(static_cast<std::size_t>(end - start + 1));
      for (int value = start; value <= end; ++value)
      {
         result.push_back(value);
      }

      return result;
   }

   std::vector<int_range<int>> parse_sites(std::vector<std::string> const& tokens)
   {
      std::set<int> unique_sites;
      for (auto const& token : tokens)
      {
         for (auto part : split_commas(token))
         {
            auto expanded = parse_site_token(std::move(part));
            unique_sites.insert(expanded.begin(), expanded.end());
         }
      }

      std::vector<int_range<int>> ranges;
      if (unique_sites.empty())
         return ranges;

      auto it = unique_sites.begin();
      while (it != unique_sites.end())
      {
         int const start = *it;
         int end = start;
         ++it;
         while (it != unique_sites.end() && *it == end + 1)
         {
            end = *it;
            ++it;
         }
         ranges.push_back(int_range<int>{start, end});
      }

      return ranges;
   }
}

namespace prog_opt = boost::program_options;

int main(int argc, char** argv)
{
   try
   {
      int Verbose = 0;
      std::string Filename;
      int MaxEigenvalues =-1;
      bool Base2 = false;
      bool ShowDensity = false;
      bool ShowDegen = false;
      bool Quiet = false;
      std::vector<std::string> SiteTokens;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "Show this help message")
         ("density-matrix,d", prog_opt::bool_switch(&ShowDensity), "Show the eigenspectrum of the density matrix")
         ("degen", prog_opt::bool_switch(&ShowDegen), "Show degeneracies in the density matrix as repeated eigenvalues (implies -d)")
         ("limit,l", prog_opt::value<int>(&MaxEigenvalues), "Limit the density matrix display to N eigenvalues (implies -d)")
         ("base2,2", prog_opt::bool_switch(&Base2), "Show the entropy using base 2 instead of base e")
         ("quiet", prog_opt::bool_switch(&Quiet), "Do not show column headings")
         ("verbose,v",  prog_opt_ext::accum_value(&Verbose), "Increase verbosity (can be used more than once)")
         ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("psi", prog_opt::value(&Filename), "psi")
         ("sites", prog_opt::value<std::vector<std::string>>(&SiteTokens)->multitoken(), "sites")
         ;

      prog_opt::positional_options_description p;
      p.add("psi", 1);
      p.add("sites", -1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);

      if (vm.count("help") || vm.count("sites") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: " << basename(argv[0]) << " [options] <psi> <sites>" << std::endl;
         std::cerr << desc << std::endl;
         std::cerr << "Specify sites using integers or inclusive ranges (e.g. 0-3)." << std::endl;
         return 1;
      }

      std::cout.precision(getenv_or_default("MP_PRECISION", 14));
      std::cerr.precision(getenv_or_default("MP_PRECISION", 14));

      if (vm.count("limit"))
         ShowDensity = true;

      if (ShowDegen)
         ShowDensity = true;

      auto SiteRanges = parse_sites(SiteTokens);
      if (SiteRanges.empty())
      {
         std::cerr << "No sites specified." << std::endl;
         return 1;
      }

      pvalue_ptr<MPWavefunction> PsiPtr = pheap::OpenPersistent(Filename, mp_pheap::CacheSize(), true);

      FiniteWavefunctionLeft Psi = PsiPtr->get<FiniteWavefunctionLeft>();

      SimpleOperator Rho = make_vacuum_simple(Psi.GetSymmetryList());
      for (auto const& range : SiteRanges)
      {
         StateComponent C = make_vacuum_state(Psi.GetSymmetryList());
         for (int site : range)
         {
            C = local_tensor_prod(C, Psi[site]);
         }
         C = prod(C, Psi.lambda(range.last_ + 1));
         Rho = tensor_prod(Rho, trace_prod(C, herm(C)));
      }

      DensityMatrix<SimpleOperator> DM(Rho);
      if (!ShowDensity)
      {
         // Just print the entropy
         std::cout << DensityEntropy(DM.begin(), DM.end(), DM.EigenSum(), Base2) << std::endl;
      }
      else
      {
         // Print the full density matrix eigenspectrum
         DM.DensityMatrixReport(std::cout, MaxEigenvalues, Base2, ShowDegen, Quiet);
         std::cout << std::endl;
      }

      pheap::Shutdown();
   }
   catch (prog_opt::error& e)
   {
      std::cerr << "Exception while processing command line options: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (pheap::PHeapCannotCreateFile& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      if (e.Why == "File exists")
         std::cerr << "Note: Use --force (-f) option to overwrite." << std::endl;
      return 1;
   }
   catch (std::exception& e)
   {
      std::cerr << "Exception: " << e.what() << std::endl;
      pheap::Cleanup();
      return 1;
   }
   catch (...)
   {
      std::cerr << "Unknown exception!" << std::endl;
      pheap::Cleanup();
      return 1;
   }
}
