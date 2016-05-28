// -*- C++ -*-
//
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// mp/mp-local-fourpoint-ex.cpp
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// ENDHEADER

#include "matrixproduct/lattice.h"
#include "matrixproduct/mpwavefunction.h"
#include "matrixproduct/centerwavefunction.h"
#include "matrixproduct/mpoperatorlist.h"
#include "quantumnumbers/all_symmetries.h"
#include "pheap/pheap.h"
#include <iostream>
#include <iomanip>
#include "mp/copyright.h"
#include "common/proccontrol.h"
#include <boost/optional.hpp>
#include <boost/none.hpp>
#include "common/terminal.h"
#include "common/environment.h"
#include <boost/program_options.hpp>

namespace prog_opt = boost::program_options;

// A structure for keeping track of all sets (i,j,k,l) where we want to evaluate
// the correlator
struct PositionMap
{
   PositionMap(int Size);

   void set_imin(int n) { imin_ = n; }
   void set_imax(int n) { imax_ = n; }
   void set_jmin(int n) { jmin_ = n; }
   void set_jmax(int n) { jmax_ = n; }
   void set_kmin(int n) { kmin_ = n; }
   void set_kmax(int n) { kmax_ = n; }
   void set_lmin(int n) { lmin_ = n; }
   void set_lmax(int n) { lmax_ = n; }

   void set_ijmin(int n) { ijmin_ = n; }
   void set_ijmax(int n) { ijmax_ = n; }
   void set_ikmin(int n) { ikmin_ = n; }
   void set_ikmax(int n) { ikmax_ = n; }
   void set_ilmin(int n) { ilmin_ = n; }
   void set_ilmax(int n) { ilmax_ = n; }
   void set_jkmin(int n) { jkmin_ = n; }
   void set_jkmax(int n) { jkmax_ = n; }
   void set_jlmin(int n) { jlmin_ = n; }
   void set_jlmax(int n) { jlmax_ = n; }
   void set_klmin(int n) { klmin_ = n; }
   void set_klmax(int n) { klmax_ = n; }

   // This must be called after setting the limits, before calling any of the functions
   // below.  Currently this resets some limits based on the constraints
   void fix();

   // returns true if the combination Op1(Loc1)*Op2(Loc2) is required   
   bool Need12(int Loc1, int Loc2) const;

   // returns true if Op1(Loc1) * Op2(n) is needed for any n > Loc2
   bool Need1Later(int Loc1, int Loc2) const;

   // returns true if Op1(Loc1)*Op2(Loc2)*Op3(n) is needed for any n > Loc3
   bool Need12Later(int Loc1, int Loc2, int Loc3) const;

   // returns true if Op4(Loc4) is required at a position Loc3 < n
   bool Need4Later(int Loc3, int Loc4) const;

   // returns true if the pair Op3(Loc3)*Op4(Loc4) is needed.
   bool Need34(int Loc3, int Loc4) const;

   bool ShouldEvaluate(int i, int j, int k, int l) const;

   // returns the lowest site number that Op1 is needed for
   int FirstSite() const;

   // returns the highest site number that Op4 is needed for
   int LastSite() const;

   int Size_;
   int imin_, imax_, jmin_, jmax_, kmin_, kmax_, lmin_, lmax_;
   int ijmin_, ijmax_, ikmin_, ikmax_, ilmin_, ilmax_;
   int jkmin_, jkmax_, jlmin_, jlmax_;
   int klmin_, klmax_;
};

PositionMap::PositionMap(int Size) 
   : Size_(Size),
     imin_(1), imax_(Size-3), jmin_(2), jmax_(Size-2), kmin_(3), kmax_(Size-1), lmin_(4), lmax_(Size),
     ijmin_(1), ijmax_(Size-3), ikmin_(2), ikmax_(Size-2), ilmin_(3), ilmax_(Size-1),
     jkmin_(1), jkmax_(Size-3), jlmin_(2), jlmax_(Size-2),
     klmin_(1), klmax_(Size-3)
{
}

bool PositionMap::Need12(int i, int j) const
{
   if (i < imin_ || i > imax_) return false;
   if (j < jmin_ || j > jmax_) return false;
   if (j-i < ijmin_ || j-i > ijmax_) return false;
   return true;
}

bool PositionMap::Need1Later(int i, int j) const
{
   if (i < imin_ || i > imax_) return false;
   if (j > jmax_) return false;
   if (j-i > ijmax_) return false;
   return true;
}

bool PositionMap::Need12Later(int i, int j, int k) const
{
   if (!Need12(i,j)) return false;
   if (k > kmax_) return false;
   if (k-i > ikmax_) return false;
   if (k-j > jkmax_) return false;
   return true;
}

bool PositionMap::Need4Later(int k, int l) const
{
   if (l < lmin_ || l > lmax_) return false;
   if (k < kmin_) return false;
   if (l-k > klmax_) return false;
   return true;
}

bool PositionMap::Need34(int k, int l) const
{
   if (l < lmin_ || l > lmax_) return false;
   if (k < kmin_ || k > kmax_) return false;
   if (l-k < klmin_ || l-k > klmax_) return false;
   return true;
}

bool PositionMap::ShouldEvaluate(int i, int j, int k, int l) const
{
   if (i < imin_ || i > imax_) return false;
   if (j < jmin_ || j > jmax_) return false;
   if (l < lmin_ || l > lmax_) return false;
   if (k < kmin_ || k > kmax_) return false;
   if (j-i < ijmin_ || j-i > ijmax_) return false;
   if (k-i < ikmin_ || k-i > ikmax_) return false;
   if (l-i < ilmin_ || l-i > ilmax_) return false;
   if (k-j < jkmin_ || k-j > jkmax_) return false;
   if (l-j < jlmin_ || l-j > jlmax_) return false;
   if (l-k < klmin_ || l-k > klmax_) return false;
   return true;
}

void PositionMap::fix()
{
   // sanity check the pairwise minima
   if (ijmin_ < 1) ijmin_ = 1;
   if (ikmin_ < 1) ikmin_ = 1;
   if (ilmin_ < 1) ilmin_ = 1;
   if (jkmin_ < 1) jkmin_ = 1;
   if (jlmin_ < 1) jlmin_ = 1;
   if (klmin_ < 1) klmin_ = 1;

#if 0
   // argh, we've got the logic here the wrong way around.  but it isn't so critical.
   // make the minimum and maximums of i,j,k,l internally consistent
   if (lmin_ > Size_) lmin_ = Size_;
   if (kmin_ > lmin_ - klmin_) kmin_ = lmin_ - klmin_;
   if (jmin_ > kmin_ - jkmin_) jmin_ = kmin_ - jkmin_;
   if (jmin_ > lmin_ - jlmin_) jmin_ = lmin_ - jlmin_;
   if (imin_ > jmin_ - ijmin_) imin_ = jmin_ - ijmin_;
   if (imin_ > kmin_ - ikmin_) imin_ = jmin_ - ikmin_;
   if (imin_ > lmin_ - ilmin_) imin_ = jmin_ - ilmin_;

   if (imax_ < 1) imax_ = 1;
   if (jmax_ < imax_ + ijmin_) jmax_ = imax_ + ijmin_;
   if (kmax_ < imax_ + ikmin_) kmax_ = imax_ + ikmin_;
   if (kmax_ < jmax_ + jkmin_) kmax_ = jmax_ + jkmin_;
   if (lmax_ < imax_ + ilmin_) lmax_ = imax_ + ilmin_;
   if (lmax_ < jmax_ + jlmin_) lmax_ = jmax_ + jlmin_;
   if (lmax_ < kmax_ + klmin_) lmax_ = kmax_ + klmin_;
#endif

   // check the pairwise maxima
   if (ijmax_ > jmax_ - imin_) ijmax_ = jmax_ - imin_;
   if (ikmax_ > kmax_ - imin_) ikmax_ = kmax_ - imin_;
   if (ilmax_ > lmax_ - imin_) ilmax_ = lmax_ - imin_;
   if (jkmax_ > kmax_ - jmin_) jkmax_ = kmax_ - jmin_;
   if (jlmax_ > lmax_ - jmin_) jlmax_ = lmax_ - jmin_;
   if (klmax_ > lmax_ - kmin_) klmax_ = lmax_ - kmin_;
}

int PositionMap::FirstSite() const
{
   return imin_;
}

int PositionMap::LastSite() const
{
   return lmax_;
}




int main(int argc, char** argv)
{
   try
   {
      bool Verbose = false;
      std::string LatticeStr, PsiStr, Op1, Op2, Op3, Op4;

      prog_opt::options_description desc("Allowed options", terminal::columns());
      desc.add_options()
         ("help", "show this help message")
	 ("imin", prog_opt::value<int>(), "minimum value for i [default 1]")
	 ("imax", prog_opt::value<int>(), "maximum value for i [default size-3]")
	 ("jmin", prog_opt::value<int>(), "minimum value for j [default 2]")
	 ("jmax", prog_opt::value<int>(), "maximum value for j [default size-2]")
	 ("kmin", prog_opt::value<int>(), "minimum value for k [default 3]")
	 ("kmax", prog_opt::value<int>(), "maximum value for k [default size-1]")
	 ("lmin", prog_opt::value<int>(), "minimum value for l [default 4]")
	 ("lmax", prog_opt::value<int>(), "maximum value for l [default size]")

	 ("ijmin", prog_opt::value<int>(), "minimum separation between i and j [default 1]")
	 ("ijmax", prog_opt::value<int>(), "maximum separation between i and j [default infinity]")
	 ("jkmin", prog_opt::value<int>(), "minimum separation between j and k [default 1]")
	 ("jkmax", prog_opt::value<int>(), "maximum separation between j and k [default infinity]")
	 ("klmin", prog_opt::value<int>(), "minimum separation between k and l [default 1]")
	 ("klmax", prog_opt::value<int>(), "maximum separation between k and l [default infinity]")

	 ("ikmin", prog_opt::value<int>(), "minimum separation between i and k [default 2]")
	 ("ikmax", prog_opt::value<int>(), "maximum separation between i and k [default infinity]")
	 ("jlmin", prog_opt::value<int>(), "minimum separation between j and l [default 2]")
	 ("jlmax", prog_opt::value<int>(), "maximum separation between j and l [default infinity]")

	 ("ilmin", prog_opt::value<int>(), "minimum separation between i and l [default 3]")
	 ("ilmax", prog_opt::value<int>(), "maximum separation between i and l [default infinity]")

	 ("middle-quantum", prog_opt::value<std::string>(), 
	  "project op1*op2 onto this quantum number (only needed in non-Abelian case)")

         ("verbose,v", prog_opt::bool_switch(&Verbose),
          "extra debug output (currently does nothing)")
	 ;

      prog_opt::options_description hidden("Hidden options");
      hidden.add_options()
         ("lattice", prog_opt::value<std::string>(&LatticeStr), "latticex")
         ("psi", prog_opt::value<std::string>(&PsiStr), "psix")
         ("op1", prog_opt::value<std::string>(&Op1), "op1x")
         ("op2", prog_opt::value<std::string>(&Op2), "op2x")
         ("op3", prog_opt::value<std::string>(&Op3), "op3x")
         ("op4", prog_opt::value<std::string>(&Op4), "op4x")
         ;

      prog_opt::positional_options_description p;
      p.add("lattice", 1);
      p.add("psi", 1);
      p.add("op1", 1);
      p.add("op2", 1);
      p.add("op3", 1);
      p.add("op4", 1);

      prog_opt::options_description opt;
      opt.add(desc).add(hidden);

      prog_opt::variables_map vm;        
      prog_opt::store(prog_opt::command_line_parser(argc, argv).
                      options(opt).positional(p).run(), vm);
      prog_opt::notify(vm);    

      if (vm.count("help") > 0 || vm.count("op4") == 0)
      {
         print_copyright(std::cerr, "tools", basename(argv[0]));
         std::cerr << "usage: mp-local-fourpoint-ex [options] <lattice> <psi> <op1> <op2> <op3> <op4>\n";
	 std::cerr << "evaluates the 4-point correlation function < op1(i) op2(j) op3(k) op4(l) >\n";
         std::cerr << desc << '\n';
         return 1;
      }

      long CacheSize = getenv_or_default("MP_CACHESIZE", 655360);
      pvalue_ptr<OperatorList> System = pheap::OpenPersistent(LatticeStr, CacheSize, true);
      pvalue_ptr<MPWavefunction> Psi1 = pheap::ImportHeap(PsiStr);

      Lattice Lat = System->GetLattice();

      CenterWavefunction Psi = *Psi1;
      Psi1 = pvalue_ptr<MPWavefunction>();

      std::cout.precision(12);


      PositionMap Position(Lat.size());

      if (vm.count("imin") > 0) Position.set_imin(vm["imin"].as<int>());
      if (vm.count("imax") > 0) Position.set_imax(vm["imax"].as<int>());
      if (vm.count("jmin") > 0) Position.set_jmin(vm["jmin"].as<int>());
      if (vm.count("jmax") > 0) Position.set_jmax(vm["jmax"].as<int>());
      if (vm.count("kmin") > 0) Position.set_kmin(vm["kmin"].as<int>());
      if (vm.count("kmax") > 0) Position.set_kmax(vm["kmax"].as<int>());
      if (vm.count("lmin") > 0) Position.set_lmin(vm["lmin"].as<int>());
      if (vm.count("lmax") > 0) Position.set_lmax(vm["lmax"].as<int>());

      if (vm.count("ijmin") > 0) Position.set_ijmin(vm["ijmin"].as<int>());
      if (vm.count("ijmax") > 0) Position.set_ijmax(vm["ijmax"].as<int>());
      if (vm.count("ikmin") > 0) Position.set_ikmin(vm["ikmin"].as<int>());
      if (vm.count("ikmax") > 0) Position.set_ikmax(vm["ikmax"].as<int>());
      if (vm.count("ilmin") > 0) Position.set_ilmin(vm["ilmin"].as<int>());
      if (vm.count("ilmax") > 0) Position.set_ilmax(vm["ilmax"].as<int>());
      if (vm.count("jkmin") > 0) Position.set_jkmin(vm["jkmin"].as<int>());
      if (vm.count("jkmax") > 0) Position.set_jkmax(vm["jkmax"].as<int>());
      if (vm.count("jlmin") > 0) Position.set_jlmin(vm["jlmin"].as<int>());
      if (vm.count("jlmax") > 0) Position.set_jlmax(vm["jlmax"].as<int>());
      if (vm.count("klmin") > 0) Position.set_klmin(vm["klmin"].as<int>());
      if (vm.count("klmax") > 0) Position.set_klmax(vm["klmax"].as<int>());

      Position.fix();

      boost::optional<QuantumNumbers::QuantumNumber> MiddleProjection = boost::none;
      if (vm.count("middle-quantum") > 0)
	 MiddleProjection = QuantumNumbers::QuantumNumber(Psi.GetSymmetryList(), 
							  vm["middle-quantum"].as<std::string>());

      // This holds the E matrices for the left system
      typedef std::map<int, MatrixOperator> OpMapType;


      // Assemble the right blocks.  Op4Position is a list of positions of operator 4.
      // Op3Position is a list of pairs (Op3Pos, Op4Pos) of positions of operators 3 and 4.
      typedef std::list<MatrixOperator> OpListType;
      typedef std::list<int> Op4PosListType;
      typedef std::list<std::pair<int, int> > Op3PosListType;
      typedef std::list<pvalue_handle<MatrixOperator> > OpStackType;
      typedef std::list<OpStackType> RightBlockType;

      RightBlockType Op3Block;
      OpListType Op4Block;
      Op4PosListType Op4Position;
      Op3PosListType Op3Position;

      DEBUG_TRACE("Constructing right block stack");
      for (int i = Position.LastSite(); i >= Position.FirstSite(); --i)
      {
	 DEBUG_TRACE("Top of right block consstruction loop")(i);
	 MPStateComponent PsiR = Psi.LookupLinear(i-1);
	 // Update Op3Block as necessary
	 RightBlockType::iterator I = Op3Block.begin();
	 while (I != Op3Block.end())
	 {
	    I->push_back(new MatrixOperator(operator_prod(PsiR, *I->back().load(), herm(PsiR))));
	    ++I;
	 }

	 // Update Op4Block
	 Op4PosListType::iterator Ip4 = Op4Position.begin();
	 OpListType::iterator I4 = Op4Block.begin();
	 SiteBlock::const_iterator I3 = Lat[i].find(Op3);
	 bool Op3ExistsHere = (I3 != Lat[i].end());
	 while (I4 != Op4Block.end())
	 {
	    if (Op3ExistsHere && Position.Need34(i, *Ip4))
	    {
	       DEBUG_TRACE("Adding op3*op4")(i)(*Ip4);
	       Op3Position.push_front(std::make_pair(i, *Ip4));
	       if (MiddleProjection)
		  Op3Block.push_front(OpStackType(1, new MatrixOperator
						  (operator_prod(SimpleOperator(I3->second),
								 PsiR, *I4, herm(PsiR), 
								 adjoint(*MiddleProjection)))));
	       else
		  Op3Block.push_front(OpStackType(1, new MatrixOperator
						  (operator_prod(SimpleOperator(I3->second),
								 PsiR, *I4, herm(PsiR)))));
	    }

	    if (Position.Need4Later(i, *Ip4))
	    {
	       *I4 = operator_prod(PsiR, *I4, herm(PsiR));
	       ++I4;
	       ++Ip4;
	    }
	    else
	    {
	       // don't need Op4 anymore
	       DEBUG_TRACE("Discarding op4")(*Ip4);
	       I4 = Op4Block.erase(I4);
	       Ip4 = Op4Position.erase(Ip4);
	    }
	 }

	 // add current site to Op4Block, as necessary
	 SiteBlock::const_iterator Si4 = Lat[i].find(Op4);
	 bool Op4ExistsHere = (Si4 != Lat[i].end());
	 if (Op4ExistsHere && Position.Need4Later(i, i))
	 {
	    DEBUG_TRACE("Adding op4")(i);
	    MatrixOperator Ident = MatrixOperator::make_identity(PsiR.Basis2());
	    Op4Position.push_back(i);
	    Op4Block.push_back(operator_prod(SimpleOperator(Si4->second),
					     PsiR, Ident, herm(PsiR)));
	 }
      }

      // left blocks that represent Op1,Op2 at various locations.
      // Op1Block is just Op1, at sites labelled by the corresponding
      // Op1Position.  Op2Block is the product of Op1 * Op2, at the
      // positions given by Op2Position.
      // we construct Op3 and Op4 as right blocks later.
      typedef std::list<int> Op1PosListType;
      typedef std::list<std::pair<int, int> > Op2PosListType;
      std::list<MatrixOperator> Op1Block, Op2Block;
      std::list<int> Op1Position;
      std::list<std::pair<int, int> > Op2Position;
      
      DEBUG_TRACE("Rotating right");

      // rotate Psi to the first needed position for the left blocks
      while (Psi.LeftSize() < Position.FirstSite())
	 Psi.RotateRight();

      for (int i = Position.FirstSite(); i < Position.LastSite(); ++i)
      {
	 DEBUG_TRACE("Top of left block construction loop")(i);
	 // update the Op2Block as necessary
	 OpListType::iterator I = Op2Block.begin();
	 Op2PosListType::iterator Ip2 = Op2Position.begin();
	 while (I != Op2Block.end())
	 {
	    // do we need to keep this entry?
	    if (Position.Need12Later(Ip2->first, Ip2->second, i))
	    {
	       *I = operator_prod(herm(Psi.Left()), *I, Psi.Left());
	       ++I;
	       ++Ip2;
	    }
	    else
	    {
	       DEBUG_TRACE("Discarding left block operator")(Ip2->first)(Ip2->second);
	       I = Op2Block.erase(I);
	       Ip2 = Op2Position.erase(Ip2);
	    }
	 }
	 CHECK(Ip2 == Op2Position.end());

	 // Update the Op1 list.
	 SiteBlock::const_iterator I2 = Lat[i].find(Op2);
	 bool Op2ExistsHere = (I2 != Lat[i].end());
	 I = Op1Block.begin();
	 Op1PosListType::iterator Ip1 = Op1Position.begin();
	 while (I != Op1Block.end())
	 {
	    // do we need the combination (*I) * Op2(i) ?
	    if (Op2ExistsHere && Position.Need12(*Ip1, i))
	    {
	       DEBUG_TRACE("Inserting op1*op2 combination")(*Ip1)(i);
	       Op2Position.push_back(std::make_pair(*Ip1, i));
	       SimpleOperator LocalOp2 = I2->second;
	       if (MiddleProjection)
		  Op2Block.push_back(operator_prod(herm(LocalOp2), herm(Psi.Left()), *I,
						   Psi.Left(), adjoint(*MiddleProjection)));
	       else
		  Op2Block.push_back(operator_prod(herm(LocalOp2), herm(Psi.Left()), *I,
						   Psi.Left()));
	    }
	 
	    // either update this operator with the local identity, or discard it
	    if (Position.Need1Later(*Ip1, i))
	    {
	       *I = operator_prod(herm(Psi.Left()), *I, Psi.Left());
	       ++I;
	       ++Ip1;
	    }
	    else
	    {
	       DEBUG_TRACE("Discarding op1")(*Ip1);
	       I = Op1Block.erase(I);
	       Ip1 = Op1Position.erase(Ip1);
	    }
	 }
	 CHECK(Ip1 == Op1Position.end());

	 // Add Op1 at the current site to Op1Block as necessary
	 SiteBlock::const_iterator I1 = Lat[i].find(Op1);
	 bool Op1ExistsHere = (I1 != Lat[i].end());
	 if (Position.Need1Later(i,i) && Op1ExistsHere)
	 {
	    MatrixOperator Ident = MatrixOperator::make_identity(Psi.Left().Basis1());
	    Op1Position.push_back(i);
	    Op1Block.push_back(operator_prod(herm(SimpleOperator(I1->second)), 
					     herm(Psi.Left()), Ident, Psi.Left()));
	 }

	 // Go through the right stack and update (ie, pop off the stack the un-needed operator)
	 RightBlockType::iterator R = Op3Block.begin();
	 Op3PosListType::iterator Rp = Op3Position.begin();
	 while (R != Op3Block.end())
	 {
	    R->pop_back();
	    if (R->empty())
	    {
	       R = Op3Block.erase(R);
	       Rp = Op3Position.erase(Rp);
	    }
	    else
	    {
	       ++R;
	       ++Rp;
	    }
	 }
	 
	 // evaluate the expectation values with the right blocks.  We do this at the
	 // point where Op2 is evaluated at the left site.  this is a somewhat
	 // arbitrary choice.
	 I = Op2Block.begin();
	 Ip2 = Op2Position.begin();
	 while (I != Op2Block.end())
	 {
	    if (Ip2->second != i)  // evaluate only at Op2(i)
	    {
	       ++I;
	       ++Ip2;
	       continue;
	    }
	    
	    R = Op3Block.begin();
	    Rp = Op3Position.begin();
	    while (R != Op3Block.end())
	    {
	       if (Position.ShouldEvaluate(Ip2->first, Ip2->second, Rp->first, Rp->second))
	       {
		  std::complex<double> Res = inner_prod(Psi.Center(),
							triple_prod(*I, Psi.Center(),
								    herm(*R->back().load())));
		  std::cout << std::setw(5) << Lat.coordinate_at_site(Ip2->first) << "   " 
			    << std::setw(5) << Lat.coordinate_at_site(Ip2->second) << "   "
			    << std::setw(5) << Lat.coordinate_at_site(Rp->first) << "   "
			    << std::setw(5) << Lat.coordinate_at_site(Rp->second) << "   "
			    << std::setw(18) << Res.real() << "   " 
			    << std::setw(18) << Res.imag() << '\n';
	       }
	       ++R;
	       ++Rp;
	    }
	    ++I;
	    ++Ip2;
	 }
         
	 Psi.RotateRight();
      }
      pheap::Shutdown();

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
