// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://physics.uq.edu.au/people/ianmcc/mptoolkit/
//
// common/niftycounter.h
//
// Copyright (C) 2016 Ian McCulloch <ianmcc@physics.uq.edu.au>
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
/* -*- C++ -*-
  niftycounter.h

  An implementation of the nifty counter used to ensure correct ordering of static data initialization.
  To use: write some functions MyInitFunc() and MyExitFunc().  Write a .h file that opens an unnamed 
  namespace (unusual!) and declares a namespace scope nifty_counter<MyInitFunc, MyExitFunc> with extern linkage.  
  This ensures that MyInitFunc() is called before any other static data members of any translation units that 
  #include your header.

  2003-10-17 Ian McCulloch: Allow functors now, with a default of DoNothing for the ExitFunc.
                            You can supply the functors to the constructor, if you want.
*/

#if !defined(MPTOOLKIT_COMMON_NIFTYCOUNTER_H)
#define MPTOOLKIT_COMMON_NIFTYCOUNTER_H

namespace NiftyCounter
{

typedef void(*InitFuncType)();

void DoNothing();

template <InitFuncType Init, InitFuncType Exit = DoNothing>
class nifty_counter
{
   public:
      nifty_counter()
      {
	 if (count++ == 0) Init();
      }

      ~nifty_counter()
      {
	 if (--count == 0) Exit();
      }

   private:
      static int count;
};

template <InitFuncType InitFunc, InitFuncType ExitFunc>
int nifty_counter<InitFunc, ExitFunc>::count = 0;

} // namespace NiftyCounter

#endif
