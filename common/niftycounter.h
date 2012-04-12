/* -*- C++ -*- $Id$
  niftycounter.h

  An implementation of the nifty counter used to ensure correct ordering of static data initialization.
  To use: write some functions MyInitFunc() and MyExitFunc().  Write a .h file that opens an unnamed 
  namespace (unusual!) and declares a namespace scope nifty_counter<MyInitFunc, MyExitFunc> with extern linkage.  
  This ensures that MyInitFunc() is called before any other static data members of any translation units that 
  #include your header.

  2003-10-17 Ian McCulloch: Allow functors now, with a default of DoNothing for the ExitFunc.
                            You can supply the functors to the constructor, if you want.
*/

#if !defined(NIFTYCOUNTER_H_HDFJKFHD8357YGUIR3746ERUYTY7HERUI)
#define NIFTYCOUNTER_H_HDFJKFHD8357YGUIR3746ERUYTY7HERUI

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
