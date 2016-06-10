// -*- C++ -*-
// ENDHEADER

#include "color-output.h"
#include <unistd.h>

bool is_cout_terminal()
{
   return isatty(1);
}

