/* -*- C++ -*- $Id$

  Attempts to determine the native data format.
  This is done by preprocessor checks only, no attempt is made at automatic
  detection.

  If we ever wanted to cross-compile, the checks would need to be modified as we are
  assuming here that the host platform is the same as the target platform.

  Currently recognized data formats:

  Linux_i386  : little endian LP32, 32 bit int/long, 64 bit long long, IEEE 32/64/80 floating point
  Tru64_Alpha : little endian LP64, 32 bit int, 64 bit long, IEEE 32/64/128 floating point

  Instructions for adding a new platform:
  1. Devise an identifier that describes the platform
  2. Edit this header file to add preprocessor checks to detect the current platform
     and #define NATIVE_DATA_FORMAT to be the identifier devised in step (1), together
     with NATIVE_DATA_FORMAT_STR (string version of NATIVE_DATA_FORMAT), and 
     a symbol IS_DATA_FORMAT_<name> defined to be 1.
  3. Construct platformformat_<name>.h and currentformat_<name>.h
     See platformformats.h and currentformat.h for what these files should contain.
  4. Edit platformformats.h and currentformat.h to include the files created by step (3).
*/

#if !defined(MCONFIG_H_RY437RWYUT786RT782TYU2366T3FEWG)
#define MCONFIG_H_RY437RWYUT786RT782TYU2366T3FEWG

#if defined(__linux__) && defined(__i386__)

#define CURRENT_DATA_FORMAT Linux_i386
#define CURRENT_DATA_FORMAT_STR "Linux_i386"
#define DATA_FORMAT_IS_Linux_i386 1

#elif defined(__unix__) && defined(__alpha) && defined(__arch64__)

#define CURRENT_DATA_FORMAT Tru64_Alpha
#define CURRENT_DATA_FORMAT_STR "Tru64_Alpha"
#define DATA_FORMAT_IS_Tru64_Alpha 1

#endif

#if !defined(CURRENT_DATA_FORMAT)
#error "Unknown data format!  Edit platformconfig.h"
#endif

#endif
