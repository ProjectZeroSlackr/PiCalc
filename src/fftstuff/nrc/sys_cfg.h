/*
** This file contains various system specific items.
**
** This include various data types, macros, etc. that you might need.
*/

#ifndef SYS_CFG_H_
#define SYS_CFG_H_

#include <float.h>

/*
** C has a very slow floor function.  With a bit of clever coding, we
** can take advantage of the x86's FPU registers and do it much faster.
** This does require a few lines of assembly code to put the FPU into
** the proper mode.  Which is why it's GNU C and x86 specific.
*/
#define USE_FAST_X86_FLOOR 1
#if !defined(__GNUC__) || !defined(i386)
//#error The USE_FAST_X86_FLOOR was set in sys_cfg.h
#undef USE_FAST_X86_FLOOR
#endif


/*
** How many digits should we put into the FFT.  The normal is to put
** 4 digits.  That provides a total multiplication size of 8 million
** digits, getting a 16 million digit result.  This takes 32m bytes.
** (Of course, allowing for how you do your trig, this might be as
** low as just 1m digits.)
**
** Or, we can just put two digits.  Right off the bat, this doubles
** the memory used, but also increases the total multiplication
** size to somewhere around one billion (1024m) digits.  Of course,
** to do this means using _lots_ more memory.  Multiplying 8m digits
** would take 64m of memory.  Going to 16m digits would take 128m of
** memory.  Going all the way to 1g digits would take a massive
** 8g of memory (which is beyond what 32 bit processors can even
** directly access.)
**
** In short, it's not really all that useful for desktop systems.  Just
** leave it at 4.  The only reason this ability is added are for those
** workstations that have gobs of memory and can do floating point super
** quick, and possibly even have custom designed FFT libraries that are
** much more tuned for the hardware than mine is.  It may or may not
** be useful, but I figured I'd add it for just in case.
#define FFT_DIGITS 2
#define FFT_DIGITS 4
*/
#define FFT_DIGITS 4

/*
** Whether to use the default convolution code in bigmul.c
** That code is for standard style FFT output.  This is sometimes
** called the 'Numerical Recipes' style, although it predates NR.
**
** If your replace FFT needs its own, then of course you shouldn't
** use the default one!
*/
#define USE_DEFAULT_CONVOLUTION 1

typedef      double            FFT_DATA_TYPE;
typedef unsigned           int UINT32;
typedef   signed short     int  INT16;
typedef   signed           int  INT32;

#endif


