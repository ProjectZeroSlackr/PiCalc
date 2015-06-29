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
** If your replacement FFT needs its own, then of course you shouldn't
** use the default one!
*/
#define USE_DEFAULT_CONVOLUTION 1

/*
** **NOTE** Although the basic 'bigmul.c' can handle 'long double'
** fft's, that does not mean the FFT itself can handle it.  This
** option is here in case your own FFT can use it.
*/

/*
** Whether to try the 'long double' version or the regular 'double'
** version.  The regular 'double' can multiply up to two 8 million digit
** numbers, getting a 16 million digit answer.  It consumes 32 mega-
** bytes.  The 'long double' version will be able to multiply at least
** two 32 million digit numbers, getting a 64 million digit answer,
** but I don't know for sure because I don't have access to a system
** with enough memory to test.  The memory consumption will be at
** least 25% more, probably 50% more, and as much as 100% more over
** the 'double' version.  Plan on it being 100% more.
**
** The 'long double' ** MAY NOT BE AVAILABLE ** on all systems.
** Many CPUs don't have 'long double' in hardware.  They have to use
** software.  And many compilers don't support it even on the platforms
** that do have it.
**
** Since this is so untested, I can't guarantee anything about it.
** I'm not even sure which 'trig' macro would be best.  Probably
** the explicit trig one, but the sin_recur might work, as might the
** 'wide_fpu' one (even though the FPU registers won't acutally be
** larger than the data, there will be some extra 'slack' that might
** allow it to work.)
**
#define LONG_DOUBLE_FFT 1
*/

#ifdef LONG_DOUBLE_FFT
#if LDBL_MANT_DIG == DBL_MANT_DIG
#error It looks like this compiler does not have 'long double'.  It looks the same as 'double'.
#endif
typedef long double            FFT_DATA_TYPE;
#else
typedef      double            FFT_DATA_TYPE;
#endif
typedef unsigned           int UINT32;
typedef   signed short     int  INT16;
typedef   signed           int  INT32;

#endif


