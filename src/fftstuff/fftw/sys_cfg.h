/*
** This file contains various system specific items.
**
** This include various data types, macros, etc. that you might need.
*/

#ifndef SYS_CFG_H_
#define SYS_CFG_H_

#include <fftw.h>

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
*/
#define USE_DEFAULT_CONVOLUTION 1

/*
** There is quite a bit of variation as to how floating point math is
** handled.
**
** Some systems, such as the x87 and 68882 have 64 bit internal FPU
** registers, where math can be operated on in precision greater than
** the 'double' we are actually using.
**
** Other systems only have 'double' hardware, and all FPU math is
** done in that precision.  Also some compilers or OSs will force
** this behavior even if your system does have wider FPU registers.
**
** Some systems do have a 'long double' data type, but it's done
** using some sort of software emulation.  In other words, it works
** if you need it, but you can't write fast code with it.
**
** Although the FFTW package does its own trig, and this macro doesn't
** effect that, it does effect the real<->complex wrapper function.
**
**
** The FFT_WIDE_FPU_RECUR is for the x87 and other systems that have
** FPU registers longer than a standard 'double'.  It does an initial
** sin() call and then trusts the compiler to keep the values in the
** wider FPU registers.  This allows it to do a simple, quick trig
** recurance for the next trig value.  When this works, such as with
** DJGPP on the x86/x87, you can multiply up to 8 million digits
** together.  When it doesn't work (such as on a pure 'double' system),
** you can only multiply 1 million digit numbers together.  Beyond that,
** the errors accumulate too much to be reliable.  With some compilers,
** it might even be a bit less.
**
** DJGPP works fine with the FFT_WIDE_FPU_RECUR.
** Watcom doesn't work with it.  It requires one from below.
**
** The FFT_TRIG_SIN is for those systems that don't fit in the above
** category.  This does explicit sin() calls every time it needs a
** new trig value.  This does take extra time.  On my 486/66 system,
** this version takes about 10% longer than the recurance.  However,
** this version is capable of working correctly on all standard
** 'double' systems.  It doesn't depend on the FPU registers or
** compilers keeping a value in greater precision.  You can multiply
** two 8 million digit numbers together with this.  (However, this
** does depend on the accuracy of the sin() function.  But if your
** system is that inaccurate, it's probably not a very good system
** and you wouldn't be able to do much numerical work anyway.)
**
** The FFT_TRIG_SIN_RECUR is a compromise between the two.  It runs
** a wee bit faster than the _SIN version because it only does the
** sin() calls every 8th time.  (Currently.  That could be changed
** by adjusting the NEXT_TRIG_POW macro to a value such as 15.)
** The error is a bit higher, but it should work.
**
** (There are other ways I could do this besides using explicit trig.
** However, this was easy to write and I don't know how many people
** would even need a solution like this, so I don't know if it would
** be worth the effort to write and test it.)
**
** There is no way to detect this at compile time.  No, the values in
** float.h aren't useful.  That only says that a 'long double' is
** available, not whether its FPU registers are that long, or whether
** the compiler will use them, or even whether the 'long double' will
** be done in hardware or software emulation.  It can be done at
** runtime, but it would involve too much work, and result in quite
** a bit of code bloat.  It's easier for you to just set this for your
** system.
**
** If you aren't sure, chose the FFT_TRIG_SIN version.  It should
** always work.  If you want to see if the faster FFT_TRIG_RECUR works,
** select it.  Then, when you run the program, give it a negative size
** (ie: pi -4m) and the program will run a number of tests to see if
** things work.
**
** If all else fails, you can always keep halving the MaxFFTLen
** in fft.h, until things work.  That will reduce the maximum length
** that the FFT will be called with.  Anything beyond that will break
** the data into chunks small enough for the FFT to work with.
** Breaking the data into smaller chunks does cost some time, but
** if all else fails, it will work.
*/
/*
#define FFT_WIDE_FPU_RECUR  1
#define FFT_TRIG_SIN        1
#define FFT_TRIG_SIN_RECUR  1
*/
#define FFT_WIDE_FPU_RECUR  1


/*
** If the FFTW package is working with 'long double' data type,
** (ie: if you built a special fftw package with fftw_real declared
** as 'long double') then you'll need to enable this.  No way to do
** it automatically.
#define LONG_DOUBLE_FFT 1
*/

#ifdef LONG_DOUBLE_FFT
typedef long double LDouble;
#define SINE(A)  sinl(A)
#define FLOOR(F) floorl(F)
/* Or whatever your system happens to need. */
#else
typedef fftw_real   LDouble;
#define SINE(A)  sin(A)
#define FLOOR(F) floor(F)
#endif

typedef              fftw_real FFT_DATA_TYPE;
typedef unsigned           int UINT32;
typedef   signed short     int  INT16;
typedef   signed           int  INT32;

#endif


