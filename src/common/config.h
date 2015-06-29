#ifndef __CONFIG_H_
#define __CONFIG_H_

/*
** One (and only one) of these *MUST* be defined.  It controls how the
** program handles the numbers and allocates the memory.
#define USE_VIRT_MEM  1
#define USE_DISK_NUM  1
*/
#define USE_VIRT_MEM  1

/*
** Whether to use virtual memory for the FFT/NTT caching etc.
** If you are using disk based numbers, it's certainly counter
** productive.  If you are pushing the limits of your virtual
** memory system, you may get better results without using
** virtual memory for the FFT/NTT caches.
*/
#define VIRTUAL_CACHE 1

/*
** DJGPP has a very, very, very poor malloc.  It won't let you allocate
** more than half of the remaining free memory.
**
** If you have 30 megs of memory free, you'd expect to be able to allocate
** at least 16 megs of it at once, but with DJGPP, you can't!  It takes the
** size requested, adds a few bytes of overhead, and then rounds that UP
** to the next power of two!  So if you had 30 megs, and you tried to just
** 8 megs, it'd add a few bytes, round that up to 16 megs, allocate 16 of
** the 30 megs free, and then give you the 8 megs you requested, and wasting
** the other 8 megs!  If you tried to malloc(16meg), DJGPP would add a few
** bytes, round that up to 32 megs, and try to allocate a full 32 megs
** of memory.  If you are using virtual memory, it'd probably succeed, but
** it'd still wast 16megs of virtual memory space.
**
** By defining the below symbol, you are telling my program that it should
** use the replacement malloc() package (based on malloc6() by D.J. Delorie)
** instead of the one normally used in DJGPP 2.7.2.1
**
** **NOTE** DJGPP 2.8.1 uses a new malloc routine (in fact, it's probably
** related to the one that I use.)  So, you might not need this.  It
** wont hurt anything to use it anyway, though.
*/
#define USE_SPECIAL_DJGPP_MALLOC 1
#ifndef __DJGPP__
/* You probably don't need it, but you can override this if you want. */
/* You might need to modify malloc6.c, for a different compiler. */
#undef USE_SPECIAL_DJGPP_MALLOC
#endif

/*
** I don't entirely trust the replacement malloc package.  I don't have
** any reason to distrust it, it's just that I don't have enough
** experience with it to completely trust it.  So, I'm allowing you
** to not use it, by not defining the symbol above.  However, DJGPP
** is still stuck with a poor malloc and something still needs to be
** done.  So, I'm allowing you to let the disk numbers allocate their
** memory via sbrk().  That's an old Unix function that most C compilers
** support.  You can't free sbrk() memory, so this only works with the
** disk numbers.
**
** Anyway, you now have a total of 3 choices.  You can use the regular
** malloc package (with or without alignment).  You can use the new
** malloc6 package (with or without alignment).  Or you can use sbrk()
** when doing disk number runs.
#define USE_SBRK
*/

#if !defined(USE_DISK_NUM)
#undef USE_SBRK
#endif


/*
** Whether to do an aligned malloc.  Some malloc() routines only allocate
** on 4 byte boundaries, even though many systems, such as the x86, get
** much better fpu performance when things are aligned on an 8 byte (or
** more) boundary.  This option will have the routines in malloc6.c
** allocate a few bytes more, and then align the result.  Of course,
** the program will always deallocate that memory with the modified free().
**
** This may NOT work with all compilers.  I tried to make it generic,
** but no guarantees.
*/
#define ALIGN_MALLOC 1

/*
** Whether to gather some timing stats during the run.
**
** The accuracy of the timings, and how much/little they impact
** the performance of the program depends upon you.
**
** Under DOS and DJGPP, I can directly read the PC's 18.2 timer
** with little performance penalty.  Obviously that wont work for
** other systems, so for them I have it set up to just call the
** time() function.  This will almost certaily cause a *SIGNIFICANT*
** performance hit, but you can modify pi.h for something a bit
** less painful if your system supplies it.
**
** The pentium version uses the built in 64 bit cycle counter.
**
#define DO_TIMINGS 1
#define DO_PENTIUM_TIMINGS 1
*/


#if defined(DO_PENTIUM_TIMINGS) && defined(__GNUC__)
/* Since it uses asm, it currently only works with DJGPP. */
typedef unsigned long long TIMER_TYPE;
#define DO_TIMINGS 1
#else
#undef DO_PENTIUM_TIMINGS
typedef time_t TIMER_TYPE;
#endif

/*
** If your computer/OS doesn't have a command line, define the line
** below.  Even though there is an option in the pi.ini, you can
** still define this to override that.
#define NO_COMMAND_LINE 1
*/


/*
** Some safety, sensibility self-checks.
*/

#if !defined(USE_DISK_NUM) && !defined(USE_VIRT_MEM)
#define USE_VIRT_MEM  1
#endif

#ifndef USE_VIRT_MEM
#undef VIRTUAL_CACHE
#endif


#endif


