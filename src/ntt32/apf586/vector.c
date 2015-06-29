#include "pi.h"
#include "ntt.h"
#include "modmath.h"

/*
** 586 specific assembly language 'Vector' routines.
*/

/*
   this file contains the assembly routines for a
   high-speed number-theoretic FFT. It can work
   modulo any number < 2^31-1; approximate timings
   on a 200MHz Pentium MMX:

      256 pts   168us
      512 pts   322us
     1024 pts   710us
     2048 pts  1560us

   I hereby place this code into the public domain. Use
   it for whatever you want, no strings attached, but
   optionally please be nice and tell me if you find it useful.

                                        Jason Papadopoulos
                                        jasonp@glue.umd.edu
                                        10/24/98

The command line I used was

gcc -o test.exe ntt586.c -O2 -m486 -fomit-frame-pointer

**********************************************************************
**********************************************************************
You can optionally go to -O3 with -fno-inline-functions; in any case,
I need ebp all the time so -fomit-frame-pointer is vital.
**********************************************************************
**********************************************************************

*/

/*
** Obviously, this file has been heavily modified, but the asm routines
** are still his, or based on his.  The PD status is still applicable.
**
** CEB
*/


/*
** The asm code needs to access these two variables.  But, various
** GCC compilers handle the exact asm name differently.  Some do
** an underscore, others don't.  So we need to explicitly name them.
**
** We could use the 'Prime' variable, instead of JPrime, but we
** don't know how that one will be named.
**
** GCC v2.8.1 has absolutely terrible data alignment problems.
** What's worse, the method GCC tells you (in GCC.I9 doc file)
** to use to fix alignment problems, doesn't work!
**
** This means you are going to have to fix this yourself!  *AND*
** this will change depending upon whether you are using disk numbers
** or virtual memory, or the Pentium timers.  Just move the alignment
** variable out of the comments.  You may need to adjust the number
** of elements in the array, too.

** The BigInt is there since that data type changes size.  This allows
** a resonably consistant alignment.  The int array is there in case
** you need to align more.
**
** **HOWEVER** there is no guarantee that the assembler or linker
** will actually put the alignment variable where it will do any
** good!!  There is absolutely no way around this.  Sorry.  I really
** wish there was a better generic solution.
*/
BigInt VectorAlign1;
int VectorAlign2[1];
ModInt JPrime asm("_JPrime");
double JSP_stuff[8] asm("_JSP_stuff");
/* contents: 3*2^62  3*2^51  n  temp1 temp2 temp3  1/n  */

void
VectorModAdd(ModInt *Num1, ModInt *Num2, size_t Len)
{
while (Len)
  {
   Num1[0]=ModAdd(Num1[0],Num2[0]);
   Num1[1]=ModAdd(Num1[1],Num2[1]);
   Num1[2]=ModAdd(Num1[2],Num2[2]);
   Num1[3]=ModAdd(Num1[3],Num2[3]);
   Num1+=4;Num2+=4;Len-=4;
  }
}

void
VectorModSub(ModInt *Num1, ModInt *Num2, size_t Len)
{
while (Len)
  {
   Num1[0]=ModSub(Num1[0],Num2[0]);
   Num1[1]=ModSub(Num1[1],Num2[1]);
   Num1[2]=ModSub(Num1[2],Num2[2]);
   Num1[3]=ModSub(Num1[3],Num2[3]);
   Num1+=4;Num2+=4;Len-=4;
  }
}

void
VectorModMul(ModInt *Num1, ModInt *Num2, size_t Len)
{
   /*  Num1[0...len] <- Num1[0...len] * Num2[0...len] mod n

      len must be a multiple of 4. Running time: 64 cycles
      every 4 multiplies. The loop itself is coded here too, to
      avoid the overhead of pushing and popping registers all
      the time (saves about 10 cycles). The final answers are
      always positive; even with sign correction, using the FPU
      is almost four times faster than using mul/div, though
      *much* more complicated.

      A general rule with gcc 2.7.2.1 is that if you ask for
      all the ALU registers, invariably something will go very
      wrong. This routine was no execption, but the register al-
      location should work fine now.     */

   asm volatile ("                    # ! means not ready yet

     push %%ebp          # gcc *never* saves ebp, even if you ask it to
     .align 4
    0:
     fildl (%0)               # a0!
     fildl (%1)               # b0!  a0!
     fildl 4(%0)              # a1!  b0!  a0!
     fildl 4(%1)              # b1!  a1!  b0!  a0
     fildl 8(%0)              # a2!  b1!  a1!  b0   a0
     fildl 8(%1)              # b2!  a2!  b1!  a1   b0   a0
     fildl 12(%0)             # a3!  b2!  a2!  b1   a1   b0   a0
     fildl 12(%1)             # b3!  a3!  b2!  a2   b1   a1   b0   a0
     fxch %%st(7)             # a0   a3!  b2!  a2   b1   a1   b0   b3!
     fmulp %%st, %%st(6)      # a3!  b2   a2   b1   a1   ab0! b3!
     fxch %%st(4)             # a1   b2   a2   b1   a3!  ab0! b3!
     fmulp %%st, %%st(3)      # b2   a2   ab1! a3   ab0! b3        stall
     fld %%st(4)              # ab0  b2   a2   ab1! a3   ab0  b3
     fmull _JSP_stuff+48      # q0!  b2   a2   ab1! a3   ab0  b3
     fld %%st(3)              # ab1  q0!  b2   a2   ab1  a3   ab0  b3
     fmull _JSP_stuff+48      # q1!  q0!  b2   a2   ab1  a3   ab0  b3
     fxch %%st(1)             # q0!  q1!  b2   a2   ab1  a3   ab0  b3
     faddl _JSP_stuff         # q0+! q1!  b2   a2   ab1  a3   ab0  b3
     fxch %%st(2)             # b2   q1!  q0+! a2   ab1  a3   ab0  b3
     fmulp %%st, %%st(3)      # q1!  q0+! ab2! ab1  a3   ab0  b3
     faddl _JSP_stuff         # q1+! q0+! ab2! ab1  a3   ab0  b3
     fxch %%st(4)             # a3   q0+! ab2! ab1  q1+! ab0  b3
     fmulp %%st, %%st(6)      # q0+  ab2! ab1  q1+! ab0  ab3!
     fsubl _JSP_stuff         # q0-! ab2  ab1  q1+! ab0  ab3!
     fld %%st(1)              # ab2  q0-! ab2  ab1  q1+  ab0  ab3!
     fmull _JSP_stuff+48      # q2!  q0-! ab2  ab1  q1+  ab0  ab3!
     fld %%st(6)              # ab3  q2!  q0-  ab2  ab1  q1+  ab0  ab3
     fmull _JSP_stuff+48      # q3!  q2!  q0-  ab2  ab1  q1+  ab0  ab3
     fxch %%st(5)             # q1+  q2!  q0-  ab2  ab1  q3!  ab0  ab3
     fsubl _JSP_stuff         # q1-! q2   q0-  ab2  ab1  q3!  ab0  ab3
     fxch %%st(2)             # q0-  q2   q1-! ab2  ab1  q3!  ab0  ab3
     fmull _JSP_stuff+16      # nq0! q2   q1-! ab2  ab1  q3!  ab0  ab3
     fxch %%st(1)             # q2   nq0! q1-! ab2  ab1  q3!  ab0  ab3
     faddl _JSP_stuff         # q2+! nq0! q1-! ab2  ab1  q3   ab0  ab3
     fxch %%st(5)             # q3   nq0! q1-! ab2  ab1  q2+! ab0  ab3
     faddl _JSP_stuff         # q3+! nq0! q1-  ab2  ab1  q2+! ab0  ab3
     fxch %%st(2)             # q1-  nq0! q3+! ab2  ab1  q2+! ab0  ab3
     fmull _JSP_stuff+16      # nq1! nq0  q3+! ab2  ab1  q2+! ab0  ab3
     fxch %%st(5)             # q2+! nq0  q3+! ab2  ab1  nq1! ab0  ab3
     fsubl _JSP_stuff         # q2-! nq0  q3+! ab2  ab1  nq1! ab0  ab3
     fxch %%st(1)             # nq0  q2-! q3+! ab2  ab1  nq1! ab0  ab3
     fsubrp %%st, %%st(6)     # q2-! q3+! ab2  ab1  nq1! R0!  ab3
     fxch %%st(1)             # q3+! q2-! ab2  ab1  nq1! R0!  ab3
     fsubl _JSP_stuff         # q3-! q2-! ab2  ab1  nq1  R0!  ab3
     fxch %%st(4)             # nq1  q2-! ab2  ab1  q3-! R0!  ab3
     fsubrp %%st, %%st(3)     # q2-  ab2  R1!  q3-! R0!  ab3
     fmull _JSP_stuff+16      # nq2! ab2  R1!  q3-! R0   ab3
     fxch %%st(4)             # R0   ab2  R1!  q3-! nq2! ab3
     faddl _JSP_stuff+8       # R00! ab2  R1!  q3-  nq2! ab3
     fxch %%st(3)             # q3-  ab2  R1!  R00! nq2! ab3
     fmull _JSP_stuff+16      # nq3! ab2  R1   R00! nq2! ab3
     fxch %%st(4)             # nq2! ab2  R1   R00! nq3! ab3
     fsubrp %%st, %%st(1)     # R2!  R1   R00! nq3! ab3
     fxch %%st(1)             # R1   R2!  R00! nq3! ab3
     faddl _JSP_stuff+8       # R11! R2!  R00  nq3! ab3
     fxch %%st(2)             # R00  R2!  R11! nq3! ab3
     fstpl _JSP_stuff+24      # R2   R11! nq3  ab3
     faddl _JSP_stuff+8       # R22! R11  nq3  ab3
     fxch %%st(3)             # ab3  R11  nq3  R22!
     fsubp %%st, %%st(2)      # R11  R3!  R22!
     fstpl _JSP_stuff+32      # R3!  R22
     faddl _JSP_stuff+8       # R33! R22
     fxch %%st(1)             # R22  R33!
     fstpl _JSP_stuff+40      # R33!

     movl _JSP_stuff+24, %%esi
     movl _JSP_stuff+24, %%edi
     sarl $31, %%esi
     movl _JSP_stuff+32, %%ebp
     sarl $31, %%ebp
     andl %2, %%esi
     addl %%esi, %%edi
     movl _JSP_stuff+32, %%esi
     andl %2, %%ebp
     movl %%edi, (%0)
     addl %%ebp, %%esi
     movl _JSP_stuff+40, %%edi
     sarl $31, %%edi
     movl %%esi, 4(%0)
     fstpl _JSP_stuff+24
     movl _JSP_stuff+40, %%esi
     movl _JSP_stuff+24, %%ebp
     sarl $31, %%ebp
     andl %2, %%edi
     addl %%edi, %%esi
     movl _JSP_stuff+24, %%edi
     andl %2, %%ebp
     movl %%esi, 8(%0)
     addl %%ebp, %%edi
     movl %%edi, 12(%0)
     nop

     addl $16, %0
     addl $16, %1
     subl $4, %3
     jnz 0b

     popl %%ebp
   ":
    : "a"(Num1), "b"(Num2), "c"(JPrime), "d"(Len)
    : "%esi", "%edi", "%ebp", "memory" );
}

void
VectorModMulC(ModInt *Num1, ModInt Num2, size_t Len)
{ModInt *Num2Ptr=&Num2; /* The code needs it to be a pointer */
   /*  Num1[0...len] <- Num1[0...len] * Num2 mod n

      len must be a multiple of 4. Running time: 64 cycles
      every 4 multiplies. The loop itself is coded here too, to
      avoid the overhead of pushing and popping registers all
      the time (saves about 10 cycles). The final answers are
      always positive; even with sign correction, using the FPU
      is almost four times faster than using mul/div, though
      *much* more complicated.

      A general rule with gcc 2.7.2.1 is that if you ask for
      all the ALU registers, invariably something will go very
      wrong. This routine was no execption, but the register al-
      location should work fine now.     */

   asm volatile ("                    # ! means not ready yet

     push %%ebp          # gcc *never* saves ebp, even if you ask it to
     .align 4
    0:
     fildl (%0)               # a0!
     fildl (%1)               # b0!  a0!
     fildl 4(%0)              # a1!  b0!  a0!
     fildl (%1)               # b1!  a1!  b0!  a0
     fildl 8(%0)              # a2!  b1!  a1!  b0   a0
     fildl (%1)               # b2!  a2!  b1!  a1   b0   a0
     fildl 12(%0)             # a3!  b2!  a2!  b1   a1   b0   a0
     fildl (%1)               # b3!  a3!  b2!  a2   b1   a1   b0   a0
     fxch %%st(7)             # a0   a3!  b2!  a2   b1   a1   b0   b3!
     fmulp %%st, %%st(6)      # a3!  b2   a2   b1   a1   ab0! b3!
     fxch %%st(4)             # a1   b2   a2   b1   a3!  ab0! b3!
     fmulp %%st, %%st(3)      # b2   a2   ab1! a3   ab0! b3        stall
     fld %%st(4)              # ab0  b2   a2   ab1! a3   ab0  b3
     fmull _JSP_stuff+48      # q0!  b2   a2   ab1! a3   ab0  b3
     fld %%st(3)              # ab1  q0!  b2   a2   ab1  a3   ab0  b3
     fmull _JSP_stuff+48      # q1!  q0!  b2   a2   ab1  a3   ab0  b3
     fxch %%st(1)             # q0!  q1!  b2   a2   ab1  a3   ab0  b3
     faddl _JSP_stuff         # q0+! q1!  b2   a2   ab1  a3   ab0  b3
     fxch %%st(2)             # b2   q1!  q0+! a2   ab1  a3   ab0  b3
     fmulp %%st, %%st(3)      # q1!  q0+! ab2! ab1  a3   ab0  b3
     faddl _JSP_stuff         # q1+! q0+! ab2! ab1  a3   ab0  b3
     fxch %%st(4)             # a3   q0+! ab2! ab1  q1+! ab0  b3
     fmulp %%st, %%st(6)      # q0+  ab2! ab1  q1+! ab0  ab3!
     fsubl _JSP_stuff         # q0-! ab2  ab1  q1+! ab0  ab3!
     fld %%st(1)              # ab2  q0-! ab2  ab1  q1+  ab0  ab3!
     fmull _JSP_stuff+48      # q2!  q0-! ab2  ab1  q1+  ab0  ab3!
     fld %%st(6)              # ab3  q2!  q0-  ab2  ab1  q1+  ab0  ab3
     fmull _JSP_stuff+48      # q3!  q2!  q0-  ab2  ab1  q1+  ab0  ab3
     fxch %%st(5)             # q1+  q2!  q0-  ab2  ab1  q3!  ab0  ab3
     fsubl _JSP_stuff         # q1-! q2   q0-  ab2  ab1  q3!  ab0  ab3
     fxch %%st(2)             # q0-  q2   q1-! ab2  ab1  q3!  ab0  ab3
     fmull _JSP_stuff+16      # nq0! q2   q1-! ab2  ab1  q3!  ab0  ab3
     fxch %%st(1)             # q2   nq0! q1-! ab2  ab1  q3!  ab0  ab3
     faddl _JSP_stuff         # q2+! nq0! q1-! ab2  ab1  q3   ab0  ab3
     fxch %%st(5)             # q3   nq0! q1-! ab2  ab1  q2+! ab0  ab3
     faddl _JSP_stuff         # q3+! nq0! q1-  ab2  ab1  q2+! ab0  ab3
     fxch %%st(2)             # q1-  nq0! q3+! ab2  ab1  q2+! ab0  ab3
     fmull _JSP_stuff+16      # nq1! nq0  q3+! ab2  ab1  q2+! ab0  ab3
     fxch %%st(5)             # q2+! nq0  q3+! ab2  ab1  nq1! ab0  ab3
     fsubl _JSP_stuff         # q2-! nq0  q3+! ab2  ab1  nq1! ab0  ab3
     fxch %%st(1)             # nq0  q2-! q3+! ab2  ab1  nq1! ab0  ab3
     fsubrp %%st, %%st(6)     # q2-! q3+! ab2  ab1  nq1! R0!  ab3
     fxch %%st(1)             # q3+! q2-! ab2  ab1  nq1! R0!  ab3
     fsubl _JSP_stuff         # q3-! q2-! ab2  ab1  nq1  R0!  ab3
     fxch %%st(4)             # nq1  q2-! ab2  ab1  q3-! R0!  ab3
     fsubrp %%st, %%st(3)     # q2-  ab2  R1!  q3-! R0!  ab3
     fmull _JSP_stuff+16      # nq2! ab2  R1!  q3-! R0   ab3
     fxch %%st(4)             # R0   ab2  R1!  q3-! nq2! ab3
     faddl _JSP_stuff+8       # R00! ab2  R1!  q3-  nq2! ab3
     fxch %%st(3)             # q3-  ab2  R1!  R00! nq2! ab3
     fmull _JSP_stuff+16      # nq3! ab2  R1   R00! nq2! ab3
     fxch %%st(4)             # nq2! ab2  R1   R00! nq3! ab3
     fsubrp %%st, %%st(1)     # R2!  R1   R00! nq3! ab3
     fxch %%st(1)             # R1   R2!  R00! nq3! ab3
     faddl _JSP_stuff+8       # R11! R2!  R00  nq3! ab3
     fxch %%st(2)             # R00  R2!  R11! nq3! ab3
     fstpl _JSP_stuff+24      # R2   R11! nq3  ab3
     faddl _JSP_stuff+8       # R22! R11  nq3  ab3
     fxch %%st(3)             # ab3  R11  nq3  R22!
     fsubp %%st, %%st(2)      # R11  R3!  R22!
     fstpl _JSP_stuff+32      # R3!  R22
     faddl _JSP_stuff+8       # R33! R22
     fxch %%st(1)             # R22  R33!
     fstpl _JSP_stuff+40      # R33!

     movl _JSP_stuff+24, %%esi
     movl _JSP_stuff+24, %%edi
     sarl $31, %%esi
     movl _JSP_stuff+32, %%ebp
     sarl $31, %%ebp
     andl %2, %%esi
     addl %%esi, %%edi
     movl _JSP_stuff+32, %%esi
     andl %2, %%ebp
     movl %%edi, (%0)
     addl %%ebp, %%esi
     movl _JSP_stuff+40, %%edi
     sarl $31, %%edi
     movl %%esi, 4(%0)
     fstpl _JSP_stuff+24
     movl _JSP_stuff+40, %%esi
     movl _JSP_stuff+24, %%ebp
     sarl $31, %%ebp
     andl %2, %%edi
     addl %%edi, %%esi
     movl _JSP_stuff+24, %%edi
     andl %2, %%ebp
     movl %%esi, 8(%0)
     addl %%ebp, %%edi
     movl %%edi, 12(%0)
     nop

     addl $16, %0
     subl $4, %3
     jnz 0b

     popl %%ebp
   ":
    : "a"(Num1), "b"(Num2Ptr), "c"(JPrime), "d"(Len)
    : "%esi", "%edi", "%ebp", "memory" );
}

void
VectorModButterfly(ModInt *Num1, ModInt *Num2, size_t Len)
{
     /* performs "len" FFT butterflies mod n on arrays a and b.

        {a,b} -> {a+b, a-b} mod n

        It's assumed that all the numbers involved are positive,
        and the code checks that they stay that way. len must be
        a multiple of 4. Running time: 35 cycles per 4 butterflies. */

   asm("
     pushl %%ebp
    .align 4
    0:
     movl (%0), %%esi
     movl (%0), %%ebp
     pushl %3
     movl (%1), %%edi
     addl %%edi, %%ebp
     subl %%edi, %%esi
     movl %%esi, %%edi
     subl %2, %%ebp
     sarl $31, %%edi
     movl %%ebp, %3
     sarl $31, %3
     andl %2, %%edi
     addl %%edi, %%esi
     andl %2, %3
     addl %3, %%ebp
     movl %%esi, (%1)

     movl %%ebp, (%0)
     movl 4(%1), %%edi
     movl 4(%0), %%esi
     movl 4(%0), %%ebp
     addl %%edi, %%ebp
     subl %%edi, %%esi
     movl %%esi, %%edi
     subl %2, %%ebp
     sarl $31, %%edi
     movl %%ebp, %3
     sarl $31, %3
     andl %2, %%edi
     addl %%edi, %%esi
     andl %2, %3
     addl %3, %%ebp
     movl %%esi, 4(%1)

     movl %%ebp, 4(%0)
     movl 8(%1), %%edi
     movl 8(%0), %%esi
     movl 8(%0), %%ebp
     addl %%edi, %%ebp
     subl %%edi, %%esi
     movl %%esi, %%edi
     subl %2, %%ebp
     sarl $31, %%edi
     movl %%ebp, %3
     sarl $31, %3
     andl %2, %%edi
     addl %%edi, %%esi
     andl %2, %3
     addl %3, %%ebp
     movl %%esi, 8(%1)

     movl %%ebp, 8(%0)
     movl 12(%1), %%edi
     movl 12(%0), %%esi
     movl 12(%0), %%ebp
     addl %%edi, %%ebp
     subl %%edi, %%esi
     movl %%esi, %%edi
     subl %2, %%ebp
     sarl $31, %%edi
     movl %%ebp, %3
     sarl $31, %3
     andl %2, %%edi
     addl %%edi, %%esi
     andl %2, %3
     addl %3, %%ebp
     movl %%esi, 12(%1)

     popl %3
     movl %%ebp, 12(%0)
     addl $16, %0
     addl $16, %1
     subl $4, %3
     jnz 0b

     popl %%ebp
    ":
     : "a"(Num1), "b"(Num2), "c"(JPrime), "d"(Len)
     : "%esi", "%edi", "%ebp", "memory" );
}

void
VectorNTT_First2(ModInt *Data, ModInt Trig, size_t Len)
/* Perform the first two passes of the DiT NTT */
{
/*
** The asm expects the Trig to be in JSP_stuff[3] so we have to humor it
*/
JSP_stuff[3]=Trig;

   /* performs the first two levels of the decimation in time NTT
      for "len" points in a row. len above must be a multiple of 4.
      Running time: 48 cycles per 4 array elements. This involves
      four butterflies but only one modular multiplication (by the
      contents of JSP_stuff[3] ) per group of 4 points. The mod step
      for each butterfly and the multiply is *very* tedious, and
      more than doubles the runtime. There are also many more stalls
      here than elsewhere, and they're unavoidable for want of more
      registers.                                       */
   asm("
     pushl %%ebp
    .align 4
    0:
     movl 8(%0), %%ebp
     movl 12(%0), %%edi
     leal (%%ebp,%%edi), %%esi
     subl %%edi, %%ebp
     movl %%ebp, %%edi
     subl %1, %%esi
     sarl $31, %%edi
     movl %%esi, %%edx
     sarl $31, %%edx
     andl %1, %%edi
     addl %%edi, %%ebp
     andl %1, %%edx
     addl %%edx, %%esi
     movl %%ebp, 12(%0)

     movl %%esi, 8(%0)
    fildl 12(%0)
     movl (%0), %%ebp
     movl 4(%0), %%edi
     leal (%%ebp,%%edi), %%esi
     subl %%edi, %%ebp
    fmull _JSP_stuff+24
     subl %1, %%esi
     movl %%ebp, %%edi
     sarl $31, %%edi
     movl %%esi, %%edx
    fld %%st
    fmull _JSP_stuff+48
     sarl $31, %%edx
     andl %1, %%edi
     addl %%edi, %%ebp
     andl %1, %%edx
    faddl _JSP_stuff
     addl %%edx, %%esi
     movl %%ebp, 4(%0)

     movl 8(%0), %%edi
    fsubl _JSP_stuff
     leal (%%esi,%%edi), %%ebp
     subl %%edi, %%esi
     movl %%esi, %%edi
     subl %1, %%ebp
    fmull _JSP_stuff+16
     sarl $31, %%edi
     movl %%ebp, %%edx
     sarl $31, %%edx
     andl %1, %%edi
    fsubrp %%st, %%st(1)
     addl %%edi, %%esi
     andl %1, %%edx
     addl %%edx, %%ebp
     movl %%esi, 8(%0)
    faddl _JSP_stuff+8
     movl %%ebp, (%0)
     movl 4(%0), %%ebp

    fstpl _JSP_stuff+32         # stall x 2
     movl _JSP_stuff+32, %%esi
     addl $16, %0
     sarl $31, %%esi
     movl _JSP_stuff+32, %%edi
     andl %1, %%esi
     addl %%esi, %%edi          # stall
     leal (%%ebp,%%edi), %%esi
     subl %%edi, %%ebp
     movl %%ebp, %%edi
     subl %1, %%esi
     sarl $31, %%edi
     movl %%esi, %%edx
     sarl $31, %%edx
     andl %1, %%edi
     addl %%edi, %%ebp
     andl %1, %%edx
     addl %%edx, %%esi
     movl %%ebp, -4(%0)

     movl %%esi, -12(%0)
     nop
     subl $4, %2
     jnz 0b

     popl %%ebp
    ":
     : "a"(Data), "b"(JPrime), "c"(Len)
     : "%edx", "%esi", "%edi", "%ebp", "memory" );
}

void
VectorNTT_Last2(ModInt *Data, ModInt Trig, size_t Len)
/* Perform the last two passes of the DiF NTT */
{
/*
** The asm expects the Trig to be in JSP_stuff[3] so we have to humor it
*/
JSP_stuff[3]=Trig;

   /* performs the first two levels of the decimation in time NTT
      for "len" points in a row. len above must be a multiple of 4.
      Running time: 48 cycles per 4 array elements. This involves
      four butterflies but only one modular multiplication (by the
      contents of JSP_stuff[3] ) per group of 4 points. The mod step
      for each butterfly and the multiply is *very* tedious, and
      more than doubles the runtime. There are also many more stalls
      here than elsewhere, and they're unavoidable for want of more
      registers.                                       */
   asm("
     pushl %%ebp
    .align 4
    0:
# do modadd/sub d1 & d3
     movl 4(%0), %%ebp  # was 8 load d1
     movl 12(%0), %%edi #       load d3
     leal (%%ebp,%%edi), %%esi
     subl %%edi, %%ebp
     movl %%ebp, %%edi
     subl %1, %%esi
     sarl $31, %%edi
     movl %%esi, %%edx
     sarl $31, %%edx
     andl %1, %%edi
     addl %%edi, %%ebp
     andl %1, %%edx
     addl %%edx, %%esi
     movl %%ebp, 12(%0) #       store d3
     movl %%esi, 4(%0)  # was 8 store d1

# interleave modmul of d3
    fildl 12(%0)       # load d3 to mul
# do modadd/sub of d0 & d2
     movl (%0), %%ebp  #       load d0
     movl 8(%0), %%edi # was 4 load d2
     leal (%%ebp,%%edi), %%esi
     subl %%edi, %%ebp
    fmull _JSP_stuff+24
     subl %1, %%esi
     movl %%ebp, %%edi
     sarl $31, %%edi
     movl %%esi, %%edx
    fld %%st
    fmull _JSP_stuff+48
     sarl $31, %%edx
     andl %1, %%edi
     addl %%edi, %%ebp
     andl %1, %%edx
    faddl _JSP_stuff
     addl %%edx, %%esi
     movl %%ebp, 8(%0) # was 4  store d2

# no need to store d0 because we are going to use it.
     movl 4(%0), %%edi # was 8  load d1
    fsubl _JSP_stuff
     leal (%%esi,%%edi), %%ebp
     subl %%edi, %%esi
     movl %%esi, %%edi
     subl %1, %%ebp
    fmull _JSP_stuff+16
     sarl $31, %%edi
     movl %%ebp, %%edx
     sarl $31, %%edx
     andl %1, %%edi
    fsubrp %%st, %%st(1)
     addl %%edi, %%esi
     andl %1, %%edx
     addl %%edx, %%ebp
     movl %%esi, 4(%0)  # was 4 store d1
    faddl _JSP_stuff+8
     movl %%ebp, (%0)   # store d0
     movl 8(%0), %%ebp  # was 4 load d2

    fstpl _JSP_stuff+32         # stall x 2
     movl _JSP_stuff+32, %%esi # load d3 (just finished modmul)
     addl $16, %0
     sarl $31, %%esi
     movl _JSP_stuff+32, %%edi
     andl %1, %%esi
     addl %%esi, %%edi          # stall
     leal (%%ebp,%%edi), %%esi
     subl %%edi, %%ebp
     movl %%ebp, %%edi
     subl %1, %%esi
     sarl $31, %%edi
     movl %%esi, %%edx
     sarl $31, %%edx
     andl %1, %%edi
     addl %%edi, %%ebp
     andl %1, %%edx
     addl %%edx, %%esi
     movl %%ebp, -4(%0) # store d3

     movl %%esi, -8(%0) # was -12 store d2
     nop
     subl $4, %2
     jnz 0b

     popl %%ebp
    ":
     : "a"(Data), "b"(JPrime), "c"(Len)
     : "%edx", "%esi", "%edi", "%ebp", "memory" );
}

void
PrepVector(UINT32 Prime, size_t Len)
/*
** Prepare for a vector of length 'Len' with prime 'Prime'.  This
** will be valid until this routine is called again, not just for
** the next vector opeartion.
*/
{
 JSP_stuff[2] = Prime;
 JSP_stuff[6] = 1.0/Prime;
 JPrime=Prime;
}

/*
** Initialize for vector mod math.  This also includes checking the
** rounding mode, which must be in 'round to nearest' mode.  If it
** isn't, then setting it will be system & compiler dependant.
**
** Even though FLT_ROUNDS can be a value (rather than a constant) and
** can change during the program run, rather than checking it in
** every VectorMod...() call, I'll assume that it's not going to change
** during the program run.
*/
void
InitVectorModMath(size_t Len)
{volatile long double LD;

/*
** This assumes that the radix is 2!!  Too much trouble to try and code it
** for any other radix.
*/
if (LDBL_MANT_DIG < 64)
  FatalError("This version needs long doubles to be at least 64 bits.  See vector.c\n");

LD=1.0;
LD+=LDBL_EPSILON;
if (LD==1.0)
  FatalError("It seems that you can't do long double math.  Sorry.  See vector.c\n");

 if (FLT_ROUNDS != 1)
   FatalError("It appears that the FPU isn't rounding to nearest.\nSee docs, version.txt, and vector.c");

 JSP_stuff[0] = 3.0*2147483648.0*2147483648.0;
 JSP_stuff[1] = 6755399441055744.0;

CHECK_ALIGNMENT(JSP_stuff);
}

void
DeInitVectorModMath(void)
{
}


