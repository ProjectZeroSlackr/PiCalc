/*---------------------------------------------------------------------*/
void single_fft( Cmplx *c, Cmplx *w, long size );
void single_ifft( Cmplx *c, Cmplx *w, long size );

static void single_fft2( Cmplx *c0, Cmplx *c1 );
static void single_splitrad_DIF( Cmplx *c, Cmplx *w, long size );
static void single_splitrad_DIT( Cmplx *c, Cmplx *w, long size );

void single_fft( Cmplx *c, Cmplx *w, long size ) {

   /* performs a single forward (DiF) split-radix FFT */
   if( size==2 ) {
      single_fft2(c, c+1);
   }
   else if ( size==4 ) {
      single_splitrad_DIF( c, w, size );
      single_fft2(c, c+1);
   }
   else {
      single_splitrad_DIF( c, w, size );
      single_fft(c, w+size/2, size/2);
      single_fft(c+size/2, w+3*size/4, size/4);
      single_fft(c+3*size/4, w+3*size/4, size/4);
   }
}

/*---------------------------------------------------------------------*/
void single_ifft( Cmplx *c, Cmplx *w, long size ) {

   /* performs a single inverse (DiT) split-radix FFT */
   if( size==2 ) {
      single_fft2(c, c+1);
   }
   else if ( size==4 ) {
      single_fft2(c, c+1);
      single_splitrad_DIT( c, w, size );
   }
/*  From some other code Jason sent me.  These funcs are in asm, not C.
   if( size==4 ) {
      rad4_DIT(c);
   }
   else if ( size==8 ) {
      rad4_DIT(c);
      rad22(c+4);
      single_splitrad_DIT( c, w, size );
   }
*/
   else {
      single_ifft(c, w+size/2, size/2);
      single_ifft(c+size/2, w+3*size/4, size/4);
      single_ifft(c+3*size/4, w+3*size/4, size/4);
      single_splitrad_DIT( c, w, size );
   }
}


/*---------------------------------------------------------------------*/

      /*-------------------------------------------------*/
      /* These are the lowest level routines. There's no */
      /* assembly language here, but the C is written so */
      /* that certain compilers produce efficient code.  */
      /* The default is for intel-based gcc, and makes   */
      /* efficient use of the Pentium FPU (the DIT code  */
      /* much more than the DIF code). Everything in the */
      /* RISC section produces efficient code for Sun cc */
      /* on an UltraSPARC, and more aggressive compilers */
      /* will do even better. BE CAREFUL; the two dif-   */
      /* ferent coding styles are much slower when used  */
      /* on the wrong architecture.                      */
      /*                                                 */
      /* Both are also a terrible mess. Abandon all hope */
      /* ye who press enter here.                        */
      /*-------------------------------------------------*/

#define RISC_DIF(q)                                                   \
      f10 = w[0].r;   f11 = w[0].i;                                   \
      f12 = w[1].r;   f13 = w[1].i;                                   \
      for(i=0;i<q;i++) {                                              \
         f0 = c0[i].r;   f1 = c0[dist+i].r;  f2=f0+f1; f6=f0-f1;      \
         f0 = c0[i].i;   f1 = c0[dist+i].i;  f3=f0+f1; f7=f0-f1;      \
         f0 = c1[i].i;   f1 = c1[dist+i].i;  f5=f0+f1; f8=f0-f1;      \
         f0 = c1[i].r;   f1 = c1[dist+i].r;  f4=f0+f1; f9=f1-f0;      \
         c0[i].r = f2;   f2=f6+f8;         c0[i].i = f3;   f3=f7+f9;  \
         c1[i].r = f4;   f4=f6-f8;         c1[i].i = f5;   f5=f7-f9;  \
         f6 = f2 * f10;        f7 = f3 * f11;                         \
         f8 = f4 * f12;        f9 = f5 * f13;                         \
         f2 = f2 * f11;        f3 = f3 * f10;    f6 = f6 - f7;        \
         f4 = f4 * f13;        f5 = f5 * f12;    f8 = f8 - f9;        \
         f2 = f2 + f3;        c0[dist+i].r = f6;                      \
         f4 = f4 + f5;        c1[dist+i].r = f8;                      \
         c0[dist+i].i = f2;                                           \
         c1[dist+i].i = f4;                                           \
      }

#define PENTIUM_DIF(q)                                                     \
        f0 = c0[q].r;      f1 = c0[dist+q].r;      f0 = f0 - f1;           \
        f1 = c0[q].i;      f2 = c0[dist+q].i;      f1 = f1 - f2;           \
        f2 = c1[q].i;      f3 = c1[dist+q].i;      f2 = f2 - f3;           \
        f3 = c1[q].r;      f4 = c1[dist+q].r;      f3 = f4 - f3;           \
        f4 = f0;           f0 = f0 + f2;         f2 = f4 - f2;             \
        f4 = f1;           f1 = f1 + f3;         f3 = f4 - f3;             \
        f4 = f0;           f5 = w[0].i; f4 = f4 * f5;    f0 = f0 * w[0].r; \
        f5 = f1;           f6 = w[0].r; f5 = f5 * f6;    f1 = f1 * w[0].i; \
        f6 = c0[q].r;      f6 = f6 + c0[dist+q].r; f4 = f4 + f5;           \
        f5 = c0[q].i;      f5 = f5 + c0[dist+q].i; f0 = f0 - f1;           \
        c0[q].r = f6;      f6 = f2;              f6 = f6 * w[1].r;         \
        c0[q].i = f5;      f5 = f3;              f5 = f5 * w[1].i;         \
        c0[dist+q].r = f0;   f2 = f2 * w[1].i;                             \
        c0[dist+q].i = f4;   f4 = c1[q].r;         f3 = f3 * w[1].r;       \
        f4 = f4 + c1[dist+q].r;                                            \
        f6 = f6 - f5;      f5 = c1[q].i;         f5 = f5 + c1[dist+q].i;   \
        f2 = f2 + f3;      c1[q].r = f4;         c1[dist+q].r = f6;        \
                           c1[q].i = f5;         c1[dist+q].i = f2;

#define RISC_DIT(q)                                                        \
      f12 = w[0].r;   f13 = w[0].i;                                        \
      f14 = w[1].r;   f15 = w[1].i;                                        \
      for(i=0;i<q;i++) {                                                   \
         f0 = c0[dist+i].r;  f1 = c0[dist+i].i;                            \
         f2 = c1[dist+i].r;  f3 = c1[dist+i].i;                            \
         f4 = f0 * f12;    f5 = f1 * f13;                                  \
         f8 = f2 * f14;    f9 = f3 * f15;                                  \
         f6 = f0 * f13;    f7 = f1 * f12;  f4 = f4 + f5;                   \
         f10 = f2 * f15;   f11 = f3 * f14; f8 = f8 + f9;                   \
         f6 = f7 - f6;     f0 = c0[i].r;  f10 = f11 - f10; f1 = c0[i].i;   \
         f5 = f4 + f8;     f2 = c1[i].r;  f9 = f4 - f8;    f3 = c1[i].i;   \
         f7 = f6 + f10;    f11 = f6 - f10;                                 \
                                                                           \
         f4 = f0 + f5;   f8 = f0 - f5;                                     \
         f6 = f1 + f7;   f10 = f1 - f7;  c0[i].r = f4;                     \
         f0 = f2 - f11;  c0[dist+i].r = f8;                                \
         f5 = f2 + f11;  c0[i].i = f6;                                     \
         f1 = f3 + f9;   c0[dist+i].i = f10;                               \
         f7 = f3 - f9;   c1[i].r = f0;  c1[dist+i].r = f5;                 \
                         c1[i].i = f1;  c1[dist+i].i = f7;                 \
      }

#define PENTIUM_DIT(q)                                                 \
      f0 = c0[dist+q].r;  f1 = w[0].r;  f0 = f0 * f1;                  \
      f1 = c0[dist+q].i;  f2 = w[0].i;  f1 = f1 * f2;                  \
      f2 = c0[dist+q].r;  f3 = w[0].i;  f2 = f2 * f3;                  \
      f3 = c0[dist+q].i;  f4 = w[0].r;  f3 = f3 * f4;                  \
      f4 = c1[dist+q].r;  f5 = w[1].r;  f4 = f4 * f5;                  \
      f5 = c1[dist+q].i;  f6 = w[1].i;  f5 = f5 * f6;                  \
      f6 = c1[dist+q].r;  f7 = w[1].i;  f6 = f6 * f7;                  \
      f7 = c1[dist+q].i;  f7 = f7 * w[1].r;                            \
      f0 = f0 + f1;     f2 = f3 - f2;   f4 = f4 + f5;   f6 = f7 - f6;  \
      f1 = f0;          f1 = f1 - f4;   f0 = f0 + f4;                  \
      f3 = f2;          f3 = f3 - f6;   f2 = f2 + f6;                  \
      f4 = c0[q].r;     f4 = f4 + f0;                                  \
      f5 = c0[q].r;     f0 = f5 - f0;                                  \
      f5 = c0[q].i;     f5 = f5 + f2;  c0[q].r = f4;                   \
      f6 = c0[q].i;     f2 = f6 - f2;  c0[dist+q].r = f0;              \
      f4 = c1[q].r;     f4 = f4 - f3;  c0[q].i = f5;                   \
      f0 = c1[q].r;     f3 = f0 + f3;  c0[dist+q].i = f2;              \
      f5 = c1[q].i;     f5 = f5 + f1;  c1[q].r = f4;                   \
      f2 = c1[q].i;     f1 = f2 - f1;  c1[dist+q].r = f3;              \
                                       c1[q].i = f5;                   \
                                       c1[dist+q].i = f1;

/*---------------------------------------------------------------------*/

static void single_fft2( Cmplx *c0, Cmplx *c1 ) {

#ifdef RISC
   double f0,f1,f2,f3,f6,f7;
#else
   double f0,f1,f2,f3,f4;
#endif

#ifdef RISC

   f0 = c0[0].r;   f1 = c1[0].r;  f2=f0+f1; f6=f0-f1;
   f0 = c0[0].i;   f1 = c1[0].i;  f3=f0+f1; f7=f0-f1;
   c0[0].r = f2;  c0[0].i = f3;   c1[0].r = f6;  c1[0].i = f7;

#else

   f0 = c0[0].r;   f1 = c1[0].r;   f0 = f0 + f1;
   f1 = c0[0].r;   f2 = c1[0].r;   f1 = f1 - f2;
   f2 = c0[0].i;   f3 = c1[0].i;   f2 = f2 + f3;
   f3 = c0[0].i;   f4 = c1[0].i;   f3 = f3 - f4;

   c0[0].r = f0;   c1[0].r = f1;
   c0[0].i = f2;   c1[0].i = f3;

#endif

}
/*---------------------------------------------------------------------*/

static void single_splitrad_DIF( Cmplx *c, Cmplx *w, long size ) {

#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13;
   long   i;
#else
   double f0,f1,f2,f3,f4,f5,f6;
#endif

   long dist = size/4;
   Cmplx *c0 = c,
        *c1 = c + dist;

   dist = 2*dist;
   size = size/4;

#ifdef RISC

   do {
      RISC_DIF(1);

      c0++;  c1++;
      w+=2;   size--;
   } while(size);

#else

   do {
        PENTIUM_DIF(0)

        c0++;  c1++;
        w+=2;   size--;
   } while(size);

#endif
}


/*---------------------------------------------------------------------*/

static void single_splitrad_DIT( Cmplx *c, Cmplx *w, long size ) {

#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
   long   i;
#else
   double f0,f1,f2,f3,f4,f5,f6,f7;
#endif

   long dist = size/4;
   Cmplx *c0 = c,
        *c1 = c + dist;

   dist = 2*dist;
   size = size/4;

#ifdef RISC

   do {
      RISC_DIT(1)

      c0++;   c1++;
      w+=2;   size--;
   } while(size);

#else

   do {
      PENTIUM_DIT(0)

      c0++;   c1++;
      w+=2;   size--;
   } while(size);

#endif

}

#undef RISC_DIF(q)
#undef RISC_DIT(q)
#undef PENTIUM_DIF(q)
#undef PENTIUM_DIT(q)


