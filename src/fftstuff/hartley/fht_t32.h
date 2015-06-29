#include <math.h>

#define SQRT12   0.70710678118654752440084436210484
#define SQRT2    1.414213562373095048801688724209698078569
#define sumdiff4(a,b,s,d)  { s=a+b; d=a-b; }
#define cmult6(c,s,c1,s1,u,v)  { u=c1*c-s1*s; v=c1*s+s1*c; }
#define sumdiff2(s,d)  { double t=s-d; s+=d; d=t; }

void
fht_t32(double *f)  // for length 32
{
//scramble(f,32);

{
double s1,c1,s2,c2;
double g0,f0,f1,g1;
// no pow 4 loop:
sumdiff4(f[0],f[1],s1,c1);
sumdiff4(f[2],f[3],s2,c2);
sumdiff4(s1,s2,f0,f1);
sumdiff4(c1,c2,g0,g1);
sumdiff4(f[4],f[5],s1,c1);
sumdiff4(f[6],f[7],s2,c2);
sumdiff2(s1,s2);
sumdiff4(f0,s1,f[0],f[4]);
sumdiff4(f1,s2,f[2],f[6]);
c1 *= SQRT2;
c2 *= SQRT2;
sumdiff4(g0,c1,f[1],f[5]);
sumdiff4(g1,c2,f[3],f[7]);
// no pow 4 loop:
sumdiff4(f[8],f[9],s1,c1);
sumdiff4(f[10],f[11],s2,c2);
sumdiff4(s1,s2,f0,f1);
sumdiff4(c1,c2,g0,g1);
sumdiff4(f[12],f[13],s1,c1);
sumdiff4(f[14],f[15],s2,c2);
sumdiff2(s1,s2);
sumdiff4(f0,s1,f[8],f[12]);
sumdiff4(f1,s2,f[10],f[14]);
c1 *= SQRT2;
c2 *= SQRT2;
sumdiff4(g0,c1,f[9],f[13]);
sumdiff4(g1,c2,f[11],f[15]);
// no pow 4 loop:
sumdiff4(f[16],f[17],s1,c1);
sumdiff4(f[18],f[19],s2,c2);
sumdiff4(s1,s2,f0,f1);
sumdiff4(c1,c2,g0,g1);
sumdiff4(f[20],f[21],s1,c1);
sumdiff4(f[22],f[23],s2,c2);
sumdiff2(s1,s2);
sumdiff4(f0,s1,f[16],f[20]);
sumdiff4(f1,s2,f[18],f[22]);
c1 *= SQRT2;
c2 *= SQRT2;
sumdiff4(g0,c1,f[17],f[21]);
sumdiff4(g1,c2,f[19],f[23]);
// no pow 4 loop:
sumdiff4(f[24],f[25],s1,c1);
sumdiff4(f[26],f[27],s2,c2);
sumdiff4(s1,s2,f0,f1);
sumdiff4(c1,c2,g0,g1);
sumdiff4(f[28],f[29],s1,c1);
sumdiff4(f[30],f[31],s2,c2);
sumdiff2(s1,s2);
sumdiff4(f0,s1,f[24],f[28]);
sumdiff4(f1,s2,f[26],f[30]);
c1 *= SQRT2;
c2 *= SQRT2;
sumdiff4(g0,c1,f[25],f[29]);
sumdiff4(g1,c2,f[27],f[31]);
}

{ // k4=32
double f0,f1,f2,f3;
// do loop:
sumdiff4(f[0],f[8],f0,f1);
sumdiff4(f[16],f[24],f2,f3);
sumdiff4(f0,f2,f[0],f[16]);
sumdiff4(f1,f3,f[8],f[24]);
sumdiff4(f[4],f[12],f0,f1);
f3 = SQRT2 * f[28];
f2 = SQRT2 * f[20];
sumdiff4(f0,f2,f[4],f[20]);
sumdiff4(f1,f3,f[12],f[28]);
}

{ // kx=4
double s1, c1;
double a,b, g0,f0,f1,g1, f2,g2,f3,g3;
double c2, s2;
c1=0.980785280403230430579242238;
s1=0.195090322016128248083788321;
c2=0.923879532511286738483136105;
s2=0.382683432365089726268081449;
// do loop II:
cmult6(s2,c2,f[9],f[15],b,a);
sumdiff4(f[1],a,f0,f1);
sumdiff4(f[7],b,g0,g1);
cmult6(s2,c2,f[25],f[31],b,a);
sumdiff4(f[17],a,f2,f3);
sumdiff4(f[23],b,g2,g3);
cmult6(s1,c1,f2,g3,b,a);
sumdiff4(f0,a,f[1],f[17]);
sumdiff4(g1,b,f[15],f[31]);
cmult6(c1,s1,g2,f3,b,a);
sumdiff4(g0,a,f[7],f[23]);
sumdiff4(f1,b,f[9],f[25]);
c1=0.923879532511286738483136105;
s1=0.38268343236508978177923268;
// do loop II:
sumdiff4(f[10],f[14],a,b);
a *= SQRT12;
b *= SQRT12;
sumdiff4(f[2],a,f0,f1);
sumdiff4(f[6],b,g0,g1);
sumdiff4(f[26],f[30],a,b);
a *= SQRT12;
b *= SQRT12;
sumdiff4(f[18],a,f2,f3);
sumdiff4(f[22],b,g2,g3);
cmult6(s1,c1,f2,g3,b,a);
sumdiff4(f0,a,f[2],f[18]);
sumdiff4(g1,b,f[14],f[30]);
cmult6(c1,s1,g2,f3,b,a);
sumdiff4(g0,a,f[6],f[22]);
sumdiff4(f1,b,f[10],f[26]);
c1=0.831469612302545235671402679;
s1=0.555570233019602177648721408;
c2=0.382683432365089837290383912;
s2=0.923879532511286627460833643;
// do loop II:
cmult6(s2,c2,f[11],f[13],b,a);
sumdiff4(f[3],a,f0,f1);
sumdiff4(f[5],b,g0,g1);
cmult6(s2,c2,f[27],f[29],b,a);
sumdiff4(f[19],a,f2,f3);
sumdiff4(f[21],b,g2,g3);
cmult6(s1,c1,f2,g3,b,a);
sumdiff4(f0,a,f[3],f[19]);
sumdiff4(g1,b,f[13],f[29]);
cmult6(c1,s1,g2,f3,b,a);
sumdiff4(g0,a,f[5],f[21]);
sumdiff4(f1,b,f[11],f[27]);
}

// opct:  #mult=63=1.96875/pt   #add=180=5.625/pt
}  // end of fht_32()
