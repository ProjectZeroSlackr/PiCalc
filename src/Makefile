# Makefile for Pi AGM 2.3.1, Carey Bloodworth's
# premiere Pi calculation program for *nix or DOS/gcc.
# Warning: this file has tabs in it.  It is required for GNU Make.
#
# This version may not produce optimal runs.  You'll need to experiment
# with the compiler options.
PI_VER=2.4

#Set this to your remove command.
RM=rm -f

#Set this to your copy command.
CP=cp

#Set this to your strip command.
STRIP=echo

MATH_LIBS=-lm

#Optimization flags. Choose the right ones for your system.

#For any gcc these options should be uncommented
#CC=gcc
#CFLAGS=-Wall -O6 -fno-exceptions -freduce-all-givs -fmove-all-movables -fomit-frame-pointer -fexpensive-optimizations -ffast-math -funroll-loops
LDFLAGS=-static -o
LD=${CC}

#For gcc on a 486 these options should be uncommented
#CFLAGS+=-malign-double -march=i486

#For gcc on a 586 these options should be uncommented
#CFLAGS+=-malign-double -mwide-multiply -march=i486

#For gcc on a 686 these options should be uncommented
#CFLAGS+=-malign-double -mwide-multiply -march=i686

#For gcc on a k6 these options should be uncommented
#CFLAGS+=-malign-double -march=k6

#For gcc on 680X0 these options should be uncommented
#CFLAGS+=-m68020 -m68881 

#For gcc on PowerPC these options should be uncommented
#CFLAGS+=-mpowerpc

#For gcc on PowerPC with libmoto these options should be uncommented
#MATH_LIBS=-lmoto -lm

#For gcc on SPARC these options should be uncommented
#CFLAGS+=-mcypress

#For gcc on Alpha these options should be uncommented
#CFLAGS+=-mno-soft-float



#For Compaq cc on Alpha these options should be uncommented
#CC=ccc
#LD=${CC}
#CFLAGS=-fast
#LDFLAGS=-non_shared -o

#For BeOS Intel these options should be uncommented
#CC=cc
#CFLAGS=-mpentium -Wall -freduce-all-givs -fmove-all-movables -fomit-frame-pointer -O3 -fno-strict-aliasing -fexpensive-optimizations -ffast-math -funroll-loops
#MATH_LIBS=
#LD=${CC}
#LDFLAGS=-static -s -o

#For BeOS PowerPC these options should be uncommented
#CC=mwcc
#CFLAGS=-O4
#LD=mwld
#LDFLAGS=-o
#MATH_LIBS=
#STRIP=echo

#For BeOS PowerPC with libmoto these options should be uncommented
#MATH_LIBS=-lmoto

#For HP-UX with HP cc this option should be uncommented
#CC=cc
#LD=${CC}
#CFLAGS=-Ae +O4 +Oaggressive +Ovectorize
#LDFLAGS=-Wl,-a, archive -Ae

#For PA-1.1 processors these options should be uncommented
#CFLAGS+= +DA1.1
#LDFLAGS+= +DA1.1 -o

#For PA-2.0 processors these options should be uncommented
#CFLAGS+= +DA2.0W
#LDFLAGS+= +DA2.0W -o

#For Tru64 with Compaq cc these options should be uncommented
#CC=cc
#LD=${CC}
#CFLAGS=-fast
#LDFLAGS=-non_shared -o

#For IRIX with SGI cc these option should be uncommented
#CC=cc
#LD=${CC}
#CFLAGS=-n32 -O3 -Ofast
#LDFLAGS=-n32 -IPA -Bstatic -o

#For Solaris with Sun cc these options should be uncommented
#CC=cc
#LD=${CC}
#CFLAGS=-xcrossfile=2 -fsimple=2 -fast -xO5
#LDFLAGS=-B static -o

#For MacOS X Server the following should be uncommented
#CC=cc
#CFLAGS=-Wall -freduce-all-givs -fmove-all-movables -fomit-frame-pointer -O3 -fno-strict-aliasing -fexpensive-optimizations -ffast-math -funroll-loops
#LD=${CC}
#LDFLAGS=-o
#MATH_LIBS=


#Leave items below this point untouched!

#Directories.
OUT_FILE=pi-agm
PREFIX=/usr/local

HDR_DIR=common

FFT_DIR=fftstuff

EXAMPLE_DIR=${FFT_DIR}/example
EXAMPLE_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${EXAMPLE_DIR}

FFTW_DIR=${FFT_DIR}/fftw
FFTW_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${FFTW_DIR}

JSPLIT_DIR=${FFT_DIR}/jsplit
JSPLIT_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${JSPLIT_DIR}

QUAD2_DIR=${FFT_DIR}/quad2
QUAD2_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${QUAD2_DIR}

HARTLEY_DIR=${FFT_DIR}/hartley
HARTLEY_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${HARTLEY_DIR}

OOURA_DIR=${FFT_DIR}/ooura
OOURA_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${OOURA_DIR}

OOURA2_DIR=${FFT_DIR}/ooura2
OOURA2_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${OOURA2_DIR}

OOURA3_DIR=${FFT_DIR}/ooura3
OOURA3_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${OOURA3_DIR}

NRC_DIR=${FFT_DIR}/nrc
NRC_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${NRC_DIR}


NTT32_DIR=ntt32

NTT486_DIR=${NTT32_DIR}/gcc486
NTT486_INC=-I ${HDR_DIR} -I ${NTT32_DIR} -I ${NTT486_DIR}

NTT586_DIR=${NTT32_DIR}/gcc586
NTT586_INC=-I ${HDR_DIR} -I ${NTT32_DIR} -I ${NTT586_DIR}

NTTLONG_DIR=${NTT32_DIR}/longlong
NTTLONG_INC=-I ${HDR_DIR} -I ${NTT32_DIR} -I ${NTTLONG_DIR}

NTTGEN_DIR=${NTT32_DIR}/generic
NTTGEN_INC=-I ${HDR_DIR} -I ${NTT32_DIR} -I ${NTTGEN_DIR}

NTTC586_DIR=${NTT32_DIR}/c586
NTTC586_INC=-I ${HDR_DIR} -I ${NTT32_DIR} -I ${NTTC586_DIR}

NTTAPF586_DIR=${NTT32_DIR}/apf586
NTTAPF586_INC=-I ${HDR_DIR} -I ${NTT32_DIR} -I ${NTTAPF586_DIR}

all:
	@echo "Pi AGM version ${PI_VER}:"
	@echo "The fft-hartley version is fastest in general but are limited to 8m digits."
	@echo "The fft-quad2 versions can be faster than the fft-hartley in a few situations."
	@echo "Usually the fft-hartley versions are the fastest though"
	@echo "The ntt versions are much slower but should work up to 1g digits!"
	@echo "Please use one of the following commands to build"
	@echo ""
	@echo "FFT versions are designed for faster runs but are limited to 8m digits."
	@echo "fft-example  -use your own custom made fft"
	@echo "fft-fftw     -use the fftw fft"
	@echo "fft-hartley  -use the hartley fft (usually the best)"
	@echo "fft-jspit    -use jasonp's split fft"
	@echo "fft-nrc      -use the numeric recipe style fft"
	@echo "fft-ooura    -use Takuyo Ooura's original fft"
	@echo "fft-ooura2   -use Takuyo Ooura's improved fft"
	@echo "fft-ooura3   -use Takuyo Ooura's even more improved fft"
	@echo "fft-quad2    -use the quad2 style fft (occasionally the best)"
	@echo ""
	@echo "NTT versions are designed for longer runs (up to 1g digits) but are slower."
	@echo "ntt-gen      -use the generic C ntt"
	@echo "ntt-long     -use a C ntt that requires a long long"
	@echo "ntt-c586     -use a C ntt that requires an 80bit FPU.  Little endian."
	@echo "ntt-apf586   -use ntt from Mikko Tommila's apfloat WON'T WORK WITH GCC-2.95.*"
	@echo "ntt-486      -use an assembly ntt that requires 486 & GCC"
	@echo "ntt-586      -use an assembly ntt that requires 586 & GCC"
	@echo ""

ntt-long:
	@echo "Compiling Pi AGM ${PI_VER} NTT LongLong Version..."
	${CC} ${NTTLONG_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${NTTLONG_INC} ${CFLAGS} -c ${NTT32_DIR}/*.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	${cleanup}

ntt-gen:
	@echo "Compiling Pi AGM ${PI_VER} NTT Generic Version..."
	${CC} ${NTTGEN_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${NTTGEN_INC} ${CFLAGS} -c ${NTT32_DIR}/*.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

ntt-c586:
	@echo "Compiling Pi AGM ${PI_VER} NTT c586 Version..."
	${CC} ${NTTC586_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${NTTC586_INC} ${CFLAGS} -c ${NTT32_DIR}/bigmul.c ${NTT32_DIR}/primes.c ${NTT32_DIR}/ntt.c ${NTTC586_DIR}/vector.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

ntt-apf586:
	@echo "Compiling Pi AGM ${PI_VER} NTT apf586 Version..."
	${CC} ${NTTAPF586_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${NTTAPF586_INC} ${CFLAGS} -c ${NTT32_DIR}/bigmul.c ${NTT32_DIR}/primes.c ${NTTAPF586_DIR}/ntt.c ${NTTAPF586_DIR}/vector.c ${NTTAPF586_DIR}/sixstep.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

ntt-486:
	@echo "Compiling Pi AGM ${PI_VER} NTT 486 Version..."
	${CC} ${NTT486_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${NTT486_INC} ${CFLAGS} -c ${NTT32_DIR}/*.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

ntt-586:
	@echo "Compiling Pi AGM ${PI_VER} NTT 586 Version..."
	${CC} ${NTT586_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${NTT586_INC} ${CFLAGS} -c ${NTT32_DIR}/bigmul.c ${NTT32_DIR}/primes.c ${NTT32_DIR}/ntt.c ${NTT586_DIR}/vector.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-example:
	@echo "Compiling Pi AGM ${PI_VER} FFT example Version..."
	${CC} ${EXAMPLE_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${EXAMPLE_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${EXAMPLE_INC} ${CFLAGS} -c ${EXAMPLE_DIR}/fft.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-fftw:
	@echo "Compiling Pi AGM ${PI_VER} FFT fftw Version..."
	${CC} ${FFTW_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${FFTW_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${FFTW_INC} ${CFLAGS} -c ${FFTW_DIR}/fft.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS} -lfftw
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-hartley:
	@echo "Compiling Pi AGM ${PI_VER} FFT hartley Version..."
	${CC} ${HARTLEY_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${HARTLEY_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${HARTLEY_INC} ${CFLAGS} -c ${HARTLEY_DIR}/fft.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-jsplit:
	@echo "Compiling Pi AGM ${PI_VER} FFT split Version..."
	${CC} ${JSPLIT_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${JSPLIT_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${JSPLIT_INC} ${CFLAGS} -c ${JSPLIT_DIR}/fft.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-ooura:
	@echo "Compiling Pi AGM ${PI_VER} FFT ooura Version..."
	${CC} ${OOURA_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${OOURA_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${OOURA_INC} ${CFLAGS} -c ${OOURA_DIR}/fft.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-ooura2:
	@echo "Compiling Pi AGM ${PI_VER} FFT ooura2 Version..."
	${CC} ${OOURA2_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${OOURA2_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${OOURA2_INC} ${CFLAGS} -c ${OOURA2_DIR}/fft.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-ooura3:
	@echo "Compiling Pi AGM ${PI_VER} FFT ooura3 Version..."
	${CC} ${OOURA3_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${OOURA3_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${OOURA3_INC} ${CFLAGS} -c ${OOURA3_DIR}/fft.c 
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-quad2:
	@echo "Compiling Pi AGM ${PI_VER} FFT quad2 Version..."
	${CC} ${QUAD2_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${QUAD2_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${QUAD2_INC} ${CFLAGS} -c ${QUAD2_DIR}/fft.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-nrc:
	@echo "Compiling Pi AGM ${PI_VER} FFT nrc Version..."
	${CC} ${NRC_INC} ${CFLAGS} -c ${HDR_DIR}/*.c
	${CC} ${NRC_INC} ${CFLAGS} -c ${FFT_DIR}/*.c
	${CC} ${NRC_INC} ${CFLAGS} -c ${NRC_DIR}/fft.c
	@echo "Linking files..."
	${LD} ${LDFLAGS} ${OUT_FILE} *.o ${MATH_LIBS}
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

clean:
	@echo "Cleaning up files generated by compiling Pi AGM ${PI_VER}"
	${RM} *.o
	${RM} pi.ini
	${RM} ${OUT_FILE}

install:
	@echo "Installing Pi AGM ${PI_VER} into ${BIN_DIR}."
	mkdir -p ${PREFIX}/share/pi-agm
	${CP} ${OUT_FILE} ${PREFIX}/share/pi-agm
	${CP} pi.ini ${PREFIX}/share/pi-agm
	ln -sf ${PREFIX}/share/pi-agm/pi-agm ${PREFIX}/bin/pi-agm

