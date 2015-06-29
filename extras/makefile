# Makefile for Pi AGM 2.3, Carey Bloodworth's
# premiere Pi calculation program for *nix or DOS/gcc.
# Warning: this file has tabs in it.  It is required for GNU Make.
#
# This version may not produce optimal runs.  You'll need to experiment
# with the compiler options.
PI_VER=2.3

#Set this to your remove command.
RM=rm -f

#Set this to your copy command.
CP=cp

#Set this to your strip command.
STRIP=strip

#For any gcc these options should be uncommented
CC=gcc
CC_OPTS=-Wall -O3 -ffast-math -funroll-loops -c
LD_OPTS=-static -s -o

#For gcc on a 486 these options should be uncommented
#CC_OPTS=-m486 -Wall -O3 -ffast-math -funroll-loops -c

#For gcc on a 586 these options should be uncommented
#CC_OPTS=-mpentium -Wall -O3 -ffast-math -funroll-loops -c

#For HP-UX with HP cc these options should be uncommented
#CC=cc
#CC_OPTS=-c -Ae +O3 +Oaggressive +Ovectorize
#LD_OPTS=-Wl,-a, archive -Ae -o

#For IRIX with SGI cc these options should be uncommented
#CC=cc
#CC_OPTS=-64 -O3 -Ofast -c
#LD_OPTS=-64 -IPA -Bstatic -o

#For Solaris with Sun cc these options should be uncommented
#CC=cc
#CC_OPTS=-fast -dalign -xO4 -c
#LD_OPTS=-Bstatic -o

#Leave items below this point untouched

OUT_FILE=piagm23
BIN_DIR=/usr/local/bin
HDR_DIR=common

FFT_DIR=fftstuff

EXAMPLE_DIR=${FFT_DIR}/example
EXAMPLE_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${EXAMPLE_DIR}

FFTW_DIR=${FFT_DIR}/fftw
FFTW_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${FFTW_DIR}

QUAD2_DIR=${FFT_DIR}/quad2
QUAD2_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${QUAD2_DIR}

HARTLEY_DIR=${FFT_DIR}/hartley
HARTLEY_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${HARTLEY_DIR}

OOURA_DIR=${FFT_DIR}/ooura
OOURA_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${OOURA_DIR}

OOURA2_DIR=${FFT_DIR}/ooura2
OOURA2_INC=-I ${HDR_DIR} -I ${FFT_DIR} -I ${OOURA2_DIR}

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
	@echo "fft-nrc      -use the numeric recipe style fft"
	@echo "fft-ooura    -use Takuyo Ooura's original fft"
	@echo "fft-ooura2   -use Takuyo Ooura's improved fft"
	@echo "fft-quad2    -use the quad2 style fft (occasionally the best)"
	@echo ""
	@echo "NTT versions are designed for longer runs (up to 1g digits) but are slower."
	@echo "ntt-gen      -use the generic C ntt"
	@echo "ntt-long     -use a C ntt that requires a long long"
	@echo "ntt-c586     -use a C ntt that requires an 80bit FPU.  Little endian."
	@echo "ntt-486      -use an assembly ntt that requires 486 & GCC"
	@echo "ntt-586      -use an assembly ntt that requires 586 & GCC"
	@echo ""


ntt-long:
	@echo "Compiling Pi AGM ${PI_VER} NTT LongLong Version..."
	${CC} ${NTTLONG_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${NTTLONG_INC} ${CC_OPTS} ${NTT32_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

ntt-gen:
	@echo "Compiling Pi AGM ${PI_VER} NTT Generic Version..."
	${CC} ${NTTGEN_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${NTTGEN_INC} ${CC_OPTS} ${NTT32_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

ntt-c586:
	@echo "Compiling Pi AGM ${PI_VER} NTT c586 Version..."
	${CC} ${NTTC586_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${NTTC586_INC} ${CC_OPTS} ${NTT32_DIR}/bigmul.c ${NTT32_DIR}/primes.c ${NTT32_DIR}/ntt.c ${NTTC586_DIR}/vector.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

ntt-486:
	@echo "Compiling Pi AGM ${PI_VER} NTT 486 Version..."
	${CC} ${NTT486_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${NTT486_INC} ${CC_OPTS} ${NTT32_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

ntt-586:
	@echo "Compiling Pi AGM ${PI_VER} NTT 586 Version..."
	${CC} ${NTT586_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${NTT586_INC} ${CC_OPTS} ${NTT32_DIR}/bigmul.c ${NTT32_DIR}/primes.c ${NTT32_DIR}/ntt.c ${NTT586_DIR}/vector.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-example:
	@echo "Compiling Pi AGM ${PI_VER} FFT example Version..."
	${CC} ${EXAMPLE_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${EXAMPLE_INC} ${CC_OPTS} ${FFT_DIR}/*.c
	${CC} ${EXAMPLE_INC} ${CC_OPTS} ${EXAMPLE_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-fftw:
	@echo "Compiling Pi AGM ${PI_VER} FFT fftw Version..."
	${CC} ${FFTW_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${FFTW_INC} ${CC_OPTS} ${FFT_DIR}/*.c
	${CC} ${FFTW_INC} ${CC_OPTS} ${FFTW_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm -lfftw
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-hartley:
	@echo "Compiling Pi AGM ${PI_VER} FFT hartley Version..."
	${CC} ${HARTLEY_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${HARTLEY_INC} ${CC_OPTS} ${FFT_DIR}/*.c
	${CC} ${HARTLEY_INC} ${CC_OPTS} ${HARTLEY_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-ooura:
	@echo "Compiling Pi AGM ${PI_VER} FFT ooura Version..."
	${CC} ${OOURA_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${OOURA_INC} ${CC_OPTS} ${FFT_DIR}/*.c
	${CC} ${OOURA_INC} ${CC_OPTS} ${OOURA_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-ooura2:
	@echo "Compiling Pi AGM ${PI_VER} FFT ooura2 Version..."
	${CC} ${OOURA2_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${OOURA2_INC} ${CC_OPTS} ${FFT_DIR}/*.c
	${CC} ${OOURA2_INC} ${CC_OPTS} ${OOURA2_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-quad2:
	@echo "Compiling Pi AGM ${PI_VER} FFT quad2 Version..."
	${CC} ${QUAD2_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${QUAD2_INC} ${CC_OPTS} ${FFT_DIR}/*.c
	${CC} ${QUAD2_INC} ${CC_OPTS} ${QUAD2_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
	@echo "Cleaning up..."
	${STRIP} ${OUT_FILE}
	${RM} *.o
	${CP} ${HDR_DIR}/pi.ini .

fft-nrc:
	@echo "Compiling Pi AGM ${PI_VER} FFT nrc Version..."
	${CC} ${NRC_INC} ${CC_OPTS} ${HDR_DIR}/*.c
	${CC} ${NRC_INC} ${CC_OPTS} ${FFT_DIR}/*.c
	${CC} ${NRC_INC} ${CC_OPTS} ${NRC_DIR}/*.c
	@echo "Linking files..."
	${CC} ${LD_OPTS} ${OUT_FILE} *.o -lm
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
	${CP} ${OUT_FILE} ${BIN_DIR}
	${CP} pi.ini ${BIN_DIR}


