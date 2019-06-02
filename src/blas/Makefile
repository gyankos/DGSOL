#
# BLAS directory.
#

DG_ARCH = linux
CC 	= gcc
FC 	 = gfortran
L_ARCH   = $(DG_ARCH)
LIB_NAME = d-$(L_ARCH).a

OPTFLAGS = -O
CFLAGS   =  $(OPTFLAGS) 
FFLAGS   =  $(OPTFLAGS)

AR       = ar -rcv

FILES = dasum.o daxpy.o dcopy.o ddot.o dgemm.o dgemv.o \
        dger.o dnrm2.o dscal.o dswap.o dsyrk.o dtrmm.o \
        dtrmv.o dtrsm.o dtrsv.o idamax.o lsame.o xerbla.o


CFLAGS  =  $(OPTFLAGS) 
FFLAGS  =  $(OPTFLAGS)

blas: $(FILES) 
	$(AR) $(LIB_NAME) $(FILES)  
	rm -f *.o

clean:
	- rm -f *.o
	- rm -f *.a

.c.o:
	$(CC) $(CFLAGS) -c $*.c
.f.o:
	$(FC) $(FFLAGS) -c $*.f


