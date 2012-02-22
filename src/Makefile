NDIM=5000
NILEVEL=1
CC = gcc
COPT = -O3 -Wall 
INCDIR = 
LIBDIR = 
RANLIB = echo
LIBS = libcrn.a -lm
CFLAGS = $(COPT) -DNDIM=$(NDIM) -DNI_LEVEL=$(NILEVEL) \
	$(INCDIR) $(LIBDIR)
LFLAGS = -framework vecLib
LIBCRN = libcrn.a
OBJS = eprintf.o blast.o xesn.o chkaccu.o sov.o secdp.o
PROGS = xpredm$(NDIM) lpredm
INSTALL_DIR = ../bin

all: $(LIBCRN) $(PROGS)
install: all
	mkdir -p $(INSTALL_DIR)
	cp $(PROGS) ../bin
xpredm$(NDIM): xpredm.o xpredmsub.o $(LIBCRN)
	$(CC) $(CFLAGS) xpredm.o xpredmsub.o $(LIBS) -o xpredm$(NDIM)
lpredm: lpredm.o xpredmsub.o $(LIBCRN)
	$(CC) $(CFLAGS) lpredm.o xpredmsub.o $(LIBS) -o lpredm
xpredm.o: xpredm.c xpredm.h $(LIBCRN)
	$(CC) -c $(CFLAGS) xpredm.c
lpredm.o: lpredm.c xpredm.h $(LIBCRN)
	$(CC) -c $(CFLAGS) lpredm.c
xpredmsub.o: xpredmsub.c xpredm.h $(LIBCRN)
	$(CC) -c $(CFLAGS) xpredmsub.c
libcrn.a: $(OBJS)
	ar r libcrn.a $(OBJS)
	$(RANLIB) libcrn.a 
blast.o: blast.h blast.c  
	$(CC) -c $(CFLAGS) blast.c
eprintf.o: eprintf.h eprintf.c  
	$(CC) -c $(CFLAGS) eprintf.c
xesn.o: xesn.h blast.h xesn.c
	$(CC) -c $(CFLAGS) xesn.c
chkaccu.o: chkaccu.h blast.h chkaccu.c
	$(CC) -c $(CFLAGS) chkaccu.c
sov.o: sov.h blast.h sov.c
	$(CC) -c $(CFLAGS) sov.c
secdp.o: secdp.c blast.h
	$(CC) -c $(CFLAGS) secdp.c
clean:
	-rm -f $(PROGS) *.o libcrn.a
realclean: clean
	(cd $(INSTALL_DIR); rm -f $(PROGS))
