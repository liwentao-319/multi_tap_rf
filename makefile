CFL=-g
F77=f77 -v
.f.o: ; $(F77) $(CFL) -c $*.f
.c.o: ; ${CC}   -c $*.c

CC = gcc -g 



INCLUDE         =-I$(NRCINC)

DEST=.


ICEPACK=jtridib.o jtinvit.o 
SIGU=sigstuff.o jfour1.o jrealft.o 
TAP=multitap.o mult_tap_spec.o hires.o adwait.o ftest.o 

sendmt: sendtap.o $(TAP) $(ICEPACK) $(SIGU)
	$(CC)  -o $(DEST)/sendmt  sendtap.o $(TAP) $(SIGU) $(ICEPACK) -lm
jtap: jtap.o jfour1.o jrealft.o
	$(CC)  -o $(DEST)/jtap jtap.o jfour1.o jrealft.o -lm


DAT=`date +%Y.%m.%d`
TARFILE=MTM-2.2_$(DAT).tar

clean:
	rm -f *.o core sendmt jtap

all:
	make sendmt jtap
tar:
	tar cvf $(TARFILE) *.[ch] makefile njeff.dat
