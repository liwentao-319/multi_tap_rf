#TAR1 = rfstack
TAR1 = rf_stack
#TAR2 = zrfstack_gcarc
CC = gcc 
CFLAGS=-g -lm 
OBJS = jfour1.o jrealft.o  jtinvit.o jtridib.o  sacio.o \
multi_tap_rf.o sigstuff.o multitap.o  

$(TAR1):$(OBJS) rf_stack.o
	$(CC)  $^ -o $@ -lm -g

$(OBJS):%.o:%.c 
	$(CC) -c -g $^ -o $@

rf_stack.o:rf_stack.c 


.PHONY : clean
all:$(TAR1)  cleanobj
clean:
	-rm  $(TAR1) 
cleanobj:
	-rm *.o 





