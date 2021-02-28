TAR1 = rfstack
TAR2 = zrfs
CC = gcc 
CFLAGS=-g -lm 
OBJS = jfour1.o jrealft.o  jtinvit.o jtridib.o  sacio.o \
multi_tap_rf.o sigstuff.o multitap.o  

$(TAR1):$(OBJS) rf_stack.o
	$(CC)  $^ -o $@ -lm -g
$(TAR2):$(OBJS) zrf.o
	$(CC)  $^ -o $@ -lm -g
$(OBJS):%.o:%.c 
	$(CC) -c -g $^ -o $@
rf_stack.o:rf_stack.c
	$(CC) -c -g $^ -o $@
zrf.o:zrf.c 
	$(CC) -c -g $^ -o $@
.PHONY : clean
all:$(TAR1) $(TAR2) cleanobj
clean:
	-rm *.o $(TAR1) $(TAR2)
cleanobj:
	-rm *.o 



