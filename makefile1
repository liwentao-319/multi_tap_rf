TAR = rfstack
CC = gcc 
CFLAGS=-g -lm 
OBJS = jfour1.o jrealft.o  jtinvit.o jtridib.o  sacio.o \
multi_tap_rf.o sigstuff.o multitap.o rf_stack.o 

all:$(TAR)
$(TAR):$(OBJS)
	$(CC)  $^ -o $@ -lm -g

$(OBJS):%.o:%.c 
	$(CC) -c -g $^ -o $@


.PHONY : clean
clean:
	-rm *.o $(TAR)




