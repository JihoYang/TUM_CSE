###########################
#                                                                  #
#             Makefile for FVM Code                   #
#                                                                  #
###########################

CC = g++
CFLAGS = -g -O3
LIBS = -lm

all: springDamper 


SRCG = damper.cpp helper.c init.c
OBJG =  $(addsuffix .o, $(basename $(SRCG)))


%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

springDamper: $(OBJG) 
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
	
clean:
	rm -rf *.o springDamper
