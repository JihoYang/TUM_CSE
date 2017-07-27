#----------------------------------------------------------------------#
#																	   #
#	Makefile for 4th order backward differentiation formula    	   	   #
#																	   #
#	Written by: Jiho Yang (MEng)									   #
#			    M.Sc. candidate, Computational Science & Engineering   #
#               Technische Universitat Munchen                         #
#																	   #
#----------------------------------------------------------------------#

#Include files
SOURCES=helper.c main.c bdf.c visual.c

#Compiler
#--------
CC = gcc

CFLAGS = -Wall -Werror -std=gnu99 -g

#Linker flags
#------------
LDFLAGS= -lm

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=sim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@

