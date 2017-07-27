# This is a makefile template


# Compiler
# --------
CC=g++
#CC=icpc

include files.mk


# Compiler flags
# -------------------------------------------------------------------------
CFLAGS=-g -O3 

# Linker flags
# ------------
LDFLAGS= -Wl,--no-as-needed -lcppunit 

INCLUDES=

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Unittesting

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 

clean:
	rm $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
 
