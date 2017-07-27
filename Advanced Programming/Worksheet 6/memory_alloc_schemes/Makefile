# This is a makefile template


# Compiler
# --------
CC=g++
#CC= g++
#CC=icpc

SOURCES=memory_alloc_schemes.cpp


# Compiler flags
# -------------------------------------------------------------------------
CFLAGS = -O0 #-fno-stack-protector
#CFLAGS=-g -O3 
#CFLAGS=-fopenmp -g -O3

# Linker flags
# ------------
LDFLAGS= 

INCLUDES= 

OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=memory_alloc_schemes

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 

clean:
	rm $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

