# Compilers & Flags
CC       = cc
CXX      = g++
CFLAGS   = -g -O3 -Wall -fno-strict-aliasing -ffp-contract=off \
           -fPIC -I"/usr/include/x86_64-linux-gnu"
LDFLAGS  = -shared -g -O3 -Wall -fno-strict-aliasing -ffp-contract=off \
           -fPIC -Wl,-shared -Wl,-z,relro -lc -lm \
           -L/usr/lib/x86_64-linux-gnu -lpari
CXXFLAGS = -g -O0 -Wall

# Targets
LIBRARY  = findpairs.gp.so
TOOLS    = fk2d printfk2d

.PHONY: all library tools clean

all: library tools

library: $(LIBRARY)

tools: $(TOOLS)

# Build the shared object library for PARI/GP
$(LIBRARY): findpairs.gp.o
	$(CC) -o $@ $(LDFLAGS) findpairs.gp.o

findpairs.gp.o: findpairs.gp.c
	$(CC) -c -o $@ $(CFLAGS) findpairs.gp.c

# Build standalone verification utilities
fk2d: fk2d.cc
	$(CXX) $(CXXFLAGS) $< -o $@

printfk2d: printfk2d.cc
	$(CXX) $(CXXFLAGS) $< -o $@

# Clean workspace
clean:
	rm -f *.o *.so $(TOOLS)
