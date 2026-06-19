# 2D Lattice Sieve Number Field Intersection Matcher

This repository contains an optimized implementation of the 
Franke-Kleinjung (FK) 2D Lattice Sieve integrated into PARI/GP 
via a custom C extension. The goal of the project is to efficiently 
search across a large pool of quadratic number fields (fundamental 
discriminants) to find pairs of fields that share the maximum number 
of identical, signed, smooth algebraic norms.

## File Structure

* README.md: This configuration and architectural guide.
* Makefile: Build script to automate compiling the shared object 
  and verification utilities.
* findpairs.gp.c: The performance-critical core. This is a 
  gp2c-transpiled file derived from the original GP logic, heavily 
  modified with a high-performance Franke-Kleinjung 2D lattice 
  sieve written in native C. It uses pthread multithreading via 
  PARI's parallel engine to process fields concurrently.
* findpairs.gp: A short GP wrapper script that configures runtime 
  defaults (stack size, thread memory allocations) and maps the 
  shared object help descriptors via addhelp().
* originalpari.gp: The baseline implementation written entirely 
  in pure PARI/GP script to discover common smooth norms.
* fk2d.cc: A standalone C++ implementation of the Franke-Kleinjung 
  2D lattice sieve algorithm used to verify correctness against 
  theoretical benchmarks.
* printfk2d.cc: A diagnostic variant of the standalone sieve that 
  prints randomized sieve hits and index mappings for debugging 
  and analysis.

## Prerequisites

To compile and execute this project, you must have the PARI/GP 
development libraries installed.

On Debian/Ubuntu-based systems run:
sudo apt-get install pari-gp libpari-dev gp2c

## Compilation

You can build the entire suite using the provided Makefile:

# Build both the PARI/GP C-extension and standalone utilities
make all

# Build only the PARI/GP shared object
make library

# Build only the verification tools
make tools

## How to Run

1. Execute the Main Sieve Search
Load the wrapper script inside the PARI/GP interactive shell:
```
$ gp
? \r findpairs.gp
```
Once loaded, use the built-in online help to check the parameters:
```
? ?findpairs
```
Execute a search across a range of fields:
```
? findpairs(1000, 1000, 300, 1000, 1000000, 1, 16, 1024)
```
Note, that you may have to increase thread memory with e.g.
```
? default(threadsize, "500M")
```
or main memory with e.g.
```
default(parisizemax,32*1024^3);
```

2. Standalone Verification Tools
To run the algorithmic verification binaries:
# Run correctness checks
./fk2d

# View sample raw sieve streams
./printfk2d

## Mathematical Summary of the Sieve
The sieve targets a region dynamically scaled by the input parameter A. 
For any given element (cx, cy) inside the lattice space, the true 
search region sweeps across:

[-A/2, A/2) x [1, 2A]

Norm calculation evaluates the algebraic resultant between the field's 
defining polynomial f(x) = ax^2 + bx + c and the coordinates' 
linear form L(x) = cy * x - cx. Intersections across distinct 
fields strictly match identical norms including signs (i.e., +N != -N).

