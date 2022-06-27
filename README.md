# Easy353

## Introduction

Easy353 is a tool for recovering Angiosperms353 gene set(AGS), which can filter and de novo assemble reads from sequencing data based on databases, helping users capture AGS accurately and effectively.

## Instructions

It is recommended to download the latest binary release (Linux or OSX) there: https://github.com/GATB/minia/releases

Otherwise, Minia may be compiled from sources as follows:

```
# get a local copy of minia source code
git clone --recursive https://github.com/GATB/minia.git

# compile the code an run a simple test on your computer
cd minia
sh INSTALL
```

## Requirements

CMake 3.10+; see http://www.cmake.org/cmake/resources/software.html

C++11 compiler; (g++ version>=4.7 (Linux), clang version>=4.3 (Mac OSX))

## User manual

Type `minia` without any arguments for usage instructions.

A more complete manual is here: https://github.com/GATB/minia/raw/master/doc/manual.pdf

## What is new ? (2018)

Minia version 1 was implementing a rather unusual way to perform the assembly: traverse the graph and attempt to jump over errors and variants. This worked rather okay but not for e.g. repeated regions with many sequencing errors. Minia version 2 also followed the same philosophy, and had major improvements coming from the integration of the GATB library (mostly speed improvements) and cascading Bloom filter. Minia version 3 uses newer techniques and has virtually nothing in common with Minia 1: there is no Bloom filter anymore (the data structure is based on unitigs produced by the BCALM software). The assembly is performed using graph simplifications that are heavily inspired by the SPAdes assembler.

## Contact

Questions, suggestions, comments, etc? You can email the developer for specific support: zzhen0302@163.com
