Library: phyclust-c
Version: 0.1-2
Date: 2010-02-26
Author: Wei-Chen Chen
Maintainer: Wei-Chen Chen <phyclust@gmail.com>
Description: Model Based Phylogenetic Clustering
License: GPL (>= 2)
URL: http://thirteen-01.stat.iastate.edu/snoweye/phyclust/

Examples of main functions are in "test/" and short commands to compile
is in "test/make.test". A quick example is "test_toy.c" and one can
type the following,

> gcc -std=gnu99 -O3 -Wall -o test_toy test_toy.c ../*.c -I../ -lm
> test_toy Tt.009.anc.015.phy 4 0 0 3 2 0 2 3 10 0

:Wei-Chen Chen
