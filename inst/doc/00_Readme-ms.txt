ms version May 29, 2007.

This directory contains "ms" source codes partially.
See the other "readme" in this directory for more details.
The pdf file in "ms.tar" (download from author's webpage) is also
useful to understand the command line options.

I turn "ms" in a library format by adding and modifing some functions, "R_*.c",
for directorly using "msdir" in R, and the standardlone version will not
be guaranteed here. Try the original "ms" for standardlone.

I modify the source codes in order to avoid enormous warnings from "gcc" and
stuck by pointer missed reallociation including output results to file
instead of stdout, and memory leaks and release.
I also modify the random seed and random number generator and replace them by
the functions of R.

Wei-Chen Chen
