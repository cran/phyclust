seq-gen versin 1.3.2 - 7 Jan, 2005

This directory contains "seq-gen" source codes partially.
See the other "README" in this directory for more details.
The html file in "Seq-Gen.v.1.3.2" (download from author's webpage) is also
useful to understand the command line options.

I turn "seq-gen" in a library format by adding and modifing some functions,
"R_*.c", for directorly using "seq-gen" in R, and the standardlone version
will not be guaranteed here. Try the original "seq-gen" for standardlone.

I modify the source codes in order to avoid enormous warnings from "gcc".
I also modify the random seed and random number generator and replace them by
the functions of R.

Known bug:
There are some leaking memory for cases without inputing ancestor
state. But, it is OK since tree_fv is not closed and reassigned at somewhere.

==6256== 352 bytes in 1 blocks are still reachable in loss record 1 of 1
==6256==    at 0x402703E: malloc (vg_replace_malloc.c:207)
==6256==    by 0x40BE9CE: (within /lib/tls/i686/cmov/libc-2.9.so)
==6256==    by 0x40BEA9B: fopen (in /lib/tls/i686/cmov/libc-2.9.so)
==6256==    by 0x80514EB: ReadFileParams (seq-gen.c:676)
==6256==    by 0x8053057: seq_gen_main (seq-gen.c:936)
==6256==    by 0x8048E2B: main (test.c:20)


Wei-Chen Chen
