# hpc_finalproject
explicit.c/implicit.c is that using explicit/implicit euler method to slove problem;
Makefile is used to make explicit/implicit; command line is "make explicit/implicit";
subscript is the file to sub in the cluster. 
"-h" means the element length, "-dt" means the time of ever step;
"-maxit" means the number of step, "-restart" means that whether reading the restart file to restart program, if " -restart 1" means restart, else do not restart;
explicit.h5/implicit.h5 is the restart data file;
ex_valg.sh/im_valg.sh is the subscript with valgrind to check the memory leak.
