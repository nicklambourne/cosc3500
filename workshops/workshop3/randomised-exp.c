/*  Copyleft Bromley brom@physics.uq.edu.au August 2018
 *  Implementing a couple of the simplest random methods in c

 *  basic compile eg. with:
 *  gcc -Wall -Wextra randomised-simple.c -o randomised-simple.exe -lm
 *  then test usage eg.:
 *  ./randomised-simple.exe 12345 10
 *  
 *  older gnu compilers might need -std=gnu99 for drand48 support
 *
 *  ************80 char length page width best for a2ps ***********************
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int main(int argc, char *argv[]) {
  unsigned int iseed;  /* input 1 random number seed */
  int nloops;          /* input 2 number of loops */

  int jloops;

  // int irandcur;
  // double drandcur;

  if (argc != 3) {
    fprintf(stderr, "Usage: %s <iseed> <nloops>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  iseed = atoi(argv[1]);
  nloops = atoi(argv[2]);

  printf("Running %s with iseed=%u and nloops=%d \n", argv[0],iseed,nloops);

  
  double x;
  long iseedlong = (long int) iseed;
  double sum = 0;

  struct drand48_data drandbuf;

  printf("method2 using drand48()");

  #pragma omp parallel private(x) shared(iseed, iseedlong)
  {
      srand48_r(iseedlong * omp_get_thread_num(), &drandbuf);
      #pragma omp for reduction(+:sum)
      for (jloops = 0; jloops < nloops; jloops++) {
        x = drand48() * 1000;
        double y = exp(-x);
        sum += y;
      }
  }

  double average = sum / nloops;
  double integral = average * 1000; 

  printf("integral exp(-x) = %f\n", integral);

  exit(EXIT_SUCCESS);
}
