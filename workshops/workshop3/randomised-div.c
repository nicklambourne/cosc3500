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

  int count = 0;
  double x, y;
  long iseedlong = (long int) iseed;
  double sum = 0;

  struct drand48_data drandbuf;

  printf("method2 using drand48()");

  #pragma omp parallel private(x, y) shared(iseed, iseedlong, count)
  {
      srand48_r(iseedlong * omp_get_thread_num(), &drandbuf);
      #pragma omp for reduction(+:sum,count)
      for (jloops = 0; jloops < nloops; jloops++) {
        drand48_r(&drandbuf, &x);
        x *= 2;
        if (x > 0.99 && x < 1.01) {
            continue;
        }
        y = exp(-x)/(x-1);
        sum += y; //fabs(y);
        count++;
      }
  }

  double average = sum / count;
  double integral = average * 2; 

  printf("integral exp(-x) = %f\n", integral);

  exit(EXIT_SUCCESS);
}
