/*

   see http://www.ibiblio.org/pub/languages/fortran/ch2-5.html

*/


#include <stdio.h>
#include <stdlib.h>

void memalloc_(long *a, long *m, long *n, long *b)
{
  long     *new;
  new = (long *)calloc(*m, *n);
  if (new == NULL)
    {
      printf("\n memalloc: not enough memory \n");
      exit(1);
    }
  else
    *b = (long)(new - a) + 1;
}
  
