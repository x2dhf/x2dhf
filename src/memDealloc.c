/*

   see http://www.ibiblio.org/pub/languages/fortran/ch2-5.html

*/


#include <stdio.h>
#include <stdlib.h>

void memdealloc_(long *a, long *m, long *n, long *b)
{
  free(a + *b - 1);
}

