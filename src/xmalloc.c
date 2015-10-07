# include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include "xmalloc.h"

void *xcalloc (size_t m, size_t n)
{
   void *ptr;

   if (!(ptr = calloc (m, n))) {
#ifdef DEBUG
      fputs ("Error in calloc()\n", stderr);
      exit (EXIT_FAILURE);
#endif
   }

   return ptr;
}

void *xmalloc (size_t n)
{
   void *ptr;

   if (!(ptr = malloc (n))) {
#ifdef DEBUG
      fputs ("Error in malloc()\n", stderr);
      exit (EXIT_FAILURE);
#endif
   }

   return ptr;
}

void *xrealloc (void *ptr, size_t n)
{
   if (!(ptr = realloc (ptr, n))) {
#ifdef DEBUG
      fputs ("Error in realloc()\n", stderr);
      exit (EXIT_FAILURE);
#endif
   }

   return ptr;
}

void xfree (void *ptr)
{
   free (ptr);
}
