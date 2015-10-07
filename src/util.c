#include <stdio.h>
#include <ctype.h>
#include <libio.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"

#ifdef NO_GU_BASEPAIRS
const int BPI[6][6] = {
   {6, 6, 6, 0, 6, 6},
   {6, 6, 1, 6, 6, 6},
   {6, 2, 6, 6, 6, 6},
   {3, 6, 6, 6, 6, 6},
   {6, 6, 6, 6, 6, 6},
   {6, 6, 6, 6, 6, 6}
};
#else
const int BPI[6][6] = {
   {6, 6, 6, 0, 6, 6},
   {6, 6, 1, 6, 6, 6},
   {6, 2, 6, 4, 6, 6},
   {3, 6, 5, 6, 6, 6},
   {6, 6, 6, 6, 6, 6},
   {6, 6, 6, 6, 6, 6}
};
#endif


int roundInt (double d)
{
   return (int) (d + .5);
}

void strcatc (char *str, char c)
{
   str[strlen (str) + 1] = 0;
   str[strlen (str)] = c;
}

char *filename (char *file)
{
   char *name;

   for (name = file; *file; ++file)
      if (*file == '/' || *file == '\\')
	 name = file + 1;

   return name;
}

unsigned char toNum (char c)
{
   c = toupper (c);
   switch (c) {
   case 'A':
   case '0':
      return 0;
   case 'C':
   case '1':
      return 1;
   case 'G':
   case '2':
      return 2;
   case 'T':
   case 'U':
   case '3':
      return 3;
   }
   return 4;
}

int seqcmp (unsigned char *seq1, unsigned char *seq2, int length)
{
   int i;

   for (i = 0; i < length; ++i)
      if (seq1[i] < seq2[i])
	 return -1;
      else if (seq1[i] > seq2[i])
	 return 1;
   return 0;
}


int min3 (int a, int b, int c)
{
   if (a <= b && a <= c)
      return a;
   if (b <= c)
      return b;
   return c;
}

int same (unsigned char *a, unsigned char *b, int len)
{
   int i;

   for (i = 1; i <= len; ++i)
      if (a[i] != b[i])
	 return 0;
   return 1;
}

void readOrDie (unsigned int num, const char *name, FILE * file, const char *format, ...)
{
   va_list arg;
   va_start (arg, format);
   if (vfscanf (file, format, arg) != num) {
#ifdef DEBUG
      fprintf (stderr, "Error: %s file is corrupt\n", name);
#endif
      exit (EXIT_FAILURE);
   }
   va_end (arg);
}
