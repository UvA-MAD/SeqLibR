# include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#include <string.h>
#include "energy.h"
#include "xmalloc.h"

#ifndef isinf
# define isinf(x) (!finite(x) && x == x)
#endif

/* functions to load energy rules
 * programs should include energy.h and link with energy.o
 */

const double R1 = .0019872;
const char BASES[5] = { 'A', 'C', 'G', 'U', 'N' };
const char BASE_PAIRS[6][4] = { "A-U", "C-G", "G-C", "U-A", "G-U", "U-G" };
unsigned char toNum (char);	/* in util.h */

#ifdef INTEGER
#define scale(d) (isinf(d) ? INFINITY : floor((d) * PRECISION + 0.5))
const ENERGY INFINITY = 999999;
#else
#define scale(d) ((d) * PRECISION)
const ENERGY INFINITY = 1.0 / 0.0;
#endif

double ion (int NA, int polymer, double naConc, double mgConc)
{
   if (NA == 0)
      return 0;
   else
      return -0.114 * log (naConc + 3.3 * sqrt (mgConc));
}

int min (int a, int b)
{
   return (a < b) ? a : b;
}

FILE *openFile (char *name, int NA, int Enthalpy)
{
   FILE *file;
   char *buffer, *fname;

   fname = xmalloc (strlen (name) + 5);
   strcpy (fname, name);
   strcat (fname, ".D");
   if (Enthalpy)
      strcat (fname, "H");
   else
      strcat (fname, "G");
   if (NA)
      strcat (fname, "D");

   file = fopen (fname, "rt");

   if (!file && getenv ("UNAFOLDDAT")) {
      buffer = xmalloc (strlen (getenv ("UNAFOLDDAT")) + strlen (fname) + 2);
      strcpy (buffer, getenv ("UNAFOLDDAT"));
      strcat (buffer, "/");
      strcat (buffer, fname);
      file = fopen (buffer, "rt");
      free (buffer);
   }

   if (!file) {
      buffer = xmalloc (strlen (PKGDATADIR) + strlen (fname) + 1);
      strcpy (buffer, PKGDATADIR);
      strcat (buffer, fname);
      if (!(file = fopen (buffer, "rt"))) {
	 perror (fname);
#ifdef DEBUG
	 exit (EXIT_FAILURE);
#endif
      }
      free (buffer);
   }
   free (fname);
   return file;
}

ENERGY combineEE (double energy, double enthalpy, double tRatio, double eCor)
{
   if (!isFinite (energy) || !isFinite (enthalpy))
      return INFINITY;
   else
      return scale (tRatio * (energy + eCor) + (1.0 - tRatio) * enthalpy);
}

void VcombineEE (double *energy, double *enthalpy, ENERGY * result, int tl, double tRatio, double eCor)
{
   int i;
   for (i = 0; i < tl; i++)
      result[i] = combineEE (energy[i], enthalpy[i], tRatio, eCor);

}


void loadStack (emat4 stackEnergies, emat4e stackEnthalpies, int NA)
{
   int i, j, ii, jj;
   FILE *gFile, *hFile;

   gFile = openFile ("stack", NA, 0);
   hFile = openFile ("stack", NA, 1);

   for (i = 0; i < 5; ++i)
      for (ii = 0; ii < 5; ++ii)
	 for (j = 0; j < 5; ++j)
	    for (jj = 0; jj < 5; ++jj)
	       if (i == 4 || j == 4 || ii == 4 || jj == 4)
		  stackEnthalpies[i][j][ii][jj] = INFINITY;
	       else {
		  readOrDie (1, "stack", gFile, "%lf", &stackEnergies[i][j][ii][jj]);
		  readOrDie (1, "stack", hFile, "%lf", &stackEnthalpies[i][j][ii][jj]);
	       }

   fclose (gFile);
   fclose (hFile);
}


void combineStack (emat4 stackEnergies, emat4e stackEnthalpies, double tRatio, double saltCorrection, ENERGY stack[5][5][5][5])
{
   int i, j, ii, jj;

   for (i = 0; i < 5; ++i)
      for (ii = 0; ii < 5; ++ii)
	 for (j = 0; j < 5; ++j)
	    for (jj = 0; jj < 5; ++jj)
	       if (i == 4 || j == 4 || ii == 4 || jj == 4)
		  stack[i][j][ii][jj] = INFINITY;
	       else
		  stack[i][j][ii][jj] =
		     combineEE (stackEnergies[i][j][ii][jj], stackEnthalpies[i][j][ii][jj], tRatio, saltCorrection);
}

void symmetryCheckStack (emat4 stack, char *which)
{
   int i, j, ii, jj;

   for (i = 0; i < 4; ++i)
      for (j = 0; j < 4; ++j)
	 for (ii = 0; ii < 4; ++ii)
	    for (jj = 0; jj < 4; ++jj)
	       if (stack[i][j][ii][jj] != stack[jj][ii][j][i])
		  fprintf (stderr,
			   "Warning: %c-%c/%c-%c stack %s is %g; %c-%c/%c-%c stack %s is %g\n",
			   BASES[i], BASES[j], BASES[ii], BASES[jj], which,
			   stack[i][j][ii][jj], BASES[jj], BASES[ii], BASES[j], BASES[i], which, stack[jj][ii][j][i]);
}

void loadDangle (emat3 dangleEnergies3, emat3e dangleEnthalpies3, emat3 dangleEnergies5, emat3e dangleEnthalpies5, int NA)
{
   int i, j, k;
   FILE *gFile, *hFile;

   gFile = openFile ("dangle", NA, 0);
   hFile = openFile ("dangle", NA, 1);
   for (i = 0; i < 5; ++i)
      for (j = 0; j < 5; ++j)
	 for (k = 0; k < 6; ++k)
	    if (i == 4 || j == 4)
	       dangleEnthalpies3[i][j][k] = INFINITY;
	    else if (k == 4)
	       dangleEnthalpies3[i][j][k] = INFINITY;
	    else if (k == 5)
	       dangleEnthalpies3[i][j][k] = INFINITY;
	    else {
	       readOrDie (1, "dangle", gFile, "%lf", &dangleEnergies3[i][j][k]);
	       readOrDie (1, "dangle", hFile, "%lf", &dangleEnthalpies3[i][j][k]);
	    }

   for (i = 0; i < 5; ++i)
      for (j = 0; j < 5; ++j)
	 for (k = 0; k < 6; ++k)
	    if (i == 4 || j == 4)
	       dangleEnthalpies5[i][j][k] = INFINITY;
	    else if (k == 4)
	       dangleEnthalpies5[i][j][k] = INFINITY;
	    else if (k == 5)
	       dangleEnthalpies5[i][j][k] = INFINITY;
	    else {
	       readOrDie (1, "dangle", gFile, "%lf", &dangleEnergies5[i][j][k]);
	       readOrDie (1, "dangle", hFile, "%lf", &dangleEnthalpies5[i][j][k]);
	    }

   fclose (gFile);
   fclose (hFile);
}

void combineDangle (emat3 dangleEnergies3, emat3 dangleEnergies5, emat3e dangleEnthalpies3, emat3e dangleEnthalpies5, double tRatio,
		    double saltCorrection, ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6])
{
   int i, j, k;

   for (i = 0; i < 5; ++i)
      for (j = 0; j < 5; ++j)
	 for (k = 0; k < 6; ++k) {
	    if (i == 4 || j == 4 || k == 4 || k == 5) {
	       dangle3[i][j][k] = INFINITY;
	       dangle5[i][j][k] = INFINITY;
	    }
	    else {
	       dangle3[i][j][k] = combineEE (dangleEnergies3[i][j][k], dangleEnthalpies3[i][j][k], tRatio, 0.5 * saltCorrection);
	       dangle5[i][j][k] = combineEE (dangleEnergies5[i][j][k], dangleEnthalpies5[i][j][k], tRatio, 0.5 * saltCorrection);
	    }
	 }
}

void loadLoop (double hairpinLoopEnergies[30], double interiorLoopEnergies[30], double bulgeLoopEnergies[30],
	       double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], int NA)
{
   int k;
   FILE *gFile, *hFile;

   gFile = openFile ("loop", NA, 0);
   hFile = openFile ("loop", NA, 1);

   for (k = 0; k < 30; ++k) {
      readOrDie (3, "loop", gFile, "%*f%lf%lf%lf", &interiorLoopEnergies[k], &bulgeLoopEnergies[k], &hairpinLoopEnergies[k]);
      readOrDie (3, "loop", hFile, "%*f%lf%lf%lf", &interiorLoopEnthalpies[k], &bulgeLoopEnthalpies[k], &hairpinLoopEnthalpies[k]);
   }

   fclose (gFile);
   fclose (hFile);
}

void combineLoop (double hairpinLoopEnergies[30], double interiorLoopEnergies[30], double bulgeLoopEnergies[30],
		  double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30],
		  double tRatio, double saltCorrection, ENERGY hairpinLoop[30], ENERGY interiorLoop[30], ENERGY bulgeLoop[30])
{
   int k;
   double sc;

   for (k = 0; k < 30; ++k) {
      sc = saltCorrection * (1.0 + 0.5 * min (k, 10));
      hairpinLoop[k] = combineEE (hairpinLoopEnergies[k], hairpinLoopEnthalpies[k], tRatio, 0.0);
      interiorLoop[k] = combineEE (interiorLoopEnergies[k], interiorLoopEnthalpies[k], tRatio, sc);
      bulgeLoop[k] = combineEE (bulgeLoopEnergies[k], bulgeLoopEnthalpies[k], tRatio, sc);
   }
}

void loadSint2 (emat_s2e sint2Energies, emat_s2g sint2Enthalpies, int NA)
{
   int b, c, i, j;
   FILE *gFile, *hFile;

   gFile = openFile ("sint2", NA, 0);
   hFile = openFile ("sint2", NA, 1);

   for (b = 0; b < 7; ++b)
      for (i = 0; i < 5; ++i)
	 for (c = 0; c < 7; ++c)
	    for (j = 0; j < 5; ++j)
	       if (b == 6 || c == 6)
		  sint2Enthalpies[b][c][i][j] = INFINITY;
	       else if (i == 4 || j == 4)
		  sint2Enthalpies[b][c][i][j] = 0;
	       else {
		  readOrDie (1, "sint2", gFile, "%lf", &sint2Energies[b][c][i][j]);
		  readOrDie (1, "sint2", hFile, "%lf", &sint2Enthalpies[b][c][i][j]);
	       }

   fclose (gFile);
   fclose (hFile);
}

void combineSint2 (emat_s2e sint2Energies, emat_s2g sint2Enthalpies, double tRatio, double saltCorrection, ENERGY sint2[7][7][5][5])
{
   int b, c, i, j;

   for (b = 0; b < 6; ++b)
      for (c = 0; c < 6; ++c)
	 for (i = 0; i < 5; ++i)
	    for (j = 0; j < 5; ++j) {
	       if (b == 6 || c == 6)
		  sint2[b][c][i][j] = INFINITY;
	       else if (i == 4 || j == 4)
		  sint2[b][c][i][j] = 0;
	       else
		  sint2[b][c][i][j] =
		     combineEE (sint2Energies[b][c][i][j], sint2Enthalpies[b][c][i][j], tRatio, 2.0 * saltCorrection);
	    }
}

void symmetryCheckSint2 (emat_s2e sint2, char *which)
{
   int b, c, i, j;

   for (b = 0; b < 6; ++b)
      for (c = 0; c < 6; ++c)
	 for (i = 0; i < 4; ++i)
	    for (j = 0; j < 4; ++j) {
	       int bb, cc;
	       bb = (b > 3) ? 9 - b : 3 - b;
	       cc = (c > 3) ? 9 - c : 3 - c;

	       if (sint2[b][c][i][j] != sint2[cc][bb][j][i])
		  fprintf (stderr,
			   "Warning: %s/%s/%c/%c sint2 %s is %g; %s/%s/%c/%c sint2 %s is %g\n",
			   BASE_PAIRS[b], BASE_PAIRS[c], BASES[i], BASES[j],
			   which, sint2[b][c][i][j], BASE_PAIRS[cc], BASE_PAIRS[bb], BASES[j], BASES[i], which,
			   sint2[cc][bb][j][i]);
	    }
}

void loadAsint1x2 (emat_as1e asint1x2Energies, emat_as1g asint1x2Enthalpies, int NA)
{
   int b, c, i, j, k;
   FILE *gFile, *hFile;

   gFile = openFile ("asint1x2", NA, 0);
   hFile = openFile ("asint1x2", NA, 1);

   for (b = 0; b < 7; ++b)
      for (k = 0; k < 5; ++k)
	 for (i = 0; i < 5; ++i)
	    for (c = 0; c < 7; ++c)
	       for (j = 0; j < 5; ++j)
		  if (b == 6 || c == 6)
		     asint1x2Enthalpies[b][c][i][j][k] = INFINITY;
		  else if (i == 4 || j == 4 || k == 4)
		     asint1x2Enthalpies[b][c][i][j][k] = 0;
		  else {
		     readOrDie (1, "asint1x2", gFile, "%lf", &asint1x2Energies[b][c][i][j][k]);
		     readOrDie (1, "asint1x2", hFile, "%lf", &asint1x2Enthalpies[b][c][i][j][k]);
		  }

   fclose (gFile);
   fclose (hFile);
}

void combineAsint1x2 (emat_as1e asint1x2Energies, emat_as1g asint1x2Enthalpies, double tRatio,
		      double saltCorrection, ENERGY asint1x2[7][7][5][5][5])
{
   int b, c, i, j, k;

   for (b = 0; b < 7; ++b)
      for (c = 0; c < 7; ++c)
	 for (i = 0; i < 5; ++i)
	    for (j = 0; j < 5; ++j)
	       for (k = 0; k < 5; ++k) {
		  if (b == 6 || c == 6)
		     asint1x2[b][c][i][j][k] = INFINITY;
		  else if (i == 4 || j == 4 || k == 4)
		     asint1x2[b][c][i][j][k] = 0;
		  else
		     asint1x2[b][c][i][j][k] =
			combineEE (asint1x2Energies[b][c][i][j][k], asint1x2Enthalpies[b][c][i][j][k], tRatio,
				   2.5 * saltCorrection);
	       }
}

void loadSint4 (emat_s4e sint4Energies, emat_s4g sint4Enthalpies, int NA)
{
   int b, c, i1, j1, i2, j2;
   FILE *gFile, *hFile;

   gFile = openFile ("sint4", NA, 0);
   hFile = openFile ("sint4", NA, 1);

   for (b = 0; b < 7; ++b)
      for (c = 0; c < 7; ++c)
	 for (i1 = 0; i1 < 5; ++i1)
	    for (j1 = 0; j1 < 5; ++j1)
	       for (i2 = 0; i2 < 5; ++i2)
		  for (j2 = 0; j2 < 5; ++j2)
		     if (b == 6 || c == 6)
			sint4Enthalpies[b][c][i1][j1][i2][j2] = INFINITY;
		     else if (i1 == 4 || j1 == 4 || i2 == 4 || j2 == 4)
			sint4Enthalpies[b][c][i1][j1][i2][j2] = 0;
		     else {
			readOrDie (1, "sint4", gFile, "%lf", &sint4Energies[b][c][i1][j1][i2][j2]);
			readOrDie (1, "sint4", hFile, "%lf", &sint4Enthalpies[b][c][i1][j1][i2][j2]);
		     }

   fclose (gFile);
   fclose (hFile);
}

void
combineSint4 (emat_s4e sint4Energies, emat_s4g sint4Enthalpies, double tRatio, double saltCorrection,
	      ENERGY sint4[7][7][5][5][5][5])
{
   int b, c, i1, j1, i2, j2;

   for (b = 0; b < 7; ++b)
      for (c = 0; c < 7; ++c)
	 for (i1 = 0; i1 < 5; ++i1)
	    for (j1 = 0; j1 < 5; ++j1)
	       for (i2 = 0; i2 < 5; ++i2)
		  for (j2 = 0; j2 < 5; ++j2) {
		     if (b == 6 || c == 6)
			sint4[b][c][i1][j1][i2][j2] = INFINITY;
		     else if (i1 == 4 || j1 == 4 || i2 == 4 || j2 == 4)
			sint4[b][c][i1][j1][i2][j2] = 0;
		     else
			sint4[b][c][i1][j1][i2][j2] =
			   combineEE (sint4Energies[b][c][i1][j1][i2][j2],
				      sint4Enthalpies[b][c][i1][j1][i2][j2], tRatio, 3.0 * saltCorrection);
		  }
}

void symmetryCheckSint4 (emat_s4e sint4, char *which)
{
   int b, c, i1, j1, i2, j2;

   for (b = 0; b < 6; ++b)
      for (c = 0; c < 6; ++c)
	 for (i1 = 0; i1 < 4; ++i1)
	    for (j1 = 0; j1 < 4; ++j1)
	       for (i2 = 0; i2 < 4; ++i2)
		  for (j2 = 0; j2 < 4; ++j2) {
		     int bb, cc;
		     bb = (b > 3) ? 9 - b : 3 - b;
		     cc = (c > 3) ? 9 - c : 3 - c;

		     if (sint4[b][c][i1][j1][i2][j2] != sint4[cc][bb][j2][i2][j1][i1])
			fprintf (stderr,
				 "Warning: %s/%s/%c/%c/%c/%c sint4 %s is %g; %s/%s/%c/%c/%c/%c sint4 %s is %g\n",
				 BASE_PAIRS[b], BASE_PAIRS[c], BASES[i1], BASES[j1],
				 BASES[i2], BASES[j2], which,
				 sint4[b][c][i1][j1][i2][j2], BASE_PAIRS[cc],
				 BASE_PAIRS[bb], BASES[j2], BASES[i2], BASES[j1], BASES[i1], which, sint4[cc][bb][j2][i2][j1][i1]);
		  }
}

void loadTstack (emat4 tsEng, emat4e tsEnth, int NA, char *fn)
{
   int i, j, ii, jj;
   FILE *gFile, *hFile;
   gFile = openFile (fn, NA, 0);
   hFile = openFile (fn, NA, 1);
   for (i = 0; i < 5; ++i)
      for (ii = 0; ii < 5; ++ii)
	 for (j = 0; j < 5; ++j)
	    for (jj = 0; jj < 5; ++jj)
	       if (i == 4 || j == 4)
		  tsEnth[i][j][ii][jj] = INFINITY;
	       else if (ii == 4 || jj == 4)
		  tsEnth[i][j][ii][jj] = 0;
	       else {
		  readOrDie (1, fn, gFile, "%lf", &tsEng[i][j][ii][jj]);
		  readOrDie (1, fn, hFile, "%lf", &tsEnth[i][j][ii][jj]);
	       }
   fclose (gFile);
   fclose (hFile);
}

void loadTstack2 (emat4 tsEng, emat4s tsEnth, int NA, char *fn)
{
   int i1, j1, i2, j2;
   FILE *gFile, *hFile;

   gFile = openFile (fn, NA, 0);
   hFile = openFile (fn, NA, 1);

   for (i1 = 0; i1 < 5; ++i1)
      for (i2 = 0; i2 < 6; ++i2)
	 for (j1 = 0; j1 < 5; ++j1)
	    for (j2 = 0; j2 < 6; ++j2)
	       if (i1 == 4 || j1 == 4)
		  tsEnth[i1][j1][i2][j2] = INFINITY;
	       else if (i2 == 5 || j2 == 5)
		  tsEnth[i1][j1][i2][j2] = INFINITY;
	       else if (i2 == 4 || j2 == 4)
		  tsEnth[i1][j1][i2][j2] = 0;
	       else {
		  readOrDie (1, "tstackm", gFile, "%lf", &tsEng[i1][j1][i2][j2]);
		  readOrDie (1, "tstackm", hFile, "%lf", &tsEnth[i1][j1][i2][j2]);
	       }

   fclose (gFile);
   fclose (hFile);
}

void combineTstack1 (emat4 tstackEnergies, emat4e tstackEnthalpies, double tRatio, ENERGY tstack[5][5][5][5])
{
   int i1, j1, i2, j2;

   for (i1 = 0; i1 < 5; ++i1)
      for (j1 = 0; j1 < 5; ++j1)
	 for (i2 = 0; i2 < 5; ++i2)
	    for (j2 = 0; j2 < 5; ++j2) {
	       if (i1 == 4 || j1 == 4)
		  tstack[i1][j1][i2][j2] = INFINITY;
	       else if (i2 == 4 || j2 == 4)
		  tstack[i1][j1][i2][j2] = 0;
	       else if (!isFinite (tstackEnergies[i1][j1][i2][j2]) || !isFinite (tstackEnthalpies[i1][j1][i2][j2]))
		  tstack[i1][j1][i2][j2] = INFINITY;
	       else
		  tstack[i1][j1][i2][j2] =
		     scale (tRatio * tstackEnergies[i1][j1][i2][j2] + (1.0 - tRatio) * tstackEnthalpies[i1][j1][i2][j2]);
	    }
}

void combineTstack2 (emat4 tstackEnergies, emat4s tstackEnthalpies, double tRatio, double saltCorrection, ENERGY tstack[5][5][6][6])
{
   int i1, j1, i2, j2;

   for (i1 = 0; i1 < 5; ++i1)
      for (j1 = 0; j1 < 5; ++j1)
	 for (i2 = 0; i2 < 6; ++i2)
	    for (j2 = 0; j2 < 6; ++j2) {
	       if (i1 == 4 || j1 == 4)
		  tstack[i1][j1][i2][j2] = INFINITY;
	       else if (i2 == 5 || j2 == 5)
		  tstack[i1][j1][i2][j2] = INFINITY;
	       else if (i2 == 4 || j2 == 4)
		  tstack[i1][j1][i2][j2] = 0;
	       else
		  tstack[i1][j1][i2][j2] =
		     combineEE (tstackEnergies[i1][j1][i2][j2], tstackEnthalpies[i1][j1][i2][j2], tRatio, saltCorrection);
	    }
}

void loadMulti (double multiEnergies[3], double multiEnthalpies[3], int NA)
{
   FILE *gFile, *hFile;

   gFile = openFile ("miscloop", NA, 0);
   hFile = openFile ("miscloop", NA, 1);

   readOrDie (3, "miscloop", gFile, "%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnergies[0], &multiEnergies[1], &multiEnergies[2]);
   readOrDie (3, "miscloop", hFile, "%*g%*g%*g%*g%*g%*g%lg%lg%lg", &multiEnthalpies[0], &multiEnthalpies[1], &multiEnthalpies[2]);

   fclose (gFile);
   fclose (hFile);
}

void combineMulti (double multiEnergies[3], double multiEnthalpies[3], double tRatio, ENERGY multi[3])
{
   VcombineEE (multiEnergies, multiEnthalpies, multi, 3, tRatio, 0.0);
}

void loadMisc (double miscEnergies[13], double miscEnthalpies[13], int NA)
{
   FILE *gFile, *hFile;

   gFile = openFile ("miscloop", NA, 0);
   hFile = openFile ("miscloop", NA, 1);

   readOrDie (13, "miscloop", gFile,
	      "%lg%lg%lg%lg%lg%lg%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg%lg",
	      &miscEnergies[12], &miscEnergies[4], &miscEnergies[0],
	      &miscEnergies[1], &miscEnergies[2], &miscEnergies[3],
	      &miscEnergies[6], &miscEnergies[8], &miscEnergies[9],
	      &miscEnergies[10], &miscEnergies[11], &miscEnergies[5], &miscEnergies[7]);
   readOrDie (13, "miscloop", hFile,
	      "%lg%lg%lg%lg%lg%lg%*g%*g%*g%*g%*g%*g%lg%lg%lg%lg%lg%lg%lg",
	      &miscEnthalpies[12], &miscEnthalpies[4], &miscEnthalpies[0],
	      &miscEnthalpies[1], &miscEnthalpies[2], &miscEnthalpies[3],
	      &miscEnthalpies[6], &miscEnthalpies[8], &miscEnthalpies[9],
	      &miscEnthalpies[10], &miscEnthalpies[11], &miscEnthalpies[5], &miscEnthalpies[7]);

   fclose (gFile);
   fclose (hFile);
}

void combineMisc (double miscEnergies[13], double miscEnthalpies[13], double tRatio, ENERGY misc[13])
{
   VcombineEE (miscEnergies, miscEnthalpies, misc, 13, tRatio, 0.0);
   misc[7] = miscEnergies[7] == 1 || miscEnthalpies[7] == 1;
}

void makeAUPenalty (ENERGY misc[13], ENERGY aup[5][5], int isPF)
{
   int i, j;

   for (i = 0; i < 5; ++i)
      for (j = 0; j < 5; ++j)
	 aup[i][j] = isPF ? 1.0 : 0.0;

   aup[0][3] = aup[3][0] = aup[2][3] = aup[3][2] = misc[6];
}

void loadTriloop (struct triloopE **triloopEnergies, struct triloopE **triloopEnthalpies, int *num, int NA)
{
   FILE *gFile, *hFile;
   int i, size;
   double energy;

   gFile = openFile ("triloop", NA, 0);

   *num = 0;
   size = 16;
   *triloopEnergies = calloc (16, sizeof (struct triloopE));

   while (fscanf (gFile, "%5s %lg", (*triloopEnergies)[*num].loop, &energy) == 2) {
      for (i = 0; i < 5; ++i)
	 (*triloopEnergies)[*num].loop[i] = toNum ((*triloopEnergies)[*num].loop[i]);
      (*triloopEnergies)[*num].energy = energy;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *triloopEnergies = realloc (*triloopEnergies, size * sizeof (struct triloopE));
      }
   }

   *triloopEnergies = realloc (*triloopEnergies, *num * sizeof (struct triloopE));

   fclose (gFile);

   hFile = openFile ("triloop", NA, 1);

   *num = 0;
   size = 16;
   *triloopEnthalpies = calloc (16, sizeof (struct triloopE));

   while (fscanf (hFile, "%5s %lg", (*triloopEnthalpies)[*num].loop, &energy) == 2) {
      for (i = 0; i < 5; ++i)
	 (*triloopEnthalpies)[*num].loop[i] = toNum ((*triloopEnthalpies)[*num].loop[i]);
      (*triloopEnthalpies)[*num].energy = energy;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *triloopEnthalpies = realloc (*triloopEnthalpies, size * sizeof (struct triloopE));
      }
   }

   *triloopEnthalpies = realloc (*triloopEnthalpies, *num * sizeof (struct triloopE));
   fclose (hFile);
}

void
combineTriloop (const struct triloopE *triloopEnergies,
		const struct triloopE *triloopEnthalpies, double tRatio, struct triloop *triloop, int num)
{
   int i;

   for (i = 0; i < num; ++i) {
      memcpy (triloop[i].loop, triloopEnergies[i].loop, 5);
      triloop[i].energy = combineEE (triloopEnergies[i].energy, triloopEnthalpies[i].energy, tRatio, 0.0);
   }
}

int triloopcmp (const void *loop1, const void *loop2)
{
   int i;
   const unsigned char *h1 = loop1;
   const struct triloop *h2 = loop2;

   for (i = 0; i < 5; ++i)
      if (h1[i] < h2->loop[i])
	 return -1;
      else if (h1[i] > h2->loop[i])
	 return 1;
   return 0;
}

void loadTloop (struct tloopE **tloopEnergies, struct tloopE **tloopEnthalpies, int *num, int NA)
{
   FILE *gFile, *hFile;
   int i, size;
   double energy;

   gFile = openFile ("tloop", NA, 0);

   *num = 0;
   size = 16;
   *tloopEnergies = calloc (16, sizeof (struct tloopE));

   while (fscanf (gFile, "%6s %lg", (*tloopEnergies)[*num].loop, &energy) == 2) {
      for (i = 0; i < 6; ++i)
	 (*tloopEnergies)[*num].loop[i] = toNum ((*tloopEnergies)[*num].loop[i]);
      (*tloopEnergies)[*num].energy = energy;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *tloopEnergies = realloc (*tloopEnergies, size * sizeof (struct tloopE));
      }
   }

   *tloopEnergies = realloc (*tloopEnergies, *num * sizeof (struct tloopE));

   fclose (gFile);

   hFile = openFile ("tloop", NA, 1);

   *num = 0;
   size = 16;
   *tloopEnthalpies = calloc (16, sizeof (struct tloopE));

   while (fscanf (hFile, "%6s %lg", (*tloopEnthalpies)[*num].loop, &energy) == 2) {
      for (i = 0; i < 6; ++i)
	 (*tloopEnthalpies)[*num].loop[i] = toNum ((*tloopEnthalpies)[*num].loop[i]);
      (*tloopEnthalpies)[*num].energy = energy;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *tloopEnthalpies = realloc (*tloopEnthalpies, size * sizeof (struct tloopE));
      }
   }

   *tloopEnthalpies = realloc (*tloopEnthalpies, *num * sizeof (struct tloopE));

   fclose (hFile);
}

void
combineTloop (const struct tloopE *tloopEnergies, const struct tloopE *tloopEnthalpies, double tRatio, struct tloop *tloop, int num)
{
   int i;

   for (i = 0; i < num; ++i) {
      memcpy (tloop[i].loop, tloopEnergies[i].loop, 6);
      tloop[i].energy = combineEE (tloopEnergies[i].energy, tloopEnthalpies[i].energy, tRatio, 0.0);
   }
}

int tloopcmp (const void *loop1, const void *loop2)
{
   int i;
   const unsigned char *h1 = loop1;
   const struct tloop *h2 = loop2;

   for (i = 0; i < 6; ++i)
      if (h1[i] < h2->loop[i])
	 return -1;
      else if (h1[i] > h2->loop[i])
	 return 1;

   return 0;
}

void loadHexaloop (struct hexaloopE **hexaloopEnergies, struct hexaloopE **hexaloopEnthalpies, int *num, int NA)
{
   FILE *gFile, *hFile;
   int i, size;
   double energy;

   gFile = openFile ("hexaloop", NA, 0);
   *num = 0;
   size = 16;
   *hexaloopEnergies = calloc (16, sizeof (struct hexaloopE));

   while (fscanf (gFile, "%8s %lg", (*hexaloopEnergies)[*num].loop, &energy) == 2) {
      for (i = 0; i < 8; ++i)
	 (*hexaloopEnergies)[*num].loop[i] = toNum ((*hexaloopEnergies)[*num].loop[i]);
      (*hexaloopEnergies)[*num].energy = energy;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *hexaloopEnergies = realloc (*hexaloopEnergies, size * sizeof (struct hexaloopE));
      }
   }

   *hexaloopEnergies = realloc (*hexaloopEnergies, *num * sizeof (struct hexaloopE));

   fclose (gFile);

   hFile = openFile ("hexaloop", NA, 1);

   *num = 0;
   size = 16;
   *hexaloopEnthalpies = calloc (16, sizeof (struct hexaloopE));

   while (fscanf (hFile, "%8s %lg", (*hexaloopEnthalpies)[*num].loop, &energy) == 2) {
      for (i = 0; i < 8; ++i)
	 (*hexaloopEnthalpies)[*num].loop[i] = toNum ((*hexaloopEnthalpies)[*num].loop[i]);
      (*hexaloopEnthalpies)[*num].energy = energy;
      ++*num;
      if (*num == size) {
	 size *= 2;
	 *hexaloopEnthalpies = realloc (*hexaloopEnthalpies, size * sizeof (struct hexaloopE));
      }
   }

   *hexaloopEnthalpies = realloc (*hexaloopEnthalpies, *num * sizeof (struct hexaloopE));

   fclose (hFile);
}

void
combineHexaloop (const struct hexaloopE *hexaloopEnergies,
		 const struct hexaloopE *hexaloopEnthalpies, double tRatio, struct hexaloop *hexaloop, int num)
{
   int i;

   for (i = 0; i < num; ++i) {
      memcpy (hexaloop[i].loop, hexaloopEnergies[i].loop, 8);
      hexaloop[i].energy = combineEE (hexaloopEnergies[i].energy, hexaloopEnthalpies[i].energy, tRatio, 0.0);
   }
}

int hexaloopcmp (const void *loop1, const void *loop2)
{
   int i;
   const unsigned char *h1 = loop1;
   const struct hexaloop *h2 = loop2;

   for (i = 0; i < 8; ++i)
      if (h1[i] < h2->loop[i])
	 return -1;
      else if (h1[i] > h2->loop[i])
	 return 1;

   return 0;
}
