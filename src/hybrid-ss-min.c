# include "config.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#if HAVE_IEEEFP_H
# include <ieeefp.h>
#endif

#include "energy.h"
#include "getopt.h"
#include "util.h"
#include "xmalloc.h"
#include "sequence.h"

/* hybrid-ss-min is developed by 
 * by Nick Markham & Michael Zuker
 * © 2006 Rensselaer Polytechnic Institute
 * Troy, NY 12180
 * compute minimum energy folding of NA sequence 
 */

#define Q(i, j) q[(g_len - 1) * (i - 1) + j - 1]
#define Qprime(i, j) qprime[(g_len - 1) * (i - 1) + j - 1]
#define QM(i, j) qm[(g_len - 1) * (i - 1) + j - 1]
#define Q5(i) q5[i]
#define Q3(j) q3[j - 1]
#define ssOK(i, j) 1

#define auPenalty(i, j) g_aup[g_seq[i]][g_seq[j]]
#define min2(a, b) ((a) < (b) ? (a) : (b))

#define IMPOSSIBLE_TEMPERATURE -300.0
#define IMPOSSIBLE_CONCENTRATION -10.0

static int g_len;
static ENERGY *q = NULL, *qprime = NULL, *qm = NULL, *q5 = NULL, *q3 = NULL;

static int g_maxLoop = 30;
static int g_prefilter1 = 2;
static int g_prefilter2 = 2;
static int g_circular = 0;
static unsigned char *g_seq;	/* [0-4] for [A,C,G,TU,N] */

static ENERGY g_dangle3[5][5][6];
static ENERGY g_dangle5[5][5][6];
static ENERGY g_stack[5][5][5][5];
static ENERGY g_hairpinLoop[30];
static ENERGY g_interiorLoop[30];
static ENERGY g_bulgeLoop[30];
static ENERGY g_sint2[7][7][5][5];
static ENERGY g_asint1x2[7][7][5][5][5];
static ENERGY g_sint4[7][7][5][5][5][5];
static ENERGY g_tstackh[5][5][5][5];
static ENERGY g_tstacki[5][5][5][5];
static ENERGY g_tstacki23[5][5][5][5];
static ENERGY g_tstackm[5][5][6][6];
static ENERGY g_tstacke[5][5][6][6];
static struct triloop *g_triloop = NULL;
static int numTriloops = 0;
static struct tloop *g_tloop = NULL;
static int numTloops = 0;
static struct hexaloop *g_hexaloop = NULL;
static int numHexaloops = 0;
static ENERGY g_multi[3];
static ENERGY g_misc[13];
static ENERGY g_aup[5][5];

static emat3 dangleEnergies3, dangleEnergies5;
static emat4 stackEnergies, tstackhEnergies, tstackiEnergies, tstacki23Energies, tstackmEnergies, tstackeEnergies;
static emat4e stackEnthalpies, tstackhEnthalpies, tstackiEnthalpies, tstacki23Enthalpies;
static emat4s tstackmEnthalpies, tstackeEnthalpies;
static emat3e dangleEnthalpies3, dangleEnthalpies5;
static emat_s2e sint2Energies;
static emat_s2g sint2Enthalpies;
static emat_as1e asint1x2Energies;
static emat_as1g asint1x2Enthalpies;
static emat_s4e sint4Energies;
static emat_s4g sint4Enthalpies;


static double interiorLoopEnergies[30], bulgeLoopEnergies[30], hairpinLoopEnergies[30];
static double interiorLoopEnthalpies[30], bulgeLoopEnthalpies[30], hairpinLoopEnthalpies[30];

static struct triloopE *triloopEnergies;
static struct triloopE *triloopEnthalpies;
static struct tloopE *tloopEnergies;
static struct tloopE *tloopEnthalpies;
static struct hexaloopE *hexaloopEnergies;
static struct hexaloopE *hexaloopEnthalpies;
static double multiEnergies[3];
static double multiEnthalpies[3];
static double miscEnergies[13];
static double miscEnthalpies[13];

static int loaded_data = -1;
static double last_temp = IMPOSSIBLE_TEMPERATURE;
static double lastNaConc = IMPOSSIBLE_CONCENTRATION;
static double lastMgConc = IMPOSSIBLE_CONCENTRATION;

ENERGY *recalloc2 (ENERGY * ptr, int n)
{
   return xrealloc (ptr, n * n * sizeof (ENERGY));
}

void prefilter ()
{
   char **in;
   int i, j, k, count;

   in = xcalloc (g_len, sizeof (char *));
   for (i = 1; i <= g_len; ++i)
      in[i - 1] = xcalloc (g_len, 1);

   for (i = 1; i <= g_len - g_prefilter2 + 1; ++i)
      for (j = g_len; j >= g_prefilter2 && j >= i; --j) {
	 count = 0;
	 for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	    if (isFinite (Qprime (i + k, j - k)))
	       ++count;
	 if (count >= g_prefilter1)
	    for (k = 0; k < g_prefilter2 && k <= (j - i) / 2; ++k)
	       ++in[i + k - 1][j - k - 1];
      }

   for (i = 1; i <= g_len; ++i) {
      for (j = g_len; j >= i; --j)
	 if (!in[i - 1][j - 1])
	    Qprime (i, j) = INFINITY;
      xfree (in[i - 1]);
   }
   xfree (in);
}

void initializeMatrices ()
{
   int i, j;

   /* Q' is initialized to +infinity iff base pair is illegal; 0 otherwise
      Q and QM are always initialized to +infinity */
   for (i = 1; i <= g_len; ++i) {
      for (j = i; j <= g_len; ++j) {
	 if (j - i < TURN + 1 || basePairIndex (g_seq[i], g_seq[j]) == 6) {
	    Q (i, j) = Qprime (i, j) = QM (i, j) = INFINITY;
	 }
	 else {
	    Q (i, j) = QM (i, j) = INFINITY;
	    Qprime (i, j) = 0;
	 }
      }
   }
   if (g_prefilter1 > 1)
      prefilter ();
}

ENERGY Ed3 (int i, int j)
{
   return ssOK (i + 1, i + 1) ? g_dangle3[g_seq[i]][g_seq[j]][g_seq[i + 1]] : INFINITY;
}

ENERGY Ed5 (int i, int j)
{
   return ssOK (j - 1, j - 1) ? g_dangle5[g_seq[i]][g_seq[j]][g_seq[j - 1]] : INFINITY;
}

ENERGY Etstacke (int i, int j)
{
   return (ssOK (i + 1, i + 1)
	   && ssOK (j - 1, j - 1)) ? g_tstacke[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] : INFINITY;
}

ENERGY Q5_1 (int i)
{
   int k;
   ENERGY min;

   min = INFINITY;
   for (k = 0; k <= i - TURN - 2; ++k)
      min = min2 (min, min2 (ssOK (1, k) ? 0 : INFINITY, Q5 (k)) + auPenalty (k + 1, i) + Qprime (k + 1, i));

   return min;
}

ENERGY Q5_2 (int i)
{
   int k;
   ENERGY min;

   min = INFINITY;
   for (k = 0; k <= i - TURN - 3; ++k)
      min = min2 (min, min2 (ssOK (1, k) ? 0 : INFINITY, Q5 (k)) + auPenalty (k + 2, i) + Ed5 (i, k + 2) + Qprime (k + 2, i));

   return min;
}

ENERGY Q5_3 (int i)
{
   int k;
   ENERGY min;

   min = INFINITY;
   for (k = 0; k <= i - TURN - 3; ++k)
      min =
	 min2 (min,
	       min2 (ssOK (1, k) ? 0 : INFINITY, Q5 (k)) + auPenalty (k + 1, i - 1) + Ed3 (i - 1, k + 1) + Qprime (k + 1, i - 1));
   return min;
}

ENERGY Q5_4 (int i)
{
   int k;
   ENERGY min;

   min = INFINITY;
   for (k = 0; k <= i - TURN - 4; ++k)
      min =
	 min2 (min,
	       min2 (ssOK (1, k) ? 0 : INFINITY, Q5 (k)) + auPenalty (k + 2, i - 1) +
	       Etstacke (i - 1, k + 2) + Qprime (k + 2, i - 1));
   return min;
}

ENERGY Q3_1 (int j)
{
   int k;
   ENERGY min;

   min = INFINITY;
   for (k = j + TURN + 2; k <= g_len + 1; ++k)
      min = min2 (min, auPenalty (j, k - 1) + Qprime (j, k - 1) + min2 (ssOK (k, g_len) ? 0 : INFINITY, Q3 (k)));
   return min;
}

ENERGY Q3_2 (int j)
{
   int k;
   ENERGY min;

   min = INFINITY;
   for (k = j + TURN + 3; k <= g_len + 1; ++k)
      min = min2 (min, auPenalty (j, k - 2) + Qprime (j, k - 2) + Ed3 (k - 2, j) + min2 (ssOK (k, g_len) ? 0 : INFINITY, Q3 (k)));

   return min;
}

ENERGY Q3_3 (int j)
{
   int k;
   ENERGY min;

   min = INFINITY;
   for (k = j + TURN + 3; k <= g_len + 1; ++k)
      min =
	 min2 (min,
	       auPenalty (j + 1, k - 1) + Qprime (j + 1, k - 1) + Ed5 (k - 1, j + 1) +
	       min2 (ssOK (k, g_len) ? 0 : INFINITY, Q3 (k)));

   return min;
}

ENERGY Q3_4 (int j)
{
   int k;
   ENERGY min;

   min = INFINITY;
   for (k = j + TURN + 4; k <= g_len + 1; ++k)
      min =
	 min2 (min,
	       auPenalty (j + 1, k - 2) + Qprime (j + 1, k - 2) + Etstacke (k - 2, j + 1) +
	       min2 (ssOK (k, g_len) ? 0 : INFINITY, Q3 (k)));

   return min;
}

ENERGY Etstackm (int i, int j)
{
   return (ssOK (i + 1, i + 1)
	   && ssOK (j - 1, j - 1)) ? g_tstackm[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] : INFINITY;
}

ENERGY Eh (int i, int j)
{
   ENERGY energy;
   int loopSize = j - i - 1;
   int k;

   if (loopSize < TURN)
      return INFINITY;

   if (i <= g_len && g_len < j)
      return INFINITY;
   else if (i > g_len) {
      i -= g_len;
      j -= g_len;
   }

   if (loopSize <= 30)
      energy = g_hairpinLoop[loopSize - 1];
   else
      energy = g_hairpinLoop[29] + g_misc[12] * log ((double) loopSize / 30);

   if (loopSize > 3)
      energy += g_tstackh[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
   else
      energy += auPenalty (i, j);

   if (loopSize == 3) {
      struct triloop *loop;
      if (numTriloops)
	 if ((loop = bsearch (g_seq + i, g_triloop, numTriloops, sizeof (struct triloop), triloopcmp)))
	    energy += loop->energy;
   }
   else if (loopSize == 4) {
      struct tloop *loop;
      if (numTloops)
	 if ((loop = bsearch (g_seq + i, g_tloop, numTloops, sizeof (struct tloop), tloopcmp)))
	    energy += loop->energy;
   }
   else if (loopSize == 6) {
      struct hexaloop *loop;
      if (numHexaloops)
	 if ((loop = bsearch (g_seq + i, g_hexaloop, numHexaloops, sizeof (struct hexaloop), hexaloopcmp)))
	    energy += loop->energy;
   }

   /* GGG */
   if (i >= 3 && g_seq[i - 2] == 2 && g_seq[i - 1] == 2 && g_seq[i] == 2 && g_seq[j] == 3)
      energy += g_misc[8];

   /* poly-C */
   if (loopSize == 3 && g_seq[i + 1] == 1 && g_seq[i + 2] == 1 && g_seq[i + 3] == 1)
      energy += g_misc[11];
   else {
      for (k = 1; k <= loopSize; ++k)
	 if (g_seq[i + k] != 1)
	    return energy;
      energy += g_misc[9] * loopSize + g_misc[10];
   }

   return energy;
}

ENERGY Es (int i, int j)
{
   if (i >= j)
      return INFINITY;
   /* fputs("Error in Es(): i isn't less than j\n", stderr); */

   if (i == g_len || j == g_len + 1)
      return INFINITY;

   if (i > g_len)
      i -= g_len;
   if (j > g_len)
      j -= g_len;

   return g_stack[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
}

ENERGY Ebi (int i, int j, int ii, int jj)
{
   int loopSize1, loopSize2;
   ENERGY loopEnergy, asPenalty;

   loopSize1 = ii - i - 1;
   loopSize2 = j - jj - 1;
   if (loopSize1 + loopSize2 > g_maxLoop)
      return INFINITY;

   if (loopSize1 == 0) {
      if (loopSize2 == 1) {
	 return g_bulgeLoop[0] + g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]];
      }
      else if (loopSize2 <= 30) {
	 return g_bulgeLoop[loopSize2 - 1] + auPenalty (i, j) + auPenalty (ii, jj);
      }
      else {
	 return g_bulgeLoop[29] + g_misc[12] * log ((double) loopSize2 / 30) + auPenalty (i, j) + auPenalty (ii, jj);
      }
   }
   else if (loopSize2 == 0) {
      if (loopSize1 == 1) {
	 return g_bulgeLoop[0] + g_stack[g_seq[i]][g_seq[j]][g_seq[ii]][g_seq[jj]];
      }
      else if (loopSize1 <= 30) {
	 return g_bulgeLoop[loopSize1 - 1] + auPenalty (i, j) + auPenalty (ii, jj);
      }
      else {
	 return g_bulgeLoop[29] + g_misc[12] * log ((double) loopSize1 / 30) + auPenalty (i, j) + auPenalty (ii, jj);
      }
   }
   else if (loopSize1 == 1 && loopSize2 == 1) {
      return g_sint2[basePairIndex (g_seq[i], g_seq[j])][basePairIndex (g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]];
   }
   else if (loopSize1 == 1 && loopSize2 == 2) {
      return
	 g_asint1x2[basePairIndex (g_seq[i], g_seq[j])][basePairIndex
							(g_seq[ii], g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[j - 2]];
   }
   else if (loopSize1 == 2 && loopSize2 == 1) {
      return
	 g_asint1x2[basePairIndex (g_seq[jj], g_seq[ii])][basePairIndex
							  (g_seq[j], g_seq[i])][g_seq[jj + 1]][g_seq[ii - 1]][g_seq[ii - 2]];
   }
   else if (loopSize1 == 2 && loopSize2 == 2) {
      return
	 g_sint4[basePairIndex (g_seq[i], g_seq[j])][basePairIndex
						     (g_seq[ii],
						      g_seq[jj])][g_seq[i + 1]][g_seq[j - 1]][g_seq[i + 2]][g_seq[j - 2]];
   }
   else if ((loopSize1 == 2 && loopSize2 == 3) || (loopSize1 == 3 && loopSize2 == 2)) {
      return g_tstacki23[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]] +
	 g_tstacki23[g_seq[jj]][g_seq[ii]][g_seq[jj + 1]][g_seq[ii - 1]];
   }
   else {
      if (loopSize1 + loopSize2 <= 30) {
	 loopEnergy = g_interiorLoop[loopSize1 + loopSize2 - 1];
      }
      else {
	 loopEnergy = g_interiorLoop[29] + g_misc[12] * log ((double) (loopSize1 + loopSize2) / 30);
      }
      if (g_misc[7] && (loopSize1 == 1 || loopSize2 == 1)) {
	 loopEnergy += g_tstacki[g_seq[i]][g_seq[j]][0][0];
	 loopEnergy += g_tstacki[g_seq[jj]][g_seq[ii]][0][0];
      }
      else {
	 loopEnergy += g_tstacki[g_seq[i]][g_seq[j]][g_seq[i + 1]][g_seq[j - 1]];
	 loopEnergy += g_tstacki[g_seq[jj]][g_seq[ii]][g_seq[jj + 1]][g_seq[ii - 1]];
      }
      asPenalty = abs (loopSize1 - loopSize2) * g_misc[min3 (4, loopSize1, loopSize2) - 1];
      if (asPenalty > g_misc[4]) {
	 asPenalty = g_misc[4];
      }
      loopEnergy += asPenalty;

      return loopEnergy;
   }
}

ENERGY QBI (int i, int j)
{
   int d, ii, jj;
   ENERGY energy = INFINITY;

   for (d = j - i - 3; d >= TURN + 1 && d >= j - i - 2 - g_maxLoop; --d)
      for (ii = i + 1; ii < j - d && ii <= g_len; ++ii) {
	 jj = d + ii;
	 if (isFinite (Qprime (ii, jj)))
	    energy = min2 (energy, Ebi (i, j, ii, jj) + Qprime (ii, jj));
      }
   return energy;
}

ENERGY min4 (ENERGY a, ENERGY b, ENERGY c, ENERGY d)
{
   if (a < b && a < c && a < d)
      return a;
   else if (b < c && b < d)
      return b;
   else if (c < d)
      return c;
   else
      return d;
}

ENERGY min5 (ENERGY a, ENERGY b, ENERGY c, ENERGY d, ENERGY e)
{
   if (a < b && a < c && a < d && a < e)
      return a;
   else if (b < c && b < d && b < e)
      return b;
   else if (c < d && c < e)
      return c;
   else if (d < e)
      return d;
   else
      return e;
}

void fillMatrices1 ()
{
   int i, j, k;

   /* start at top left, fill each column bottom->top
      when Q' is +infinity, don't consider it */
   for (j = 2; j <= g_len; ++j)
      for (i = j - TURN - 1; i >= (g_circular && j > g_len / 2 ? j - g_len / 2 : 1); --i) {
	 ENERGY au;
	 au = auPenalty (i, j);

	 if (g_circular) {
	    if (j - i > g_len / 2)
	       continue;
	    if (i > g_len / 2) {
	       Q (i, j) = Q (i - g_len / 2, j - g_len / 2);
	       Qprime (i, j) = Qprime (i - g_len / 2, j - g_len / 2);
	       QM (i, j) = QM (i - g_len / 2, j - g_len / 2);
	       continue;
	    }
	 }

	 if (isFinite (Qprime (i, j))) {
	    Qprime (i, j) = min4 (Eh (i, j),
				  Es (i, j) + Qprime (i + 1, j - 1), QBI (i, j), g_multi[0] + g_multi[2] + au + QM (i + 1, j - 1));

	    if (j > 2)
	       Qprime (i, j) = min2 (Qprime (i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed5 (i, j) + QM (i + 1, j - 2));
	    if (i < g_len - 1)
	       Qprime (i, j) = min2 (Qprime (i, j), g_multi[0] + g_multi[1] + g_multi[2] + au + Ed3 (i, j) + QM (i + 2, j - 1));
	    if (j > 2 && i < g_len - 1)
	       Qprime (i, j) =
		  min2 (Qprime (i, j), g_multi[0] + 2.0 * g_multi[1] + g_multi[2] + au + Etstackm (i, j) + QM (i + 2, j - 2));

	 }

	 QM (i, j) = INFINITY;
	 for (k = i + TURN + 1; k <= j - TURN - 2; ++k)
	    QM (i, j) = min2 (QM (i, j), Q (i, k) + Q (k + 1, j));

	 Q (i, j) = min4 (ssOK (i, i) ? g_multi[1] + Q (i + 1, j) : INFINITY,
			  ssOK (j, j) ? g_multi[1] + Q (i, j - 1) : INFINITY, g_multi[2] + au + Qprime (i, j), QM (i, j));

	 Q (i, j) = min2 (Q (i, j), g_multi[1] + g_multi[2] + auPenalty (i + 1, j) + Ed5 (j, i + 1) + Qprime (i + 1, j));
	 Q (i, j) = min2 (Q (i, j), g_multi[1] + g_multi[2] + auPenalty (i, j - 1) + Ed3 (j - 1, i) + Qprime (i, j - 1));
	 if (i < j - TURN - 2) {
	    Q (i, j) =
	       min2 (Q (i, j),
		     2.0 * g_multi[1] + g_multi[2] + auPenalty (i + 1, j - 1) + Etstackm (j - 1, i + 1) + Qprime (i + 1, j - 1));
	 }

      }
}

void computeQ53 ()
{
   int i, j;

   Q5 (0) = Q5 (1) = INFINITY;
   Q3 (g_len + 1) = Q3 (g_len) = INFINITY;

   for (i = 2; i <= g_len; ++i)
      Q5 (i) = min5 (ssOK (i, i) ? Q5 (i - 1) : INFINITY, Q5_1 (i), Q5_2 (i), Q5_3 (i), Q5_4 (i));

   for (j = g_len - 1; j >= 1; --j)
      Q3 (j) = min5 (ssOK (j, j) ? Q3 (j + 1) : INFINITY, Q3_1 (j), Q3_2 (j), Q3_3 (j), Q3_4 (j));

}

int LoadMatrices (int NA)
{
   loadStack (stackEnergies, stackEnthalpies, NA);
   loadDangle (dangleEnergies3, dangleEnthalpies3, dangleEnergies5, dangleEnthalpies5, NA);
   loadLoop (hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies,
	     hairpinLoopEnthalpies, interiorLoopEnthalpies, bulgeLoopEnthalpies, NA);
   loadSint2 (sint2Energies, sint2Enthalpies, NA);
   loadAsint1x2 (asint1x2Energies, asint1x2Enthalpies, NA);
   loadSint4 (sint4Energies, sint4Enthalpies, NA);
   loadTstack (tstackhEnergies, tstackhEnthalpies, NA, "tstackh");
   loadTstack (tstackiEnergies, tstackiEnthalpies, NA, "tstacki");
   loadTstack (tstacki23Energies, tstacki23Enthalpies, NA, "tstacki23");
   loadTstack2 (tstackmEnergies, tstackmEnthalpies, NA, "tstackm");
   loadTstack2 (tstackeEnergies, tstackeEnthalpies, NA, "tstacke");
   loadTriloop (&triloopEnergies, &triloopEnthalpies, &numTriloops, NA);
   loadTloop (&tloopEnergies, &tloopEnthalpies, &numTloops, NA);
   loadHexaloop (&hexaloopEnergies, &hexaloopEnthalpies, &numHexaloops, NA);
   loadMulti (multiEnergies, multiEnthalpies, NA);
   loadMisc (miscEnergies, miscEnthalpies, NA);

   symmetryCheckStack (stackEnergies, "energy");
   symmetryCheckSint2 (sint2Energies, "energy");
   symmetryCheckSint4 (sint4Energies, "energy");
   g_triloop = (struct triloop *) xrealloc (g_triloop, numTriloops * sizeof (struct triloop));
   g_tloop = (struct tloop *) xrealloc (g_tloop, numTloops * sizeof (struct tloop));
   g_hexaloop = (struct hexaloop *) xrealloc (g_hexaloop, numHexaloops * sizeof (struct hexaloop));
   last_temp = IMPOSSIBLE_TEMPERATURE;	/* invalidate Concentration/Temp cache */
   return NA;
}

int PrepareMatrices (double temp, double naConc, double mgConc, int NA)
{
   double tRatio, saltCorrection;
   int polymer = 0;
   saltCorrection = ion (NA, polymer, naConc, mgConc);
   tRatio = (temp + 273.15) / 310.15;
   combineStack (stackEnergies, stackEnthalpies, tRatio, saltCorrection, g_stack);
   combineDangle (dangleEnergies3, dangleEnergies5, dangleEnthalpies3, dangleEnthalpies5, tRatio, saltCorrection, g_dangle3,
		  g_dangle5);
   combineLoop (hairpinLoopEnergies, interiorLoopEnergies, bulgeLoopEnergies, hairpinLoopEnthalpies, interiorLoopEnthalpies,
		bulgeLoopEnthalpies, tRatio, saltCorrection, g_hairpinLoop, g_interiorLoop, g_bulgeLoop);
   combineSint2 (sint2Energies, sint2Enthalpies, tRatio, saltCorrection, g_sint2);
   combineAsint1x2 (asint1x2Energies, asint1x2Enthalpies, tRatio, saltCorrection, g_asint1x2);
   combineSint4 (sint4Energies, sint4Enthalpies, tRatio, saltCorrection, g_sint4);
   combineTstack1 (tstackiEnergies, tstackiEnthalpies, tRatio, g_tstacki);
   combineTstack1 (tstacki23Energies, tstacki23Enthalpies, tRatio, g_tstacki23);
   combineTstack1 (tstackhEnergies, tstackhEnthalpies, tRatio, g_tstackh);
   combineTstack2 (tstackmEnergies, tstackmEnthalpies, tRatio, saltCorrection, g_tstackm);
   combineTstack2 (tstackeEnergies, tstackeEnthalpies, tRatio, saltCorrection, g_tstacke);
   combineMulti (multiEnergies, multiEnthalpies, tRatio, g_multi);
   combineMisc (miscEnergies, miscEnthalpies, tRatio, g_misc);
   combineTriloop (triloopEnergies, triloopEnthalpies, tRatio, g_triloop, numTriloops);
   combineTloop (tloopEnergies, tloopEnthalpies, tRatio, g_tloop, numTloops);
   combineHexaloop (hexaloopEnergies, hexaloopEnthalpies, tRatio, g_hexaloop, numHexaloops);
   makeAUPenalty (g_misc, g_aup, 0);
   last_temp = temp;
   lastNaConc = naConc;
   lastMgConc = mgConc;
   return 1;
}

double hybrid_main (char *cseq, double temp, double naConc, double mgConc, int acid_type)
{
   int i, bestI, bestJ;
   ENERGY mfe;

   if (acid_type == ACID_TYPE_RNA && (naConc != 1.0 || mgConc != 0.0)) {
      fputs ("Warning: salt concentrations ignored for RNA\n", stderr);
   }

   if (g_maxLoop < 0) {
      g_maxLoop = 999999999;
   }

   /* read free energies and entropies */
   if (loaded_data != acid_type)
      loaded_data = LoadMatrices (acid_type);

   g_len = strlen (cseq);

   /* convert sequence to numbers for speed */
   g_seq = xrealloc (g_seq, g_circular ? 2 * g_len + 2 : g_len + 2);
   for (i = 1; i <= g_len; ++i)
      g_seq[i] = toNum (cseq[i - 1]);
   if (g_circular) {
      for (i = 1; i <= g_len; ++i)
	 g_seq[g_len + i] = toNum (cseq[i - 1]);
      g_seq[0] = g_seq[2 * g_len + 1] = 5;
   }
   else
      g_seq[0] = g_seq[g_len + 1] = 5;

   if (g_circular)
      g_len *= 2;

   q = xrealloc (q, g_len * g_len * sizeof (ENERGY));
   qprime = xrealloc (qprime, g_len * g_len * sizeof (ENERGY));
   qm = xrealloc (qm, g_len * g_len * sizeof (ENERGY));
   q5 = xrealloc (q5, (g_len + 1) * sizeof (ENERGY));
   q3 = xrealloc (q3, (g_len + 1) * sizeof (ENERGY));

   if ((temp != last_temp) || (lastNaConc != naConc) || (lastMgConc != mgConc)) {
      PrepareMatrices (temp, naConc, mgConc, acid_type);
   }
   bestI = bestJ = 0;

   initializeMatrices ();
   fillMatrices1 ();
   if (g_circular) {
      int i, j;
      mfe = INFINITY;
      for (i = 1; i < g_len / 2; ++i)
	 for (j = i + 1; j <= g_len / 2; ++j)
	    if (Qprime (i, j) + Qprime (j, i + g_len / 2) < mfe) {
	       bestI = i;
	       bestJ = j;
	       mfe = Qprime (i, j) + Qprime (j, i + g_len / 2);
	    }
   }
   else {
      computeQ53 ();
      if (fabs (Q5 (g_len) - Q3 (1)) > 1e-12)
	 fprintf (stderr, "Warning: Q5(n) = %g but Q3(1) = %g.\n", (double) Q5 (g_len) / PRECISION, (double) Q3 (1) / PRECISION);
      mfe = Q5 (g_len);
   }
   return (double) mfe / PRECISION;
}
