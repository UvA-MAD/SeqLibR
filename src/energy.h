#ifndef ENERGY_H
#define ENERGY_H

/* include this file to use functions in energy.c
 * then link with energy.o
 */

#ifndef ENERGY
# define ENERGY double
#endif

#ifndef PRECISION
# define PRECISION 1
#endif

#ifdef INTEGER
# define isFinite(x) (x < INFINITY / 2)
#else
# define isFinite(x) finite(x)
#endif

#ifdef INFINITY
# undef INFINITY
#endif

typedef double emat3[4][4][4];
typedef double emat3e[5][5][6];
typedef double emat4[4][4][4][4];
typedef double emat4e[5][5][5][5];
typedef double emat4s[5][5][6][6];
typedef double emat_s2e[6][6][4][4];
typedef double emat_s2g[7][7][5][5];
typedef double emat_as1e[6][6][4][4][4];
typedef double emat_as1g[7][7][5][5][5];
typedef double emat_s4e[6][6][4][4][4][4];
typedef double emat_s4g[7][7][5][5][5][5];


extern const ENERGY INFINITY;
extern const double R1;
extern const char BASES[5];
extern const char BASE_PAIRS[6][4];

double ion (int NA, int polymer, double naConc, double mgConc);

void loadStack (emat4 stackEnergies, emat4e stackEnthalpies, int NA);

void combineStack (emat4 stackEnergies, emat4e stackEnthalpies, double tRatio, double saltCorrection, ENERGY stack[5][5][5][5]);

void symmetryCheckStack (emat4 stack, char *which);

void loadDangle (emat3 dangleEnergies3, emat3e dangleEnthalpies3, emat3 dangleEnergies5, emat3e dangleEnthalpies5, int NA);

void combineDangle (emat3 dangleEnergies3, emat3 dangleEnergies5,
		    emat3e dangleEnthalpies3,
		    emat3e dangleEnthalpies5, double tRatio,
		    double saltCorrection, ENERGY dangle3[5][5][6], ENERGY dangle5[5][5][6]);

void loadLoop (double hairpinLoopEnergies[30],
	       double interiorLoopEnergies[30], double bulgeLoopEnergies[30],
	       double hairpinLoopEnthalpies[30], double interiorLoopEnthalpies[30], double bulgeLoopEnthalpies[30], int NA);

void combineLoop (double hairpinLoopEnergies[30],
		  double interiorLoopEnergies[30],
		  double bulgeLoopEnergies[30],
		  double hairpinLoopEnthalpies[30],
		  double interiorLoopEnthalpies[30],
		  double bulgeLoopEnthalpies[30], double tRatio,
		  double saltCorrection, ENERGY hairpinLoop[30], ENERGY interiorLoop[30], ENERGY bulgeLoop[30]);

void loadSint2 (emat_s2e sint2Energies, emat_s2g sint2Enthalpies, int NA);

void combineSint2 (emat_s2e sint2Energies, emat_s2g sint2Enthalpies, double tRatio, double saltCorrection,
		   ENERGY sint2[7][7][5][5]);

void symmetryCheckSint2 (emat_s2e sint2, char *which);

void loadAsint1x2 (emat_as1e asint1x2Energies, emat_as1g asint1x2Enthalpies, int NA);

void combineAsint1x2 (emat_as1e asint1x2Energies, emat_as1g asint1x2Enthalpies, double tRatio,
		      double saltCorrection, ENERGY asint1x2[7][7][5][5][5]);

void loadSint4 (emat_s4e sint4Energies, emat_s4g sint4Enthalpies, int NA);

void combineSint4 (emat_s4e sint4Energies, emat_s4g sint4Enthalpies, double tRatio,
		   double saltCorrection, ENERGY sint4[7][7][5][5][5][5]);

void symmetryCheckSint4 (emat_s4e sint4, char *which);

void loadTstack (emat4 tsEng, emat4e tsEnth, int NA, char *fn);

void loadTstack2 (emat4 tstackmEnergies, emat4s tstackmEnthalpies, int NA, char *fn);

void combineTstack1 (emat4 tstackEnergies, emat4e tstackEnthalpies, double tRatio, ENERGY tstack[5][5][5][5]);

void combineTstack2 (emat4 tstackEnergies, emat4s tstackEnthalpies,
		     double tRatio, double saltCorrection, ENERGY tstack[5][5][6][6]);

void loadMulti (double multiEnergies[3], double multiEnthalpies[3], int NA);

void combineMulti (double multiEnergies[3], double multiEnthalpies[3], double tRatio, ENERGY multi[3]);

void loadMisc (double miscEnergies[13], double miscEnthalpies[13], int NA);

void combineMisc (double miscEnergies[13], double miscEnthalpies[13], double tRatio, ENERGY misc[13]);

void makeAUPenalty (ENERGY misc[13], ENERGY aup[5][5], int isPF);

struct triloopE
{
   char loop[5];
   double energy;
};

struct triloop
{
   char loop[5];
   ENERGY energy;
};

void loadTriloop (struct triloopE **, struct triloopE **, int *, int);

void combineTriloop (const struct triloopE *, const struct triloopE *, double, struct triloop *, int);

int triloopcmp (const void *, const void *);

struct tloopE
{
   char loop[6];
   double energy;
};

struct tloop
{
   char loop[6];
   ENERGY energy;
};

void loadTloop (struct tloopE **, struct tloopE **, int *, int);

void combineTloop (const struct tloopE *, const struct tloopE *, double, struct tloop *, int);

int tloopcmp (const void *, const void *);

struct hexaloopE
{
   char loop[8];
   double energy;
};

struct hexaloop
{
   char loop[8];
   ENERGY energy;
};

void loadHexaloop (struct hexaloopE **, struct hexaloopE **, int *, int);

void combineHexaloop (const struct hexaloopE *, const struct hexaloopE *, double, struct hexaloop *, int);

int hexaloopcmp (const void *, const void *);

void readOrDie (unsigned int num, const char *name, FILE * file, const char *format, ...);

#endif
