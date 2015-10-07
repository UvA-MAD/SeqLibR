#ifndef UTIL_H
#define UTIL_H

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define TURN 3			/* minimum size of hairpin loop */

/* #define NO_GU_BASEPAIRS */

extern const int BPI[6][6];
#define basePairIndex(a, b) BPI[a][b]

int roundInt (double d);
void strcatc (char *str, char c);
char *filename (char *file);
unsigned char toNum (char c);
int seqcmp (unsigned char *seq1, unsigned char *seq2, int length);
int min3 (int a, int b, int c);
int same (unsigned char *a, unsigned char *b, int len);
void readOrDie (unsigned int num, const char *name, FILE * file, const char *format, ...);

#endif
