#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "sequence.h"
#include "hybrid-ss-min.h"

#define MINV(a,b)   (((a)<(b))?(a):(b))
#define MAXV(a,b)   (((a)>(b))?(a):(b))
#define ABSV(a)   (((a)>(0))?(a):(-(a)))
#define GET_NUC_2B(s,p) (((s)->elms[(p)>>2] >> stab2b[(p)&3])&3)

#define WARNING_STRING_LENGTH 2048


static char nuc_tab_enum[16] = { 'A', 'C', 'G', 'T', 'N', 'R', 'S', 'V', 'W', 'Y', 'H', 'K', 'D', 'B', 'M', '-' };
nuc_tp comp_tab[16] = { T, G, C, A, N, Y, S, B, W, R, D, M, H, V, K, gap };
nuc_tp code4b_to_nuc[16] = { gap, A, C, M, G, R, S, V, T, W, Y, H, K, D, B, N };
nuc_tp code2b_to_nuc[4] = { A, C, G, T };
static unsigned char nuc_code_2b[16] = { 0, 1, 2, 3, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };
static unsigned char nuc_code_4b[16] = { 1, 2, 4, 8, 15, 5, 6, 7, 9, 10, 11, 12, 13, 14, 3, 0 };

/* enum nucleotides_2b {As,Cs,Gs,Ts};
typedef enum nucleotides_2b nuc_2b_tp; */
static int mask_2b[4] = { 0x3f, 0xcf, 0xf3, 0xfc };
static char nuc_tab_2b[4] = { 'a', 'c', 'g', 't' };
static char cap_nuc_tab_2b[4] = { 'A', 'C', 'G', 'T' };

/*enum nucleotides_4b {le,Ae,Ce,Me,Ge,Re,Se,Ve,Te,We,Ye,He,Ke,De,Be,Ne};
typedef enum nucleotides_ext nuc_4b_tp; */
static int mask_4b[2] = { 0x0f, 0xf0 };
static char nuc_tab_4b[16] = { '-', 'a', 'c', 'm', 'g', 'r', 's', 'v', 't', 'w', 'y', 'h', 'k', 'd', 'b', 'n' };
static char cap_nuc_tab_4b[16] = { '-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N' };

static void (*ext_warn) (char *message) = NULL;

typedef struct restriction_iterator
{
   sequence_tp *sequence;
   sequence_tp *restriction_sequence;
   sequence_tp *restriction_sequence_rc;
   int cut_pos;
   int end;
   int current_location;
} restriction_iterator;



void sequence_catch_warning (void (*my_warn) (char *))
{
   ext_warn = my_warn;
}


static void seq_warning (char *w)
{
   if (ext_warn)
      ext_warn (w);
#ifdef DEBUG
   else
      fprintf (stderr, "%s\n", w);
#endif
}


static void handle_error ()
{
#ifdef DEBUG
   fprintf (stderr, "Encountered error bye\n");
   exit (-1);
#endif
}

void sequence_free (sequence_tp * s)
{
   free (s->elms);
   free (s);
}

static void sequence_extend (sequence_tp * s)
{
/* independent of type */
   unsigned char *new;
   int ns;
   ns = MAXV (s->elm_alloced_bytes + 20, (int) (s->elm_alloced_bytes * 1.2));
/*	fprintf(stderr,"Extend to %d\n",ns); */
   new = (unsigned char *) malloc (ns);
   if (s->elm_alloced_bytes) {
      memcpy (new, s->elms, s->elm_alloced_bytes);
      free (s->elms);
   }
   s->elms = new;
   s->elm_alloced_bytes = ns;
}

char nuc_char (nuc_tp n)
{
   return nuc_tab_enum[n];
}

void nucleotide_print (nuc_tp n, FILE * df)
{
   fprintf (df, "%c", nuc_tab_enum[n]);
}

static char *sequence_chars_4b (sequence_tp * s, char *res, char *nuctable)
{
   int j, pos = 0;
   unsigned char qrt;
   if (!s) {
      res[pos] = '\0';
      return res;
   }
   for (j = 0; j < s->len / 2; j++) {
      qrt = s->elms[j];
      res[pos++] = nuctable[qrt >> 4];
      res[pos++] = nuctable[qrt & 15];
   }
   qrt = s->elms[s->len / 2];
   if (s->len % 2)
      res[pos++] = nuctable[qrt >> 4];
   res[pos] = '\0';
   return res;
}
static char *sequence_chars_2b (sequence_tp * s, char *res, char *nuctable)
{
   int i, j, pos = 0;
   unsigned char qrt;
   if (!s) {
      res[pos] = '\0';
      return res;
   }
   for (j = 0; j < s->len / 4; j++) {
      qrt = s->elms[j];
      res[pos++] = nuctable[qrt >> 6];
      res[pos++] = nuctable[(qrt >> 4) & 3];
      res[pos++] = nuctable[(qrt >> 2) & 3];
      res[pos++] = nuctable[qrt & 3];
   }
   qrt = s->elms[s->len / 2];
   i = s->len % 4;
   if (i)
      res[pos++] = nuctable[qrt >> 6];
   if (i > 1)
      res[pos++] = nuctable[(qrt >> 4) & 3];
   if (i > 2)
      res[pos++] = nuctable[(qrt >> 2) & 3];
   res[pos] = '\0';

   return res;
}

char *sequence_chars (sequence_tp * s, char *res)
{
   int t;
   t = s->type;

   if (t == SEQ_TYPE_2BIT)
      return sequence_chars_2b (s, res, cap_nuc_tab_2b);
   else
      return sequence_chars_4b (s, res, cap_nuc_tab_4b);
}


static void sequence4b_print_tab (sequence_tp * s, FILE * df, int ll, char *nuctable)
{
   int j, i, cc = 0;
   unsigned char qrt;
   if (!s)
      return;
   for (j = 0; j < s->len / 2; j++) {
      if (ll && (cc >= ll)) {
	 fprintf (df, "\n");
	 cc = 0;
      }
      qrt = s->elms[j];
      fprintf (df, "%c%c", nuctable[qrt >> 4], nuctable[qrt & 15]);
      cc += 2;
   }
   qrt = s->elms[s->len / 2];
   i = s->len % 2;
   if (i && ll && (cc >= ll))
      fprintf (df, "\n");
   if (i)
      fprintf (df, "%c", nuctable[qrt >> 4]);
   if (i > 1)
      fprintf (df, "%c", nuctable[(qrt >> 4) & 15]);
}

static void sequence_print_tab (sequence_tp * s, FILE * df, int ll, char *nuctable)
{
   int j, i, cc = 0;
   unsigned char qrt;
   if (!s)
      return;
   for (j = 0; j < s->len / 4; j++) {
      if (ll && (cc >= ll)) {
	 fprintf (df, "\n");
	 cc = 0;
      }
      qrt = s->elms[j];
      fprintf (df, "%c%c%c%c", nuctable[qrt >> 6], nuctable[(qrt >> 4) & 3], nuctable[(qrt >> 2) & 3], nuctable[qrt & 3]);
      cc += 4;
   }
   qrt = s->elms[s->len / 4];
   i = s->len % 4;
   if (i && ll && (cc >= ll))
      fprintf (df, "\n");
   if (i)
      fprintf (df, "%c", nuctable[qrt >> 6]);
   if (i > 1)
      fprintf (df, "%c", nuctable[(qrt >> 4) & 3]);
   if (i > 2)
      fprintf (df, "%c", nuctable[(qrt >> 2) & 3]);
}



static int GetSequenceType (sequence_tp * s)
{
   if ((s->type == SEQ_TYPE_2BIT) || (s->type == SEQ_TYPE_4BIT))
      return s->type;
#ifdef DEBUG
   fprintf (stderr, "CODING-ERROR in sequence.c got sequence with invalid type\n");
#endif
   handle_error ();
   return SEQ_TYPE_2BIT;
}

void sequence_print (sequence_tp * s, FILE * df, int ll, int lettertype)
{
   int t;
   t = GetSequenceType (s);
   if (lettertype == SEQ_PRINT_CAPS) {
      if (t == SEQ_TYPE_2BIT)
	 sequence_print_tab (s, df, ll, cap_nuc_tab_2b);
      else
	 sequence4b_print_tab (s, df, ll, cap_nuc_tab_4b);
   }
   else {
      if (t == SEQ_TYPE_2BIT)
	 sequence_print_tab (s, df, ll, nuc_tab_2b);
      else
	 sequence4b_print_tab (s, df, ll, nuc_tab_4b);
   }
}


void sequence_write_fasta (FILE * df, sequence_tp * s, char *name, int lettertype)
{
   if (name)
      fprintf (df, ">%s\n", name);
   else
      fprintf (df, ">sequence\n");
   sequence_print (s, df, 72, lettertype);
   fprintf (df, "\n");
}

void sequence_write_ts (FILE * df, sequence_tp * s, char *name, int lettertype)
{
   if (name)
      fprintf (df, "%s\t", name);
   else
      fprintf (df, "sequence\t");
   sequence_print (s, df, 0, lettertype);
   fprintf (df, "\n");
}

static unsigned char stab2b[4] = { 6, 4, 2, 0 };
static unsigned char stab4b[2] = { 4, 0 };

nuc_tp sequence_get_nuc (sequence_tp * s, int loc)
{
   int t;
   t = GetSequenceType (s);
   if (t == SEQ_TYPE_2BIT)
      return (s->elms[loc >> 2] >> stab2b[loc & 3]) & 3;
//      fprintf(stderr,"get nuc %d elm = %d cont = %x shift = %d val = %d \n",loc,loc>>1,s->elms[loc >>1],stab4b[loc&1],(s->elms[loc >>1]  >> stab4b[loc&1]) & 15 );
   return code4b_to_nuc[(s->elms[loc >> 1] >> stab4b[loc & 1]) & 15];
}

static int IsNuc (nuc_tp nuc, int t)
{
   if ((nuc < 0) || (nuc > 15) || ((t == SEQ_TYPE_2BIT) && (nuc > 3)))
      return 0;
   return 1;
}
static int CheckNuc (nuc_tp nuc, int t)
{
   int ok;
   ok = IsNuc (nuc, t);
   if (!ok) {
#ifdef DEBUG
      fprintf (stderr, "Illegal nucleotide %d encountered in sequence.c \n", nuc);
#endif
      handle_error ();
   }
   return ok;
}


void sequence_append_nuc (sequence_tp * s, nuc_tp nuc)
{
   int loc, t;
   loc = s->len;
   t = GetSequenceType (s);
   if (!CheckNuc (nuc, t))
      return;
   if (t == SEQ_TYPE_2BIT) {
      if ((loc >> 2) >= s->elm_alloced_bytes)
	 sequence_extend (s);
      s->elms[loc >> 2] &= mask_2b[loc & 3];
      s->elms[loc >> 2] |= nuc << (2 * (3 - (loc & 3)));
      s->len++;
   }
   else {
      if ((loc >> 1) >= (s->elm_alloced_bytes - 1))
	 sequence_extend (s);
      s->elms[loc >> 1] &= mask_4b[loc % 2];
      s->elms[loc >> 1] |= nuc_code_4b[nuc] << (4 * (1 - (loc & 1)));
      s->len++;
   }
}


void sequence_prepend_nuc (sequence_tp * s, nuc_tp nuc)
{
   unsigned char rm, rms;
   int i, t;
   t = GetSequenceType (s);
   if (!CheckNuc (nuc, t))
      return;
   if (t == SEQ_TYPE_2BIT) {
      if ((s->len >> 2) >= s->elm_alloced_bytes)
	 sequence_extend (s);
      rm = nuc;			/* using the 1-1 mapping */
      for (i = 0; i <= s->len >> 2; i++) {
	 rms = (s->elms[i]) & 3;
	 s->elms[i] = (s->elms[i] >> 2) | (rm << 6);
	 rm = rms;
      }
      s->len++;
   }
   else {
      if ((s->len >> 1) >= s->elm_alloced_bytes)
	 sequence_extend (s);
      rm = nuc_code_4b[nuc];
      for (i = 0; i <= s->len >> 1; i++) {
	 rms = (s->elms[i]) & 15;
	 s->elms[i] = (s->elms[i] >> 4) | (rm << 4);
	 rm = rms;
      }
      s->len++;
   }
}

void sequence_set_nuc (sequence_tp * s, int loc, nuc_tp nuc)
{
   int t;
   t = GetSequenceType (s);
   if ((loc < 0) || (loc >= s->len)) {
#ifdef DEBUG
      fprintf (stderr, "loc = %d len = %d\n", loc, s->len);
#endif
      seq_warning ("Not setting nucleotide outside sequence");
      handle_error ();
   }
   if (t == SEQ_TYPE_2BIT) {
      if (CheckNuc (nuc, t)) {
	 s->elms[loc >> 2] &= mask_2b[loc % 4];
	 s->elms[loc >> 2] |= nuc_code_2b[nuc] << (2 * (3 - (loc % 4)));
      }
   }
   else {
      s->elms[loc >> 1] &= mask_4b[loc % 2];
      s->elms[loc >> 1] |= nuc_code_4b[nuc] << 4 * (1 - (loc % 2));
   }
}

long int sequence_get_index (sequence_tp * s, int loc, int ml)
{
   int t, bits = 2;
   long int i, r = 0;
   t = GetSequenceType (s);
   if (t == SEQ_TYPE_4BIT)
      bits = 4;
   for (i = 0; i < ml; i++)
      r = (r << bits) | sequence_get_nuc (s, loc + i);
   return r;
}


int sequence_delete (sequence_tp * s)
{
   if (!s) {
      seq_warning ("Trying to free NULL-sequence");
      return -1;
   }
   if (s->elm_alloced_bytes)
      free (s->elms);
   free (s);
   return 1;
}

#define INITIAL_ALLOCATION 40

sequence_tp *sequence_new (int type)
{
   sequence_tp *n;
   n = (sequence_tp *) malloc (sizeof (sequence_tp));
   n->len = 0;
   n->type = type;
   n->elms = (unsigned char *) malloc (sizeof (char) * INITIAL_ALLOCATION);
   n->elm_alloced_bytes = INITIAL_ALLOCATION;
   return n;
}

int sequence_length (sequence_tp * s)
{
   return s->len;
}


sequence_tp *sequence_complement (sequence_tp * s)
{
   int i;
   sequence_tp *n;
   n = (sequence_tp *) malloc (sizeof (sequence_tp));
   n->elms = malloc (s->elm_alloced_bytes);
   n->elm_alloced_bytes = s->elm_alloced_bytes;
   n->type = GetSequenceType (s);
   n->len = 0;
   for (i = 0; i < s->len; i++) {
      sequence_append_nuc (n, comp_tab[sequence_get_nuc (s, i)]);
   }
   return n;
}


/* 
* Number of passes needed to generate the sequence as an oligo on a NimbleGen array  
* Input: sequence
* Output: number of cycles needed to generate the sequence
*   In case of an non-determinate nucleotide in the sequence -1 is returned;
*/
nuc_tp NimblegenNucOrder[4] = { A, C, G, T };

int sequence_Nimblegen_passes (sequence_tp * s)
{
   int i, passes = 0, cidx = 0;
   nuc_tp nn;
   for (i = s->len - 1; i >= 0; i--) {
      nn = sequence_get_nuc (s, i);
      if ((nn != A) && (nn != C) && (nn != G) && (nn != T))
	 return -1;
      while (nn != NimblegenNucOrder[cidx]) {
	 cidx = (cidx + 1) % 4;
	 passes++;
      }
      cidx = (cidx + 1) % 4;
      passes++;
   }
   return passes;
}


sequence_tp *sequence_reverse_complement (sequence_tp * s)
{
   int i;
   sequence_tp *n;
   n = (sequence_tp *) malloc (sizeof (sequence_tp));
   n->type = GetSequenceType (s);
   n->elms = malloc (s->elm_alloced_bytes);
   n->elm_alloced_bytes = s->elm_alloced_bytes;
   n->type = GetSequenceType (s);
   n->len = 0;
   for (i = 0; i < s->len; i++) {
      sequence_append_nuc (n, comp_tab[sequence_get_nuc (s, s->len - 1 - i)]);
   }
   return n;
}

long int sequence_content (sequence_tp * s, nuc_tp n)
{
   long int i, cs = 0;
   for (i = 0; i < s->len; i++) {
      if (n == sequence_get_nuc (s, i))
	 cs++;
   }
   return cs;
}

double sequence_cg_count (sequence_tp * s)
{
   int i;
   double cgs = 0;
   char n;
   for (i = 0; i < s->len; i++) {
      n = sequence_get_nuc (s, i);
      switch (n) {
      case C:
      case G:
      case S:
	 cgs += 1.0;
	 break;
      case A:
      case T:
      case W:
	 break;
      case B:
      case V:
	 cgs += 0.666666667;
	 break;
      case K:
      case M:
      case N:
      case R:
      case Y:
      case gap:
	 cgs += 0.5;
	 break;
      case D:
      case H:
	 cgs += 0.3333;
	 break;
      }
   }
   return cgs;
}

double sequence_cg_content (sequence_tp * s)
{
   return sequence_cg_count (s) / s->len;
}


/* Table 1 (Unified) /2 deltaH and deltaS NN parameters in 1 M NaCl  
	Santa Lucia in Proc.Natl.Acad Sci. USA 95 (1998)
	 	deltaG	deltaH	deltaS
   AA/TT	-1.00	 -7.9	-22.2   AA || TT
   AT/TA	-0.88	 -7.2	-20.4   AT
   TA/AT	-0.58	 -7.2	-21.3   TA
   CA/GT	-1.45	 -8.5	-22.7   CA || TG
   GT/CA	-1.44	 -8.4	-22.4   GT || AC
   CT/GA	-1.28	 -7.8	-21.4   CT || AG
   GA/CT	-1.30	 -8.2	-22.2   GA || TC
   CG/GC	-2.17	-10.6	-27.2   CG  
   GC/CG	-2.24	 -9.8	-24.4   GC
   GG/CC	-1.84	 -8.0	-19.9   GG || CC
*/

#define GAS_CONSTANT_R 1.987	/* [ cal/K.mol ] */
#define KELVIN_ZERO  -273.15

double dHsantaLuciaNN[4][4] =	/*  [ cal / mol ] */
{ {-7900, -8400, -7800, -7200},
{-8500, -8000, -10600, -7800},
{-8200, -9800, -8000, -8400},
{-7200, -8200, -8500, -7900}
};

double dHsantaLuciaInit[4] = { 2300, 100, 100, 2300 };

double dSsantaLuciaNN[4][4] =	/*  [ cal / K . mol ]  */
{ {-22.2, -22.4, -21.4, -20.4},
{-22.7, -19.9, -27.2, -21.4},
{-22.2, -24.4, -19.9, -22.4},
{-21.3, -22.2, -22.7, -22.2}
};

double dSsantaLuciaInit[4] = { 4.1, -2.8, -2.8, 4.1 };

double dGsantaLuciaNN[4][4] =	/* [ cal / mol ] */
{ {-1000, -1440, -1280, -880},
{-1450, -1840, -2170, -1280},
{-1300, -2240, -1840, -1440},
{-580, -1300, -1450, -1000}
};

double dGsantaLuciaInit[4] = { 1030, 980, 980, 1030 };


double sequence_melt_nn_SantaLucia (sequence_tp * s, double Ct)
{
   int i = 0;
   double dH, dS;
   double rv;
   char n1, n2;
   if (s->type != SEQ_TYPE_2BIT) {
      seq_warning
	 ("sequence_melt_nn_SantaLucia: can not determine melting temperature on sequences with undetermined nucleotides\n");
      return -1e10;
   }
   n1 = GET_NUC_2B (s, 0);
   dS = dSsantaLuciaInit[(int) n1];
   dH = dHsantaLuciaInit[(int) n1];
   for (i = 1; i < s->len; i++) {
      n2 = GET_NUC_2B (s, i);
      dS += dSsantaLuciaNN[(int) n1][(int) n2];
      dH += dHsantaLuciaNN[(int) n1][(int) n2];
      n1 = n2;
   }
   dS += dSsantaLuciaInit[(int) n1];
   dH += dHsantaLuciaInit[(int) n1];
   rv = dH / (dS + GAS_CONSTANT_R * log (Ct)) + KELVIN_ZERO;
   return rv;
}


double sequence_dg_nearest_neighbor (sequence_tp * s, double nnTable[4][4], double InitTable[4])
{
   int i = 0;
   double dG;
   char n1, n2;
   if (s->type != SEQ_TYPE_2BIT) {
      seq_warning
	 ("sequence_dg_nearest_neighbor: can not determine melting temperature on sequences with undetermined nucleotides\n");
      return 1e10;
   }
   if (!nnTable)
      nnTable = dGsantaLuciaNN;
   if (!InitTable)
      InitTable = dGsantaLuciaInit;
   n1 = GET_NUC_2B (s, 0);
   dG = InitTable[(int) n1];
   for (i = 1; i < s->len; i++) {
      n2 = GET_NUC_2B (s, i);
      dG += nnTable[(int) n1][(int) n2];
      n1 = n2;
   }
   dG += InitTable[(int) n1];
   return dG;
}


// Calculate longes stretch of  self dimerization
double sequence_dimerization_longest_match (sequence_tp * s, sequence_tp * q, int pos, int nwindow)
{
   int i, stq, sts, eq, es, ovl, seq_p, query_p, l, lmax = 0, w_s, w_e, mpos = 0;
   if (s->type != SEQ_TYPE_2BIT) {
      seq_warning
	 ("sequence_dimerization_longest_match: can not determine melting temperature on sequences with undetermined nucleotides\n");
      return 1e10;
   }
   stq = MAXV (0, -pos);
   eq = MINV (q->len, s->len - pos);
   sts = MAXV (0, pos);
   es = MINV (s->len, pos + q->len);
   ovl = es - sts;
   query_p = stq;
   seq_p = sts;
   l = 0;
   w_s = ovl / 2 - nwindow / 2;
   w_e = ovl / 2 + nwindow / 2;

   for (i = 0; i < ovl; i++) {
      if ((i <= w_s) || (i >= w_e)) {
	 if (GET_NUC_2B (q, query_p) == ((GET_NUC_2B (s, seq_p))))
	    l++;
	 else {
	    if (l > lmax) {
	       lmax = l;
	       mpos = i;
	    }
	    l = 0;
	 }
      }
      else {
	 l = 0;
      }
      query_p++;
      seq_p++;
   }
#if 0
   if (lmax > 15) {
      fprintf (stdout, "pos = %d\n", pos);
      seq_p = sts;
      query_p = stq;
      for (i = 0; i < -pos; i++)
	 fprintf (stdout, " ");
      sequence_print (s, stdout, 0);
      printf ("\n");
      for (i = 0; i < pos; i++)
	 fprintf (stdout, " ");
      sequence_print (q, stdout, 0);
      printf ("\n");
      for (i = 0; i < ABSV (pos); i++)
	 printf (" ");

      for (i = 0; i < ovl; i++) {
	 if ((i <= w_s) || (i >= w_e)) {
	    if (GET_NUC (q, query_p) == (GET_NUC (s, seq_p)))
	       printf ("%c", nuc_tab[GET_NUC (q, query_p)]);
	    else
	       printf (".");
	 }
	 else
	    printf ("X");
	 query_p++;
	 seq_p++;
      }
      fprintf (stdout, "\n len = %d AT pos %d \n", lmax, mpos);

   }
#endif
   return (double) lmax;
}

double sequence_self_dimerization_score (sequence_tp * s, int nwindow)
{
   int i;
   double v, vmax = 0;
   sequence_tp *query;
   query = sequence_reverse_complement (s);
   for (i = -query->len + 1 + nwindow; i < s->len - nwindow; i++) {
      v = sequence_dimerization_longest_match (s, query, i, nwindow);
      if (v > vmax)
	 vmax = v;
   }
   return vmax;
}

sequence_tp *sequence_copy (sequence_tp * s)
{
   sequence_tp *n;
   n = (sequence_tp *) malloc (sizeof (sequence_tp));
   n->len = s->len;
   n->type = GetSequenceType (s);
   n->elms = malloc (s->elm_alloced_bytes);
   n->elm_alloced_bytes = s->elm_alloced_bytes;
   memcpy (n->elms, s->elms, s->elm_alloced_bytes);
   return n;
}

static int append_character_to_sequence (char rc, sequence_tp * s)
{
   int t, added_nuc = 1;
   t = GetSequenceType (s);
   switch (rc) {
   case 'a':
   case 'A':
      sequence_append_nuc (s, A);
      break;
   case 'c':
   case 'C':
      sequence_append_nuc (s, C);
      break;
   case 'g':
   case 'G':
      sequence_append_nuc (s, G);
      break;
   case 't':
   case 'T':
      sequence_append_nuc (s, T);
      break;
   default:
      added_nuc = 0;
   }
   if ((t == SEQ_TYPE_2BIT) || added_nuc)
      return added_nuc;
   added_nuc = 1;
   switch (rc) {
   case 'n':
   case 'N':
      sequence_append_nuc (s, N);
      break;
   case 'm':
   case 'M':
      sequence_append_nuc (s, M);
      break;
   case 'r':
   case 'R':
      sequence_append_nuc (s, R);
      break;
   case 's':
   case 'S':
      sequence_append_nuc (s, S);
      break;
   case 'v':
   case 'V':
      sequence_append_nuc (s, V);
      break;
   case 'w':
   case 'W':
      sequence_append_nuc (s, W);
      break;
   case 'y':
   case 'Y':
      sequence_append_nuc (s, Y);
      break;
   case 'h':
   case 'H':
      sequence_append_nuc (s, H);
      break;
   case 'k':
   case 'K':
      sequence_append_nuc (s, K);
      break;
   case 'd':
   case 'D':
      sequence_append_nuc (s, D);
      break;
   case 'b':
   case 'B':
      sequence_append_nuc (s, B);
      break;
   case '-':
      sequence_append_nuc (s, gap);
      break;
   default:
      added_nuc = 0;
   }
   return added_nuc;
}


sequence_tp *sequence_read (FILE * df, int seqtype)
{
   char rc, ws[WARNING_STRING_LENGTH];
   int spos = 0, end_seq = 0;
   sequence_tp *s;
   s = sequence_new (seqtype);
   while (!end_seq && ((rc = fgetc (df)) != EOF)) {
      spos++;
      if (!append_character_to_sequence (rc, s))
	 switch (rc) {
	 case '\n':
	 case '\r':
	 case '*':
	    if (s->len > 0)
	       end_seq = 1;
	    break;
	 case ' ':
	 case '\t':
	    break;
	 default:
	    sprintf (ws, "Skipping unrecognized character in sequence [%c] at position %d  ", rc, spos);
	    seq_warning (ws);
	 }

   }
   if (s->len == 0) {
      sequence_delete (s);
      return NULL;

   }
   return s;
}

#define FNA_STATE_START       0
#define FNA_STATE_COMMENT     1
#define FNA_STATE_SEQUENCE    2
#define FNA_STATE_STOP        4

sequence_tp *sequence_read_fasta (FILE * df, char *comment, int seqtype)
{
   int com_len = 0, state = FNA_STATE_START, spos;
   char rc, ws[WARNING_STRING_LENGTH];
   comment[com_len] = 0;
   sequence_tp *s;
   while (state == FNA_STATE_START) {
      rc = getc (df);
      switch (rc) {
      case '>':
	 state = FNA_STATE_COMMENT;
	 break;
      case '\n':
      case '\t':
      case ' ':
      case '\r':
	 break;
      case EOF:
	 state = FNA_STATE_STOP;
	 break;
      default:
	 state = FNA_STATE_SEQUENCE;
      }
   }
   if (state == FNA_STATE_STOP)
      return NULL;
   s = sequence_new (seqtype);
   spos = 0;
   while (state != FNA_STATE_STOP) {
      switch (state) {
      case FNA_STATE_COMMENT:
	 switch (rc) {
	 case EOF:
	    state = FNA_STATE_STOP;
	    break;
	 case '\n':
	 case '\r':
	    state = FNA_STATE_SEQUENCE;
	    comment[com_len] = 0;
	    break;
	 default:
	    if ((com_len < SEQUENCE_MAX_COMMENT_LEN - 1) && ((com_len > 0) || rc != '>'))
	       comment[com_len++] = rc;
	 }
	 break;
      case FNA_STATE_SEQUENCE:
	 spos++;
//                      fprintf(stderr,"%d got char %d : %c \n",spos,rc,rc);
	 if (!append_character_to_sequence (rc, s))
	    switch (rc) {
	    case '>':
	       state = FNA_STATE_STOP;
	       ungetc (rc, df);
	       break;

	    case EOF:
	       state = FNA_STATE_STOP;
	       break;
	    case '*':
	       state = FNA_STATE_STOP;
	       break;
	    case '\n':
	    case '\r':
	    case ' ':
	    case '\t':
	       spos--;
	       break;
	    default:
	       sprintf (ws, "While reading sequence: skipping unrecognized character [%c] in sequence %s at position %d", rc, comment, spos);
	       seq_warning (ws);
	    }
	 break;
      default:
	 seq_warning ("Unknown state in fasta-read");
      }
      if (!(state == FNA_STATE_STOP))
	 rc = getc (df);
   }
   comment[com_len] = 0;
   return s;
}

sequence_tp *sequence_from_string (const char *seq)
{
   sequence_tp *s;
   char ws[WARNING_STRING_LENGTH];
   int rc = 0;
   s = sequence_new (SEQ_TYPE_4BIT);
   while (*seq) {

      if (!append_character_to_sequence (*seq, s)) {
	 sprintf (ws, "Get sequence: skipping unrecognized character [%c] at position %d", rc, *seq);
	 seq_warning (ws);
      }
      rc++;
      seq++;
   }
   return s;
}

sequence_tp *sequence_read_ts (FILE * df, char *comment, int seqtype, int id_col, int seq_col, char sep)
{
   int have_id = 0, have_seq = 0, cur_col = 1;
   int com_len = 0, state = FNA_STATE_START;
   char rc, ws[120];
   comment[com_len] = 0;
   sequence_tp *s;
   s = sequence_new (seqtype);
   if (id_col == seq_col)
      return NULL;
   state = FNA_STATE_START;
   if (id_col == 1)
      state = FNA_STATE_COMMENT;
   if (seq_col == 1)
      state = FNA_STATE_SEQUENCE;

   while (state != FNA_STATE_STOP) {
      rc = getc (df);
      //      fprintf(stderr,"Got [%c]\n",rc);
      if ((rc == EOF) || (rc == '\n') || (rc == '\r')) {
	 state = FNA_STATE_STOP;
      }
      else if (rc == sep) {
	 cur_col++;
	 if (cur_col == id_col)
	    state = FNA_STATE_COMMENT;
	 else if (cur_col == seq_col)
	    state = FNA_STATE_SEQUENCE;
	 else
	    state = FNA_STATE_START;
//                      fprintf(stderr,"Got separator %d %d\n",cur_col,state);
      }
      else {
	 switch (state) {
	 case FNA_STATE_COMMENT:
	    //                      fprintf(stderr,"Adding to comment\n");
	    if (com_len < SEQUENCE_MAX_COMMENT_LEN) {
	       comment[com_len++] = rc;
	       have_id = 1;
	    }
	    break;
	 case FNA_STATE_SEQUENCE:
	    //                      fprintf(stderr,"Adding to sequence\n");
	    if (!append_character_to_sequence (rc, s)) {
	       sprintf (ws, "Skipping unrecognized character in sequence [%c] ", rc);
	       seq_warning (ws);
	    }
	    else {
	       have_seq = 1;
	    }
	    break;
	 }
      }
   }
   comment[com_len] = 0;
   if (have_seq && ((id_col == 0) || (have_id)))
      return s;
   else
      free (s);
   return NULL;
}

int sequence_set_nucs (sequence_tp * s, int loc_dest, int loc_src, int len_src, unsigned char *sequence)
{
   int i, end, byte;
   nuc_tp nucv;
   if ((loc_dest < 0) || (loc_dest + len_src > s->len)) {
      seq_warning ("Cannot write nuc outside sequence\n");
      return -1;
   }
   i = loc_src;
   end = loc_src + len_src;
   while (i < end) {
      byte = i >> 2;
      nucv = ((*(sequence + byte)) >> (2 * (3 - (i & 3)))) & 3;
      sequence_set_nuc (s, loc_dest++, nucv);
      i++;
   }
   return 1;

}

sequence_tp *sequence_select_part (sequence_tp * s, int start, int end)
{

   sequence_tp *n;
   int isize, i, t;
   if ((start < 1) || (start > s->len) || (end < 1) || (end > s->len))
      return NULL;
   n = (sequence_tp *) malloc (sizeof (sequence_tp));
   n->len = 0;
   t = GetSequenceType (s);
   if (t == SEQ_TYPE_2BIT)
      isize = abs (end - start) / 4 + 1;
   else
      isize = abs (end - start) / 2 + 1;
   n->type = t;
   n->elms = malloc (isize);
   n->elm_alloced_bytes = isize;
   if (end >= start) {
      for (i = start; i <= end; i++)
	 sequence_append_nuc (n, sequence_get_nuc (s, i - 1));
   }
   else {
      for (i = start; i >= end; i--) {
	 sequence_append_nuc (n, comp_tab[sequence_get_nuc (s, i - 1)]);
      }
   }
   return n;
}



long int sequence_find_sequence (sequence_tp * h, sequence_tp * n, long int start)
{
   long int tp, sp;
   tp = start - 1;
   if (tp < 0)
      return -1;

   while (tp < sequence_length (h) - sequence_length (n) + 1) {
      sp = 0;
      while (sp < sequence_length (n) && (sequence_get_nuc (h, tp + sp) == sequence_get_nuc (n, sp)))
	 sp++;
      if (sp == sequence_length (n)) {
	 return tp + 1;
      }
      tp++;
   }
   return -1;
}

void *sequence_get_restriction_iterator (sequence_tp * s, sequence_tp * restr, int cut_pos)
{
   restriction_iterator *ri;
   ri = malloc (sizeof (restriction_iterator));
   ri->restriction_sequence = sequence_copy (restr);
   ri->restriction_sequence_rc = sequence_reverse_complement (restr);
   ri->cut_pos = cut_pos;
   ri->current_location = 1;
   ri->end = 0;
   ri->sequence = s;
   return ri;

}


sequence_tp *sequence_get_restriction_next (void *i)
{
   long int id1, id2, old_loc;
   restriction_iterator *ri;
   ri = (restriction_iterator *) i;
   if (ri->end)
      return NULL;
   old_loc = ri->current_location;
   id1 = sequence_find_sequence (ri->sequence, ri->restriction_sequence, old_loc);
   id2 = sequence_find_sequence (ri->sequence, ri->restriction_sequence_rc, old_loc);
   if ((id1 == -1) && (id2 == -1)) {
      ri->end = 1;
      return sequence_select_part (ri->sequence, old_loc, sequence_length (ri->sequence));
   }
   else if (id1 == -1)
      ri->current_location = id2 + sequence_length (ri->restriction_sequence) - ri->cut_pos;
   else if (id2 == -1)
      ri->current_location = id1 + ri->cut_pos;
   else if (id1 < id2)
      ri->current_location = id1 + ri->cut_pos;
   else
      ri->current_location = id2 + sequence_length (ri->restriction_sequence) - ri->cut_pos;
   return sequence_select_part (ri->sequence, old_loc, ri->current_location);
}

void sequence_free_restriction_iterator (void *i)
{
   free (i);
}



void sequence_append (sequence_tp * s, sequence_tp * t)
{
   int i, l;
   l = sequence_length (t);
   for (i = 0; i < l; i++)
      sequence_append_nuc (s, sequence_get_nuc (t, i));
}

void sequence_append_nucs (sequence_tp * s, int loc_src, int len_src, unsigned char *sequence)
{
   unsigned char n;
   int i, end, byte;
   nuc_tp nucv;
   i = loc_src;
   end = loc_src + len_src;
   while (i < end) {
      byte = i >> 2;
      nucv = (*(sequence + byte)) << (((3 - (i & 3)) << 2));
      sequence_append_nuc (s, nucv);
      n++;
   }
}

void sequence_set_len (sequence_tp * s, int len)
{
   unsigned char *new;
   int ns, t;
   t = GetSequenceType (s);
   if (t == SEQ_TYPE_2BIT)
      ns = (len + 3) >> 2;
   else
      ns = (len + 3) >> 1;
   if (ns > s->elm_alloced_bytes) {
      new = malloc (ns);
      if (s->elm_alloced_bytes) {
	 memcpy (new, s->elms, s->elm_alloced_bytes);
	 free (s->elms);
      }
      s->elms = new;
      s->elm_alloced_bytes = ns;
   }
   s->len = len;
}

int sequence_conv_to_acgt_only (sequence_tp * s)
{
   int t, i;
   sequence_tp *tns;
   nuc_tp n;
   t = GetSequenceType (s);
   if (t == SEQ_TYPE_2BIT)
      return 1;
   tns = sequence_new (SEQ_TYPE_2BIT);
   for (i = 0; i < s->len; i++) {
      n = sequence_get_nuc (s, i);
      if (IsNuc (n, SEQ_TYPE_2BIT)) {
	 sequence_append_nuc (tns, n);
      }
      else {
	 sequence_free (tns);
	 return 0;
      }
   }
   s->type = SEQ_TYPE_2BIT;
   s->elm_alloced_bytes = tns->elm_alloced_bytes;
   free (s->elms);
   s->elms = tns->elms;
   free (tns);
   return 1;
}



void sequence_set_nucs_long (sequence_tp * s, int loc, int ilen, long int sequence)
{
   unsigned char n;
   int i;
   for (i = ilen - 1; i >= 0; i--) {
      n = sequence & 3;
      sequence_set_nuc (s, loc + i, n);
      sequence = sequence >> 2;
   }
}

/* creates a string of characters for a sequence */
char *sequence_to_string (sequence_tp * s)
{
   int i;
   char *st;
   st = (char *) malloc ((1 + sequence_length (s)) * sizeof (char));
   for (i = 0; i < s->len; i++) {
      st[i] = nuc_tab_enum[sequence_get_nuc (s, i)];
      /*fprintf(stdout,".. %c\n",st[i]); */
   }
   st[i] = (char) 0;
   return st;
}

double sequence_hybrid_ss_min (sequence_tp * s, double folding_temp, double naConc, double mgConc, int acid_type)
{

   char *cseq;
   cseq = sequence_to_string (s);
   double dGhybrid = hybrid_main (cseq, folding_temp, naConc, mgConc, acid_type);	/* calls in hybrid-ss-min.c */
   free (cseq);
   return dGhybrid;
}
