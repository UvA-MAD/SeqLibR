#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <stdio.h>
#define SEQUENCE_MAX_COMMENT_LEN 1024
#define SEQ_TYPE_2BIT 0		/* ACGT only */
#define SEQ_TYPE_4BIT 1		/* all fasta  characters */
#define SEQ_PRINT_LOWER 0	/* sequence_print lettertype */
#define SEQ_PRINT_CAPS  1

#define ACID_TYPE_RNA   0
#define ACID_TYPE_DNA   1

enum nucleotides
{ A, C, G, T, N, R, S, V, W, Y, H, K, D, B, M, gap };
typedef enum nucleotides nuc_tp;

typedef struct sequence_tp
{
   int len;
   int type;
   unsigned char *elms;
   int elm_alloced_bytes;
} sequence_tp;


void sequence_free (sequence_tp *);
void sequence_print (sequence_tp * s, FILE *, int ll, int lettertype);
char *sequence_chars (sequence_tp * s, char *str);
void nucleotide_print (nuc_tp n, FILE *);
char nuc_char (nuc_tp n);

void sequence_set_nuc (sequence_tp * s, int loc, nuc_tp nuc);
nuc_tp sequence_get_nuc (sequence_tp * s, int loc);

void sequence_append_nuc (sequence_tp * s, nuc_tp nuc);
void sequence_append (sequence_tp * s, sequence_tp * t);
void sequence_prepend_nuc (sequence_tp * s, nuc_tp nuc);

int sequence_set_nucs (sequence_tp * s, int loc_dest, int loc_src, int len_src, unsigned char *sequence);
void sequence_append_nucs (sequence_tp * s, int loc_src, int len_src, unsigned char *sequence);

long int sequence_get_index (sequence_tp * s, int loc, int ml);
void sequence_set_nucs_long (sequence_tp * s, int loc, int ilen, long int sequence);
double sequence_self_dimerization_score (sequence_tp * s, int nwindow);
void sequence_set_len (sequence_tp * s, int len);
int sequence_delete (sequence_tp * s);
int sequence_length (sequence_tp * s);
long int sequence_content (sequence_tp * s, nuc_tp n);
double sequence_cg_content (sequence_tp * s);
double sequence_cg_count (sequence_tp * s);
double sequence_melt_nn_SantaLucia (sequence_tp * s, double Ct);
int sequence_Nimblegen_passes (sequence_tp * s);
double sequence_hybrid_ss_min (sequence_tp * s, double folding_temp, double naConc, double mgConc, int acid_type);
int sequence_conv_to_acgt_only (sequence_tp * s);

sequence_tp *sequence_read (FILE * df, int seqtype);
sequence_tp *sequence_read_fasta (FILE * df, char *comment, int seqtype);
sequence_tp *sequence_read_ts (FILE * df, char *comment, int seqtype, int id_col, int seq_col, char sep);
sequence_tp *sequence_new (int);
sequence_tp *sequence_copy (sequence_tp * s);
sequence_tp *sequence_complement (sequence_tp * s);
sequence_tp *sequence_reverse_complement (sequence_tp * s);
sequence_tp *sequence_select_part (sequence_tp * s, int start, int end);
sequence_tp *sequence_from_string (const char *);

char *sequence_to_string (sequence_tp * s);

void sequence_write_fasta (FILE * df, sequence_tp * s, char *name, int lettertype);
void sequence_write_ts (FILE * df, sequence_tp * s, char *name, int lettertype);

void *sequence_get_restriction_iterator (sequence_tp * s, sequence_tp * restr, int cut_pos);
sequence_tp *sequence_get_restriction_next (void *i);
void sequence_free_restriction_iterator (void *i);

long int sequence_find_sequence (sequence_tp * h, sequence_tp * n, long int start);


void sequence_catch_warning (void (*my_warn) (char *));

#endif
