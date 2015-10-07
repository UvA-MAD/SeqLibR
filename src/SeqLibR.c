#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include "sequence.h"

#define DEBUG 0

typedef struct oligoDesignParms {
	int     oligo_len;
	int     NimbleGenMaxPass;
	int     nb_probes; 
	int     side_bias;
	double  tm_min;
	double  tm_max;
	double  cg_min;
	double  cg_max;
	sequence_tp **forbidden_seqs; 
	} oligoDesignParms;


SEXP seqlib_tm_santa_lucia(SEXP sequences,SEXP ct)
{
	int vlen,i,wg=0;
	sequence_tp*ms;
	double Ct;
	SEXP res;
	if(!isString(sequences))
	   error("sequence must have character type");
	if (!isReal(ct) || length(ct) != 1)
	   error("ct value must be single real");
	vlen = length(sequences);
	Ct = REAL(ct)[0];
	PROTECT(res = NEW_NUMERIC(vlen));
	for (i=0; i< vlen; i++)
	{
		ms = sequence_from_string(CHAR(STRING_ELT(sequences,i)));
		if (sequence_conv_to_acgt_only(ms))
			NUMERIC_POINTER(res)[i] = sequence_melt_nn_SantaLucia(ms,Ct);
		else 
		{
			NUMERIC_POINTER(res)[i] = NA_REAL;
			if (!wg) 
			{
				warning("Non-determined nucleotides in sequences");
				wg = 1;
			}
		}
		free(ms);
	}
	UNPROTECT(1);
	return res;
}




static SEXP reftab;

static int my_blast_elm_cmp(const void *e1,const void *e2)
{
	int v1,v2;
	v1 = *(int*)e1;
	v2 = *(int*)e2;
	return INTEGER(reftab)[v1] - INTEGER(reftab)[v2];
}

static void sort_cut_list(SEXP seqno,SEXP qstart,SEXP qstop,int *sorttab)
{
	int i_start,i_end,nb_cuts;
	i_start=0;
	i_end=0;
	nb_cuts = length(seqno);
	reftab = qstart;
	while (i_start < nb_cuts) {
		while ( (i_end < nb_cuts) && (INTEGER(seqno)[i_start] == INTEGER(seqno)[i_end]) )
			i_end++;
		qsort(sorttab+i_start,(i_end-i_start),sizeof(int),my_blast_elm_cmp);
		i_start = i_end;
	}
}

SEXP sequence_clean(SEXP seqs,SEXP seqno,SEXP qstart,SEXP qstop,SEXP min_len)
{
	sequence_tp  **results,*ms;
	int  i,nb_cuts,nb_seqs,cidx,nb_res,qs,qe,cp,ml,maxlen = 0,seq_name_len,pcnt,*sorttab;
	char *seqstring,**rnames;
	const char *seq_name;
	SEXP res,names;
#if DEBUG
	fprintf(stderr,"Cleaning\n");
#endif
	if(!isString(seqs))
	   error("sequence must have character type");
	if (!isInteger(seqno))
	   error("seqno value must be integer");
	if (!isInteger(qstart))
	   error("start-values must be integer");
	if (!isInteger(qstop))
	   error("start-values must be integer");
	if (!isInteger(min_len) || (length(min_len) != 1))
	   error("start-values must be a single integer");
	nb_cuts = length(seqno);
	if ((length(qstop) != nb_cuts) || (length(qstart) != nb_cuts))
	   error("non-equal number of sequence references and start/stop-positions.");
	nb_seqs = length(seqs); 
	results = (sequence_tp**)malloc(sizeof(sequence_tp*)*nb_seqs*3);
	rnames  = (char**)malloc(sizeof(char*)*nb_seqs*3);
	nb_res = 0;
	cidx = 0;
	ml = INTEGER(min_len)[0];
	for (i=0; i<nb_cuts; i++)
		INTEGER(seqno)[i]--;
#if DEBUG
	fprintf(stderr,"Sorting \n");	
#endif
	sorttab = (int*)malloc(sizeof(int)*nb_cuts);
	for (i=0; i<nb_cuts; i++) sorttab[i] = i;
	
	sort_cut_list(seqno,qstart,qstop,sorttab);
#if DEBUG
	fprintf(stderr,"Done sorting \n");	
	
#endif
	for (i=0; i<nb_seqs; i++)
	{
		pcnt=1;
		seq_name = CHAR(STRING_ELT(getAttrib(seqs,R_NamesSymbol), i));
		seq_name_len = strlen(seq_name);
		ms = sequence_from_string(CHAR(STRING_ELT(seqs,i)));
#if DEBUG
		fprintf(stderr,"Doing sequence %d/%d: (1-%d)\n",i,nb_seqs,sequence_length(ms));
#endif
		if ((cidx < nb_cuts) && (INTEGER(seqno)[cidx] == i))
		{
			cp = 1;
			while ((INTEGER(seqno)[cidx] == i) && (cidx < nb_cuts))
			{
				qs = INTEGER(qstart)[sorttab[cidx]];
				qe = INTEGER(qstop)[sorttab[cidx]];
				if ((qs -cp) > ml) 
				{
					results[nb_res] = sequence_select_part(ms,cp,qs-1);
					rnames[nb_res] = (char*)malloc(sizeof(char)*(seq_name_len+5));
					sprintf(rnames[nb_res],"%s_p%d",seq_name,pcnt++);
#if DEBUG
					fprintf(stderr,"Add result %d : %d (%d-%d)\n",nb_res,i,cp,qs-1);
#endif
					nb_res++;
				}
				if ((qe+1) > cp) cp = qe+1;
				cidx++;
			}
			if ((sequence_length(ms)-cp) >= ml)
			{
				results[nb_res] = sequence_select_part(ms,cp,sequence_length(ms));
				rnames[nb_res] = (char*)malloc(sizeof(char)*(seq_name_len+5));
				sprintf(rnames[nb_res],"%s_p%d",seq_name,pcnt);
#if DEBUG
				fprintf(stderr,"Add tail %d : %d (%d-%d)\n",nb_res,i,cp,sequence_length(ms));
#endif
				nb_res++;
			}
			sequence_free(ms);
		}
		else
		{
			results[nb_res] = ms;
#if DEBUG
			fprintf(stderr,"Not cutting sequence %d\n",i);
#endif
			rnames[nb_res] = (char*)malloc(sizeof(char)*(seq_name_len+5));
			sprintf(rnames[nb_res],"%s",seq_name);
			nb_res++;
		}
	}
	for (i=0; i<nb_res; i++)
		if (sequence_length(results[i]) > maxlen) { maxlen = sequence_length(results[i]);}	
	maxlen+=10;
	PROTECT(res = allocVector(STRSXP,nb_res));
	PROTECT(names = allocVector(VECSXP, nb_res));
	seqstring= (char*)malloc(sizeof(char)*(maxlen));
	for (i=0; i< nb_res; i++)
	{
#if DEBUG
		fprintf(stderr,"setting sequence %d \n",i);
#endif
		SET_STRING_ELT(res,i,mkChar(sequence_chars(results[i],seqstring)));
		SET_VECTOR_ELT(names,i,mkChar(rnames[i]));
	}
	free(seqstring);
	free(sorttab);
	for (i=0; i<nb_res; i++)
	{
		free(rnames[i]);
		sequence_free(results[i]);
	}
	free(results);
	free(rnames);
	setAttrib(res, R_NamesSymbol, names);
	UNPROTECT(2);
#if DEBUG
	fprintf(stderr,"Done cleaning\n");
#endif
	return res;
}

SEXP seqlib_nimblegen_cycles(SEXP sequences)
{
	int vlen,i;
	sequence_tp*ms;
	SEXP res;
	if(!isString(sequences))
	   error("sequence must have character type");
	vlen = length(sequences);
	PROTECT(res = NEW_NUMERIC(vlen));
	for (i=0; i< vlen; i++)
	{
		ms = sequence_from_string(CHAR(STRING_ELT(sequences,i)));
		NUMERIC_POINTER(res)[i] = sequence_Nimblegen_passes(ms);
		free(ms);
	}
	UNPROTECT(1);
	return res;
}


SEXP seqlib_cg_count(SEXP sequences)
{
	int vlen,i;
	sequence_tp*ms;
	SEXP res;
	if(!isString(sequences))
	   error("sequence must have character type");
	vlen = length(sequences);
	PROTECT(res = NEW_NUMERIC(vlen));
	for (i=0; i< vlen; i++)
	{
		ms = sequence_from_string(CHAR(STRING_ELT(sequences,i)));
		NUMERIC_POINTER(res)[i] = sequence_cg_count(ms);
		free(ms);
	}
	UNPROTECT(1);
	return res;
}

SEXP seqlib_dg_selfhyb(SEXP sequences,SEXP t,SEXP naC,SEXP mgC,SEXP acidtype)
{
	int vlen,i,atype;
	double ct,nac,mgc,rv;
	sequence_tp *ms;
	SEXP res;
	if(!isString(sequences))
	   error("sequence must have character type");
	vlen = length(sequences);
	PROTECT(res = NEW_NUMERIC(vlen));

	if (!isReal(t) || length(t) != 1)
              error("t value must be single real");
        ct = REAL(t)[0];
	if (!isReal(naC) || length(naC) != 1)
              error("naC value must be single real");
        nac = REAL(naC)[0];
	if (!isReal(mgC) || length(mgC) != 1)
              error("mgC value must be single real");
        mgc = REAL(mgC)[0];
	if (!isInteger(acidtype) || length(acidtype) != 1)
              error("acidtype value must be single int");
        atype = INTEGER(acidtype)[0];
	for (i=0; i< vlen; i++)
	{
		ms = sequence_from_string(CHAR(STRING_ELT(sequences,i)));
		rv = sequence_hybrid_ss_min(ms,ct,nac,mgc,atype);
		NUMERIC_POINTER(res)[i] = rv;
		free(ms);
	}
	UNPROTECT(1);
	return res;
}


SEXP seqlib_cg_content(SEXP sequences)
{
	int vlen,i;
	sequence_tp*ms;
	SEXP res;
	if(!isString(sequences))
	   error("sequence must have character type");
	vlen = length(sequences);
	PROTECT(res = NEW_NUMERIC(vlen));
	for (i=0; i< vlen; i++)
	{
		ms = sequence_from_string(CHAR(STRING_ELT(sequences,i)));
		NUMERIC_POINTER(res)[i] = sequence_cg_content(ms);
		free(ms);
	}
	UNPROTECT(1);
	return res;
}

SEXP seqlib_rev_comp(SEXP sequences)
{
	int vlen,i;
	sequence_tp *ms, *rms;
	char *seqstring;
	SEXP res;
	if(!isString(sequences))
	   error("sequence must have character type");
	vlen = length(sequences);
	PROTECT(res = allocVector(STRSXP,vlen));
	for (i=0; i< vlen; i++)
	{
		ms = sequence_from_string(CHAR(STRING_ELT(sequences,i)));
		rms = sequence_reverse_complement(ms);
		seqstring= (char*)malloc(sizeof(char)*(2+sequence_length(rms)));
		SET_STRING_ELT(res,i,mkChar(sequence_chars(rms,seqstring)));
		free(seqstring);
		free(rms);
		free(ms);
	}
	UNPROTECT(1);
	return res;
}

SEXP seqlib_gen_tile(SEXP Rsequence,SEXP Rlen,SEXP Rstep,SEXP Rcircular)
{
	sequence_tp *ms,*oligo;
	int cnt,pos,len,olen=60,step=1,circular=0,i=0;
	SEXP res;
	char *seqstring;
	if(!isString(Rsequence) || length(Rsequence) != 1)
	   error("gen_tile: sequence is not a single string");
	if (!isInteger(Rlen) || length(Rlen) != 1)
	   error("gen_tile: len value must be single int");
	if (!isInteger(Rstep) || length(Rstep) != 1)
	   error("gen_tile: step value must be single int");
	if (!isInteger(Rcircular) || length(Rcircular) != 1)
	   error("gen_tile: circular value must be single int");
	olen = INTEGER(Rlen)[0];
	step = INTEGER(Rstep)[0];
	circular = INTEGER(Rcircular)[0];
	if (olen < 1)
	   error("gen_tile: length value must be 1 or higher");
	if (step < 1)
	   error("gen_tile: step value must be 1 or higher");

	ms = sequence_from_string(CHAR(STRING_ELT(Rsequence,i)));
	len = sequence_length(ms);
	if (len < olen)
	{
	    sequence_free(ms);
	    error("gen_tile: tile len must be smaller than sequence length");
	}
	if (circular) 
		cnt = 1+ (len-1)/step;
	else
		cnt = 1+(len-olen)/step;
	pos = 1;
	seqstring= (char*)malloc(sizeof(char)*olen);
	PROTECT(res = allocVector(STRSXP,cnt));
	while (pos+olen-1 <= len)
	{
		oligo = sequence_select_part(ms,pos,pos+olen-1);
		SET_STRING_ELT(res,i,mkChar(sequence_chars(oligo,seqstring)));
		sequence_free(oligo);
		pos = pos+step;
		i++;
	}
	if (circular)
	{
		while (pos <= len)
		{
			oligo = sequence_select_part(ms,pos,len);
			sequence_append(oligo,sequence_select_part(ms,1,olen-(len-pos)-1));
			SET_STRING_ELT(res,i,mkChar(sequence_chars(oligo,seqstring)));
			sequence_free(oligo);
			pos = pos+step;
			i++;
		}
	}
	free(seqstring);
	UNPROTECT(1);
	return res;
}

void warning_stub(char *warnstring)
{
	warning(warnstring);
}

SEXP seqlib_read_fasta(SEXP filename)
{
	FILE *df;
	int i,cnt,sl,maxlen=0;
	sequence_tp *ds;
	SEXP res,names;
	char comment[SEQUENCE_MAX_COMMENT_LEN], *seqstring;
	if(!isString(filename) || length(filename) != 1)
	   error("filename is not a single string");

	df = fopen(CHAR(STRING_ELT(filename, 0)),"r");
	if (!df) 
	{
	   	error("can not open file");
	}
	cnt =0;
	sequence_catch_warning(warning_stub);
	while ((ds=sequence_read_fasta(df,comment,SEQ_TYPE_4BIT)))
	{
		sl = sequence_length(ds);
		if (sl > maxlen)
			maxlen = sl;
		sequence_free(ds);
		cnt++;
	}
	if (cnt)
	{
		rewind(df);
		seqstring= (char*)malloc(sizeof(char)*(maxlen+1));
		PROTECT(res = allocVector(STRSXP,cnt));
		PROTECT(names = allocVector(VECSXP, cnt));
		for (i=0; i< cnt; i++)
		{
			ds=sequence_read_fasta(df,comment,SEQ_TYPE_4BIT);
			SET_STRING_ELT(res,i,mkChar(sequence_chars(ds,seqstring)));
			SET_VECTOR_ELT(names,i,mkChar(comment));
		}
		free(seqstring);
		setAttrib(res, R_NamesSymbol, names);
		UNPROTECT(2);
	}
	else
	{
	   	error("can not open file");
	}
	fclose(df);
	return res;
}

static oligoDesignParms *NewOligoDesignParmRec()
{
	oligoDesignParms *nrec;
	nrec = (oligoDesignParms*)malloc(sizeof(oligoDesignParms));
	nrec->oligo_len = 60;
	nrec->nb_probes = 1;
	nrec->side_bias = 2;
	nrec->tm_min = 0;
	nrec->tm_max = 100;
	nrec->cg_min = 0.0;
	nrec->cg_max = 1.0;
	nrec->NimbleGenMaxPass = 0;
	nrec->forbidden_seqs  = NULL; 
}

static int ContainsForbidden(sequence_tp *s,sequence_tp **fseqs)
{
	int i=0;
	if (!fseqs)
		return 0;
	while (fseqs[i])
	{
		if (sequence_find_sequence(s,fseqs[i],1) != -1)
			return 1;
		i++;
	}
	return 0;
}

static int CheckProbe(sequence_tp *pr,oligoDesignParms *parms)
{
	int Nimpass,forbidden;
	double cgc,tm;
	if (!sequence_conv_to_acgt_only(pr))
		return 0;
	cgc = sequence_cg_content(pr);
	tm = sequence_melt_nn_SantaLucia(pr,1e-12);
	Nimpass = sequence_Nimblegen_passes(pr);
	forbidden = ContainsForbidden(pr,parms->forbidden_seqs);
	if ((cgc >= parms->cg_min) && (cgc <= parms->cg_max) && (tm >=parms->tm_min) && (tm <=parms->tm_max) && ((parms->NimbleGenMaxPass==0) || (Nimpass <= parms->NimbleGenMaxPass)) && !forbidden)
		return 1;
	return 0;
}

static int GenerateProbes(sequence_tp *gen,char *name,oligoDesignParms *parms)
{
	sequence_tp*sel;
	int nbpos,*is_ok,nb_ok=0,nb_ok3p=0,nb_ok5p=0,i,is,pcnt=0,mm;
	double pos,step;
	nbpos = sequence_length(gen)-parms->oligo_len+1;
	is_ok = malloc((nbpos+1)*sizeof(int));
	for (i = 1; i <= nbpos; i++)
	{
		sel = sequence_select_part(gen,i,i+(parms->oligo_len-1));
		if (CheckProbe(sequence_select_part(gen,i,i+(parms->oligo_len-1)),parms))
		{
			is_ok[nb_ok] = i;
			nb_ok++;
			if (i <= nbpos/2) nb_ok5p++; else nb_ok3p++;
		}
	}
	if ((parms->nb_probes == 0) || (nb_ok <= parms->nb_probes)  || ((parms->side_bias == 0) && (nb_ok5p <= parms->nb_probes)) || ((parms->side_bias == 1) && (nb_ok3p <= parms->nb_probes)) )
	{
		if (parms->nb_probes == 0)
			mm = nb_ok;
		else
			mm = parms->nb_probes;
		if (parms->side_bias==1) i =nb_ok-1; else i = 0;
		if (parms->side_bias==1) is = -1;  else is = 1;
		while ((i<mm) && (i >= 0) && (i <nb_ok))
		{
/*			PrintProbe(sequence_select_part(gen,is_ok[i],is_ok[i]+(parms->oligo_len-1)),name,is_ok[i]); */
			i += is; pcnt++;
		}
	}
	else  /* More than maximum probes possible generate subset */
	{
		if (parms->side_bias == 2) /* uniformly distributed probes */
		{
			step = (double)nb_ok/parms->nb_probes;
			pos = 0;
		}
		if (parms->side_bias == 0)  /* 5' probes */
		{
			step = (double)nb_ok5p/parms->nb_probes;
			pos = 0;
		}
		if (parms->side_bias == 1)  /* 3' probes */
		{
			step = -(double)nb_ok3p/parms->nb_probes;
			pos = nb_ok-0.1;
		}
		while ((pcnt < parms->nb_probes) && ((int)pos >= 0) && ((int)pos < nb_ok))
		{
			i = is_ok[(int)pos],
/*			PrintProbe(sequence_select_part(gen,i,i+(parms->oligo_len-1)),name,i); */
			pos += step;
			pcnt++;
		}

	}
	free(is_ok);
	return pcnt;
}

SEXP seqlib_generate_oligos(SEXP seqs, SEXP parameters)
{
	const char *pname;
	int i,j;
	oligoDesignParms *parms;
	parms = NewOligoDesignParmRec();
	SEXP res,elm;
	for (i=0; i< length(parameters); i++)
	{
		pname = CHAR(STRING_ELT(getAttrib(parameters,R_NamesSymbol), i));
		elm = VECTOR_ELT(parameters,i);
		if ((strcmp(pname,"t.melt")==0) && isReal(elm) && (length(elm) == 2))
		{
			parms->tm_min = REAL(elm)[0];
			parms->tm_max = REAL(elm)[1];
		}
		else if ((strcmp(pname,"cg.content")==0) && isReal(elm) && (length(elm) == 2))
		{
			parms->cg_min = REAL(elm)[0];
			parms->cg_max = REAL(elm)[1];
		}
		else if ((strcmp(pname,"oligo.len")==0) && isInteger(elm) && (length(elm) == 1))
		{
			parms->oligo_len = (INTEGER)(elm)[0];
		}
		else if ((strcmp(pname,"nimblegen.passes")==0) && isInteger(elm) && (length(elm) == 1))
		{
			parms->NimbleGenMaxPass = (INTEGER)(elm)[0];
		}
		else if ((strcmp(pname,"forbidden.sequences") == 0) && isString(elm) && (length(elm) >= 1))
		{
			parms->forbidden_seqs = (sequence_tp**)malloc(sizeof(sequence_tp*)*(length(elm)+1));
			for (j=0; j<length(elm); j++)
				parms->forbidden_seqs[j] = sequence_from_string(CHAR(STRING_ELT(elm,j)));
			parms->forbidden_seqs[length(elm)] = NULL;
		}
#ifdef DEBUG
		else
		{
		fprintf(stderr,"Skipping  unknown parameter %s : \n",pname);
		if (isString(elm))
		{
			for (j=0; j<length(elm); j++)
				fprintf(stderr,"String Value[%d] = %s\n",j,CHAR(STRING_ELT(elm, j)));
		}
		if (isReal(elm))
		{
			for (j=0; j<length(elm); j++)
				fprintf(stderr," Real Value[%d] = %f\n",j,REAL(elm)[j]);
		}
		if (isInteger(elm))
		{
			for (j=0; j<length(elm); j++)
				fprintf(stderr,"Int Value[%d] = %d\n",j,INTEGER(elm)[j]);
		}
		if (isLogical(elm))
		{
			for (j=0; j<length(elm); j++)
				fprintf(stderr,"Logical Value[%d] = %d\n",j,LOGICAL(elm)[j]);
		}
		}
#endif
		
	}
#ifdef DEBUG
		fprintf(stderr," ----Generate oligos ----- \n");
		fprintf(stderr,"cg	 : %f -%f \n",parms->cg_min,parms->cg_max);
		fprintf(stderr,"tm	 : %f -%f \n",parms->tm_min,parms->tm_max);
		fprintf(stderr,"oligo-len  : %d \n",parms->oligo_len);
#endif
	PROTECT(res = mkString("Done gen olis"));
	UNPROTECT(1);
	return res;
	
}
