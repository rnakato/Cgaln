/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "seq.h"
#include "BA.h"
#include <math.h>

static void Colony_method(outputalignment *, int **, unsigned short **, TYPE_SEEDT *, TYPE_BLKT **);
static void copy_to_BlAlign(BlAlign *, struct colony *);
static void filter_cln(struct colonyset *);

void BA(outputalignment *aln, int **table_value, unsigned short **table_num, int strand){
  switch(strand){
  case FORWARD:
    read_table();
    Colony_method(aln, table_value, table_num, seedtable_a, blktable_a);
    break;
  case REVERSE:
    read_table_revcom();
    Colony_method(aln, table_value, table_num, seedtable_a_revcom, blktable_a_revcom);
    break;
  }

  if(!reverse || strand==REVERSE){
    free(seedtable_b);
    free(poistable_a);
    free(poistable_b);
  }
  return;
}

static struct colonyset *colonyset_new(){
  struct colonyset *colonyset = (struct colonyset *)my_malloc(sizeof(struct colonyset), "colonyset");
  colonyset->num = 0;
  colonyset->nummax = CLNMAX_DEFAULT;
  colonyset->cln = (struct colony *)my_malloc(sizeof(struct colony)*colonyset->nummax, "colonyset->cln");
  return colonyset;
}

static void colonyset_delete(struct colonyset *colonyset){
  free(colonyset->cln);
  free(colonyset);
}

static void Colony_method(outputalignment *aln, int **table_value, unsigned short **table_num, TYPE_SEEDT *seedtable, TYPE_BLKT **blktable){
  int i;
  int M_a = idata_a->blocknum[idata_a->cnt];  /* block number of genomeA */
  int M_b = idata_b->blocknum[idata_b->cnt];  /* block number of genomeB */
  double *simscore = (double *)my_malloc(M_a*sizeof(double), "simscore");
  struct column *pre_clm = (struct column *)my_malloc(sizeof(struct column)*(M_a +1), "pre_clm");
  struct colonyset *colonyset = colonyset_new();

  for(i=0; i<M_a+1; i++){ 
    pre_clm[i].id = -1; 
    pre_clm[i].score = 0;
  }
  for(i=0; i<M_b; i++){
    Score_Blocksim(table_value[i], table_num[i], simscore, seedtable, blktable, 0, M_a);
    Find_colony(colonyset, simscore, pre_clm, i);
  }
  free(simscore);
  free(pre_clm);

  aln->blarraynum = SIZE_BL_DEFAULT;
  aln->bl = (BlAlign *)my_malloc(sizeof(BlAlign)*aln->blarraynum, "aln->bl");
  aln->blnum = 0;

  if(filtercln) filter_cln(colonyset);  /* if -fc is on */

  /* global alignment in the candidates of colonies  */
  for(i=0; i<colonyset->num; i++){
    if(colonyset->cln[i].num == -1) continue;
    if(colonyset->cln[i].m_value > par.xdrop){
      Align_colony(&(colonyset->cln[i]), table_value, table_num, seedtable, blktable);
      copy_to_BlAlign(&(aln->bl[aln->blnum]), &(colonyset->cln[i]));
      aln->blnum++;
      if(aln->blnum >= aln->blarraynum){
	aln->blarraynum += SIZE_BL_DEFAULT;
	aln->bl = (BlAlign *)my_realloc(aln->bl, sizeof(BlAlign)*aln->blarraynum, "aln->bl");
      }
      free(colonyset->cln[i].clnarray);
    }
  }

  colonyset_delete(colonyset);
  free(seedtable);
  blktable_delete(blktable, size[par.kmer_ba]);
  return;
}

void Score_Blocksim(int *table_value, unsigned short *table_num, double *simscore, TYPE_SEEDT *seedtable, TYPE_BLKT **blktable, int start, int end){
  int i,j, num, kmin, M_num = -1, seedvalue;
  int geta = par.xdrop / (ALLOW_RAND_CELL +1);
  int a, b, nk_a, nk_b;
  for(i=start; i<end; i++) simscore[i]=0;

  for(i=0; table_value[i] != -1; i++){
    seedvalue = table_value[i];
    if(!seedtable[seedvalue] || !seedtable_b[seedvalue]) continue;
    nk_a = seedtable[seedvalue];
    nk_b = seedtable_b[seedvalue];
    if(nk_a > par.seedthre_a || nk_b > par.seedthre_b) continue;

    for(j=0; j<nk_a; j++){
      M_num = blktable[seedvalue][j];
      if(M_num < start) continue;
      if(M_num >= end) break;

      num=1;
      while(j+1 < nk_a){
	if(blktable[seedvalue][j+1] == M_num){
	  j++; num++;
	}else break;
      }
      kmin = min(table_num[i], num);
      if(kmin < POISSON_MAX){
	a = nk_a*POISSON_MAX + kmin;
	b = nk_b*POISSON_MAX + kmin;
       	if(kmin==num) simscore[M_num] -= poistable_a[a].first + poistable_b[b].second;
	else          simscore[M_num] -= poistable_b[b].first + poistable_a[a].second;
      }
    }
  }
  for(i=start; i<end; i++) simscore[i] = simscore[i]/(double)TEMP_TERM - geta;
}

static void copy_to_BlAlign(BlAlign *bl, struct colony *cln){
  int i;
  bl->ini_i = cln->ini_i;
  bl->m_i = cln->m_i;
  bl->num = cln->num;
  
  int num = cln->num * 2;
  bl->cln = (struct clnarray *)my_malloc(num * sizeof(struct clnarray), "bl->cln");
  for(i=0; i<num; i++) bl->cln[i] = cln->clnarray[i];

  return;
}

static int compare_ini_i(const void *_a, const void *_b){
  struct colony *a = (struct colony *)_a;
  struct colony *b = (struct colony *)_b;
  if(a->ini_i < b->ini_i)      return (-1);
  else if(a->ini_i > b->ini_i) return (1);
  else return (0);
}

static int compare_ini_j(const void *_a, const void *_b){
  struct colony *a = (struct colony *)_a;
  struct colony *b = (struct colony *)_b;
  if(a->ini_j < b->ini_j)      return (-1);
  else if(a->ini_j > b->ini_j) return (1);
  else return (0);
}

static void compare_cln(struct colony *cln, int i, int j, int size, int *array, int *arraynum, int seq){
  int k, on;
  TYPE_BLKT end, start;
  if(!seq){ start=cln[j].ini_i; end=cln[i].m_i;}
  else{     start=cln[j].ini_j; end=cln[i].m_j;}

  if(end < start + CONS_FC) return;

  if(end >= start) goto filter;
  else{
    on=0;
    for(k=0; k<(*arraynum); k++){
      if(array[k] == end){ 
	on=1; 
	break;
      }
    }
    if(on) goto filter;
    else array[(*arraynum)++] = end;
  }
  return;

 filter:
  if(cln[i].m_value > cln[j].m_value){
    cln[j].num = -1;
    if(j!=size-1) compare_cln(cln, i, j+1, size, array, arraynum, seq);
  }else{
    cln[i].num = -1;
  }
}

static int filter(struct colony *cln, int size, int(*compare)(const void *, const void *), int seq){
  int i, number=0, arraynum=0;
  int *array=(int *)my_calloc(size, sizeof(int), "array_filter");

  qsort(cln, size, sizeof(struct colony), compare);
  for(i=0; i<size-1; i++) if(cln[i].num != -1) compare_cln(cln, i, i+1, size, array, &arraynum, seq);
  for(i=0; i<size;   i++) if(cln[i].num != -1 && i != number) cln[number++] = cln[i];
  
  free(array);
  return number;
}

static void filter_cln(struct colonyset *colonyset){
  /* genomeA */
  int number = filter(colonyset->cln, colonyset->num, compare_ini_i, 0);
  /* genomeB */
  int number2 = filter(colonyset->cln, number, compare_ini_j, 1);
  colonyset->num = number2;
  return;
}

