/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "NA.h"

/* flags */
struct range{
  int bef, af;
  int blklow, blkup;
};

static void align_bl(BlAlign *, Fasta *, Fasta *, int);
static void FindHSP(BlAlign *, struct hsplist *, Fasta *, Fasta *);
static void extendhsp(struct hsplist *, Fasta *, Fasta *, int, int);

void NA(FILE *IN, Fasta *fst1, Fasta *fst2, outputalignment *aln, int strand, int *state){
  int i;
  read_multifasta(IN, fst1, strand, state);
  for(i=0; i<aln->blnum; i++){
    if(!strand) fprintf(par.OUT, "#####{ (forward) genomeA-fasta%d\tgenomeB-fasta%d\tBl #%d\n", idata_a->cnt+1, idata_b->cnt+1, i+1);
    else        fprintf(par.OUT, "#####{ (reverse) genomeA-fasta%d_revcom\tgenomeB-fasta%d\tBl #%d\n", idata_a->cnt+1, idata_b->cnt+1, i+1);

    align_bl(&(aln->bl[i]), fst1, fst2, strand);
    free(aln->bl[i].cln);

    fprintf(par.OUT, "#####}\n\n\n");
  }
  free(aln->bl);
  free(fst1->head);
  free(fst1->body);
}

/* nucleotide-level alignment within the region of bl */
static void align_bl(BlAlign *bl, Fasta *fst1, Fasta *fst2, int strand){
  struct region region;
  subst_region(&region, bl->cln[bl->num-1].a * par.blocksize, 
	       bl->cln[bl->num-1].b * par.blocksize,
	       min((bl->cln[0].a +1)*par.blocksize, fst1->length) -1, 
	       min((bl->cln[0].b +1)*par.blocksize, fst2->length) -1);

  struct hsplist *hsplist = hsplist_new(par.kmer_na, PENALTY_WEIGHT);
  FindHSP(bl, hsplist, fst1, fst2);
  chain_and_IA(hsplist, &region, fst1, fst2, strand, 0, FIRST);
  return;
}

static void fix_length_a(BlAlign *bl, int *ip, struct range *range, int *upper, int *lower){
  int i = *ip;
  range->blklow = bl->cln[i].a;
  int start_b = bl->cln[i].b;

  if(i != bl->num-1 && bl->cln[i+1].a != range->blklow) range->bef = 1;
  else range->bef = 0;

  while(bl->cln[i].b == start_b){
    i--;
    if(i == -1) break;
  }
  range->blkup = bl->cln[i+1].a;

  if(i != -1 && bl->cln[i+1].a != bl->cln[i].a) range->af = 1;
  else range->af = 0;

  *upper = (range->blkup+1) * par.blocksize -1;
  *lower = range->blklow * par.blocksize;

  *ip = i;
  return;
}

static void FindHSP(BlAlign *bl, struct hsplist *hsplist, Fasta *fst1, Fasta *fst2){
  int i,j, seedvalue, end = par.blocksize;
  int x_start = bl->ini_i * par.blocksize;
  int x_end = (bl->m_i +1) * par.blocksize - 1;
  if(bl->m_i==idata_a->blocknum[idata_a->cnt]-1) x_end = fst1->length - window[par.kmer_na];
  int posi_a, posi_b, upper, lower;
  struct range range;
  int *seedindex=NULL, *tupletable=NULL;
  int tablesize = par.blocksize*2-1;
  int mask=1;

  maketable4region(fst1, &seedindex, &tupletable, (x_end-x_start +1), x_start, par.kmer_na, mask);

  struct diag_table *diag_table = diag_table_new(tablesize, x_start);

  i = bl->num-1;
  while(i>=0){
    if(bl->cln[i].b == idata_b->blocknum[idata_b->cnt]-1) end = fst2->length%par.blocksize - window[par.kmer_na];
    posi_b = bl->cln[i].b * par.blocksize;
    fix_length_a(bl, &i, &range, &upper, &lower);
    for(j=0; j<end; j++, posi_b++, upper += range.af, lower -= range.bef){
      seedvalue = SeedScoring_masked(fst2, posi_b, par.kmer_na);
      if(seedvalue < 0) continue;

      for(posi_a = seedindex[seedvalue];; posi_a = tupletable[posi_a]){
	if(posi_a == -1) break;
	if(posi_a + x_start > upper) continue;
	else if(posi_a + x_start >= lower) renewHSP(hsplist, diag_table, fst1, fst2, tablesize, posi_a + x_start, posi_b, window[par.kmer_na], x_end);
	else break;
      }
    }
  }
  free(tupletable);
  free(seedindex);

  scan_diagtable(hsplist, diag_table, tablesize, window[par.kmer_na]);
  return;
}

static int checkfasta(Fasta *fst1, Fasta *fst2, int x, int x_e, int diag){
  int i;
  int xs = x_e;
  int xe = x;
  int ys = x_e-diag;
  int ye = x-diag;
  for(i=xs; i<xe; i++) if(fst1->body[i]==4) return 0;
  for(i=ys; i<ye; i++) if(fst2->body[i]==4) return 0;
  return 1;
}

static int add2hsplist(struct hsplist *hsplist, struct diag_table *table, Fasta *fst1, Fasta *fst2, int d, int x, int con){
  int on=0;
  if(table[d].S != con){
    hspcopy_from_diagtable(hsplist, table, d);
    extendhsp(hsplist, fst1, fst2, table[d].m_end, x);
    on=1;
  }
  table[d].len = -1;
  return on;
}

void renewHSP(struct hsplist *hsplist, struct diag_table *table, Fasta *fst1, Fasta *fst2, int tablesize, int x, int y, int seedlength, int x_end){
  int d;
  int end  = x + seedlength;
  int diag = x - y;
  int con  = seedlength * HSPCONS;

  if(diag<0) d = diag%tablesize + tablesize-1;
  else       d = diag%tablesize;

  if(table[d].len == -1) goto finish;

  if(table[d].diag == diag){
    if(x - table[d].x_e < HSP_UNITE_THRE && checkfasta(fst1, fst2, x, table[d].x_e, diag)){
      table[d].len += end - table[d].x_e;
      table[d].S    = table[d].len *HSPCONS;
      table[d].x_e  = end;
      return;
    }else{
      if(add2hsplist(hsplist, table, fst1, fst2, d, x,     con)) table[d].m_end = table[d].x_e;
    }
  }else  add2hsplist(hsplist, table, fst1, fst2, d, x_end, con);

 finish:
  table[d].x_e  = end;
  table[d].sy   = y;
  table[d].len  = seedlength;
  table[d].S    = con;
  table[d].diag = diag;
  return;
}

void hspcopy_from_diagtable(struct hsplist *hsplist, struct diag_table *diag_table, int i){
  hsplist->hsp[hsplist->num].sx  = diag_table[i].x_e - diag_table[i].len;
  hsplist->hsp[hsplist->num].sy  = diag_table[i].sy;
  hsplist->hsp[hsplist->num].len = diag_table[i].len;
  hsplist->hsp[hsplist->num].S   = diag_table[i].S;
  hsplist->hsp[hsplist->num].link = -1;
  hsplist->num++;
  if(hsplist->num >= hsplist->nummax){
    hsplist->nummax += HSP_SIZEMAX;
    hsplist->hsp = (struct hsp *)my_realloc(hsplist->hsp, sizeof(struct hsp)*(hsplist->nummax), "hsp");
  }
}

void scan_diagtable(struct hsplist *hsplist, struct diag_table *diag_table, int tablesize, int seedlength){
  int i;
  for(i=0; i<tablesize; i++){
    if(diag_table[i].len != -1 && diag_table[i].S > seedlength * HSPCONS) 
      hspcopy_from_diagtable(hsplist, diag_table, i);
  }
  free(diag_table);
  for(i=0; i<hsplist->num; i++) hsplist->hsp[i].H = hsplist->hsp[i].S;

  return;
}

static void rescorehsp(struct hsplist *hsplist, Fasta *fst1, Fasta *fst2){
  int i;
  char x, y;
  struct hsp *hsp = &(hsplist->hsp[hsplist->num-1]);
  hsp->S = 0;

  for(i=0; i<hsp->len; i++){
    x = fst1->body[hsp->sx + i];
    y = fst2->body[hsp->sy + i];
    if(x == -1 || y == -1){ /* none */ }
    else if(x == y) hsp->S += HSPCONS;
    else hsp->S -= HSPCONS;
  }
}

static int check_basepair(char ref_basex, char ref_basey, int *score, int *scoremax, int *max_x, int *max_y, int x, int y){
  int flag=0;
  char basex = define_base_nomask(ref_basex);
  char basey = define_base_nomask(ref_basey);
  if(basex == -1 || basey == -1){     // gap
    flag=1;
  }else if(basex == 4 || basey == 4){ // separator
    flag=0;
  }else if(basex == basey){           // match
    *score += HSPCONS;
    if(*score > *scoremax){
      *scoremax = *score;
      *max_x = x;
      *max_y = y;
    }
    flag=1;
  }else{                              // mismatch
    *score -= HSPCONS;
    if(*score < *scoremax - X_EXTENDHSP) flag=0;
    else flag=1;
  }
  return flag;
}

static void extendhsp(struct hsplist *hsplist, Fasta *fst1, Fasta *fst2, int x_start, int x_end){
  int score, scoremax, max_x, max_y, flag=0;
  struct hsp *hsp = &(hsplist->hsp[hsplist->num-1]);
  int sx = hsp->sx -1;
  int sy = hsp->sy -1;
  int ex = hsp->sx + hsp->len;
  int ey = hsp->sy + hsp->len;
  rescorehsp(hsplist, fst1, fst2);

  /* upstream */
  max_x = hsp->sx, max_y = hsp->sy;
  score=0, scoremax=0;
  while(sx >= x_start && sy >= 0){
    flag = check_basepair(fst1->body[sx], fst2->body[sy], &score, &scoremax, &max_x, &max_y, sx, sy);
    if(!flag) break; else{ sx--; sy--;}
  }
  hsp->sx = max_x;
  hsp->sy = max_y;
  hsp->S += scoremax;

  /* downstream */
  max_x = ex; max_y = ey;
  score=0; scoremax=0;
  while(ex < x_end && ey < fst2->length){
    flag = check_basepair(fst1->body[ex], fst2->body[ey], &score, &scoremax, &max_x, &max_y, ex, ey);
    if(!flag) break; else{ ex++; ey++;}
  }
  hsp->len = max_x - hsp->sx +1;
  hsp->S += scoremax;
  return;
}

struct hsplist *hsplist_new(int seedweight, int penaweight){
  struct hsplist *hsplist = (struct hsplist *)my_malloc(sizeof(struct hsplist), "hsplist");
  hsplist->nummax=HSP_SIZEMAX;
  hsplist->num = 0;
  hsplist->hsp = (struct hsp *)my_malloc(hsplist->nummax* sizeof(struct hsp), "hsp");
  hsplist->last_i = -1;
  hsplist->seedweight = seedweight;
  hsplist->penaweight = penaweight;
  return hsplist;
}

struct diag_table *diag_table_new(int tablesize, int x){
  int i;
  struct diag_table *diag_table = (struct diag_table *)my_malloc(tablesize * sizeof(struct diag_table), "diag_table");
  for(i=0; i<tablesize; i++){
    diag_table[i].len = -1;
    diag_table[i].m_end = x; 
  }
  return diag_table;
}

void subst_region(struct region *region, int x1, int y1, int x2, int y2){
  region->x1 = x1;
  region->y1 = y1;
  region->x2 = x2;
  region->y2 = y2;
  region->lenx = x2-x1+1;
  region->leny = y2-y1+1;
}
