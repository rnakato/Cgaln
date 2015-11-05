/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "seq.h"
#include "NA.h"

enum{BETWEENHSP, UPSTREAM, DOWNSTREAM};

static void align_region(struct region *, Fasta *, Fasta *, int, int);
static void align_betweenhsp(struct hsplist *, int, struct region *, Fasta *, Fasta *, int, int);
static int compare_hsp(const void *, const void *);

static int Penalty(struct hsplist *hsplist, int i, int j){
  struct hsp *hsp = hsplist->hsp;
  int diff = 10 + abs((hsp[i].sx - hsp[i].sy) - (hsp[j].sx - hsp[j].sy)) / hsplist->penaweight;
  int a = hsp[j].sx + hsp[j].len - hsp[i].sx;
  int b = hsp[j].sy + hsp[j].len - hsp[i].sy;
  int s = max(a,b) - HSP_CHAIN;

  if((a > 0 || b > 0) && s > 0){
    if(!consistent) return(diff + s);
    else return(10000000);  /* if option "-cs" is on */
  }else{
    return diff;
  }
}

int ChainHSP(struct hsplist *hsplist){
  int i,j,h;
  int maxi=0, maxh=0;
  struct hsp *hsp = hsplist->hsp;

  for(i=1; i<hsplist->num; i++){
    for(j=0; j<i; j++){
      h = hsp[j].H + hsp[i].S - Penalty(hsplist, i, j);
      if(h > hsp[i].H){
	hsp[i].H = h;
	hsp[i].link = j;
      }
    }
    if(hsp[i].H > maxh){
      maxh = hsp[i].H;
      maxi = i;
    }
  }
  hsplist->last_i = maxi;
  return maxi;
}

static void FindHSP_ShorterSeed(struct hsplist *hsplist, Fasta *fst1, Fasta *fst2, struct region *region, int sr){
  int i, seedvalue, posi, seedweight = hsplist->seedweight;
  int x1 = region->x1;
  int x2 = min(region->x2, fst1->length -window[seedweight]);
  int y1 = region->y1;
  int y2 = min(region->y2, fst2->length -window[seedweight]);
  int tablesize = region->lenx + region->leny-1;
  int *seedindex=NULL, *tupletable=NULL;
  int mask=0; /* do not mask lower case */

  if(!sr) maketable4region(fst1, &seedindex, &tupletable, region->lenx, x1, seedweight, mask);
  else    maketable4region(fst1, &seedindex, &tupletable, region->lenx,  0, seedweight, mask);

  struct diag_table *diag_table = diag_table_new(tablesize, x1);
  for(i=y1; i<y2; i++){
    seedvalue = SeedScoring_notmasked(fst2, i, seedweight);
    if(seedvalue < 0) continue;
    for(posi = seedindex[seedvalue];; posi = tupletable[posi]){
      if(posi == -1) break;
      if(!sr)renewHSP(hsplist, diag_table, fst1, fst2, tablesize, posi+x1, i, window[seedweight], x2-1);
      else   renewHSP(hsplist, diag_table, fst1, fst2, tablesize, posi,    i, window[seedweight], x2-1);
    }
  }
  free(tupletable);
  free(seedindex);
  
  scan_diagtable(hsplist, diag_table, tablesize, window[seedweight]);
  return;
}

void chain_and_IA(struct hsplist *hsplist, struct region *region, Fasta *fst1, Fasta *fst2, int strand, int sr, int repeat){
  int i;
  if(!hsplist->num) goto finish;
  qsort(hsplist->hsp, hsplist->num, sizeof(struct hsp), compare_hsp);

  /*-- if option "-nc" is on --*/
  if(!opt.chaining && repeat==FIRST){
    for(i=0; i<hsplist->num; i++) print_HSP(hsplist->hsp, fst1, fst2, i, strand, 0);
    goto finish;
  }

  i = ChainHSP(hsplist);
  while(i != -1){
    if(!sr){
      print_HSP(hsplist->hsp, fst1, fst2, i, strand, 0);
      if(iterative) align_betweenhsp(hsplist, i, region, fst1, fst2, strand, repeat);
    }else{
      print_HSP(hsplist->hsp, fst1, fst2, i, strand, region->x1);
      // align_betweenhsp(hsplist, i, region, fst1, fst2, strand, repeat);
    }
    i = hsplist->hsp[i].link;
  }

 finish:
  free(hsplist->hsp);
  free(hsplist);
}

static void align_bl_ShorterSeed(struct region *region, Fasta *fst1, Fasta *fst2, int strand, int sr){
  int seedweight = SEEDWEIGHT_2; /* shorter seed */
  if(min(region->lenx, region->leny) < window[seedweight]) return;

  struct hsplist *hsplist = hsplist_new(seedweight, PENALTY_WEIGHT2);
  FindHSP_ShorterSeed(hsplist, fst1, fst2, region, sr);
  chain_and_IA(hsplist, region, fst1, fst2, strand, sr, SECOND);
  return;
}

static void convert_fst2revcom(Fasta *fst, Fasta *fst_temp, int start, int len){
  int i;
  fst_temp->body = my_malloc(sizeof(char)*(len+1), "fst_temp->body");  
  for(i=0; i<len; i++) fst_temp->body[i] = convert_complement(fst->body[start + len-1-i]);
  fst_temp->body[len] = '\0';
  fst_temp->length = len;
  return;
}

static void do_short_reverse_alignment(struct region *region, Fasta *fst1, Fasta *fst2, int strand){
  Fasta *fst_temp = fasta_new();
  convert_fst2revcom(fst1, fst_temp, region->x1, region->lenx);
  fprintf(par.OUT, "###{ align (short_reverse) between HSP\n");
  align_bl_ShorterSeed(region, fst_temp, fst2, strand, 1);
  fprintf(par.OUT, "###}\n\n");

  if(fst_temp->body) free(fst_temp->body);
  free(fst_temp);
}

static void align_region(struct region *region, Fasta *fst1, Fasta *fst2, int strand, int regtype){
  if(region->lenx <= 0 || region->leny <= 0) return;
  if(!opt.otype && printregion) print_region(region, strand);

  if(short_reverse) do_short_reverse_alignment(region, fst1, fst2, strand);

  fprintf(par.OUT, "###{ align between HSP\n");

  /* use globalDP for small region */
  if(DP4smallregion(region, fst1, fst2, strand, IA_THRE1)) goto finish;
  
  /* extend hsp with gapped DP */
  if(ext_gappeddp){
    if(regtype != UPSTREAM){
      if(opt.debug) fprintf(par.OUT, "#DP_xdrop1\n");
      DP_xdrop(region, fst1, fst2, strand, FORWARD);
    }
    if(regtype != DOWNSTREAM){
      if(opt.debug) fprintf(par.OUT, "#DP_xdrop2\n");
      DP_xdrop(region, fst1, fst2, strand, REVERSE);
    }
  }
  if(!opt.otype && printregion) print_region(region, strand);
  if(region->lenx > par.blocksize*2 && region->leny > par.blocksize*2) goto finish;

  /* align rest region */
  if(DP4smallregion(region, fst1, fst2, strand, IA_THRE1)) goto finish;
  else{
    if(opt.debug) fprintf(par.OUT, "#ShorterSeed\n");
    align_bl_ShorterSeed(region, fst1, fst2, strand, 0);
  }

 finish:
  fprintf(par.OUT, "###}\n\n");
  return;
}

static void check_and_DP(struct region *region, Fasta *fst1, Fasta *fst2, int strand){
  if(region->lenx < IA_THRE2 && region->leny < IA_THRE2){
    if(opt.debug) fprintf(par.OUT, "#GlobalDP2\n");
    DP(region, fst1, fst2, strand);
  }
}

static void align_betweenhsp(struct hsplist *hsplist, int i, struct region *region_base, Fasta *fst1, Fasta *fst2, int strand, int repeat){
  struct region region;
  struct hsp *hsp = hsplist->hsp;
  int next=hsp[i].link, last=hsplist->last_i;

  if(hsp[i].link != -1){
    /* between neighboring HSPs */
    subst_region(&region, hsp[next].sx + hsp[next].len,
		 hsp[next].sy + hsp[next].len,
		 hsp[i].sx-1,
		 hsp[i].sy-1);
    if(repeat==FIRST) align_region(&region, fst1, fst2, strand, BETWEENHSP);
    else check_and_DP(&region, fst1, fst2, strand);
  }else{
    /* between region_start and first HSP */
    subst_region(&region, region_base->x1,
		 region_base->y1,
		 hsp[i].sx-1,
		 hsp[i].sy-1);
    if(repeat==FIRST) align_region(&region, fst1, fst2, strand, UPSTREAM);
    else check_and_DP(&region, fst1, fst2, strand);

    /* between last HSP and region_end */
    subst_region(&region, hsp[last].sx + hsp[last].len,
		 hsp[last].sy + hsp[last].len,
		 region_base->x2,
		 region_base->y2);
    if(repeat==FIRST) align_region(&region, fst1, fst2, strand, DOWNSTREAM);
    else check_and_DP(&region, fst1, fst2, strand);
  }
}

void maketable4region(Fasta *fst, int **seedindex, int **tupletable, int length, int start, int seedweight, int mask){
  int i, seedvalue;
  *seedindex  = (int *)my_malloc(size[seedweight] * sizeof(int), "seedindex4region");
  *tupletable = (int *)my_malloc(length           * sizeof(int), "tupletable4region");
  for(i=0; i<size[seedweight]; i++) (*seedindex)[i] = -1;
  for(i=0; i<length; i++) (*tupletable)[i] = -1;

  int len = length-window[seedweight];
  for(i=0; i<=len; i++){
    if(mask) seedvalue = SeedScoring_masked(fst, start +i, seedweight);
    else     seedvalue = SeedScoring_notmasked(fst, start +i, seedweight);
    if(seedvalue != -1 && seedvalue != -2){
      (*tupletable)[i] = (*seedindex)[seedvalue];
      (*seedindex)[seedvalue] = i;
    }
  }
  return;
}

int DP4smallregion(struct region *region, Fasta *fst1, Fasta *fst2, int strand, int thre){
  int on=0;
  if(region->lenx < thre || region->leny < thre){
    DP(region, fst1, fst2, strand);
    subst_region(region, region->x2+1, region->y1, region->x2, region->y2);
    on=1;
  }
  return on;
}

static int compare_hsp(const void *_hsp0, const void *_hsp1){
  struct hsp *hsp0 = (struct hsp *)_hsp0;
  struct hsp *hsp1 = (struct hsp *)_hsp1;

  if(hsp0->sx < hsp1->sx) return -1;
  else if(hsp0->sx > hsp1->sx) return 1;
  else return 0;
}
