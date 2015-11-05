/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "cgaln.h"

void output_BAresult(){
  int i,j;
  int geta_a=0, geta_b=0;
  for(i=0; i<idata_a->cnt; i++) geta_a += idata_a->blocknum[i];
  for(i=0; i<idata_b->cnt; i++) geta_b += idata_b->blocknum[i];

  for(i=0; i<aln_for.blnum; i++){
    fprintf(par.OUT, "\n#colony %d:\n", i);
    for(j=0; j<aln_for.bl[i].num; j++) fprintf(par.OUT, "%d\t%d\n", aln_for.bl[i].cln[j].a + geta_a, aln_for.bl[i].cln[j].b + geta_b);
  }
  if(reverse){
    for(i=0; i<aln_rev.blnum; i++){
      fprintf(par.OUT, "\n#colony_rev %d:\n", i);
      for(j=0; j<aln_rev.bl[i].num; j++) fprintf(par.OUT, "%d\t%d\n", geta_a + idata_a->blocknum[idata_a->cnt] - (aln_rev.bl[i].cln[j].a)-1, aln_rev.bl[i].cln[j].b + geta_b);
    }
  }

  /* free bl */
  bl_delete(aln_for.bl, aln_for.blnum);
  if(reverse) bl_delete(aln_rev.bl, aln_rev.blnum);
  return;
}

static void print_desc(int i, struct hsp *hsp, int strand){
  if(!strand) fprintf(par.OUT, "#HSP number: %d, length: %d, score: %d, score_cumulative: %d\n",         i, hsp->len, hsp->S, hsp->H);
  else        fprintf(par.OUT, "#HSP(revcom) number: %d, length: %d, score: %d, score_cumulative: %d\n", i, hsp->len, hsp->S, hsp->H);
}

static void print_fastahead(int i, int start, int end, struct hsp *hsp, int strand, int genome){
  if(!genome){
    if(!strand) fprintf(par.OUT, ">A_fst%d:%d-%d:", idata_b->cnt+1, start, end);
    else        fprintf(par.OUT, ">A_fst%d_revcom:%d-%d:", idata_b->cnt+1, start, end);
  }else         fprintf(par.OUT, ">B_fst%d:%d-%d:", idata_b->cnt+1, start, end);
  
  if(hsp) fprintf(par.OUT, "HSP number %d:score %d:score_cumulative %d\n", i, hsp->S, hsp->H);
  else    fprintf(par.OUT, "DP\n");
}
static void func_print(int x1, int x2, int y1, int y2, int strand, Fasta *fst1, Fasta *fst2, int i, struct hsp *hsp, int sx){
  int j;
  int x_rev1 = idata_a->length_each[idata_a->cnt] -x1 +1;
  int x_rev2 = idata_a->length_each[idata_a->cnt] -x2 +1;
  long geta_a=0, geta_b=0;
  for(j=0; j<idata_a->cnt; j++) geta_a += idata_a->length_each[j];
  for(j=0; j<idata_b->cnt; j++) geta_b += idata_b->length_each[j];

  switch(opt.otype){
  case 0: /* gnuplot */
    if(hsp) print_desc(i, hsp, strand);
    if(!strand){
      fprintf(par.OUT, "%ld\t%ld\n",   x1 + geta_a, y1 + geta_b);
      fprintf(par.OUT, "%ld\t%ld\n\n", x2 + geta_a, y2 + geta_b);
    }else{
      fprintf(par.OUT, "%ld\t%ld\n",   x_rev1 + geta_a, y1 + geta_b);
      fprintf(par.OUT, "%ld\t%ld\n\n", x_rev2 + geta_a, y2 + geta_b);
    }
    break;
  case 1: /* Ens */
    if(hsp) print_desc(i, hsp, strand);
    if(!strand) fprintf(par.OUT, "A_fst%d:%d-%d, B_fst%d:%d-%d\n\n", idata_a->cnt+1, x1,     x2,     idata_b->cnt+1, y1, y2);
    else        fprintf(par.OUT, "A_fst%d:%d-%d, B_fst%d:%d-%d\n\n", idata_a->cnt+1, x_rev1, x_rev2, idata_b->cnt+1, y1, y2);
    break;
  case 2: /* fasta */
    x1 -= sx; x2 -= sx;  /* for short reverse (in other case, sx==0) */
    if(!strand){
      print_fastahead(i, x1, x2, hsp, strand, 0);
      for(j=x1; j<=x2; j++) fprintf(par.OUT, "%c", convert_num2base(fst1->body[j]));
      fprintf(par.OUT, "\n");
      print_fastahead(i, y1, y2, hsp, strand, 1);
      for(j=y1; j<=y2; j++) fprintf(par.OUT, "%c", convert_num2base(fst2->body[j]));
      fprintf(par.OUT, "\n\n");
    }else{
      print_fastahead(i, x_rev1, x_rev2, hsp, strand, 0);
      for(j=x1; j<=x2; j++) fprintf(par.OUT, "%c", convert_num2base(fst1->body[j]));
      fprintf(par.OUT, "\n");
      print_fastahead(i, y1, y2, hsp, strand, 1);
      for(j=y1; j<=y2; j++) fprintf(par.OUT, "%c", convert_num2base(fst2->body[j]));
      fprintf(par.OUT, "\n\n");
    }
    break;
  default: break;
  }
  return;
}

void print_HSP(struct hsp *hsp, Fasta *fst1, Fasta *fst2, int i, int strand, int sx){
  int x1 = sx + hsp[i].sx;
  int x2 = sx + hsp[i].sx + hsp[i].len -1;
  int y1 = hsp[i].sy;
  int y2 = hsp[i].sy + hsp[i].len -1;
  func_print(x1, x2, y1, y2, strand, fst1, fst2, i, &(hsp[i]), sx);
  return;
}

void print_DPresult(int ex, int ey, int len, Fasta *fst1, Fasta *fst2, int strand){
  int x1 = ex-len;
  int x2 = ex-1;
  int y1 = ey-len;
  int y2 = ey-1;
  func_print(x1, x2, y1, y2, strand, fst1, fst2, 0, NULL, 0);
  return;
}

void print_region(struct region *region, int strand){
  int i;
  long x1, x2, y1, y2;
  long geta_a=0, geta_b=0;
  for(i=0; i<idata_a->cnt; i++) geta_a += idata_a->length_each[i];
  for(i=0; i<idata_b->cnt; i++) geta_b += idata_b->length_each[i];

  if(!strand){
    x1 = region->x1 + geta_a;
    x2 = region->x2 + geta_a;
    y1 = region->y1 + geta_b;
    y2 = region->y2 + geta_b;
  } else{
    x1 = idata_a->length_each[idata_a->cnt] - region->x1 -1 + geta_a;
    x2 = idata_a->length_each[idata_a->cnt] - region->x2 -1 + geta_a;
    y1 = region->y1 + geta_b;
    y2 = region->y2 + geta_b;
  }
  fprintf(par.OUT, "\n");
  fprintf(par.OUT, "%ld\t%ld\n", x1, y1);
  fprintf(par.OUT, "%ld\t%ld\n", x2, y1);
  fprintf(par.OUT, "\n");
  fprintf(par.OUT, "%ld\t%ld\n", x1, y1);
  fprintf(par.OUT, "%ld\t%ld\n", x1, y2);
  fprintf(par.OUT, "\n");
  fprintf(par.OUT, "%ld\t%ld\n", x2, y1);
  fprintf(par.OUT, "%ld\t%ld\n", x2, y2);
  fprintf(par.OUT, "\n");
  fprintf(par.OUT, "%ld\t%ld\n", x1, y2);
  fprintf(par.OUT, "%ld\t%ld\n", x2, y2);
  fprintf(par.OUT, "\n");
}

void output_fastaboundary(){
  int i;
  long geta_a=0, geta_b=0, total_a=0, total_b=0;

  if(block){
    for(i=0; i<idata_a->fstnum; i++) total_a += idata_a->blocknum[i];
    for(i=0; i<idata_b->fstnum; i++) total_b += idata_b->blocknum[i];
    for(i=0; i<idata_a->fstnum; i++){
      geta_a += idata_a->blocknum[i];
      fprintf(par.OUT, "\n\n%ld\t%d\n", geta_a, 0);
      fprintf(par.OUT, "%ld\t%ld\n\n", geta_a, total_b);
    }
    for(i=0; i<idata_b->fstnum; i++){
      geta_b += idata_b->blocknum[i];
      fprintf(par.OUT, "\n\n%d\t%ld\n",0, geta_b);
      fprintf(par.OUT, "%ld\t%ld\n\n", total_a, geta_b);
    }
  }else if(!opt.otype){
    for(i=0; i<idata_a->fstnum; i++) total_a += idata_a->length_each[i];
    for(i=0; i<idata_b->fstnum; i++) total_b += idata_b->length_each[i];
    for(i=0; i<idata_a->fstnum; i++){
      geta_a += idata_a->length_each[i];
      fprintf(par.OUT, "\n\n%ld\t%d\n", geta_a, 0);
      fprintf(par.OUT, "%ld\t%ld\n\n", geta_a, total_b);
    }
    for(i=0; i<idata_b->fstnum; i++){
      geta_b += idata_b->length_each[i];
      fprintf(par.OUT, "\n\n%d\t%ld\n",0, geta_b);
      fprintf(par.OUT, "%ld\t%ld\n\n", total_a, geta_b);
    }
  }
}
