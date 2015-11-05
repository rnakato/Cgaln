#include "util.h"

Fasta *fasta_new(){
  Fasta *fst = my_malloc(sizeof(Fasta), "fst");
  fst->head=NULL;
  fst->body=NULL;
  fst->length=0;
  return fst;
}

inputdata *idata_new(){
  inputdata *idata = my_malloc(sizeof(inputdata), "idata");
  idata->blocknum=NULL;
  idata->ksum = my_malloc(sizeof(int)*2, "ksum");
  idata->length_each=NULL;
  idata->length_total=0;
  idata->fstnum=0;
  idata->cnt=0;
  return idata;
}

void idata_delete(inputdata *idata){
  int i;
  if(idata->blocknum) free(idata->blocknum);
  if(idata->length_each) free(idata->length_each);
  for(i=0; i<2; i++){
    if(idata->ksum[i]) free(idata->ksum[i]);
  }
  free(idata->ksum);
  free(idata);
  idata=NULL;
  return;
}

void bl_delete(BlAlign *bl, int num){
  int i;
  for(i=0; i<num; i++) free(bl[i].cln);
  free(bl);
  bl=NULL;
}

void table_b_delete(int **table_value, unsigned short **table_num, int num){
  int i;
  for(i=0; i<num; i++){
    free(table_value[i]);
    free(table_num[i]);
  }
  free(table_value);
  free(table_num);
  table_value=NULL; table_num=NULL;
}

void blktable_delete(TYPE_BLKT **table, int num){
  int i;
  for(i=0; i<num; i++){
    if(table[i]) free(table[i]);
  }
  free(table);
  table=NULL;
}
