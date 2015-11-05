/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "seq.h"
#include "NA.h"
#include "dp.h"

static int HoxD70[4][4]={
  {  91, -114,  -31, -123},
  {-114,  100, -125,  -31},
  { -31, -125,  100, -114},
  {-123,  -31, -114,   91}
};
static char *str_DPtype[3]={ "GLOBAL", "SEMIGLOBAL", "LOCAL"};

static void Scoring_each_cell(struct DP_m **, int, int, int, int, char, char, int);
static void Traceback(struct DP_m **, Fasta *, Fasta *, int, int, int, int, int, int, int);
static struct DP_m **matrix_new(int, int, int);
static void matrix_delete(struct DP_m **, int);

static void renew_scoremax(struct DP_m **matrix, int i, int j, int *scoremax, int *maxi, int *maxj){
  if(*scoremax < matrix[i][j].s3){
    *scoremax = matrix[i][j].s3;
    *maxi = i; *maxj = j;
  }
}

void DP(struct region *region, Fasta *fst1, Fasta *fst2, int strand){
  int i,j;
  int xsize = region->lenx +1;
  int ysize = region->leny +1;
  int maxi=-1, maxj=-1, scoremax=0;
  if(region->lenx <=0 || region->leny <=0) return;
  int type;

  if(abs(xsize-ysize) >SEMIGLOBAL_THRE) type = SEMIGLOBAL; else type = GLOBAL;
  if(opt.debug) fprintf(par.OUT, "#DP: %s\n", str_DPtype[type]);
  
  struct DP_m **matrix = matrix_new(xsize, ysize, type);

  for(i=0; i+1<xsize; i++){
    for(j=0; j+1<ysize; j++){
      Scoring_each_cell(matrix, i, j, xsize, ysize, fst1->body[region->x1 +i], fst2->body[region->y1 +j], type);

      if(type==LOCAL){
	matrix[i+1][j+1].s3 = max4(matrix[i+1][j+1].s0, matrix[i+1][j+1].s1, matrix[i+1][j+1].s2, 0);
	renew_scoremax(matrix, i+1, j+1, &scoremax, &maxi, &maxj);
      }else matrix[i+1][j+1].s3 = max( max( matrix[i+1][j+1].s0, matrix[i+1][j+1].s1), matrix[i+1][j+1].s2);
    }
  }

  /*----- traceback -----*/
  if(type==LOCAL){ if(scoremax) Traceback(matrix, fst1, fst2, maxi,    maxj,    region->x1, region->y1, strand, FORWARD, type);}
  else{                         Traceback(matrix, fst1, fst2, xsize-1, ysize-1, region->x1, region->y1, strand, FORWARD, type);}

  matrix_delete(matrix, xsize);
}

static void Scoring_each_cell(struct DP_m **matrix, int i, int j, int xsize, int ysize, char ref_x, char ref_y, int type){
  int x = define_base_nomask(ref_x);
  int y = define_base_nomask(ref_y);
#ifdef USEHOXD
  int u = GAP_INSERT_HOXD;
  int v = GAP_EXTEND_HOXD;
#else
  int u = GAP_INSERT;
  int v = GAP_EXTEND;
#endif

#ifdef USEHOXD
  matrix[i+1][j+1].s0 = matrix[i][j].s3 + HoxD70[x][y];
#else
  if(x==y) matrix[i+1][j+1].s0 = matrix[i][j].s3 + S_SAME;
  else     matrix[i+1][j+1].s0 = matrix[i][j].s3 + S_DIFFERENT;
#endif

  if(type == SEMIGLOBAL){
    if(i!=xsize-2) matrix[i+1][j+1].s1 = max(matrix[i+1][j].s3 -u, matrix[i+1][j].s1) -v;
    else           matrix[i+1][j+1].s1 = max(matrix[i+1][j].s3,    matrix[i+1][j].s1);
    if(j!=ysize-2) matrix[i+1][j+1].s2 = max(matrix[i][j+1].s3 -u, matrix[i][j+1].s2) -v;
    else           matrix[i+1][j+1].s2 = max(matrix[i][j+1].s3,    matrix[i][j+1].s2);
  }else{
    matrix[i+1][j+1].s1 = max(matrix[i+1][j].s3 -u, matrix[i+1][j].s1) -v;
    matrix[i+1][j+1].s2 = max(matrix[i][j+1].s3 -u, matrix[i][j+1].s2) -v;
  }

  return;
}

void DP_xdrop(struct region *region, Fasta *fst1, Fasta *fst2, int strand, int forback){
  int i,j;
  int x1 = region->x1;
  int y1 = region->y1;
  int x2 = region->x2;
  int y2 = region->y2;
  int xsize = min(region->lenx, SIZE_DPX) +1;
  int ysize = min(region->leny, SIZE_DPX) +1;
  int maxi=0, maxj=0, scoremax=0, maxtemp;
  int type = GLOBAL;

  /* use globalDP for small region */
  if(DP4smallregion(region, fst1, fst2, strand, IA_THRE1)) return;

  struct DP_m **matrix = matrix_new(xsize, ysize, type);

  for(i=0; i+1<xsize; i++){
    for(j=0; j+1<ysize; j++){
      if(forback==FORWARD) Scoring_each_cell(matrix, i, j, xsize, ysize, fst1->body[x1+i], fst2->body[y1+j], type);
      else                 Scoring_each_cell(matrix, i, j, xsize, ysize, fst1->body[x2-i], fst2->body[y2-j], type);
      /* MAX */
      matrix[i+1][j+1].s3 = max( max( matrix[i+1][j+1].s0, matrix[i+1][j+1].s1), matrix[i+1][j+1].s2);
      renew_scoremax(matrix, i+1, j+1, &scoremax, &maxi, &maxj);
    }
  }

  if(forback==FORWARD) subst_region(region, maxi +x1+1, maxj +y1+1, region->x2, region->y2);
  else                 subst_region(region, region->x1, region->y1, x2 - maxi-2, y2 - maxj-2);

  if(scoremax){
    maxtemp=0;
    for(i=0; i<xsize; i++) if(matrix[i][ysize-1].s3 >maxtemp) maxtemp = matrix[i][ysize-1].s3;
    for(i=0; i<ysize; i++) if(matrix[xsize-1][i].s3 >maxtemp) maxtemp = matrix[xsize-1][i].s3;
    if(maxtemp > scoremax - X_EXTENDHSP2) DP_xdrop(region, fst1, fst2, strand, forback);
  }

  /*----- traceback -----*/
  if(forback==FORWARD) Traceback(matrix, fst1, fst2, maxi, maxj, x1, y1, strand, forback, type);
  else                 Traceback(matrix, fst1, fst2, maxi, maxj, x2, y2, strand, forback, type);

  matrix_delete(matrix, xsize);
  return;
}

static void Print_and_reset(Fasta *fst1, Fasta *fst2, int x, int y, int *len, int strand, int forback){
  if(*len){
    if(forback==FORWARD){ print_DPresult(x, y, *len, fst1, fst2, strand);}
    else{print_DPresult(x + (*len) -1, y +(*len) -1, *len, fst1, fst2, strand);}
    *len=0;
  }
}

static void Traceback(struct DP_m **matrix, Fasta *fst1, Fasta *fst2, int i, int j, int startx, int starty, int strand, int forback, int type){
  int x=-1, y=-1, len=0;
  do{
    if(matrix[i][j].s3 == matrix[i][j].s0){
      if(!len){ 
	if(forback==FORWARD){ x=startx +i; y=starty +j;}
	else{ x = startx -i; y = starty -j;}
	len=1;
      }else len++;
      i--; j--;
    }else if(matrix[i][j].s3 == matrix[i][j].s1){
      Print_and_reset(fst1, fst2, x, y, &len, strand, forback);
      j--;
    }else if(matrix[i][j].s3 == matrix[i][j].s2){
      Print_and_reset(fst1, fst2, x, y, &len, strand, forback);
      i--;
    }
    if(type == LOCAL && matrix[i][j].s3 <=0) break;
  }while(i>0 || j>0);

  Print_and_reset(fst1, fst2, x, y, &len, strand, forback);
  return;
}

static void setscore(struct DP_m **matrix, int x, int y, int value){
  matrix[x][y].s0 = value;
  matrix[x][y].s1 = value;
  matrix[x][y].s2 = value;
  matrix[x][y].s3 = value;
}

static void init_matrix(struct DP_m **matrix, int xsize, int ysize, int type){
  int i;
  int u,v;

#ifdef USEHOXD
    u = GAP_INSERT_HOXD;
    v = GAP_EXTEND_HOXD;
#else
    u = GAP_INSERT;
    v = GAP_EXTEND;
#endif

  for(i=1; i<ysize-1; i++){
    matrix[0][i].s0 = MINUS_INF;
    matrix[0][i].s2 = MINUS_INF;
    if(type == GLOBAL){
      matrix[0][i].s1 = -u - v*i;
      matrix[0][i].s3 = matrix[0][i].s1;
    }else if(type == SEMIGLOBAL){
      matrix[0][i].s1 = 0;
      matrix[0][i].s3 = 0;
    }else{ /* LOCAL*/
      matrix[0][i].s1 = MINUS_INF;
      matrix[0][i].s3 = 0;
    }
  }
  for(i=1; i<xsize-1; i++){
    matrix[i][0].s0 = MINUS_INF;
    matrix[i][0].s1 = MINUS_INF;
    if(type == GLOBAL){
      matrix[i][0].s2 = -u - v*i;
      matrix[i][0].s3 = matrix[i][0].s2;
    }else if(type == SEMIGLOBAL){
      matrix[i][0].s2 = 0;
      matrix[i][0].s3 = 0;
    }else{ /* LOCAL*/
      matrix[i][0].s2 = MINUS_INF;
      matrix[i][0].s3 = 0;
    }
  }
  setscore(matrix, 0, 0, 0);
  setscore(matrix, xsize-1, 0, MINUS_INF);
  setscore(matrix, 0, ysize-1, MINUS_INF);  
  return;
}

static struct DP_m **matrix_new(int xsize, int ysize, int type){
  int i;
  struct DP_m **matrix = (struct DP_m **)my_malloc(sizeof(struct DP_m *)* xsize, "DPmatrix");

  for(i=0; i<xsize; i++) matrix[i] = (struct DP_m *)my_malloc(sizeof(struct DP_m )* ysize, "base_line");
  init_matrix(matrix, xsize, ysize, type);
  return matrix;
}

static void matrix_delete(struct DP_m **matrix, int xsize){
  int i;
  for(i=0; i<xsize; i++) free(matrix[i]);
  free(matrix);
}
