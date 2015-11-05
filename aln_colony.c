/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#include "seq.h"
#include "BA.h"

struct swg{double s0, s1, s2, s3;};

static void renewDP(int, int, int, struct colonyset *, struct column *);
static void Make_colony(struct colonyset *, int, int, double);
static void Update_colony(struct colonyset *, int, int, int, double);
static void traceback(struct colony *, struct swg **, int, int, int, int);

void Find_colony(struct colonyset *colonyset, double *simscore, struct column *pre_clm, int M_b){
  int i, M_a=idata_a->blocknum[idata_a->cnt];
  double s0, s1, s2;
  struct column c_pre, c_cur;
  c_pre.id = -1; c_pre.score = 0;
  c_cur.id = -1; c_cur.score = 0;

  for(i=0; i<M_a; i++){
    s0 = pre_clm[i].score   + simscore[i];
    s1 = pre_clm[i+1].score + simscore[i] - par.shiftp;
    s2 = c_pre.score        + simscore[i] - par.shiftp;
    c_cur.score = max4(s0, s1, s2, 0);
    
    if(c_cur.score == s0)      renewDP(pre_clm[i].id,   i, M_b, colonyset, &(c_cur));
    else if(c_cur.score == s1) renewDP(pre_clm[i+1].id, i, M_b, colonyset, &(c_cur));
    else if(c_cur.score == s2) renewDP(c_pre.id,        i, M_b, colonyset, &(c_cur));
    else c_cur.id = -1;
    pre_clm[i] = c_pre;
    c_pre = c_cur;
  }

  pre_clm[M_a] = c_cur;
  return;
}

static void renewDP(int pre_id, int M_a, int M_b, struct colonyset *colonyset, struct column *c_cur){
  if(pre_id != -1){
    /* X-DROP OFF */
    if(c_cur->score < colonyset->cln[pre_id].m_value - par.xdrop){
      c_cur->score = 0;
      c_cur->id = -1;
    }else{
      c_cur->id = pre_id;
      Update_colony(colonyset, pre_id, M_a, M_b, c_cur->score);
    }
  }else{
    c_cur->id = colonyset->num;
    Make_colony(colonyset, M_a, M_b, c_cur->score);
  }
}


static void Make_colony(struct colonyset *colonyset, int m, int n, double score){
  int num = colonyset->num;
  colonyset->cln[num].ini_i = m;
  colonyset->cln[num].ini_j = n;
  colonyset->cln[num].m_value = score;
  colonyset->cln[num].m_i = m;
  colonyset->cln[num].m_j = n;
  colonyset->cln[num].i_j_plus = m-n;
  colonyset->cln[num].i_j_minus = m-n;
  colonyset->cln[num].num =0;
  colonyset->num++;
  if(colonyset->num >= colonyset->nummax){
    colonyset->nummax += CLNMAX_DEFAULT;
    colonyset->cln = (struct colony *)my_realloc(colonyset->cln, sizeof(struct colony)*colonyset->nummax, "colonyset->cln");
  }
  return;
}

static void Update_colony(struct colonyset *colonyset, int id, int m, int n, double value){
  if(colonyset->cln[id].m_value < value){
    colonyset->cln[id].m_value = value;
    colonyset->cln[id].m_i = m;
    colonyset->cln[id].m_j = n;
  }
  if(m-n > colonyset->cln[id].i_j_plus)  colonyset->cln[id].i_j_plus = m-n;
  if(m-n < colonyset->cln[id].i_j_minus) colonyset->cln[id].i_j_minus = m-n;
}

void Align_colony(struct colony *cln, int **table_value, unsigned short **table_num, TYPE_SEEDT *seedtable, TYPE_BLKT **blktable){
  int i,j, left, right;
  double s;
  int shift_p, shift_m;
  int width_i, width_j;
  double *simscore = (double *)my_malloc(idata_a->blocknum[idata_a->cnt]*sizeof(double), "simscore2");
  struct swg **swg=NULL;

  /*--- preparation ---*/
  shift_p = cln->i_j_plus +1;
  shift_m = cln->i_j_minus -1;
  width_i = (cln->m_i - cln->ini_i + 1) +1;
  width_j = (cln->m_j - cln->ini_j + 1) +1;

  swg = (struct swg **)my_malloc(width_i * sizeof(struct swg *), "swg");
  for(i=0; i<width_i; i++) swg[i] = (struct swg *)my_calloc(width_j, sizeof(struct swg), "swgline");

  /*--- scoring ---*/
  int geta_i = cln->ini_i -1;
  int geta_j = cln->ini_j -1;
  for(j=1; j<width_j; j++){
    left  = max((shift_m + j + geta_j), cln->ini_i);
    right = min((shift_p + j + geta_j), cln->m_i) +1;
    Score_Blocksim(table_value[j+geta_j], table_num[j+geta_j], simscore, seedtable, blktable, left, right);

    left -= geta_i;
    right -= geta_i;
    for(i=left; i<right; i++){
      s = simscore[i+geta_i];
      swg[i][j].s0 = swg[i-1][j-1].s3 + s;
      swg[i][j].s1 = swg[i-1][j  ].s3 + s - par.shiftp;
      swg[i][j].s2 = swg[i  ][j-1].s3 + s - par.shiftp;
      swg[i][j].s3 = max4(swg[i][j].s0, swg[i][j].s1, swg[i][j].s2, 0);
    }
  }
  free(simscore);

  /*----- traceback -----*/  
  traceback(cln, swg, width_i, width_j, geta_i, geta_j);

  for(i=0; i<width_i; i++) free(swg[i]);
  free(swg);
}

static void init_clnarray(struct colony *cln, int width_i, int width_j){
  int i;
  int size = width_i + width_j +1;
  cln->clnarray = (struct clnarray *)my_malloc(size * sizeof(struct clnarray), "cln->clnarray");
  for(i=0; i<size; i++) cln->clnarray[i].a = -1;
  cln->num=0;
}

static void traceback(struct colony *cln, struct swg **swg, int width_i, int width_j, int geta_i, int geta_j){
  int i_back = width_i-1;
  int j_back = width_j-1;

  init_clnarray(cln, width_i, width_j);

  while(1){
    cln->clnarray[cln->num].a = i_back + geta_i;
    cln->clnarray[cln->num].b = j_back + geta_j;
    cln->num++;

    if(swg[i_back][j_back].s3 == swg[i_back][j_back].s0){
      i_back--; j_back--;
      if(i_back < 1 || j_back < 1) break;
    }
    else if(swg[i_back][j_back].s3 == swg[i_back][j_back].s1){
      i_back--;
      if(i_back < 1) break;
    }
    else if(swg[i_back][j_back].s3 == swg[i_back][j_back].s2){
      j_back--;
      if(j_back < 1) break;
    }else{
      break;
    }
  }
}
