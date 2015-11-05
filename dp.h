/* Copyright(c)  Ryuichiro Nakato <rnakato@iam.u-tokyo.ac.jp>
 * This file is a part of Cgaln-1.*.* sources.
 */
#ifndef DP_H_
#define DP_H_

#define S_SAME 2
#define S_DIFFERENT -2
#define GAP_INSERT 4
#define GAP_EXTEND 2
#define GAP_INSERT_HOXD 400 
#define GAP_EXTEND_HOXD 30
#define SIZE_DPX 200
#define X_EXTENDHSP2 10
#define MINUS_INF -100000000
#define SEMIGLOBAL_THRE 50

#define USEHOXD 1

enum{GLOBAL, SEMIGLOBAL, LOCAL};

struct DP_m{int s0, s1, s2, s3;};

#endif /*DP_H_*/
