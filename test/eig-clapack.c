#include <stdio.h>
#include <math.h>
#include <time.h>

// GSL 
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics.h>

#ifndef MKL
#  include <f2c.h>
#  include <clapack.h>
#else
#  define integer long int
#  define doublereal double
#endif

// CLAPACK を使う．
//  see http://www.netlib.org/clapack/clapack.h

// ?geev : simple driver for eigenvalues/vectors
//         see http://www.netlib.org/lapack/lug/node32.html

integer eigenvalues( integer n, doublereal *a, doublereal *wr, doublereal *wi ) {
  /* LAPACK の _dgeev() を使って固有値（だけ）を求める */

  integer n3 = 3 * n;
  integer info; 

  doublereal *vl = (doublereal *)calloc(sizeof(doublereal), n * n);
  doublereal *vr = (doublereal *)calloc(sizeof(doublereal), n * n);
  doublereal *work = (doublereal *)calloc(sizeof(doublereal), n3);

  (void) dgeev_(
		/* char *jobvl */      "N",  /* "N" なので左固有ベクトルを計算しない */
		/* char *jobvr */      "N",  /* "N" なので右固有ベクトルを計算しない */ 
		/* integer *n */       &n,   /* 正方行列の次数 */
		/* doublereal *a, */   a,    /* A */
		/* integer *lda, */    &n,   /* A 用の作業領域 */
		/* doublereal *wr, */  wr,   /* 固有値の実部 */
		/* doublereal *wi, */  wi,   /* 固有値の虚部 */
		/* doublereal *vl, */  vl,   /* 左固有値 */
		/* integer *ldvl, */   &n,   /* 左固有値の作業用 */
		/* doublereal *vr, */  vr,   /* 右固有値 */
		/* integer *ldvr, */   &n,   /* 左固有値の作業用 */
		/* doublereal *work, */ work, /* 作業用 */
		/* integer *lwork, */  &n3,   /* 作業用の行列の次元 */
		/* integer *info */    &info);

  free( work ); 
  free( vr ); 
  free( vl ); 

  printf( "info = %d\n", info ); 
  return info;
}

integer eigenvalues_rightvectors( integer n, doublereal *a, doublereal *wr, doublereal *wi, doublereal *vr ) {
  /* LAPACK の _dgeev() を使って固有値と右固有ベクトルを求める */
  /* A * v(j) = lambda(j) * v(j), v(j) is the right eigen vector */

  integer n4 = 4 * n;
  integer info; 

  doublereal *vl = (doublereal *)calloc(sizeof(doublereal), n * n);
  doublereal *work = (doublereal *)calloc(sizeof(doublereal), n4);


  (void) dgeev_(
		/* char *jobvl */      "N",  /* "N" なので左固有ベクトルを計算しない */
		/* char *jobvr */      "V",  /* "V" なので右固有ベクトルを計算する */ 
		/* integer *n */       &n,   /* 正方行列の次数 */
		/* doublereal *a, */   a,    /* A */
		/* integer *lda, */    &n,   /* A 用の作業領域 */
		/* doublereal *wr, */  wr,   /* 固有値の実部 */
		/* doublereal *wi, */  wi,   /* 固有値の虚部 */
		/* doublereal *vl, */  vl,   /* 左固有値 */
		/* integer *ldvl, */   &n,   /* 左固有値の作業用 */
		/* doublereal *vr, */  vr,   /* 右固有値 */
		/* integer *ldvr, */   &n,   /* 左固有値の作業用 */
		/* doublereal *work, */ work, /* 作業用 */
		/* integer *lwork, */  &n4,   /* 作業用の行列の次元 */
		/* integer *info */    &info);

  free( work ); 
  free( vl ); 

  printf( "info = %d\n", info ); 
  return info;
}

int main( integer argc, char **argv ) {
  int i; 
  integer n = 2000; /* 正方行列の次数 */
  doublereal *a  = (doublereal *)calloc(sizeof(doublereal), n * n); 
  doublereal *wr = (doublereal *)calloc(sizeof(doublereal), n); 
  doublereal *wi = (doublereal *)calloc(sizeof(doublereal), n);  
  doublereal *vr = (doublereal *)calloc(sizeof(doublereal), n * n); 

  // ここからは，GSL の機能を使って，乱数を配列 a に格納
  gsl_rng_env_setup();
  gsl_rng_type *T = (gsl_rng_type *)gsl_rng_default;
  /* 乱数発生器 */
  gsl_rng *r = gsl_rng_alloc(T);
  /* システムクロックを使って乱数の初期値を設定 */
  gsl_rng_set (r, time (NULL));

  for(i = 0; i < (n * n); i++) {
    a[i] =  gsl_rng_uniform(r);  
  }
  gsl_rng_free(r);

  /* クロック開始 */
  printf( "start, \n" );
  clock_t c = clock();

  // 固有値. 実部が wr, 虚部が wi
  // eigenvalues( n, a, wr, wi ); 
  eigenvalues_rightvectors( n, a, wr, wi, vr ); 

  /* （オプション）確認出力 */
  // for(i = 0; i < n; i++) {
  //   printf("%5d %15.7e %15.7e\n", i + 1, *(wr + i), *(wi + i));
  // }

  /* クロック終了 */
  printf( "done, elapsed time = %f [sec]\n", ( (double)clock() - (double)c ) / CLOCKS_PER_SEC );

  free(wi);
  free(wr);
  free(a);

  return 0;
}