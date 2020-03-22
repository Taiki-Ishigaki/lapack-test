#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

    integer n3 =  3 * n;
    integer info; 

    doublereal *vl = (doublereal *)malloc(sizeof(doublereal)* n * n );
    doublereal *vr = (doublereal *)malloc(sizeof(doublereal)* n * n );
    doublereal *work = (doublereal *)malloc(sizeof(doublereal)* n3);

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

    printf( "info = %ld\n", info ); 
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

    printf( "info = %ld\n", info ); 
    return info;
}

int isDeclimited(char c, char del){
    if(c == del){
        return 1;
    }else{
        return 0;
    }
}

int isNumber(char c){
    if( ('0' <= c && c <= '9') || (c == '+' || c == '-') ){
        return 1;
    }else{
        return 0;
    }
}

#define BUF_SIZE 256

double dReadDataElem(FILE* fp){
    char tmp, token[BUF_SIZE];
    int i = 0;
    token[0] = '\0';
    tmp = fgetc(fp);
    if( isNumber(tmp) ){
        token[i++] = tmp;
    }
    while( !feof( fp ) ){
        if( isDeclimited(tmp = fgetc(fp), ' ') ) break;
        token[i++] = tmp;
        if(i > BUF_SIZE){
            printf("too long");
            token[0] = '\0';
            return 0;
        } 
    }
    token[i] = '\0';
    return atof(token);
}

double* ReadMat(const char *fname, double* mat, int row, int col){
    int i, j;
    FILE* fp = fopen(fname, "rb"); 
	if (fp == NULL) { 
		printf("file open error!\n");
        mat = NULL;
		return NULL;
	}
    for( i=0; i<row; i++ ){
        for( j=0; j<col; j++ ){
            double tmp = dReadDataElem( fp );
            // printf("%lf ", tmp);
            mat[i + row*j] = tmp;
        }
        // printf("\n");
    }
    return mat;
}

int main( int argc, char **argv ) {
    int i, j; 
    integer n = 20; /* 正方行列の次数 */

    double sum = 0;
    int test_time = 1;
    for(j = 0; j < test_time; j++){
        doublereal *a  = (doublereal *)malloc(sizeof(doublereal)* n * n ); 
        doublereal *wr = (doublereal *)malloc(sizeof(doublereal)* n ); 
        doublereal *wi = (doublereal *)malloc(sizeof(doublereal)* n );  
        // doublereal *vr = (doublereal *)malloc(sizeof(doublereal)* n * n );

        ReadMat("x.dat", a, n, n);
        /* クロック開始 */
        printf( "start, \n" );
        clock_t c = clock();

        // 固有値. 実部が wr, 虚部が wi
        eigenvalues( n, a, wr, wi ); 
        // eigenvalues_rightvectors( n, a, wr, wi, vr ); 

        /* （オプション）確認出力 */
        for(i = 0; i < n; i++) {
          printf("%5d %15.7e %15.7e\n", i + 1, *(wr + i), *(wi + i));
        }

        /* クロック終了 */
        double time = ( (double)clock() - (double)c ) / CLOCKS_PER_SEC;
        printf( "done, elapsed time = %f [sec]\n", time );
        sum += time;

        free(wi);
        free(wr);
        free(a);
    }
    printf( "elapsed average time = %f [sec]\n", sum / test_time );

    return 0;
}