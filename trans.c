/* 
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);
void transpose32(int M, int N, int A[N][M], int B[M][N]);
void transpose64(int M, int N, int A[N][M], int B[M][N]);
void transposeExtra(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N]){
  if (N == 32) {
      transpose32(M, N, A, B);
  } else if (N == 64) {
      transpose64(M, N, A, B);
  } else {
      transposeExtra(M, N, A, B);
  }
}


/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/*
 * transpose32 - Transposes 32x32 matrix
 */
char transpose32_desc[] = "Transpose a 32x32 matrix";
void transpose32(int M, int N, int A[N][M], int B[M][N]){
  for (int i = 0; i < N; i += 8) {
      for (int j = 0; j < M; j += 8) {
          for (int k = i; k < i + 8; k++) {
              for (int l = j; l < j + 8; l += 8) {
                for (int z = 0; z < 8; z++) {
                  B[l+z][k] = A[k][l+z];
                }
              }
          }
      }
  }
}

/*
 * transpose64 - Transposes 64x64 matrix
 */
char transpose64_desc[] = "Transpose a 64x64 matrix";
void transpose64(int M, int N, int A[N][M], int B[M][N]){
  for (int i = 0; i < N; i += 8) {
      for (int j = 0; j < M; j += 8) {
          for (int m = i; m < i + 4; m++) {
              B[j][m] = A[m][j];
              B[j + 1][m] = A[m][j + 1];
              B[j + 2][m] = A[m][j + 2];
              B[j + 3][m] = A[m][j + 3];
              B[j][m + 4] = A[m][j + 4];
              B[j + 1][m + 4] = A[m][j + 5];
              B[j + 2][m + 4] = A[m][j + 6];
              B[j + 3][m + 4] = A[m][j + 7];
          }
          for (int n = j + 4; n < j + 8; n++) {
              int val1 = B[n - 4][i + 4];
              int val2 = B[n - 4][i + 5];
              int val3 = B[n - 4][i + 6];
              int val4 = B[n - 4][i + 7];

              B[n - 4][i + 4] = A[i + 4][n - 4];
              B[n - 4][i + 5] = A[i + 5][n - 4];
              B[n - 4][i + 6] = A[i + 6][n - 4];
              B[n - 4][i + 7] = A[i + 7][n - 4];
              
              B[n][i] = val1;
              B[n][i + 1] = val2;
              B[n][i + 2] = val3;
              B[n][i + 3] = val4;
              B[n][i + 4] = A[i + 4][n];
              B[n][i + 5] = A[i + 5][n];
              B[n][i + 6] = A[i + 6][n];
              B[n][i + 7] = A[i + 7][n];
          }
      }
  }
}

/*
 * transposeExtra - Transposes non 32x32 or 64x64 matrix
 */
char transposeExtra_desc[] = "Transpose any matrix that isn't 32x32 or 64x64";
void transposeExtra(int M, int N, int A[N][M], int B[M][N]){
  for (int i = 0; i < N; i += 23) {
      for (int j = 0;  j < M; j += 23) {
          for (int k = i; k < i + 23; k++) {
            if (k >= N) {
              break;
            }
              for (int l = j; l < j + 23; l++) {
                if (l >= M) {
                  break;
                }
                  B[l][k] = A[k][l];
              }
          }
      }
  }
}

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 
    registerTransFunction(transpose32, transpose32_desc);
    registerTransFunction(transpose64, transpose64_desc);
    registerTransFunction(transposeExtra, transposeExtra_desc);

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}