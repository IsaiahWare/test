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
void transpose_32(int M, int N, int A[N][M], int B[M][N]);
void transpose_64(int M, int N, int A[N][M], int B[M][N]);
void transpose_other(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N]){

	// Determine variation of transpose function to run based on dimensions of input matrix
	switch(N) {
    	case 32:	// 32x32 matrix
    		transpose_32(M, N, A, B);
        	break;
    	case 64:	// 64x64 matrix
    		transpose_64(M, N, A, B);
      		break;
      	default:	// all matrices not 32x32 or 64x64
        	transpose_other(M, N, A, B);
        	break;
    }
}


/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/*
 * transpose_32 - Matrix transposition function optimized for a 32x32 matrix
 					Fairly optimal solution - # of misses is ~290
 * Params: 
 *	M - Number of columns in the matrix
 *	N -	Number of rows in the matrix
 *	A[N][M] - Input matrix - the one transposed from
 *	B[N][M] - Input matrix - the one transposed to
 * Returns: void
 */
char transpose_32_desc[] = "Transpose a 32x32 matrix";
void transpose_32(int M, int N, int A[N][M], int B[M][N]){

	int val = 0, pos = 0;
	// Iterates through each column and row
	for (int j = 0; j < N; j += 8) { 
		for (int i = 0; i < N; i += 8) {
			for (int n = i; n < i + 8; n++) {
				for (int m = j; m < j + 8; m++) {
					if (n != m) {
						B[m][n] = A[n][m];
					} else {
						val = A[n][m];
            pos = n;
					}
				}
				// If row and column are same, element is defined as a diagonal and our temporarily saved element is assigned
				if (i == j) {
					B[pos][pos] = val;
				}
			}	
		}
	}
}

/*
 * transpose_64 - Matrix transposition function optimized for a 64x64 matrix
 					Not the most optimal solution - couldn't determine a way to reduce # of misses below ~1750 =(
 * Params: 
 *	M - Number of columns in the matrix
 *	N -	Number of rows in the matrix
 *	A[N][M] - Input matrix - the one transposed from
 *	B[N][M] - Input matrix - the one transposed to
 * Returns: void
 */
char transpose_64_desc[] = "Transpose a 64x64 matrix";
void transpose_64(int M, int N, int A[N][M], int B[M][N]){
int i, j, k, l, t1, t2, t3, t4, t5, t6, t7, t8;
  for (i = 0; i < N; i += 8) {
      for (j = 0; j < M; j += 8) {
          for (k = i; k < i + 4; k++) {
              t1 = A[k][j];
              t2 = A[k][j + 1];
              t3 = A[k][j + 2];
              t4 = A[k][j + 3];
              t5 = A[k][j + 4];
              t6 = A[k][j + 5];
              t7 = A[k][j + 6];
              t8 = A[k][j + 7];

              B[j][k] = A[k][j];
              B[j + 1][k] = A[k][j + 1];
              B[j + 2][k] = A[k][j + 2];
              B[j + 3][k] = A[k][j + 3];
              B[j][k + 4] = A[k][j + 4];
              B[j + 1][k + 4] = A[k][j + 5];
              B[j + 2][k + 4] = A[k][j + 6];
              B[j + 3][k + 4] = A[k][j + 7];
          }
          for (l = j + 4; l < j + 8; l++) {

              t1 = B[l - 4][i + 4];
              t2 = B[l - 4][i + 5];
              t3 = B[l - 4][i + 6];
              t4 = B[l - 4][i + 7];
              t5 = A[i + 4][l - 4];
              t6 = A[i + 5][l - 4];
              t7 = A[i + 6][l - 4];
              t8 = A[i + 7][l - 4];

              B[l - 4][i + 4] = t5;
              B[l - 4][i + 5] = t6;
              B[l - 4][i + 6] = t7;
              B[l - 4][i + 7] = t8;

              B[l][i] = t1;
              B[l][i + 1] = t2;
              B[l][i + 2] = t3;
              B[l][i + 3] = t4;

              B[l][i + 4] = A[i + 4][l];
              B[l][i + 5] = A[i + 5][l];
              B[l][i + 6] = A[i + 6][l];
              B[l][i + 7] = A[i + 7][l];
          }
      }
  }
}

/*
 * transpose_other - Matrix transposition function optimized for matricies that aren't 32x32 or 64x64
 					Appears to be a fairly optimal solution - # of misses for a 61x67 matrix is ~1800
 * Params: 
 *	M - Number of columns in the matrix
 *	N -	Number of rows in the matrix
 *	A[N][M] - Input matrix - the one transposed from
 *	B[N][M] - Input matrix - the one transposed to
 * Returns: void
 */
char transpose_other_desc[] = "Transpose any matrix that isn't 32x32 or 64x64";
void transpose_other(int M, int N, int A[N][M], int B[M][N]){

	int n, m; 		// Indecies for rows and columns in matrix
	int row, col;	// Track current row and column in matrix
	int d_val = 0;	// Hold value of diagonal element found in matrix (detailed in below code)
	int diag = 0;	// Hold position of diagonal element found in matrix (detailed in below code)

	// Iterates through each column and row
	for (col = 0; col < M; col += 16) {
		for (row = 0; row < N; row += 16) {

			// For each row and column after current one, until end of matrix
			for (n = row; (n < row + 16) && (n < N); n++) {
				for (m = col; (m < col + 16) && (m < M); m++) {

					// If row and column number do not match, transposition will occur
					if (n != m) {
						B[m][n] = A[n][m];
					// Else, row and column number are same and element in matrix is defined as a diagonal
					} else {

						// Assign diagonal element to a temporary variable
						// This saves an individual cache miss on each run through the matrix where the columns and rows still match up
						diag = n;
						d_val = A[n][m];
					}
				}
				// If row and column are same, element is defined as a diagonal and our temporarily saved element is assigned
				if (row == col) {
					B[diag][diag] = d_val;
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

    // Register custom transpose functions
    registerTransFunction(transpose_32, transpose_32_desc);
    registerTransFunction(transpose_64, transpose_64_desc);
    registerTransFunction(transpose_other, transpose_other_desc);

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