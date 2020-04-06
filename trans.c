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

        int i, j, k, l, t1, t2, t3, t4, t5, t6, t7, t8;
        for (i = 0; i < N; i += 8) {
            for (j = 0; j < M; j += 8) {
                for (k = i; k < i + 8; k++) {
                    for (l = j; l < j + 8; l += 8) {
                        t1 = A[k][l];
                        t2 = A[k][l + 1];
                        t3 = A[k][l + 2];
                        t4 = A[k][l + 3];
                        t5 = A[k][l + 4];
                        t6 = A[k][l + 5];
                        t7 = A[k][l + 6];
                        t8 = A[k][l + 7];
                        B[l][k] = t1;
                        B[l + 1][k] = t2;
                        B[l + 2][k] = t3;
                        B[l + 3][k] = t4;
                        B[l + 4][k] = t5;
                        B[l + 5][k] = t6;
                        B[l + 6][k] = t7;
                        B[l + 7][k] = t8;
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
              int val1 = A[i + 4][n - 4];
              int val2 = A[i + 5][n - 4];
              int val3 = A[i + 6][n - 4];
              int val4 = A[i + 7][n - 4];
              int val5 = B[n - 4][i + 4];

              B[n - 4][i + 4] = val1;
              B[n - 4][i + 5] = val2;
              B[n - 4][i + 6] = val3;
              B[n - 4][i + 7] = val4;

              for (int z = 0; z < 8; z++) {
                if (z < 4) {
                  B[n][i+z] = val5;
                } else {
                  B[n][i+z] = A[i+z][n];
                  B[n][i+z] = A[i+z][n];
                  B[n][i+z] = A[i+z][n];
                  B[n][i+z] = A[i+z][n];
                }
              }
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