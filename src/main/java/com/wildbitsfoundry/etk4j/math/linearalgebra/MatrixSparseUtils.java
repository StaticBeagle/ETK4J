package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Arrays;

public class MatrixSparseUtils {

    /**
     * Applies the forward column and inverse row permutation specified by the two vector to the input matrix
     * and save the results in the output matrix. output[permRow[j],permCol[i]] = input[j,i]
     *
     * @param permRowInv (Input) Inverse row permutation vector. Null is the same as passing in identity.
     * @param input      (Input) Matrix which is to be permuted
     * @param permCol    (Input) Column permutation vector. Null is the same as passing in identity.
     * @param output     (Output) Matrix which has the permutation stored in it. Is reshaped.
     */
    static void permute(int[] permRowInv, MatrixSparse input, int[] permCol,
                        MatrixSparse output) {
        if (permRowInv != null && input.rows > permRowInv.length)
            throw new IllegalArgumentException("rowInv permutation vector must have at least as many elements as input has columns");
        if (permCol != null && input.cols > permCol.length)
            throw new IllegalArgumentException("permCol permutation vector must have at least as many elements as input has rows");

        output.reshape(input.rows, input.cols, input.nz_length);
        output.indicesSorted = false;
        output.nz_length = input.nz_length;

        int N = input.cols;

        // traverse through in order for the output columns
        int outputNZ = 0;
        for (int i = 0; i < N; i++) {
            int inputCol = permCol != null ? permCol[i] : i; // column of input to source from
            int inputNZ = input.col_idx[inputCol];
            int total = input.col_idx[inputCol + 1] - inputNZ; // total nz in this column

            output.col_idx[i + 1] = output.col_idx[i] + total;

            for (int j = 0; j < total; j++) {
                int row = input.nz_rows[inputNZ];
                output.nz_rows[outputNZ] = permRowInv != null ? permRowInv[row] : row;
                output.nz_values[outputNZ++] = input.nz_values[inputNZ++];
            }
        }
    }

    static int solveColB(MatrixSparse G, boolean lower, MatrixSparse B, int colB, double[] x, int[] pinv, IGrowArray g_xi, int[] w) {
        int X_rows = G.cols;
        int[] xi = adjust(g_xi, X_rows);
        int top = searchNzRowsInX(G, B, colB, pinv, xi, w);

        int idxB0;
        for (idxB0 = top; idxB0 < X_rows; ++idxB0) {
            x[xi[idxB0]] = 0.0;
        }

        idxB0 = B.col_idx[colB];
        int idxB1 = B.col_idx[colB + 1];

        int px;
        for (px = idxB0; px < idxB1; ++px) {
            x[B.nz_rows[px]] = B.nz_values[px];
        }

        for (px = top; px < X_rows; ++px) {
            int j = xi[px];
            int J = pinv != null ? pinv[j] : j;
            if (J >= 0) {
                int p;
                int q;
                if (lower) {
                    x[j] /= G.nz_values[G.col_idx[J]];
                    p = G.col_idx[J] + 1;
                    q = G.col_idx[J + 1];
                } else {
                    x[j] /= G.nz_values[G.col_idx[J + 1] - 1];
                    p = G.col_idx[J];
                    q = G.col_idx[J + 1] - 1;
                }

                while (p < q) {
                    int var10001 = G.nz_rows[p];
                    x[var10001] -= G.nz_values[p] * x[j];
                    ++p;
                }
            }
        }

        return top;
    }

    private static int searchNzRowsInX(MatrixSparse G, MatrixSparse B, int colB, int[] pinv, int[] xi, int[] w) {
        int X_rows = G.cols;
        if (xi.length < X_rows) {
            throw new IllegalArgumentException("xi must be at least G.numCols=" + G.cols);
        } else if (w.length < 2 * X_rows) {
            throw new IllegalArgumentException("w must be at least 2*G.numCols in length (2*number of rows in X) and first N elements must be zero");
        } else {
            int idx0 = B.col_idx[colB];
            int idx1 = B.col_idx[colB + 1];
            int top = X_rows;

            int i;
            for (i = idx0; i < idx1; ++i) {
                int rowB = B.nz_rows[i];
                if (rowB < X_rows && w[rowB] == 0) {
                    top = searchNzRowsInX_DFS(rowB, G, top, pinv, xi, w);
                }
            }

            for (i = top; i < X_rows; ++i) {
                w[xi[i]] = 0;
            }

            return top;
        }
    }

    private static int searchNzRowsInX_DFS(int rowB, MatrixSparse G, int top, int[] pinv, int[] xi, int[] w) {
        int N = G.cols;  // first N elements in w is the length of X
        int head = 0; // put the selected row into the FILO stack
        xi[head] = rowB; // use the head of xi to store where the stack it's searching. The tail is where
        // the graph ordered list of rows in B is stored.
        while (head >= 0) {
            // the column in G being examined
            int G_col = xi[head];
            int G_col_new = pinv != null ? pinv[G_col] : G_col;
            if (w[G_col] == 0) {
                w[G_col] = 1;
                // mark which child in the loop below it's examining
                w[N + head] = G_col_new < 0 || G_col_new >= N ? 0 : G.col_idx[G_col_new];
            }

            // See if there are any children which have yet to be examined
            boolean done = true;

            // The Right side after || is used to handle tall matrices. There will be no nodes matching
            int idx0 = w[N + head];
            int idx1 = G_col_new < 0 || G_col_new >= N ? 0 : G.col_idx[G_col_new + 1];

            for (int j = idx0; j < idx1; j++) {
                int jrow = G.nz_rows[j];
                if (jrow < N && w[jrow] == 0) {
                    w[N + head] = j + 1; // mark that it has processed up to this point
                    xi[++head] = jrow;
                    done = false;
                    break;          // It's a DFS so break and continue down
                }
            }

            if (done) {
                head--;
                xi[--top] = G_col;
            }
        }
        return top;
    }

    static int[] adjust(IGrowArray gwork, int desired, int zeroToM) {
        int[] w = adjust(gwork, desired);
        Arrays.fill(w, 0, zeroToM, 0);
        return w;
    }

    private static int[] adjust(IGrowArray gwork, int desired) {
        if (gwork == null) {
            gwork = new IGrowArray();
        }

        gwork.reshape(desired);
        return gwork.data;
    }


    static double[] adjust(DGrowArray gwork, int desired) {
        if (gwork == null) {
            gwork = new DGrowArray();
        }

        gwork.reshape(desired);
        return gwork.data;
    }

    static void permuteInv(int[] perm, double[] input, double[] output, int N) {
        for (int k = 0; k < N; ++k) {
            output[perm[k]] = input[k];
        }
    }
}
