package com.wildbitsfoundry.etk4j.math.linearalgebra;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.ColumnCounts.adjust;
import static com.wildbitsfoundry.etk4j.math.linearalgebra.LUDecompositionSparse.solveColB;

class TriangularSystemSolver {

    private TriangularSystemSolver() {}

    /*
     * Computes the solution to the triangular system.
     *
     * @param G (Input) Lower or upper triangular matrix. diagonal elements must be non-zero. Not modified.
     * @param lower true for lower triangular and false for upper
     * @param B (Input) Matrix. Not modified.
     * @param X (Output) Solution
     * @param pinv (Input, Optional) Permutation vector. Maps col j to G. Null if no pivots.
     * @param g_x (Optional) Storage for workspace.
     * @param g_xi (Optional) Storage for workspace.
     * @param g_w (Optional) Storage for workspace.
     */
    static void solve(MatrixSparse G, boolean lower,
                      MatrixSparse B, MatrixSparse X,
                      int[] pinv,
                      DGrowArray g_x, IGrowArray g_xi, IGrowArray g_w) {
        double[] x = adjust(g_x, G.rows);
        if (g_xi == null) g_xi = new IGrowArray();
        int[] xi = adjust(g_xi, G.rows);
        int[] w = adjust(g_w, G.cols * 2, G.cols);

        X.nz_length = 0;
        X.col_idx[0] = 0;
        X.indicesSorted = false;

        for (int colB = 0; colB < B.cols; colB++) {
            int top = solveColB(G, lower, B, colB, x, pinv, g_xi, w);

            int nz_count = X.rows - top;
            if (X.nz_values.length < X.nz_length + nz_count) {
                X.growMaxLength(X.nz_length * 2 + nz_count, true);
            }

            for (int p = top; p < X.rows; p++, X.nz_length++) {
                X.nz_rows[X.nz_length] = xi[p];
                X.nz_values[X.nz_length] = x[xi[p]];
            }
            X.col_idx[colB + 1] = X.nz_length;
        }
    }

    /**
     * Solution to a sparse transposed triangular system with sparse B and sparse X
     *
     * <p>G<sup>T</sup>*X = B</p>
     *
     * @param G    (Input) Lower or upper triangular matrix. diagonal elements must be non-zero. Not modified.
     * @param B    (Input) Matrix. Not modified.
     * @param X    (Output) Solution
     * @param pinv (Input, Optional) Permutation vector. Maps col j to G. Null if no pivots.
     * @param g_x  (Optional) Storage for workspace.
     * @param g_xi (Optional) Storage for workspace.
     * @param g_w  (Optional) Storage for workspace.
     */
    static void solveTran(MatrixSparse G,
                                  MatrixSparse B, MatrixSparse X,
                                  int[] pinv,
                                  DGrowArray g_x, IGrowArray g_xi, IGrowArray g_w) {
        double[] x = adjust(g_x, G.rows);

        X.zero();
        X.indicesSorted = false;

        // storage for the index of non-zero rows in X
        int[] xi = adjust(g_xi, G.rows);
        // Used to mark nodes as non-zero or not. Fill with zero initially
        int[] w = adjust(g_w, G.cols, G.cols); // Dense fill makes adds O(N) to runtime

        for (int colB = 0; colB < B.cols; colB++) {
            int idx0 = B.col_idx[colB];
            int idx1 = B.col_idx[colB + 1];

            // Sparse copy into X and mark elements are non-zero
            int X_nz_count = 0;
            for (int i = idx0; i < idx1; i++) {
                int row = B.nz_rows[i];
                x[row] = B.nz_values[i];
                w[row] = 1;
                xi[X_nz_count++] = row;
            }

            if (true) {
                for (int col = G.rows - 1; col >= 0; col--) {
                    X_nz_count = solveTranColumn(G, x, xi, w, pinv, X_nz_count, col);
                }
            } else {
                for (int col = 0; col < G.rows; col++) {
                    X_nz_count = solveTranColumn(G, x, xi, w, pinv, X_nz_count, col);
                }
            }

            // set everything back to zero for the next column
            if (colB + 1 < B.cols) {
                for (int i = 0; i < X_nz_count; i++) {
                    w[xi[i]] = 0;
                }
            }

            // Copy results into X
            if (X.nz_values.length < X.nz_length + X_nz_count) {
                X.growMaxLength(X.nz_length * 2 + X_nz_count, true);
            }
            for (int p = 0; p < X_nz_count; p++, X.nz_length++) {
                X.nz_rows[X.nz_length] = xi[p];
                X.nz_values[X.nz_length] = x[xi[p]];
            }
            X.col_idx[colB + 1] = X.nz_length;
        }
    }

    static int solveTranColumn(MatrixSparse G, double[] x, int[] xi, int[] w,
                                       int[] pinv, int x_nz_count, int col) {
        int idxG0 = G.col_idx[col];
        int idxG1 = G.col_idx[col + 1];

        int indexDiagonal = -1;
        double total = 0;
        for (int j = idxG0; j < idxG1; j++) {
            int J = pinv != null ? pinv[j] : j;
            int row = G.nz_rows[J];

            if (row == col) {
                // order matters and this operation needs to be done last
                indexDiagonal = j;
            } else if (w[row] == 1) {
                // L'[ col , row]*x[row]
                total += G.nz_values[J] * x[row];
            }
        }
        if (w[col] == 1) {
            x[col] = (x[col] - total) / G.nz_values[indexDiagonal];
        } else if (total != 0) {
            // This element in B was zero. Mark it as non-zero and add to list
            w[col] = 1;
            x[col] = -total / G.nz_values[indexDiagonal];
            xi[x_nz_count++] = col;
        }
        return x_nz_count;
    }

    static void solveL(MatrixSparse L, double[] x) {
        int N = L.cols;
        int idx0 = L.col_idx[0];

        for (int col = 0; col < N; ++col) {
            int idx1 = L.col_idx[col + 1];
            double x_j = x[col] /= L.nz_values[idx0];

            for (int i = idx0 + 1; i < idx1; ++i) {
                int row = L.nz_rows[i];
                x[row] -= L.nz_values[i] * x_j;
            }

            idx0 = idx1;
        }
    }

    static void solveU(MatrixSparse U, double[] x) {
        int N = U.cols;
        int idx1 = U.col_idx[N];

        for (int col = N - 1; col >= 0; --col) {
            int idx0 = U.col_idx[col];
            double x_j = x[col] /= U.nz_values[idx1 - 1];

            for (int i = idx0; i < idx1 - 1; ++i) {
                int row = U.nz_rows[i];
                x[row] -= U.nz_values[i] * x_j;
            }

            idx1 = idx0;
        }
    }
}
