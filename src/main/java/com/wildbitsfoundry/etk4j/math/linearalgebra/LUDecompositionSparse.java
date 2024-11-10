package com.wildbitsfoundry.etk4j.math.linearalgebra;

import java.util.Arrays;

public class LUDecompositionSparse extends LUDecomposition<MatrixSparse> {

    public LUDecompositionSparse(MatrixSparse matrix) {
        super(matrix);
        this.initialize(matrix);
        this.performLU(matrix);
    }

    private final IGrowArray gw = new IGrowArray();
    private final MatrixSparse L = new MatrixSparse(0, 0, 0);
    private final MatrixSparse U = new MatrixSparse(0, 0, 0);
    private double[] x = new double[0];

    private void initialize(MatrixSparse A) {
        int m = A.rows;
        int n = A.cols;
        int o = Math.min(m, n);
        this.L.reshape(m, m, 4 * A.nz_length + o);
        this.L.nz_length = 0;
        this.U.reshape(m, n, 4 * A.nz_length + o);
        this.U.nz_length = 0;
        this.singular = false;
        if (this.pinv.length != m) {
            this.pinv = new int[m];
            this.x = new double[m];
        }

        for (int i = 0; i < m; ++i) {
            this.pinv[i] = -1;
            this.L.col_idx[i] = 0;
        }

    }

    private final IGrowArray gxi = new IGrowArray();
    private int[] pinv = new int[0];

    private boolean singular;

    private boolean performLU(MatrixSparse A) {
        int m = A.rows;
        int n = A.cols;
        int[] w = adjust(this.gw, m * 2, m);

        int k;
        for (k = 0; k < n; ++k) {
            this.L.col_idx[k] = this.L.nz_length;
            this.U.col_idx[k] = this.U.nz_length;
            if (this.L.nz_length + n > this.L.nz_values.length) {
                this.L.growMaxLength(2 * this.L.nz_values.length + n, true);
            }

            if (this.U.nz_length + n > this.U.nz_values.length) {
                this.U.growMaxLength(2 * this.U.nz_values.length + n, true);
            }

            int col = k;
            int top = solveColB(this.L, true, A, col, this.x, this.pinv, this.gxi, w);
            int[] xi = this.gxi.data;
            int ipiv = -1;
            double a = -1.7976931348623157E308;

            for (int p = top; p < n; ++p) {
                int i = xi[p];
                if (this.pinv[i] < 0) {
                    double t;
                    if ((t = Math.abs(this.x[i])) > a) {
                        a = t;
                        ipiv = i;
                    }
                } else {
                    this.U.nz_rows[this.U.nz_length] = this.pinv[i];
                    this.U.nz_values[this.U.nz_length++] = this.x[i];
                }
            }

            if (ipiv == -1 || a <= 0.0) {
                this.singular = true;
                return false;
            }

            double pivot = this.x[ipiv];
            this.U.nz_rows[this.U.nz_length] = k;
            this.U.nz_values[this.U.nz_length++] = pivot;
            this.pinv[ipiv] = k;
            this.L.nz_rows[this.L.nz_length] = ipiv;
            this.L.nz_values[this.L.nz_length++] = 1.0;

            for (int p = top; p < n; ++p) {
                int i = xi[p];
                if (this.pinv[i] < 0) {
                    this.L.nz_rows[this.L.nz_length] = i;
                    this.L.nz_values[this.L.nz_length++] = this.x[i] / pivot;
                }

                this.x[i] = 0.0;
            }
        }

        this.L.col_idx[n] = this.L.nz_length;
        this.U.col_idx[n] = this.U.nz_length;

        for (k = 0; k < this.L.nz_length; ++k) {
            this.L.nz_rows[k] = this.pinv[this.L.nz_rows[k]];
        }

        return true;
    }

    private final DGrowArray gx = new DGrowArray();
    private final DGrowArray gb = new DGrowArray();

    public int[] getPinv() {
        return this.pinv;
    }

    public MatrixSparse getL() {
        return this.L;
    }

    public MatrixSparse getU() {
        return this.U;
    }

    @Override
    public boolean isSingular() {
        return singular;
    }

    @Override
    public MatrixSparse solve(double[] b) {
        if (b.length != this.rows) {
            int var10002 = b.length;
            throw new IllegalArgumentException("Unexpected number of rows in B based on shape of A. Found=" + var10002 + " Expected=" + this.rows);
        }

        double[] x = adjust(this.gx, b.length);
        double[] rhs = adjust(this.gb, b.length);
        MatrixSparse L = this.L;
        MatrixSparse U = this.U;

        System.arraycopy(b, 0, rhs, 0, b.length);

        permuteInv(pinv, rhs, x, b.length);

        solveL(L, x);
        solveU(U, x);

        double[][] xmat = new double[b.length][1];
        for(int i = 0; i < b.length; i++) {
            xmat[i][0] = x[i];
        }
        return MatrixSparse.from2dArray(xmat, 0);
    }

    private static int solveColB(MatrixSparse G, boolean lower, MatrixSparse B, int colB, double[] x, int[] pinv, IGrowArray g_xi, int[] w) {
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
        int N = G.cols;
        int head = 0;
        xi[head] = rowB;

        while (head >= 0) {
            int G_col = xi[head];
            int G_col_new = pinv != null ? pinv[G_col] : G_col;
            if (w[G_col] == 0) {
                w[G_col] = 1;
                w[N + head] = G_col_new >= 0 && G_col_new < N ? G.col_idx[G_col_new] : 0;
            }

            boolean done = true;
            int idx0 = w[N + head];
            int idx1 = G_col_new >= 0 && G_col_new < N ? G.col_idx[G_col_new + 1] : 0;

            for (int j = idx0; j < idx1; ++j) {
                int jrow = G.nz_rows[j];
                if (jrow < N && w[jrow] == 0) {
                    w[N + head] = j + 1;
                    ++head;
                    xi[head] = jrow;
                    done = false;
                    break;
                }
            }

            if (done) {
                --head;
                --top;
                xi[top] = G_col;
            }
        }

        return top;
    }

    private static void solveL(MatrixSparse L, double[] x) {
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

    private static void solveU(MatrixSparse U, double[] x) {
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

    private static int[] adjust(IGrowArray gwork, int desired, int zeroToM) {
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


    private static double[] adjust(DGrowArray gwork, int desired) {
        if (gwork == null) {
            gwork = new DGrowArray();
        }

        gwork.reshape(desired);
        return gwork.data;
    }

    private static void permuteInv(int[] perm, double[] input, double[] output, int N) {
        for (int k = 0; k < N; ++k) {
            output[perm[k]] = input[k];
        }
    }
}
