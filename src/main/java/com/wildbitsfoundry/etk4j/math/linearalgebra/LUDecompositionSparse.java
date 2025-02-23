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
        int[] w = MatrixSparseUtils.adjust(this.gw, m * 2, m);

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
            int top = MatrixSparseUtils.solveColB(this.L, true, A, col, this.x, this.pinv, this.gxi, w);
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
        if (this.isSingular()) {
            throw new RuntimeException("Matrix is singular.");
        }

        double[] x = MatrixSparseUtils.adjust(this.gx, b.length);
        double[] rhs = MatrixSparseUtils.adjust(this.gb, b.length);
        MatrixSparse L = this.L;
        MatrixSparse U = this.U;

        System.arraycopy(b, 0, rhs, 0, b.length);

        MatrixSparseUtils.permuteInv(pinv, rhs, x, b.length);

        TriangularSystemSolver.solveL(L, x);
        TriangularSystemSolver.solveU(U, x);

        double[][] xmat = new double[b.length][1];
        for (int i = 0; i < b.length; i++) {
            xmat[i][0] = x[i];
        }
        return MatrixSparse.from2DArray(xmat, 0);
    }

    public MatrixSparse solve(MatrixSparse B) {
        if (B.getRowCount() != this.rows) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (this.isSingular()) {
            throw new RuntimeException("Matrix is singular.");
        }

        MatrixSparse X = new MatrixSparse(0, 0, 0);
        X.reshape(cols, B.cols, X.rows);

        MatrixSparse tmp = new MatrixSparse(1, 1, 1);
        X.reshape(cols, B.cols, X.rows);

        MatrixSparse L = this.L;
        MatrixSparse U = this.U;

        MatrixSparse Bp = new MatrixSparse(1, 1, 1);

        // these are row pivots
        Bp.reshape(B.rows, B.cols, B.nz_length);
        int[] Pinv = this.getPinv();
        MatrixSparseUtils.permute(Pinv, B, null, Bp);
        DGrowArray gx = this.gx;
        IGrowArray gw = this.gw;
        IGrowArray gw1 = this.gxi;

        tmp.reshape(L.rows, B.cols, 1);

        TriangularSystemSolver.solve(L, true, Bp, tmp, null, gx, gw, gw1);
        TriangularSystemSolver.solve(U, false, tmp, X, null, gx, gw, gw1);
        return X;
    }
}
