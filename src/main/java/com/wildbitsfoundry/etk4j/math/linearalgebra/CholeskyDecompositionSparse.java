package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.ColumnCounts.adjust;
import static com.wildbitsfoundry.etk4j.math.linearalgebra.QrStructuralCounts.eliminationTree;
import static com.wildbitsfoundry.etk4j.math.linearalgebra.QrStructuralCounts.postorder;

public class CholeskyDecompositionSparse extends CholeskyDecomposition<MatrixSparse> {

    private int N;
    private boolean isSPD;

    // storage for decomposition
    MatrixSparse L = new MatrixSparse(1, 1, 0);

    // workspace storage
    IGrowArray gw = new IGrowArray(1);
    IGrowArray gs = new IGrowArray(1);
    DGrowArray gx = new DGrowArray(1);
    int[] parent = new int[1];
    int[] post = new int[1];
    int[] counts = new int[1];
    ColumnCounts columnCounter = new ColumnCounts(false);

    // true if it has successfully decomposed a matrix
    private boolean decomposed = false;
    // if true then the structure is locked and won't be computed again
    private boolean locked = false;

    public CholeskyDecompositionSparse(MatrixSparse matrix) {
        super(matrix);
        decompose(matrix);
    }

    private boolean decompose(MatrixSparse orig) {
        if (orig.cols != orig.rows) {
            throw new NonSquareMatrixException("Must be a square matrix");
        }

        if (!locked || !decomposed)
            performSymbolic(orig);

        performDecomposition(orig);
        if (isSPD) {
            decomposed = true;
            return true;
        } else {
            return false;
        }
    }

    private void performSymbolic(MatrixSparse A) {
        init(A.cols);

        eliminationTree(A, false, parent, gw);
        postorder(parent, N, post, gw);
        columnCounter.process(A, parent, post, counts);
        L.reshape(A.rows, A.cols, 0);
        L.histogramToStructure(counts);
    }

    private void init(int N) {
        this.N = N;
        if (parent.length < N) {
            parent = new int[N];
            post = new int[N];
            counts = new int[N];
            gw.reshape(3 * N);
        }
    }

    private void performDecomposition(MatrixSparse A) {
        int[] c = adjust(gw, N);
        int[] s = adjust(gs, N);
        double[] x = adjust(gx, N);

        System.arraycopy(L.col_idx, 0, c, 0, N);

        for (int k = 0; k < N; k++) {
            //----  Nonzero pattern of L(k,:)
            int top = searchNzRowsElim(A, k, parent, s, c);

            // x(0:k) is now zero
            x[k] = 0;
            int idx0 = A.col_idx[k];
            int idx1 = A.col_idx[k + 1];

            // x = full(triu(C(:,k)))
            for (int p = idx0; p < idx1; p++) {
                if (A.nz_rows[p] <= k) {
                    x[A.nz_rows[p]] = A.nz_values[p];
                }
            }
            double d = x[k]; // d = C(k,k)
            x[k] = 0; // clear x for k+1 iteration

            //---- Triangular Solve
            for (; top < N; top++) {
                int i = s[top];
                double lki = x[i] / L.nz_values[L.col_idx[i]]; // L(k,i) = x(i) / L(i,i)
                x[i] = 0;
                for (int p = L.col_idx[i] + 1; p < c[i]; p++) {
                    x[L.nz_rows[p]] -= L.nz_values[p] * lki;
                }
                d -= lki * lki; // d = d - L(k,i)**L(k,i)
                int p = c[i]++;
                L.nz_rows[p] = k;     // store L(k,i) in column i
                L.nz_values[p] = lki;
            }

            //----- Compute L(k,k)
            if (d <= 0) {
                // it's not positive definite
                this.isSPD = false;
            }
            int p = c[k]++;
            L.nz_rows[p] = k;
            L.nz_values[p] = Math.sqrt(d);
        }

        this.isSPD = true;
    }

    public boolean isSPD() {
        return isSPD;
    }

    /**
     * <p>Given an elimination tree compute the non-zero elements in the specified row of L given the
     * symmetric A matrix. This is in general much faster than general purpose algorithms</p>
     *
     * <p>Functionally equivalent to cs_ereach() in csparse</p>
     *
     * @param A      Symmetric matrix.
     * @param k      Row in A being processed.
     * @param parent elimination tree.
     * @param s      (Output) s[top:(n-1)] = pattern of L[k,:]. Must have length A.cols
     * @param w      workspace array used internally. All elements must be &ge; 0 on input. Must be of size A.cols
     * @return Returns the index of the first element in the xi list. Also known as top.
     */
    private static int searchNzRowsElim(MatrixSparse A, int k, int[] parent, int[] s, int[] w) {
        int top = A.cols;

        // Traversing through the column in A is the same as the row in A since it's symmetric
        int idx0 = A.col_idx[k], idx1 = A.col_idx[k + 1];

        w[k] = -w[k] - 2;  // makr node k as visited
        for (int p = idx0; p < idx1; p++) {
            int i = A.nz_rows[p];   // A[k,i] is not zero

            if (i > k) // only consider upper triangular part of A
                continue;

            // move up the elimination tree
            int len = 0;
            for (; w[i] >= 0; i = parent[i]) {
                s[len++] = i; // L[k,i] is not zero
                w[i] = -w[i] - 2; // mark i as being visited
            }
            while (len > 0) {
                s[--top] = s[--len];
            }
        }

        // unmark all nodes
        for (int p = top; p < A.cols; p++) {
            w[s[p]] = -w[s[p]] - 2;
        }
        w[k] = -w[k] - 2;
        return top;
    }

    public Complex computeDeterminant() {
        double value = 1;
        for (int i = 0; i < N; i++) {
            value *= L.nz_values[L.col_idx[i]];
        }
        return new Complex(value * value, 0);
    }

    public DGrowArray getGx() {
        return gx;
    }

    public MatrixSparse getL() {
        return L;
    }

    public MatrixSparse getR() { return L.transpose(); }

    public MatrixSparse solve(MatrixSparse B) {
        MatrixSparse X = new MatrixSparse(1, 1, 1);
        X.reshape(cols, B.cols, X.rows);

        IGrowArray gw1 = gw;

        MatrixSparse L = this.L;

        MatrixSparse tmp = new MatrixSparse(1, 1, 1);
        tmp.reshape(L.rows, B.cols, 1);

        TriangularSystemSolver.solve(L, true, B, tmp, null, new DGrowArray(), new IGrowArray(), gw1);
        TriangularSystemSolver.solveTran(L, tmp, X, null, new DGrowArray(), new IGrowArray(), gw1);
        return X;
    }

    public MatrixSparse solve(double[] b) {
        double[][] matrix = new double[b.length][1];
        for(int i = 0; i < b.length; i++) {
            matrix[i][0] = b[i];
        }
        MatrixSparse B = MatrixSparse.from2DArray(matrix);
        return solve(B);
    }
}