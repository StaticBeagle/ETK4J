package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.util.RefValue;

import java.util.Arrays;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.ColumnCounts.adjust;
import static com.wildbitsfoundry.etk4j.math.linearalgebra.LUDecompositionSparse.solveColB;

public class QRDecompositionSparse extends QRDecomposition<MatrixSparse> {
    int m, n, m2;

    // storage for Householder vectors
    MatrixSparse H = new MatrixSparse(1, 1, 0);
    // Storage for R matrix in QR
    MatrixSparse R = new MatrixSparse(1, 1, 0);
    // storage for beta in (I - beta*v*v')
    double[] beta = new double[0];
    RefValue.RefDouble Beta = new RefValue.RefDouble(0.0); // used to get return value from a function

    // local workspace
    double[] x = new double[0];

    QrStructuralCounts structure = new QrStructuralCounts();
    int[] structureP = new int[0];
    IGrowArray gwork = new IGrowArray();
    DGrowArray gx = new DGrowArray();

    // if true that means a singular matrix was detected
    boolean singular;

    // true if it has successfully decomposed a matrix
    private boolean decomposed = false;
    // if true then the structure is locked and won't be computed again
    private boolean locked = false;

    public QRDecompositionSparse(MatrixSparse A) {
        super(A);

        // use the same work space to reduce the overall memory foot print
        this.structure.setGwork(gwork);
        this.decompose(A);
    }


    public boolean decompose(MatrixSparse A) {


        if (!decomposed || !locked) {
            // compute the structure of V and R
            if (!structure.process(A))
                return false;

            // Initialize data structured used in the decomposition
            initializeDecomposition(A);
        }
        // perform the decomposition
        performDecomposition(A);

        decomposed = true;
        return true;
    }

    private void performDecomposition(MatrixSparse A) {
        int[] w = gwork.data;
        int[] parent = structure.getParent();
        int[] leftmost = structure.getLeftMost();
        // permutation that was done to ensure all rows have non-zero elements
        int[] pinv_structure = structure.getPinv();
        int s = m2;

        // clear mark nodes. See addRowsInAInToC
        Arrays.fill(w, 0, m2, -1);
        Arrays.fill(x, 0, m2, 0);

        // the counts from structure are actually an upper limit. the actual counts can be lower
        R.nz_length = 0;
        H.nz_length = 0;

        // compute V and R
        for (int k = 0; k < n; k++) {
            R.col_idx[k] = R.nz_length;
            int p1 = H.col_idx[k] = H.nz_length;
            w[k] = k;
            H.nz_rows[H.nz_length++] = k;                       // Add V(k,k) to V's pattern
            int top = n;

            int idx0 = A.col_idx[k];
            int idx1 = A.col_idx[k + 1];

            for (int p = idx0; p < idx1; p++) {
                int i = leftmost[A.nz_rows[p]];
                int len;
                for (len = 0; w[i] != k; i = parent[i]) {
                    w[s + len++] = i;
                    w[i] = k;
                }
                while (len > 0) {
                    w[s + --top] = w[s + --len];
                }
                i = pinv_structure[A.nz_rows[p]];
                x[i] = A.nz_values[p];
                if (i > k && w[i] < k) {
                    H.nz_rows[H.nz_length++] = i;
                    w[i] = k;
                }
            }
            // apply previously computed Householder vectors to the current columns
            for (int p = top; p < n; p++) {
                int i = w[s + p];
                applyHouseholder(H, i, beta[i], x);
                R.nz_rows[R.nz_length] = i;
                R.nz_values[R.nz_length++] = x[i];
                x[i] = 0;
                if (parent[i] == k) {
                    addRowsInAInToC(H, i, H, k, w);
                }
            }
            for (int p = p1; p < H.nz_length; p++) {
                H.nz_values[p] = x[H.nz_rows[p]];
                x[H.nz_rows[p]] = 0;
            }
            R.nz_rows[R.nz_length] = k;
            double max = findMax(H.nz_values, p1, H.nz_length - p1);
            if (max == 0.0) {
                singular = true;
                R.nz_values[R.nz_length] = 0;
                beta[k] = 0;
            } else {
                R.nz_values[R.nz_length] = computeHouseholder(H.nz_values, p1, H.nz_length, max, Beta);
                beta[k] = Beta.getValue();
            }
            R.nz_length++;
        }
        R.col_idx[n] = R.nz_length;
        H.col_idx[n] = H.nz_length;
    }

    public static double findMax(double[] u, int startU, int length) {
        double max = -1;

        int index = startU;
        int stopIndex = startU + length;
        for (; index < stopIndex; index++) {
            double val = u[index];
            val = (val < 0.0) ? -val : val;
            if (val > max)
                max = val;
        }

        return max;
    }

    private void initializeDecomposition(MatrixSparse A) {
        this.singular = false;
        this.m2 = structure.getFicticousRowCount();
        this.m = A.rows;
        this.n = A.cols;

        if (beta.length < n) {
            beta = new double[n];
        }
        if (x.length < m2) {
            x = new double[m2];
            structureP = new int[m2];
        }

        H.reshape(m2, n, structure.nz_in_V);
        R.reshape(m2, n, structure.nz_in_R);
    }

    public MatrixSparse getQ() {
        return getQ(false);
    }

    public MatrixSparse getQEconomy() {
        return getQ(true);
    }

    private MatrixSparse getQ(boolean compact) {
        MatrixSparse Q = new MatrixSparse(1, 1, 0);

        if (compact)
            Q.reshape(H.rows, n, 0);
        else
            Q.reshape(H.rows, m, 0);
        MatrixSparse I = identity(H.rows, Q.cols);

        for (int i = H.cols - 1; i >= 0; i--) {
            rank1UpdateMultR(H, i, beta[i], I, Q, gwork, gx);
            I.setTo(Q);
        }

        // Apply P transpose to Q
        permutationInverse(structure.pinv, structureP, H.rows);
        permuteRowInv(structureP, Q, I);

        // Remove fictitious rows
        if (H.rows > m)
            extractRows(I, 0, m, Q);
        else
            Q.setTo(I);

        return Q;
    }

    public static MatrixSparse identity(int length) {
        return identity(length, length);
    }

    public static MatrixSparse identity(int numRows, int numCols) {
        int min = Math.min(numRows, numCols);
        MatrixSparse A = new MatrixSparse(numRows, numCols, min);
        setIdentity(A);
        return A;
    }

    public static void setIdentity(MatrixSparse A) {
        int min = Math.min(A.rows, A.cols);
        A.growMaxLength(min, false);
        A.nz_length = min;

        Arrays.fill(A.nz_values, 0, min, 1);
        for (int i = 1; i <= min; i++) {
            A.col_idx[i] = i;
            A.nz_rows[i - 1] = i - 1;
        }
        for (int i = min + 1; i <= A.cols; i++) {
            A.col_idx[i] = min;
        }
    }

    /*
     * Creates a submatrix by extracting the specified rows from A. rows = {row0 %le; i %le; row1}.
     *
     * @param A (Input) matrix
     * @param row0 First row. Inclusive
     * @param row1 Last row+1.
     * @param out (Output, Option) Storage for output matrix
     * @return The submatrix
     */
    public static MatrixSparse extractRows(MatrixSparse A, int row0, int row1,
                                           MatrixSparse out) {

        if (out == null)
            out = new MatrixSparse(1, 1, 1);

        out.reshape(row1 - row0, A.cols, A.nz_length);
//        out.col_idx[0] = 0;
//        out.nz_length = 0;

        for (int col = 0; col < A.cols; col++) {
            int idx0 = A.col_idx[col];
            int idx1 = A.col_idx[col + 1];

            for (int i = idx0; i < idx1; i++) {
                int row = A.nz_rows[i];
                if (row >= row0 && row < row1) {
                    out.nz_values[out.nz_length] = A.nz_values[i];
                    out.nz_rows[out.nz_length++] = row - row0;
                }
            }
            out.col_idx[col + 1] = out.nz_length;
        }

        return out;
    }

    /*
     * Computes the inverse permutation vector
     *
     * @param original Original permutation vector
     * @param inverse It's inverse
     */
    public static void permutationInverse(int[] original, int[] inverse, int length) {
        for (int i = 0; i < length; i++) {
            inverse[original[i]] = i;
        }
    }

    /*
     * Applies the row permutation specified by the vector to the input matrix and save the results
     * in the output matrix. output[perm[j],:] = input[j,:]
     *
     * @param permInv (Input) Inverse permutation vector. Specifies new order of the rows.
     * @param input (Input) Matrix which is to be permuted
     * @param output (Output) Matrix which has the permutation stored in it. Is reshaped.
     */
    public static void permuteRowInv(int[] permInv, MatrixSparse input, MatrixSparse output) {
        if (input.rows > permInv.length)
            throw new IllegalArgumentException("permutation vector must have at least as many elements as input has rows");

        output.reshape(input.rows, input.cols, input.nz_length);
        output.nz_length = input.nz_length;
        output.indicesSorted = false;

        System.arraycopy(input.nz_values, 0, output.nz_values, 0, input.nz_length);
        System.arraycopy(input.col_idx, 0, output.col_idx, 0, input.cols + 1);

        int idx0 = 0;
        for (int i = 0; i < input.cols; i++) {
            int idx1 = output.col_idx[i + 1];

            for (int j = idx0; j < idx1; j++) {
                output.nz_rows[j] = permInv[input.nz_rows[j]];
            }
            idx0 = idx1;
        }
    }

    public MatrixSparse getREconomy() {
        MatrixSparse R = new MatrixSparse(0, 0, 0);

        R.setTo(this.R);
        if (m > n) {
            // there should be only zeros past row n
            R.rows = n;
        } else if (n > m && H.rows != m) {
            MatrixSparse tmp = new MatrixSparse(m, n, 0);
            extractRows(R, 0, m, tmp);
            R.setTo(tmp);
        }
        return R;
    }

    public MatrixSparse getH() {
        return H;
    }

    public MatrixSparse getR() {
        return R;
    }

    public double[] getBeta() {
        return beta;
    }

    public double getBeta(int index) {
        if (index >= n)
            throw new IllegalArgumentException("index is out of bounds");
        return beta[index];
    }


    public boolean isSingular() {
        return singular;
    }

    /*
     * <p>Applies a sparse Householder vector to a dense vector.</p>
     * <pre>
     *     x = x - v*(beta*(v'*x))</pre>
     *
     * <P>NOTE: This is the same as cs_happly() in csparse</P>
     *
     * @param V (Input) Matrix containing the Householder
     * @param colV Column in V with the Householder vector
     * @param beta scalar
     * @param x (Input and Output) vector that the Householder is applied to. Modified.
     */
    public static void applyHouseholder(MatrixSparse V, int colV, double beta,
                                        double[] x) {
        int idx0 = V.col_idx[colV];
        int idx1 = V.col_idx[colV + 1];

        // Compute tau = v'*x
        double tau = 0;
        for (int p = idx0; p < idx1; p++) {
            tau += V.nz_values[p] * x[V.nz_rows[p]];
        }
        tau *= beta;

        // x = x - v*tau
        for (int p = idx0; p < idx1; p++) {
            x[V.nz_rows[p]] -= V.nz_values[p] * tau;
        }
    }

    /*
     * <p>
     * Performs a rank-1 update operation on the submatrix specified by V with the multiply on the right.<br>
     * <br>
     * C = (I - &gamma;*v*v<sup>T</sup>)*A<br>
     * </p>
     * <p>
     * The order that matrix multiplies are performed has been carefully selected
     * to minimize the number of operations.
     * </p>
     *
     * <p>
     * Before this can become a truly generic operation the submatrix specification needs
     * to be made more generic.
     * </p>
     */
    static void rank1UpdateMultR(MatrixSparse V, int colV, double gamma,
                                 MatrixSparse A, MatrixSparse C,
                                 IGrowArray gw, DGrowArray gx) {
        if (V.rows != A.rows)
            throw new IllegalArgumentException("Number of rows in V and A must match");

        C.nz_length = 0;
        C.rows = V.rows;
        C.cols = 0;

        for (int i = 0; i < A.cols; i++) {
            double tau = dotInnerColumns(V, colV, A, i, gw, gx);
            addColAppend(1.0, A, i, -gamma * tau, V, colV, C, gw);
        }
    }

    /*
     * Adds the results of adding a column in A and B as a new column in C.<br>
     * C(:,end+1) = &alpha;*A(:,colA) + &beta;*B(:,colB)
     *
     * @param alpha scalar
     * @param A matrix
     * @param colA column in A
     * @param beta scalar
     * @param B matrix
     * @param colB column in B
     * @param C Column in C
     * @param gw workspace
     */
    static void addColAppend(double alpha, MatrixSparse A, int colA, double beta, MatrixSparse B, int colB,
                             MatrixSparse C, IGrowArray gw) {
        if (A.rows != B.rows || A.rows != C.rows)
            throw new IllegalArgumentException("Number of rows in A, B, and C do not match");

        int idxA0 = A.col_idx[colA];
        int idxA1 = A.col_idx[colA + 1];
        int idxB0 = B.col_idx[colB];
        int idxB1 = B.col_idx[colB + 1];

        C.growMaxColumns(++C.cols, true);
        C.growMaxLength(C.nz_length + idxA1 - idxA0 + idxB1 - idxB0, true);

        int[] w = adjust(gw, A.rows);
        Arrays.fill(w, 0, A.rows, -1);

        for (int i = idxA0; i < idxA1; i++) {
            int row = A.nz_rows[i];
            C.nz_rows[C.nz_length] = row;
            C.nz_values[C.nz_length] = alpha * A.nz_values[i];
            w[row] = C.nz_length++;
        }

        for (int i = idxB0; i < idxB1; i++) {
            int row = B.nz_rows[i];
            if (w[row] != -1) {
                C.nz_values[w[row]] += beta * B.nz_values[i];
            } else {
                C.nz_values[C.nz_length] = beta * B.nz_values[i];
                C.nz_rows[C.nz_length++] = row;
            }
        }
        C.col_idx[C.cols] = C.nz_length;
    }


    /*
     * Performs the operation x = x + A(:,i)*alpha
     *
     * <p>NOTE: This is the same as cs_scatter() in csparse.</p>
     */
    public static void multAddColA(MatrixSparse A, int colA,
                                   double alpha,
                                   MatrixSparse C, int mark,
                                   double[] x, int[] w) {
        int idxA0 = A.col_idx[colA];
        int idxA1 = A.col_idx[colA + 1];

        for (int j = idxA0; j < idxA1; j++) {
            int row = A.nz_rows[j];

            if (w[row] < mark) {
                if (C.nz_length >= C.nz_rows.length) {
                    C.growMaxLength(C.nz_length * 2 + 1, true);
                }

                w[row] = mark;
                C.nz_rows[C.nz_length] = row;
                C.col_idx[mark] = ++C.nz_length;
                x[row] = A.nz_values[j] * alpha;
            } else {
                x[row] += A.nz_values[j] * alpha;
            }
        }
    }

    /*
     * Adds rows to C[*,colC] that are in A[*,colA] as long as they are marked in w. This is used to grow C
     * and colC must be the last filled in column in C.
     *
     * <p>NOTE: This is the same as cs_scatter if x is null.</p>
     *
     * @param A Matrix
     * @param colA The column in A that is being examined
     * @param C Matrix
     * @param colC Column in C that rows in A are being added to.
     * @param w An array used to indicate if a row in A should be added to C. if w[i] < colC AND i is a row
     * in A[*,colA] then it will be added.
     */
    public static void addRowsInAInToC(MatrixSparse A, int colA,
                                       MatrixSparse C, int colC,
                                       int[] w) {
        int idxA0 = A.col_idx[colA];
        int idxA1 = A.col_idx[colA + 1];

        for (int j = idxA0; j < idxA1; j++) {
            int row = A.nz_rows[j];

            if (w[row] < colC) {
                if (C.nz_length >= C.nz_rows.length) {
                    C.growMaxLength(C.nz_length * 2 + 1, true);
                }

                w[row] = colC;
                C.nz_rows[C.nz_length++] = row;
            }
        }
        C.col_idx[colC + 1] = C.nz_length;
    }

    /*
     * Computes the inner product of two column vectors taken from the input matrices.
     *
     * <p>dot = A(:,colA)'*B(:,colB)</p>
     *
     *
     * @param A Matrix
     * @param colA Column in A
     * @param B Matrix
     * @param colB Column in B
     * @return Dot product
     */
    static double dotInnerColumns(MatrixSparse A, int colA, MatrixSparse B, int colB,
                                  IGrowArray gw, DGrowArray gx) {
        if (A.rows != B.rows)
            throw new IllegalArgumentException("Number of rows must match.");

        int[] w = adjust(gw, A.rows);
        Arrays.fill(w, 0, A.rows, -1);
        double[] x = adjust(gx, A.rows);

        int length = 0;

        int idx0 = A.col_idx[colA];
        int idx1 = A.col_idx[colA + 1];
        for (int i = idx0; i < idx1; i++) {
            int row = A.nz_rows[i];
            x[length] = A.nz_values[i];
            w[row] = length++;
        }

        double dot = 0;

        idx0 = B.col_idx[colB];
        idx1 = B.col_idx[colB + 1];
        for (int i = idx0; i < idx1; i++) {
            int row = B.nz_rows[i];
            if (w[row] != -1) {
                dot += x[w[row]] * B.nz_values[i];
            }
        }

        return dot;
    }

    /*
     * Creates a householder reflection.
     *
     * (I-gamma*v*v')*x = tau*e1
     *
     * <p>NOTE: Same as cs_house in csparse</p>
     *
     * @param x (Input) Vector x (Output) Vector v. Modified.
     * @param xStart First index in X that is to be processed
     * @param xEnd Last + 1 index in x that is to be processed.
     * @param gamma (Output) Storage for computed gamma
     * @return variable tau
     */
    public static double computeHouseholder(double[] x, int xStart, int xEnd, double max, RefValue.RefDouble gamma) {
        double tau = 0;
        for (int i = xStart; i < xEnd; i++) {
            double val = x[i] /= max;
            tau += val * val;
        }
        tau = Math.sqrt(tau);
        if (x[xStart] < 0) {
            tau = -tau;
        }
        double u_0 = x[xStart] + tau;
        gamma.setValue(u_0 / tau);
        x[xStart] = 1;
        for (int i = xStart + 1; i < xEnd; i++) {
            x[i] /= u_0;
        }

        return -tau * max;
    }

    public MatrixSparse solve(double[] b) {
        MatrixSparse B = new MatrixSparse(b.length, 1, b.length);
        for(int i = 0; i < b.length; i++) {
            B.unsafeSet(i, 0, b[i]);
        }
        return solve(B);
    }

    public MatrixSparse solve(MatrixSparse B) {
        MatrixSparse X = new MatrixSparse(0, 0, 0);
        X.reshape(cols, B.cols, X.rows);

        IGrowArray gw1 = gwork;
        MatrixSparse tmp = new MatrixSparse(1, 1, 1);
        IGrowArray gw = new IGrowArray();

        // Don't modify the input
        tmp.setTo(B);
        B = tmp;
        MatrixSparse B_tmp = B.createLike();
        MatrixSparse swap;

        // Apply permutation to B
        int[] pinv = this.structure.getPinv();
        permuteRowInv(pinv, B, B_tmp);
        swap = B_tmp;
        B_tmp = B;
        B = swap;

        // Apply house holders to B
        MatrixSparse V = this.getH();
        for (int i = 0; i < cols; i++) {
            rank1UpdateMultR(V, i, beta[i], B, B_tmp, gw, gx);
            swap = B_tmp;
            B_tmp = B;
            B = swap;
        }

        // Solve for X
        TriangularSystemSolver.solve(R, false, B, X, null, gx, gw, gw1);
		return X;
    }

}
