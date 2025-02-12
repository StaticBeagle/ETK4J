package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import java.util.Arrays;

import static com.wildbitsfoundry.etk4j.math.linearalgebra.ComplexTriangularSolverDense.backSubstitutionSolve;
import static com.wildbitsfoundry.etk4j.util.ComplexArrays.zeros;

/***
 * References <a href="https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition/">Mathworks blog</a>
 */
public class ComplexQRDecompositionDense extends ComplexQRDecomposition<ComplexMatrixDense> {
    protected Complex[] _data;

    private final int m;
    private final int n;


    private final ComplexMatrixDense H;
    private final ComplexMatrixDense R;

    public ComplexQRDecompositionDense(ComplexMatrixDense matrix) {
        super(matrix);

        this.m = matrix.getColumnCount();
        this.n = matrix.getRowCount();
        ComplexMatrixDense R = matrix.copy();
        ComplexMatrixDense U = new ComplexMatrixDense(zeros(m, n));
        for (int i = 0; i < Math.min(m, n); i++) {
            Complex[] rCol = Arrays.copyOfRange(getColumn(R, i), i, m);
            Complex[] u = houseGen(rCol);
            setColumn(U, u, i, i);
            Complex[][] subR = getSubMatrix(R, i, m, i, n);
            Complex[][] H = calculateReflector(u, subR);
            setSubMatrix(R, H, i, m, i, n);
            for (int j = i + 1; j < m; j++) {
                R.unsafeSet(j, i, new Complex());
            }
        }
        this.R = R;
        this.H = U;
    }

    /*
     * ------------------------ Public Methods ------------------------
     */

    /**
     * Is the matrix full rank?
     *
     * @return true if R, and hence A, has full rank.
     */
    public boolean isFullRank() {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j && R.unsafeGet(i, j).abs() == 0) {
                    return false;
                }
            }
        }
        return true;
    }


    /**
     * Return the Householder vectors
     *
     * @return Lower trapezoidal matrix whose columns define the reflections
     */
    public ComplexMatrixDense getH() {
        return H.copy();
    }

    /**
     * Return the upper triangular factor
     *
     * @return R
     */

    public ComplexMatrixDense getR() {
        return R.copy();
    }

    /**
     * Generate and return the unitary orthogonal factor
     *
     * @return Q
     */
    public ComplexMatrixDense getQ() {
        ComplexMatrixDense I = ComplexMatrixDense.Factory.identity(H.getRowCount(), H.getColumnCount());
        return houseApply(H, I);
    }

    public ComplexMatrixDense QmultiplyX(ComplexMatrixDense X) {
        return houseApply(H, X);
    }

    /**
     * Generate and return the conjugate transpose of the orthogonal factor
     *
     * @return transpose(Q)
     */
    public ComplexMatrixDense getQH() {
        ComplexMatrixDense I = ComplexMatrixDense.Factory.identity(H.getRowCount(), H.getColumnCount());
        return houseApplyTranspose(H, I);
    }

    /**
     * Least squares solution of A*X = B
     *
     * @param B A Matrix with as many rows as A and any number of columns.
     * @return X that minimizes the two norm of Q*R*X-B.
     * @throws IllegalArgumentException Matrix row dimensions must agree.
     * @throws RuntimeException         Matrix is rank deficient.
     */

    public ComplexMatrixDense solve(ComplexMatrixDense B) {
        if (B.getRowCount() != rows) {
            throw new IllegalArgumentException("Matrix row dimensions must agree.");
        }
        if (!this.isFullRank()) {
            throw new RuntimeException("Matrix is rank deficient.");
        }
        // Compute Y = transpose(Q) * B
        ComplexMatrixDense Y = houseApplyTranspose(H, B);

        // Solve R * X = Y;
        // Back Substitution
        return backSubstitutionSolve(R, Y);
    }

    /**
     * Least squares solution of A*X = B
     *
     * @param B A Matrix with as many rows as A and any number of columns.
     * @return X that minimizes the two norm of Q*R*X-B.
     * @throws IllegalArgumentException Matrix row dimensions must agree.
     * @throws RuntimeException         Matrix is rank deficient.
     */

    public ComplexMatrixDense solve(MatrixDense B) {
        return solve(ComplexMatrixDense.fromRealMatrix(B));
    }

    private static Complex sig(Complex u) {
        return u.sign().add(u.equals(new Complex()) ? 1 : 0);
    }

    private static Complex[] houseGen(Complex[] x) {
        double nu = ComplexArrays.norm(x);
        if (nu != 0) {
            Complex[] u = ComplexArrays.divideElementWise(x, nu);
            u[0] = u[0].add(sig(u[0]));
            return ComplexArrays.divideElementWise(u, Math.sqrt(u[0].abs()));
        } else {
            Complex[] u = ComplexArrays.deepCopy(x);
            u[1] = Complex.fromReal(Math.sqrt(2));
            return u;
        }
    }

    private static Complex[][] calculateReflector(Complex[] u, Complex[][] x) {
        Complex[] uH = Arrays.stream(u).map(Complex::conj).toArray(Complex[]::new);
        Complex[] uHx = zeros(Math.max(x.length, x[0].length));
        for (int i = 0; i < x[0].length; i++) {
            for (int j = 0; j < x.length; j++) {
                // (uH * x)
                uHx[i].addEquals(uH[j].multiply(x[j][i]));
            }
        }
        Complex[][] uuHx = zeros(x.length, x[0].length);
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < x[0].length; j++) {
                uuHx[i][j] = u[i].multiply(uHx[j]);
            }
        }

        Complex[][] H = new Complex[x.length][x[0].length];
        for (int i = 0; i < x.length; i++) {
            H[i] = ComplexArrays.subtractElementWise(x[i], uuHx[i]);
        }
        return H;
    }

    private static Complex[] getColumn(ComplexMatrixDense A, int col) {
        Complex[] result = new Complex[A.getRowCount()];
        for (int i = 0; i < result.length; i++) {
            result[i] = A.unsafeGet(i, col);
        }
        return result;
    }

    private static void setColumn(ComplexMatrixDense A, Complex[] values, int col, int startingRow) {
        for (int i = 0; i < A.getRowCount() - startingRow; i++) {
            A.set(i + startingRow, col, values[i]);
        }
    }

    private static Complex[][] getSubMatrix(ComplexMatrixDense A, int row0, int row1, int col0, int col1) {
        Complex[][] result = new Complex[row1 - row0][col1 - col0];
        for (int i = row0; i < row1; i++) {
            for (int j = col0; j < col1; j++) {
                result[i - row0][j - col0] = A.get(i, j);
            }
        }
        return result;
    }

    private static void setSubMatrix(ComplexMatrixDense A, Complex[][] values, int row0, int row1, int col0, int col1) {
        for (int i = row0; i < row1; i++) {
            for (int j = col0; j < col1; j++) {
                A.set(i, j, values[i - row0][j - col0]);
            }
        }
    }

    /*
     % Apply Householder reflections.
    % Z = house_apply(U,X), with U from house_qr
    % computes Q*X without actually computing Q.
     */
    private static ComplexMatrixDense houseApply(ComplexMatrixDense U, ComplexMatrixDense X) {
        ComplexMatrixDense Z = X.copy();
        int n = U.getColumnCount();
        Complex[][] H = new Complex[X.getRowCount()][X.getColumnCount()];
        for (int i = 0; i < Z.getRowCount(); i++) {
            for (int j = 0; j < Z.getColumnCount(); j++) {
                H[i][j] = Z.unsafeGet(i, j);
            }
        }
        for (int i = n - 1; i >= 0; i--) {
            H = calculateReflector(getColumn(U, i), H);
        }
        return new ComplexMatrixDense(H);
    }

    /*
     % Apply Householder transposed reflections.
    % Z = house_apply(U,X), with U from house_qr
    % computes Q'*X without actually computing Q'.
     */
    private static ComplexMatrixDense houseApplyTranspose(ComplexMatrixDense U, ComplexMatrixDense X) {
        ComplexMatrixDense Z = X.copy();
        int n = U.getColumnCount();
        Complex[][] H = new Complex[X.getRowCount()][X.getColumnCount()];
        for (int i = 0; i < Z.getRowCount(); i++) {
            for (int j = 0; j < Z.getColumnCount(); j++) {
                H[i][j] = Z.unsafeGet(i, j);
            }
        }
        for (int i = 0; i < n; i++) {
            H = calculateReflector(getColumn(U, i), H);
        }
        return new ComplexMatrixDense(H);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rows * cols; ++i) {
            if (i > 0 && i % cols == 0) {
                sb.append(System.lineSeparator());
            }
            sb.append(_data[i]).append(" ");
        }
        return sb.toString();
    }
}
