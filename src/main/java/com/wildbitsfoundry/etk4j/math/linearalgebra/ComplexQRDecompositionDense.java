package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.util.ComplexArrays;

import java.util.Arrays;

/***
 * References https://blogs.mathworks.com/cleve/2016/10/03/householder-reflections-and-the-qr-decomposition/
 */
public class ComplexQRDecompositionDense extends ComplexQRDecomposition<ComplexMatrixDense> {
    protected Complex[] _data;

    private int m = 0;
    private int n = 0;

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

    private static Complex[][] zeros(int m, int n) {
        Complex[][] zeros = new Complex[m][n];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                zeros[i][j] = new Complex();
            }
        }
        return zeros;
    }

    private static Complex[] zeros(int dim) {
        Complex[] zeros = new Complex[dim];
        for (int i = 0; i < dim; i++) {
            zeros[i] = new Complex();
        }
        return zeros;
    }

    private static Complex[][] calculateReflector(Complex[] u, Complex[][] x) {
        Complex[] uH = Arrays.stream(u).map(c -> c.conj()).toArray(Complex[]::new);
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
        for (int i = 0; i < A.getRowCount() - startingRow; i++) { // TODO could be endRow - startingRow if we pass a param
            A.set(i + startingRow, col, values[i]);
        }
    }
    private static void setColumn(ComplexMatrixDense A, Complex[] values, int col) {
        setColumn(A, values, col, 0);
    }

    private static void setColumn(ComplexMatrixDense A, Complex[] values, int row0, int row1, int col) {
        for (int i = row0; i < row1; i++) {
            A.set(i, col, values[i]);
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
        // TODO check for row0 == row1
    }

    /*
     % Apply Householder reflections.
    % Z = house_apply(U,X), with U from house_qr
    % computes Q*X without actually computing Q.
     */
    public static ComplexMatrixDense houseApply(ComplexMatrixDense U, ComplexMatrixDense X) {
        ComplexMatrixDense Z = X.copy();
        int n = U.getColumnCount();
        Complex[][] H = new Complex[X.getRowCount()][X.getColumnCount()];
        // TODO temporary 2D array method
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
    public static ComplexMatrixDense houseApplyTranspose(ComplexMatrixDense U, ComplexMatrixDense X) {
        ComplexMatrixDense Z = X.copy();
        int n = U.getColumnCount();
        Complex[][] H = new Complex[X.getRowCount()][X.getColumnCount()];
        // TODO temporary 2D array method
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

    public static ComplexMatrixDense backSubstitutionSolve(ComplexMatrixDense R, ComplexMatrixDense B) {
        int n = R.getRowCount();
        int m = B.getColumnCount();

        Complex[][] X = zeros(n, m);
        for (int col = 0; col < m; col++) {
            for (int i = n - 1; i >= 0; i--) {
                X[i][col] = B.get(i, col);
                for (int j = i + 1; j < n; j++) {
                    X[i][col].subtractEquals(R.get(i, j).multiply(X[j][col]));
                }
                X[i][col].divideEquals(R.get(i, i));
            }
        }
        return new ComplexMatrixDense(X);
    }

    private ComplexMatrixDense U;
    private ComplexMatrixDense R;

    public ComplexQRDecompositionDense(ComplexMatrixDense matrix) {
        super(matrix);

        this.m = matrix.getColumnCount();
        this.n = matrix.getRowCount();
        ComplexMatrixDense R = matrix.copy();
        ComplexMatrixDense U = new ComplexMatrixDense(zeros(m, n));
        for (int i = 0; i < Math.min(m, n); i++) {
            Complex[] rCol = Arrays.copyOfRange(getColumn(R, i), i, m); // TODO get subColumn
            Complex[] u = houseGen(rCol);
            setColumn(U, u, i, i); // TODO set subColumn
            Complex[][] subR = getSubMatrix(R, i, m, i, n);
            Complex[][] H = calculateReflector(u, subR);
            setSubMatrix(R, H, i, m, i, n);
            for (int j = i + 1; j < m; j++) {
                R.unsafeSet(j, i, new Complex());
            }
        }
        this.R = R;
        this.U = U;
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
        return U;
    }

    /**
     * Return the upper triangular factor
     *
     * @return R
     */

    public ComplexMatrixDense getR() {
        return R;
    }

    /**
     * Generate and return the unitary orthogonal factor
     *
     * @return Q
     */
    public ComplexMatrixDense getQ() {
        ComplexMatrixDense I = ComplexMatrixDense.Factory.identity(U.getRowCount(), U.getColumnCount());
        return houseApply(U, I);
    }

    public ComplexMatrixDense QmultiplyX(ComplexMatrixDense X) {
        return houseApply(U, X);
    }

    /**
     * Generate and return the conjugate transpose of the orthogonal factor
     *
     * @return transpose(Q)
     */
    public ComplexMatrixDense getQH() {
        ComplexMatrixDense I = ComplexMatrixDense.Factory.identity(U.getRowCount(), U.getColumnCount());
        return houseApplyTranspose(U, I);
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
        ComplexMatrixDense Y = houseApplyTranspose(U, B);

        // Solve R * X = Y;
        // Back Substitution
        return backSubstitutionSolve(R, Y);
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < rows * cols; ++i) {
            if (i > 0 && i % cols == 0) {
                sb.append(System.lineSeparator());
            }
            sb.append(String.format("%.4f", _data[i])).append(" ");
        }
        return sb.toString();
    }
}
