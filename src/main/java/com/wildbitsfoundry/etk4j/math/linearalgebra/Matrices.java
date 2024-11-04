package com.wildbitsfoundry.etk4j.math.linearalgebra;

public final class Matrices {
    private Matrices() {
    }

    /**
     * Solve an {@code LDU Tri-diagonal system. A Tri-diagonal system is a matrix that only has entries
     * in its sub-diagonal, its diagonal, and its super-diagonal e.g:
     * <pre>
     *     1 1 0
     *     1 1 1
     *     0 1 1
     * </pre>
     * @param lower The sub-diagonal.
     * @param diagonal The diagonal.
     * @param upper The super-diagonal.
     * @param b The solution vector.
     * @return The solution to {@code LDUX = b}
     */
    public static void solveLDUTridiagonalSystem(double[] lower, double[] diagonal, double[] upper, double[] b) {
        final int m = b.length;
        b[0] = b[0] / diagonal[0];
        // Forward Substitution
        for (int i = 0; i < m - 1; ++i) {
            upper[i] = upper[i] / diagonal[i];
            diagonal[i + 1] = diagonal[i + 1] - lower[i] * upper[i];
            b[i + 1] = (b[i + 1] - lower[i] * b[i]) / diagonal[i + 1];
        }
        // Backwards Substitution
        for (int i = m - 2; i >= 0; --i) {
            b[i] = b[i] - upper[i] * b[i + 1];
        }
    }

    /**
     * Solve an {@code LDL<sup>T</sup>} Tri-diagonal system. A Tri-diagonal system is a matrix that only has entries
     * in its sub-diagonal, its diagonal, and its super-diagonal e.g:
     * <pre>
     *     1 1 0
     *     1 1 1
     *     0 1 1
     * </pre>
     * @param lower The sub-diagonal.
     * @param diagonal The diagonal.
     * @param b The solution vector.
     * @return The solution to {@code LDL<sup>T</sup>X = b}
     */
    public static double[] solveLDLtTridiagonalSystem(double[] lower, double[] diagonal, double[] b) {
        final int length = diagonal.length;
        for (int i = 0; i < length - 1; ++i) {
            double ui = lower[i];
            lower[i] /= diagonal[i];
            diagonal[i + 1] -= ui * lower[i];
            b[i + 1] -= lower[i] * b[i];
        }
        b[length - 1] /= diagonal[length - 1];
        for (int i = length - 2; i >= 0; --i) {
            b[i] = b[i] / diagonal[i] - lower[i] * b[i + 1];
        }
        return b;
    }

    /**
     * Forward Substitution Solve.
     * @param L The Lower triangular Matrix.
     * @param b The RHS (Solution) Matrix
     * @return The solution of {@code LX = b}.
     */
    public static MatrixDense forwardSubstitutionSolve(MatrixDense L, MatrixDense b) {
        final int nx = b.getColumnCount();
        final int m = L.getRowCount();
        final int n = L.getColumnCount();
        double[] t = L.getArray();
        double[] X = b.getArrayCopy();

        for (int j = 0; j < nx; ++j) {
            for (int i = 0; i < m; ++i) {
                for (int k = 0; k < i; ++k) {
                    X[i * nx + j] -= X[k * nx + j] * t[i * n + k];
                }
                X[i * nx + j] /= t[i * n + i];
            }
        }
        return new MatrixDense(X, m, nx);
    }
}
