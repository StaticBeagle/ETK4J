package com.wildbitsfoundry.etk4j.math.linearalgebra;

public final class Matrices {
    private Matrices() {
    }

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
     * @param U The upper triangular Matrix.
     * @param B The RHS (Solution) Matrix
     * @return The solution of {@code UX = B}.
     */
    public static Matrix fwdSubsSolve(Matrix U, Matrix B) {
        final int nx = B.getColumnCount();
        final int m = U.getRowCount();
        final int n = U.getColumnCount();
        double[] t = U.getArray();
        double[] X = B.getArrayCopy();

        for (int j = 0; j < nx; ++j) {
            for (int i = 0; i < m; ++i) {
                for (int k = 0; k < i; ++k) {
                    X[i * nx + j] -= X[k * nx + j] * t[i * n + k];
                }
                X[i * nx + j] /= t[i * n + i];
            }
        }
        return new Matrix(X, m, nx);
    }
}
