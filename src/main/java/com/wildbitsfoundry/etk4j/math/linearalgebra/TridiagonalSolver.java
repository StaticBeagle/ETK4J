package com.wildbitsfoundry.etk4j.math.linearalgebra;

public final class TridiagonalSolver {

    private TridiagonalSolver() {}

    /**
     * Solve an {@code} LDU Tri-diagonal system. A Tri-diagonal system is a matrix that only has entries
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
     */
    public static void solveLDLtTridiagonalSystem(double[] lower, double[] diagonal, double[] b) {
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
    }
}
