package com.wildbitsfoundry.etk4j.math.linearalgebra;

/**
 * This class implements the Gaussian Elimination method to solve a linear system of equations
 */
public final class GaussianEliminationSolver {

    private GaussianEliminationSolver() {}

    /**
     * Solves a linear system of equations {@code A*x = b} using Gauss method with a tolerance of 1e-10.
     * @param A The coefficient matrix
     * @param b The solution vector
     * @return x = A<sup>-1</sup>*b
     */
    public static double[] solve(MatrixDense A, double[] b) {
        return solve(A, b, 1e-10);
    }

    /**
     * Solves a linear system of equations {@code A*x = b} using Gauss method to a specified tolerance.
     * @param A The coefficient matrix
     * @param b The solution vector
     * @param tol The tolerance used to stop the algorithm
     * @return x = A<sup>-1</sup>*b
     */
    public static double[] solve(MatrixDense A, double[] b, double tol) {
        int n = b.length;

        for (int p = 0; p < n; p++) {
            // Find pivot row and swap
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (Math.abs(A.unsafeGet(i, p)) > Math.abs(A.unsafeGet(max, p))) {
                    max = i;
                }
            }
            double[] temp = A.getRow(p);
            A.setRow(p, A.getRow(max));
            A.setRow(max, temp);
            double t = b[p]; b[p] = b[max]; b[max] = t;

            // Check for singular or nearly singular matrix
            if (Math.abs(A.unsafeGet(p, p)) <= tol) {
                throw new ArithmeticException("Matrix is singular or nearly singular");
            }

            // Pivot within A and b
            for (int i = p + 1; i < n; i++) {
                double alpha = A.unsafeGet(i, p) / A.unsafeGet(p, p);
                b[i] -= alpha * b[p];
                for (int j = p; j < n; j++) {
                    A.unsafeSet(i, j, A.unsafeGet(i, j) - alpha * A.unsafeGet(p, j));
                }
            }
        }

        // Back substitution
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A.unsafeGet(i, j) * x[j];
            }
            x[i] = (b[i] - sum) / A.unsafeGet(i, i);
        }
        return x;
    }
}
