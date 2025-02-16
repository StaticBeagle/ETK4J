package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.math.complex.Complex;

/**
 * This class implements the Gaussian Elimination method to solve a linear system of equations
 */
public final class ComplexGaussianEliminationSolver {

    private ComplexGaussianEliminationSolver() {}

    /**
     * Solves a linear system of equations {@code A*x = b} using Gauss method with a tolerance of 1e-10.
     * @param A The coefficient matrix
     * @param b The solution vector
     * @return x = A<sup>-1</sup>*b
     */
    public static Complex[] solve(ComplexMatrixDense A, Complex[] b) {
        return solve(A, b, 1e-10);
    }

    /**
     * Solves a linear system of equations {@code A*x = b} using Gauss method to a specified tolerance.
     * @param A The coefficient matrix
     * @param b The solution vector
     * @param tol The tolerance used to stop the algorithm
     * @return x = A<sup>-1</sup>*b
     */
    public static Complex[] solve(ComplexMatrixDense A, Complex[] b, double tol) {
        int n = b.length;

        for (int p = 0; p < n; p++) {
            // Find pivot row and swap
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (A.unsafeGet(i, p).abs() > A.unsafeGet(max, p).abs()) {
                    max = i;
                }
            }
            Complex[] temp = A.getRow(p);
            A.setRow(p, A.getRow(max));
            A.setRow(max, temp);
            Complex t = b[p]; b[p] = b[max]; b[max] = t;

            // Check for singular or nearly singular matrix
            if (A.unsafeGet(p, p).abs() <= tol) {
                throw new ArithmeticException("Matrix is singular or nearly singular");
            }

            // Pivot within A and b
            for (int i = p + 1; i < n; i++) {
                Complex alpha = A.unsafeGet(i, p).divide(A.unsafeGet(p, p));
                b[i].subtractEquals(alpha.multiply(b[p]));
                for (int j = p; j < n; j++) {
                    A.unsafeSet(i, j, A.unsafeGet(i, j).subtract(alpha.multiply(A.unsafeGet(p, j))));
                }
            }
        }

        // Back substitution
        Complex[] x = new Complex[n];
        for (int i = n - 1; i >= 0; i--) {
            Complex sum = new Complex();
            for (int j = i + 1; j < n; j++) {
                sum.addEquals(A.unsafeGet(i, j).multiply(x[j]));
            }
            x[i] = b[i].subtract(sum).divide(A.unsafeGet(i, i));
        }
        return x;
    }
}
