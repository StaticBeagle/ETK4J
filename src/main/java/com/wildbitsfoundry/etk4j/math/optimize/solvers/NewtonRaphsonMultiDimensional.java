package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

public class NewtonRaphsonMultiDimensional {
    @FunctionalInterface
    public interface Function {
        double[] evaluate(double[] x); // Represents the system of equations F(x)
    }

    @FunctionalInterface
    interface Jacobian {
        double[][] evaluate(double[] x); // Represents the Jacobian matrix J(x)
    }

    static double[] solve(Function func, Jacobian jacobian, double[] x0, double tol, int maxIter) {
        double[] x = Arrays.copyOf(x0, x0.length); // Current estimate
        double[] f = func.evaluate(x); // Evaluate F(x)
        int n = x.length;
        for (int iter = 0; iter < maxIter; iter++) {
            double normF = norm(f);
            System.out.println("Iteration " + iter + ": ||F(x)|| = " + normF);
            // Check convergence
            if (normF < tol) {
                return x;
            }
            // Compute the Jacobian matrix at current x
            double[][] J = jacobian != null ? jacobian.evaluate(x) : JacobianFiniteDifferences.computeJacobian(func, x, 1e-6);
            // Solve for the Newton step: J(x) * dx = -F(x)
            double[] dx = solveLinearSystem(J, scaleVector(f, -1));
            // Perform a line search to find an optimal step size double
            // alpha0 = initial step size
            // c1 Armijo conditions
            // c2 Step size reduction factor
            double alpha = lineSearch(func, x, dx, f, 1.0, 1e-4, 0.9, 100);
            // Update x
            // x = x + alpha * dx;
            for (int i = 0; i < n; i++) {
                x[i] += alpha * dx[i];
            }
            // Recompute F(x)
            f = func.evaluate(x);
        }
        throw new RuntimeException("Newton's method did not converge after " + maxIter + " iterations.");
    }

    static double lineSearch(Function func, double[] x, double[] dx, double[] f, double alpha0, double c1, double c2, int maxIter) {
        double alpha = alpha0;
        double f0Norm = norm(f); // Initial norm of F(x) double
        double gradFdx = dotProduct(f, dx); // Directional derivative at current point
        for (int i = 0; i < maxIter; i++) { // Compute new x along the search direction
            double[] xNew = addVectors(x, scaleVector(dx, alpha));
            double[] fNew = func.evaluate(xNew);
            double fNewNorm = norm(fNew); // Norm of F(x + alpha * dx) // Check Armijo condition (sufficient decrease)
            if (fNewNorm <= f0Norm + c1 * alpha * gradFdx) { // Check Wolfe curvature condition
                // Gradients for nonlinear systems often use residual as proxy
                if (Math.abs(dotProduct(fNew, dx)) <= c2 * Math.abs(gradFdx)) {
                    return alpha; // Wolfe conditions satisfied
                }
            } // Reduce step size
            alpha *= 0.5;
        }
        return alpha; // Return last alpha if conditions not satisfied
    }

    // Compute the Euclidean norm of a vector
    static double norm(double[] v) {
        return Math.sqrt(Arrays.stream(v).map(vi -> vi * vi).sum());
    } // Compute the dot product of two vectors

    static double dotProduct(double[] a, double[] b) {
        double result = 0.0;
        for (int i = 0; i < a.length; i++) {
            result += a[i] * b[i];
        }
        return result;
    } // Add two vectors

    static double[] addVectors(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] + b[i];
        }
        return result;
    } // Scale a vector by a scalar

    static double[] scaleVector(double[] v, double scalar) {
        return Arrays.stream(v).map(vi -> vi * scalar).toArray();
    } // Solve a linear system Ax = b using Gaussian elimination

    static double[] solveLinearSystem(double[][] A, double[] b) {
        int n = b.length;
        double[] x = Arrays.copyOf(b, n); // Gaussian elimination
        for (int i = 0; i < n; i++) {
            for (int k = i + 1; k < n; k++) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < n; j++) {
                    A[k][j] -= factor * A[i][j];
                }
                x[k] -= factor * x[i];
            }
        } // Back substitution
        for (int i = n - 1; i >= 0; i--) {
            x[i] /= A[i][i];
            for (int j = 0; j < i; j++) {
                x[j] -= A[j][i] * x[i];
            }
        }
        return x;
    }

    public static void main(String[] args) {
        // Define the system of equations
        Function func = (x) -> new double[]{x[0] * x[0] + x[1] * x[1] - 1, x[0] * x[0] - x[1]}; // f1(x, y) = x^2 + y^2 - 1 x[0] * x[0] - x[1] // f2(x, y) = x^2 - y };
        // Define the Jacobian matrix
        Jacobian jacobian = (x) -> new double[][]{{2 * x[0], 2 * x[1]},
                // df1/dx, df1/dy
                {2 * x[0], -1}
                // df2/dx, df2/dy
        };
        func = (double[] x) -> new double[] {x[0] + x[1] - 3 * x[2] + x[3] - 2, -5 * x[0] + 3 * x[1] - 4 * x[2] + x[3], x[0] + 2 * x[2] - x[3] - 1, x[0] + 2 * x[1] - 12};
        //Initial guess
        double[] x0 = {1, 5, 5, 10};
        // Solve using the Newton-Raphson method
        double tol = 1e-6;
        int maxIter = 100;
        double[] solution = solve(func, null, x0, tol, maxIter);
        // [0.7861513777574233, 0.6180339887498949]
        // [0.7861513778721655, 0.6180339887498943]
        System.out.println("Solution: " + Arrays.toString(solution));
    }

    public static class JacobianFiniteDifferences { // Function interface for multi-dimensional functions


        /**
         * Compute the Jacobian matrix using finite differences. * * @param func The multi-dimensional function to evaluate. * @param x The point at which the Jacobian is computed. * @param epsilon Small perturbation for finite differences (default 1e-8). * @return The Jacobian matrix as a 2D array.
         */
        public static double[][] computeJacobian(Function func, double[] x, double epsilon) {
            int n = x.length; // Number of variables
            double[] f0 = func.evaluate(x); // Evaluate the function at x
            int m = f0.length; // Number of equations
            double[][] jacobian = new double[m][n];
            for (int j = 0; j < n; j++) {
                double[] xPerturbed = Arrays.copyOf(x, n);
                xPerturbed[j] += epsilon; // Perturb the j-th variable
                double[] fPerturbed = func.evaluate(xPerturbed); // Evaluate f(x + epsilon)
                for (int i = 0; i < m; i++) {
                    jacobian[i][j] = (fPerturbed[i] - f0[i]) / epsilon; // Compute partial derivative
                }
            }
            return jacobian;
        }
    }
}