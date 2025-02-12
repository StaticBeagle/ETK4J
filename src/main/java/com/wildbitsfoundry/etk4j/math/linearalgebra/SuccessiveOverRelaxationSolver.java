package com.wildbitsfoundry.etk4j.math.linearalgebra;

/**
 * This class implements the SOR (Successive-Over-Relaxation) iterative method to solve a system of equations
 */
public class SuccessiveOverRelaxationSolver {
    private final double w;
    private final double[] b;
    private final MatrixDense A;
    private int iterationLimit = 100;
    private double tol = 1e-9;
    private double[] x0 = null;

    public SuccessiveOverRelaxationSolver(MatrixDense A, double[] b, double w) {
        this.A = A;
        this.b = b;
        this.w = w;
    }

    /**
     * Sets the maximum number of iterations allowed
     * @param iterationLimit the maximum number of iterations allowed
     * @return {@code this}
     */
    public SuccessiveOverRelaxationSolver iterationLimit(int iterationLimit) {
        this.iterationLimit = iterationLimit;
        return this;
    }

    /**
     * Sets absolute tolerance that determines when to stop the algorithm
     * @param tolerance the absolute tolerance
     * @return {@code this}
     */
    public SuccessiveOverRelaxationSolver tolerance(double tolerance) {
        this.tol = tolerance;
        return this;
    }

    /**
     * Sets the initial guess to the solution of the system of equations
     * @param x0 the initial guess of the solutions. All zeros if not set.
     * @return {@code this}
     */
    public SuccessiveOverRelaxationSolver initialGuess(double[] x0) {
        this.x0 = x0;
        return this;
    }

    /**
     * Solves A * x = b
     * @return An {@link IterativeSolverResults} with the solutions of the system A * x  = b
     */
    public IterativeSolverResults<double[]> solve() {
        double[] x = x0 == null ? new double[b.length] : x0;
        int n = b.length;
        int it = 0;
        double dxmax = Double.NaN;
        while (it++ < iterationLimit) {
            dxmax = 0;
            for (int i = 0; i < n; i++) {
                double residual = b[i];
                for (int j = 0; j < n; j++) {
                    residual -= A.unsafeGet(i, j) * x[j];
                }
                if (Math.abs(residual) > dxmax) {
                    dxmax = Math.abs(residual);
                }
                x[i] += w * residual / A.unsafeGet(i, i);
            }
            if (dxmax < tol) {
                IterativeSolverResults<double[]> solverResults = new IterativeSolverResults<>();
                solverResults.setSolverStatus("Converged");
                solverResults.setHasConverged(true);
                solverResults.setError(dxmax);
                solverResults.setValue(x);
                solverResults.setNumberOfIterations(it);
                return solverResults;
            }
        }
        IterativeSolverResults<double[]> solverResults = new IterativeSolverResults<>();
        solverResults.setSolverStatus("Maximum number of iterations exceeded");
        solverResults.setHasConverged(false);
        solverResults.setError(dxmax);
        solverResults.setValue(x);
        solverResults.setNumberOfIterations(it);
        return solverResults;
    }
}
