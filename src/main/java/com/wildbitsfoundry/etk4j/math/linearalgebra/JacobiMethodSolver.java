package com.wildbitsfoundry.etk4j.math.linearalgebra;

import com.wildbitsfoundry.etk4j.util.DoubleArrays;
// TODO javadocs
public class JacobiMethodSolver {

    private double[] b;
    private MatrixDense A;
    private int iterationLimit = 100;
    private double tol = 1e-9;
    private double[] x0 = null;

    public JacobiMethodSolver(MatrixDense A, double[] b) {
        this.A = A;
        this.b = b;
    }

    public JacobiMethodSolver iterationLimit(int iterationLimit) {
        this.iterationLimit = iterationLimit;
        return this;
    }

    public JacobiMethodSolver tolerance(double tolerance) {
        this.tol = tolerance;
        return this;
    }

    public JacobiMethodSolver initialGuess(double[] x0) {
        this.x0 = x0;
        return this;
    }

    public IterativeSolverResults<double[]> solve() {
        double[] x = x0 == null ? new double[b.length] : x0;
        int n = A.getRowCount();
        double[] xNew = new double[n];
        int k;
        double error = Double.NaN;
        for (k = 0; k < iterationLimit; k++) {
            for (int i = 0; i < n; i++) {
                double sum = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        sum += A.unsafeGet(i, j) * x[j];
                    }
                }
                xNew[i] = (b[i] - sum) / A.unsafeGet(i, i);
            }
            // Check for convergence
            error = DoubleArrays.norm2(DoubleArrays.subtractElementWise(xNew, x));
            if(error < tol) {
                IterativeSolverResults<double[]> solverResults = new IterativeSolverResults<>();
                solverResults.setSolverStatus("Converged");
                solverResults.setHasConverged(true);
                solverResults.setError(error);
                solverResults.setValue(xNew);
                solverResults.setNumberOfIterations(k + 1);
                return solverResults;
            }
            System.arraycopy(xNew, 0, x, 0, xNew.length);
        }
        IterativeSolverResults<double[]> solverResults = new IterativeSolverResults<>();
        solverResults.setSolverStatus("Maximum number of iterations exceeded");
        solverResults.setHasConverged(false);
        solverResults.setError(error);
        solverResults.setValue(x);
        solverResults.setNumberOfIterations(k + 1);
        return solverResults;

    }
    public static void main(String[] args) {
        double[][] A = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
        double[] b = {1, 0, 1};
        double[] x = {0, 0, 0}; // Initial guess
        int maxIterations = 100;
        double tolerance = 0.0001;
        double[] solution = new JacobiMethodSolver(new MatrixDense(A), b)
                .initialGuess(x).iterationLimit(maxIterations).tolerance(tolerance).solve().getValue();
        System.out.println("Solution:");
        for (double val : solution) {
            System.out.println(val);
        }
    }
}
