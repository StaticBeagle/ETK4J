package com.wildbitsfoundry.etk4j.math.optimize.solver;

import com.wildbitsfoundry.etk4j.math.calculus.JacobianCalculation5PointStencilStrategy;
import com.wildbitsfoundry.etk4j.math.calculus.JacobianCalculationStrategy;
import com.wildbitsfoundry.etk4j.math.function.MultivariateFunction;
import com.wildbitsfoundry.etk4j.math.linearalgebra.GaussianEliminationSolver;
import com.wildbitsfoundry.etk4j.math.linearalgebra.LUDecompositionDense;
import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;
import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

public class NewtonRaphsonMultivariate {

    private final MultivariateFunction[] functions;
    private MultivariateFunction[][] jacobian;
    private final double[] x0;
    private int maxNumberOfIterations = 100;
    private double tol = 1e-6;
    private double alpha0 = 1.0;
    private double c1 = 1e-4;
    private double c2 = 0.9;
    private int lineSearchMaxNumberOfIterations = 100;
    private double h = 1e-6;
    private JacobianCalculationStrategy jacobianCalculationStrategy = new JacobianCalculation5PointStencilStrategy();

    public NewtonRaphsonMultivariate(MultivariateFunction[] functions, double[] x0) {
        this.functions = functions;
        this.x0 = x0;
    }

    public NewtonRaphsonMultivariate jacobian(MultivariateFunction[][] jacobian) {
        this.jacobian = jacobian;
        return this;
    }

    public NewtonRaphsonMultivariate iterationLimit(int iterationLimit) {
        this.maxNumberOfIterations = iterationLimit;
        return this;
    }

    public NewtonRaphsonMultivariate tolerance(double tol) {
        this.tol = tol;
        return this;
    }

    public NewtonRaphsonMultivariate lineSearchInitialStepSize(double alpha0) {
        this.alpha0 = alpha0;
        return this;
    }

    public NewtonRaphsonMultivariate lineSearchArmijoParameter(double c1) {
        this.c1 = c1;
        return this;
    }

    public NewtonRaphsonMultivariate lineSearchStepSizeReductionFactor(double c2) {
        this.c2 = c2;
        return this;
    }

    public NewtonRaphsonMultivariate lineSearchIterationLimit(int iterationLimit) {
        this.lineSearchMaxNumberOfIterations = iterationLimit;
        return this;
    }

    public NewtonRaphsonMultivariate differentiationStepSize(double h) {
        this.h = h;
        return this;
    }

    public NewtonRaphsonMultivariate setJacobianCalculationStrategy(JacobianCalculationStrategy jacobianCalculationStrategy) {
        this.jacobianCalculationStrategy = jacobianCalculationStrategy;
        return this;
    }

    public SolverResults<double[]> solve() {
        double[] x = Arrays.copyOf(x0, x0.length); // Current estimate
        double[] f = Arrays.stream(functions).mapToDouble(fn -> fn.evaluateAt(x)).toArray(); // Evaluate F(x)
        int n = x.length;
        double normF = Double.NaN;
        for (int iter = 0; iter < maxNumberOfIterations; iter++) {
            normF = DoubleArrays.norm(f);
            // Check convergence
            if (normF <= tol) {
                SolverResults<double[]> solverResults = new SolverResults<>();
                solverResults.setSolverStatus("Converged");
                solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
                solverResults.setHasConverged(true);
                solverResults.setError(normF);
                solverResults.setValue(x);
                solverResults.setNumberOfIterations(iter);
                return solverResults;
            }
            // Compute the Jacobian matrix at current x
            double[][] J = jacobian != null ? evaluateJacobian(jacobian, x) : jacobianCalculationStrategy.calculateJacobian(functions, x, h);
            // Solve for the Newton step: J(x) * dx = -F(x)
            double[] dx = GaussianEliminationSolver.solve(J, DoubleArrays.multiplyElementWise(f, -1));
            // Perform a line search to find an optimal step size double
            // alpha0 = initial step size
            // c1 Armijo conditions
            // c2 Step size reduction factor
            double alpha = lineSearch(functions, x, dx, f, alpha0, c1, c2, lineSearchMaxNumberOfIterations);
            // Update x
            // x = x + alpha * dx;
            for (int i = 0; i < n; i++) {
                x[i] += alpha * dx[i];
            }
            // Recompute F(x)
            f = Arrays.stream(functions).mapToDouble(fn -> fn.evaluateAt(x)).toArray();
        }
        SolverResults<double[]> solverResults = new SolverResults<>();
        solverResults.setSolverStatus("Maximum number of iterations exceeded");
        solverResults.setOptimizerStatusType(OptimizerStatusType.MAXIMUM_NUMBER_OF_ITERATIONS_EXCEEDED);
        solverResults.setHasConverged(true);
        solverResults.setError(normF);
        solverResults.setValue(x);
        solverResults.setNumberOfIterations(maxNumberOfIterations);
        return solverResults;
    }

    private static double lineSearch(MultivariateFunction[] functions, double[] x, double[] dx, double[] f, double alpha0, double c1, double c2, int maxIter) {
        double alpha = alpha0;
        double f0Norm = DoubleArrays.norm(f); // Initial norm of F(x) double
        double gradFdx = DoubleArrays.dot(f, dx); // Directional derivative at current point
        for (int i = 0; i < maxIter; i++) { // Compute new x along the search direction
            double[] xNew = DoubleArrays.addElementWise(x, DoubleArrays.multiplyElementWise(dx, alpha));
            double[] fNew = Arrays.stream(functions).mapToDouble(fn -> fn.evaluateAt(xNew)).toArray();
            double fNewNorm = DoubleArrays.norm(fNew); // Norm of F(x + alpha * dx) // Check Armijo condition (sufficient decrease)
            if (fNewNorm <= f0Norm + c1 * alpha * gradFdx) { // Check Wolfe curvature condition
                // Gradients for nonlinear systems often use residual as proxy
                if (Math.abs(DoubleArrays.dot(fNew, dx)) <= c2 * Math.abs(gradFdx)) {
                    return alpha; // Wolfe conditions satisfied
                }
            } // Reduce step size
            alpha *= 0.5;
        }
        return alpha; // Return last alpha if conditions not satisfied
    }

    private static double[][] evaluateJacobian(MultivariateFunction[][] jacobian, double[] x) {
        double[][] J = new double[jacobian.length][jacobian[0].length];
        for(int i = 0; i < jacobian.length; i++) {
            for(int j = 0; j < jacobian[0].length; j++) {
                J[i][j] = jacobian[i][j].evaluateAt(x);
            }
        }
        return J;
    }
}