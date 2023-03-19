package com.wildbitsfoundry.etk4j.math.optimize.solvers;


import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;

/**
 * The {@code Ridder}'s class contains the Brent's root finding algorithm.
 * @see <a href="https://en.wikipedia.org/wiki/Ridders%27_method#:~:text=In%20numerical%20analysis%2C%20Ridders%27%20method%20is%20a%20root-finding,function.%20The%20method%20is%20due%20to%20C.%20Ridders.">Ridder's method</a>
 */
public class Ridder {

    private int maxNumberOfIterations = 100;
    private double absTol = 1e-9;
    private double relTol = 4.0 * ConstantsETK.DOUBLE_EPS;
    private double a;
    private double b;

    private UnivariateFunction function;


    /**
     * Constructs an instance of the {@code Ridder}'s root finding algorithm.
     * @param function The function whose root is to be found.
     * @param a The lower end of the bracketing interval [a, b].
     * @param b The upper end of the bracketing interval [a, b].
     */
    public Ridder(UnivariateFunction function, double a, double b) {
        this.function = function;
        this.a = a;
        this.b = b;
        if (b <= a) {
            throw new IllegalArgumentException("The upper bracket (b) must be greater than the lower bracket (a).");
        }
    }

    /**
     * Maximum number of iterations.
     * @param limit The maximum number of iterations allowed.
     */
    public Ridder iterationLimit(int limit) {
        maxNumberOfIterations = limit;
        return this;
    }

    /**
     * Absolute tolerance.
     * @param tol The maximum allowed absolute tolerance.
     */
    public Ridder absTolerance(double tol) {
        absTol = tol;
        return this;
    }

    /**
     * Relative tolerance.
     * @param tol The maximum allowed relative tolerance. This value must be bigger than 4 * {@link ConstantsETK#DOUBLE_EPS}.
     * @return
     */
    public Ridder relTolerance(double tol) {
        if (tol < 4.0 * ConstantsETK.DOUBLE_EPS) {
            throw new IllegalArgumentException("Relative tolerance cannot be smaller than 4 * ConstantsETK.DOUBLE_EPS");
        }
        relTol = tol;
        return this;
    }

    /**
     * Find the root.
     * @return The {@link SolverResults} containing the root and other solver results.
     */
    SolverResults<Double> solve() {
        double fa, fb, fc, fx, c, s, dx, x = 0.0, xold = 0.0, error = 0.0;
        int currentIteration = 0;
        fa = function.evaluateAt(a);
        fb = function.evaluateAt(b);

        if (fa * fb > 0.0) {
            error = Math.abs(Double.NaN);
            SolverResults<Double> solverResults = new SolverResults<>();
            solverResults.setSolverStatus("Root is not bracketed.");
            solverResults.setOptimizerStatusType(OptimizerStatusType.ROOT_IS_NOT_BRACKETED);
            solverResults.setHasConverged(false);
            solverResults.setError(error);
            solverResults.setValue(b);
            solverResults.setNumberOfIterations(currentIteration);
            return solverResults;
        }
        if (fa == 0) {
            error = 0.0;
            SolverResults<Double> solverResults = new SolverResults<>();
            solverResults.setSolverStatus("Converged");
            solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
            solverResults.setHasConverged(true);
            solverResults.setError(error);
            solverResults.setValue(xold);
            solverResults.setNumberOfIterations(currentIteration);
            return solverResults;
        }
        if (fb == 0) {
            error = 0.0;
            SolverResults<Double> solverResults = new SolverResults<>();
            solverResults.setSolverStatus("Converged");
            solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
            solverResults.setHasConverged(true);
            solverResults.setError(error);
            solverResults.setValue(x);
            solverResults.setNumberOfIterations(currentIteration);
            return solverResults;
        }

        for (currentIteration = 1; currentIteration <= maxNumberOfIterations; currentIteration++) {
            c = 0.5 * (a + b);
            fc = function.evaluateAt(c);
            s = Math.sqrt(fc * fc - fa * fb);
            if (s == 0.0) {
                error = Math.abs(x - xold);
                SolverResults<Double> solverResults = new SolverResults<>();
                solverResults.setSolverStatus("Converged");
                solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
                solverResults.setHasConverged(true);
                solverResults.setError(error);
                solverResults.setValue(x);
                solverResults.setNumberOfIterations(currentIteration);
                return solverResults;
            }
            dx = (c - a) * fc / s;
            if (fa - fb < 0.0) {
                dx = -dx;
            }
            x = c + dx;
            fx = function.evaluateAt(x);
            if (currentIteration > 1) { // Check if we are done
                if (MathETK.isClose(x, xold, absTol, relTol)) {
                    error = Math.abs(x - xold);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
                    solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
                    solverResults.setHasConverged(true);
                    solverResults.setError(error);
                    solverResults.setValue(x);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                }
            }
            error = Math.abs(x - xold);
            xold = x; // Update the root
            // Re-bracket and continue
            if (fc * fx > 0.0) {
                if (fa * fx < 0.0) {
                    b = x;
                    fb = fx;
                } else {
                    a = x;
                    fa = fx;
                }
            } else {
                a = c;
                b = x;
                fa = fc;
                fb = fx;
            }

        }
        SolverResults<Double> solverResults = new SolverResults<>();
        solverResults.setSolverStatus("Maximum number of iterations exceeded");
        solverResults.setOptimizerStatusType(OptimizerStatusType.MAXIMUM_NUMBER_OF_ITERATIONS_EXCEEDED);
        solverResults.setHasConverged(false);
        solverResults.setError(error);
        solverResults.setValue(x);
        solverResults.setNumberOfIterations(currentIteration);
        return solverResults;
    }
}
