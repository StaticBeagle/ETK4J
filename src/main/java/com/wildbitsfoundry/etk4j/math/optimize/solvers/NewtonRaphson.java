package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;

/*
Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
All rights reserved. see https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
 */

/**
 * The {@code NewtonRaphson}'s class contains the Newton-Raphson root finding algorithm.
 *
 * @see <a href="https://en.wikipedia.org/wiki/Newton%27s_method">Newton's method</a>
 */
public final class NewtonRaphson {

    private int maxNumberOfIterations = 100;
    private double absTol = 1e-9;
    private double relTol = 0;
    private double x0 = 0;
    private Double x1 = null;

    private UnivariateFunction function;
    private UnivariateFunction derivative;
    private UnivariateFunction secondDerivative;

    /**
     * Constructs an instance of the {@code NewtonRaphson}'s root finding algorithm.
     *
     * @param function     The function whose root is to be found.
     * @param initialGuess The initial guess of where the root might be.
     */
    public NewtonRaphson(UnivariateFunction function, double initialGuess) {
        this.function = function;
        x0 = initialGuess;
    }

    /**
     * Second initial guess.
     *
     * @param secondInitialGuess The second guess of there the where the root might be.
     * @return {@code this} instance
     */
    public NewtonRaphson secondInitialGuess(double secondInitialGuess) {
        this.x1 = secondInitialGuess;
        return this;
    }

    /**
     * Derivative of the function.
     *
     * @param derivative The first derivative of the function. If this argument is not provided, the Secant method
     *                   is used to find the root of the function.
     * @return {@code this} instance
     */
    public NewtonRaphson derivative(UnivariateFunction derivative) {
        this.derivative = derivative;
        return this;
    }

    /**
     * Second derivative of the function.
     *
     * @param secondDerivative The second derivative of the function. If this argument is provided, Halley's method
     *                         is used to find the root of the function.
     * @return {@code this} instance
     */
    public NewtonRaphson secondDerivative(UnivariateFunction secondDerivative) {
        this.secondDerivative = secondDerivative;
        return this;
    }

    /**
     * Maximum number of iterations.
     *
     * @param limit The maximum number of iterations allowed.
     * @return {@code this} instance
     */
    public NewtonRaphson iterationLimit(int limit) {
        maxNumberOfIterations = limit;
        return this;
    }

    /**
     * Absolute tolerance.
     *
     * @param tol The maximum allowed absolute tolerance.
     * @return {@code this} instance
     */
    public NewtonRaphson absTolerance(double tol) {
        absTol = tol;
        return this;
    }

    /**
     * Relative tolerance.
     *
     * @param tol The maximum allowed relative tolerance.
     * @return {@code this} instance
     */
    public NewtonRaphson relTolerance(double tol) {
        relTol = tol;
        return this;
    }

    /**
     * Find the root.
     *
     * @return The {@link SolverResults} containing the root and other solver results.
     */
    public SolverResults<Double> solve() {
        int currentIteration = 0;
        UnivariateFunction func = this.function;

        double xCurrent = x0;
        double xFinal = 0.0;
        if (derivative != null) {
            while (currentIteration++ < maxNumberOfIterations) {
                double funcValue = func.evaluateAt(xCurrent);
                if (funcValue == 0) {
                    double error = Math.abs(xFinal - xCurrent);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
                    solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
                    solverResults.setHasConverged(true);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                }
                double funcDerivativeValue = derivative.evaluateAt(xCurrent);
                if (funcDerivativeValue == 0) {
                    double error = Math.abs(xFinal - xCurrent);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Derivative was zero");
                    solverResults.setOptimizerStatusType(OptimizerStatusType.DERIVATIVE_WAS_ZERO);
                    solverResults.setHasConverged(false);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                }
                // Newton Step
                double newtonStep = funcValue / funcDerivativeValue;
                if (secondDerivative != null) {
                    // Use Halley's method
                    double funcSecondDerivativeValue = secondDerivative.evaluateAt(xCurrent);
                    double adj = newtonStep * funcSecondDerivativeValue / funcDerivativeValue / 2.0;
                    if (Math.abs(adj) < 1) {
                        newtonStep /= 1.0 - adj;
                    }
                }
                xFinal = xCurrent - newtonStep;
                if (MathETK.isClose(xFinal, xCurrent, absTol, relTol)) {
                    double error = Math.abs(xFinal - xCurrent);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
                    solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
                    solverResults.setHasConverged(true);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                }
                xCurrent = xFinal;
            }
        } else {
            // Secant method
            if (x1 != null) {
                if (x1 == x0) {
                    double error = Math.abs(xFinal - xCurrent);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Second initial guess was equal to first initial guess");
                    solverResults.setOptimizerStatusType(OptimizerStatusType.SECOND_INITIAL_GUESS_WAS_EQUAL_TO_FIRST_INITIAL_GUESS);
                    solverResults.setHasConverged(false);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                }
                xFinal = x1;
            } else {
                double eps = 1e-4;
                xFinal = xCurrent * (1 + eps);
                xFinal += xFinal >= 0 ? eps : -eps;
            }
            double q0 = func.evaluateAt(xCurrent);
            double q1 = func.evaluateAt(xFinal);
            if (Math.abs(q1) < Math.abs(q0)) {
                double swap = xCurrent;
                xCurrent = xFinal;
                xFinal = swap;
                swap = q0;
                q0 = q1;
                q1 = swap;
            }
            while (currentIteration++ < maxNumberOfIterations) {
                double x = 0;
                if (q1 == q0) {
                    if (xFinal != xCurrent) {
                        double error = Math.abs(xFinal - xCurrent);
                        SolverResults<Double> solverResults = new SolverResults<>();
                        solverResults.setSolverStatus("Tolerance was reached");
                        solverResults.setOptimizerStatusType(OptimizerStatusType.TOLERANCE_WAS_REACHED);
                        solverResults.setHasConverged(false);
                        solverResults.setError(error);
                        solverResults.setValue(xFinal);
                        solverResults.setNumberOfIterations(currentIteration);
                        return solverResults;
                    }
                    x = (xFinal + xCurrent) / 2.0;
                    double error = Math.abs(xFinal - xCurrent);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
                    solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
                    solverResults.setHasConverged(true);
                    solverResults.setError(error);
                    solverResults.setValue(x);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                } else {
                    if (Math.abs(q1) > Math.abs(q0)) {
                        x = (-q0 / q1 * xFinal + xCurrent) / (1 - q0 / q1);
                    } else {
                        x = (-q1 / q0 * xCurrent + xFinal) / (1 - q1 / q0);
                    }
                }
                if (MathETK.isClose(x, xFinal, absTol, relTol)) {
                    double error = Math.abs(xFinal - xCurrent);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
                    solverResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
                    solverResults.setHasConverged(true);
                    solverResults.setError(error);
                    solverResults.setValue(x);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                }
                xCurrent = xFinal;
                q0 = q1;
                xFinal = x;
                q1 = func.evaluateAt(xFinal);
            }
        }
        double error = Math.abs(xFinal - xCurrent);
        SolverResults<Double> solverResults = new SolverResults<>();
        solverResults.setSolverStatus("Maximum number of iterations exceeded");
        solverResults.setOptimizerStatusType(OptimizerStatusType.MAXIMUM_NUMBER_OF_ITERATIONS_EXCEEDED);
        solverResults.setHasConverged(false);
        solverResults.setError(error);
        solverResults.setValue(xFinal);
        solverResults.setNumberOfIterations(currentIteration);
        return solverResults;
    }
}