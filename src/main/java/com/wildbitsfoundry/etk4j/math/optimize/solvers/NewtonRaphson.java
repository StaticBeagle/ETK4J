package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

/*
Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
All rights reserved. see https://github.com/StaticBeagle/ETK4J/blob/master/COPYING.
 */
public final class NewtonRaphson {

    private int maxNumberOfIterations = 100;
    private double absTol = 1e-9;
    private double relTol = 0;
    private double x0 = 0;
    private Double x1 = null;

    private UnivariateFunction func;
    private UnivariateFunction derivative;
    private UnivariateFunction secondDerivative;


    public NewtonRaphson(UnivariateFunction func, double initialGuess) {
        this.func = func;
        x0 = initialGuess;
    }

    public NewtonRaphson secondInitialGuess(double secondInitialGuess) {
        this.x1 = secondInitialGuess;
        return this;
    }

    public NewtonRaphson derivative(UnivariateFunction derivative) {
        this.derivative = derivative;
        return this;
    }

    public NewtonRaphson secondDerivative(UnivariateFunction secondDerivative) {
        this.secondDerivative = secondDerivative;
        return this;
    }

    public NewtonRaphson iterationLimit(int limit) {
        maxNumberOfIterations = limit;
        return this;
    }

    public NewtonRaphson absTolerance(double tol) {
        absTol = tol;
        return this;
    }

    public NewtonRaphson relTolerance(double tol) {
        relTol = tol;
        return this;
    }

    public SolverResults<Double> solve() {
        int currentIteration = 0;
        UnivariateFunction func = this.func;

        double xCurrent = x0;
        double xFinal = 0.0;
        if (derivative != null) {
            while (currentIteration++ < maxNumberOfIterations) {
                double funcValue = func.evaluateAt(xCurrent);
                if (funcValue == 0) {
                    double error = Math.abs(xFinal - xCurrent);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
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
            if(Math.abs(q1) < Math.abs(q0)) {
                double swap = xCurrent;
                xCurrent = xFinal;
                xFinal = swap;
                swap = q0;
                q0 = q1;
                q1 = swap;
            }
            while (currentIteration++ < maxNumberOfIterations) {
                double x = 0;
                if(q1 == q0) {
                    if(xFinal != xCurrent) {
                        double error = Math.abs(xFinal - xCurrent);
                        SolverResults<Double> solverResults = new SolverResults<>();
                        solverResults.setSolverStatus("Tolerance was reached");
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
                    solverResults.setHasConverged(true);
                    solverResults.setError(error);
                    solverResults.setValue(x);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                } else {
                    if(Math.abs(q1) > Math.abs(q0)) {
                        x = (-q0 / q1 * xFinal + xCurrent) / (1 - q0 / q1);
                    } else {
                        x = (-q1 / q0 * xCurrent + xFinal) / (1 - q1 / q0);
                    }
                }
                if(MathETK.isClose(x, xFinal, absTol, relTol)) {
                    double error = Math.abs(xFinal - xCurrent);
                    SolverResults<Double> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
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
        solverResults.setHasConverged(false);
        solverResults.setError(error);
        solverResults.setValue(xFinal);
        solverResults.setNumberOfIterations(currentIteration);
        return solverResults;
    }
}