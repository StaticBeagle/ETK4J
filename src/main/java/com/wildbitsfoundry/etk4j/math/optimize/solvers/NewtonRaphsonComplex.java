package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.ComplexUnivariateFunction;

/*
Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
All rights reserved. see https://github.com/StaticBeagle/ETK4J/blob/master/COPYING.
 */
public final class NewtonRaphsonComplex {
    private int maxNumberOfIterations = 100;
    private double absTol = 1e-9;
    private double relTol = 0;
    private Complex x0 = null;
    private Complex x1 = null;

    private ComplexUnivariateFunction func;
    private ComplexUnivariateFunction derivative;
    private ComplexUnivariateFunction secondDerivative;


    public NewtonRaphsonComplex(ComplexUnivariateFunction func, Complex initialGuess) {
        this.func = func;
        x0 = initialGuess;
    }

    public NewtonRaphsonComplex secondInitialGuess(Complex secondInitialGuess) {
        this.x1 = secondInitialGuess;
        return this;
    }

    public NewtonRaphsonComplex derivative(ComplexUnivariateFunction derivative) {
        this.derivative = derivative;
        return this;
    }

    public NewtonRaphsonComplex secondDerivative(ComplexUnivariateFunction secondDerivative) {
        this.secondDerivative = secondDerivative;
        return this;
    }

    public NewtonRaphsonComplex iterationLimit(int limit) {
        maxNumberOfIterations = limit;
        return this;
    }

    public NewtonRaphsonComplex absTolerance(double tol) {
        absTol = tol;
        return this;
    }

    public NewtonRaphsonComplex relTolerance(double tol) {
        relTol = tol;
        return this;
    }

    public SolverResults<Complex> solve() {
        int currentIteration = 0;
        ComplexUnivariateFunction func = this.func;

        Complex xCurrent = x0;
        Complex xFinal = new Complex();
        if (derivative != null) {
            while (currentIteration++ < maxNumberOfIterations) {
                Complex funcValue = func.evaluateAt(xCurrent);
                if (funcValue == new Complex()) {
                    double error = xFinal.subtract(xCurrent).abs();
                    SolverResults<Complex> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
                    solverResults.setHasConverged(true);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                }
                Complex funcDerivativeValue = derivative.evaluateAt(xCurrent);
                if (funcDerivativeValue == new Complex()) {
                    double error = xFinal.subtract(xCurrent).abs();
                    SolverResults<Complex> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Derivative was zero");
                    solverResults.setHasConverged(false);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                }
                // Newton Step
                Complex newtonStep = funcValue.divide(funcDerivativeValue);
                if (secondDerivative != null) {
                    // Use Halley's method
                    Complex funcSecondDerivativeValue = secondDerivative.evaluateAt(xCurrent);
                    Complex adj = newtonStep.multiply(funcSecondDerivativeValue).divide(funcDerivativeValue).divide(2);
                    if (adj.abs() < 1) {
                        newtonStep.divideEquals(adj.uminus().add(1.0));
                    }
                }
                xFinal = xCurrent.subtract(newtonStep);
                if (MathETK.isClose(xFinal, xCurrent, absTol, relTol)) {
                    double error = xFinal.subtract(xCurrent).abs();
                    SolverResults<Complex> solverResults = new SolverResults<>();
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
                if (x1.equals(x0)) {
                    double error = xFinal.subtract(xCurrent).abs();
                    SolverResults<Complex> solverResults = new SolverResults<>();
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
                xFinal = xCurrent.multiply(1 + eps);
                xFinal.addEquals(xFinal.compareTo(new Complex()) >= 0 ? eps : -eps);
            }
            Complex q0 = func.evaluateAt(xCurrent);
            Complex q1 = func.evaluateAt(xFinal);
            if (q1.abs() < q0.abs()) {
                Complex swap = xCurrent;
                xCurrent = xFinal;
                xFinal = swap;
                swap = q0;
                q0 = q1;
                q1 = swap;
            }
            while (currentIteration++ < maxNumberOfIterations) {
                Complex x = new Complex();
                if (q1.equals(q0)) {
                    if (!xFinal.equals(xCurrent)) {
                        double error = xFinal.subtract(xCurrent).abs();
                        SolverResults<Complex> solverResults = new SolverResults<>();
                        solverResults.setSolverStatus("Tolerance was reached");
                        solverResults.setHasConverged(false);
                        solverResults.setError(error);
                        solverResults.setValue(xFinal);
                        solverResults.setNumberOfIterations(currentIteration);
                        return solverResults;
                    }
                    x = xFinal.add(xCurrent).divide(2.0);
                    double error = xFinal.subtract(xCurrent).abs();
                    SolverResults<Complex> solverResults = new SolverResults<>();
                    solverResults.setSolverStatus("Converged");
                    solverResults.setHasConverged(true);
                    solverResults.setError(error);
                    solverResults.setValue(x);
                    solverResults.setNumberOfIterations(currentIteration);
                    return solverResults;
                } else {
                    if (q1.abs() > q0.abs()) {
                        x = q0.uminus().divide(q1).multiply(xFinal).add(xCurrent).divide(q0.divide(q1).uminus().add(1.0));
                    } else {
                        x = q1.uminus().divide(q0).multiply(xCurrent).add(xFinal).divide(q1.divide(q0).uminus().add(1.0));
                    }
                }
                if (MathETK.isClose(x, xFinal, absTol, relTol)) {
                    double error = xFinal.subtract(xCurrent).abs();
                    SolverResults<Complex> solverResults = new SolverResults<>();
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
        double error = xFinal.subtract(xCurrent).abs();
        SolverResults<Complex> solverResults = new SolverResults<>();
        solverResults.setSolverStatus("Maximum number of iterations exceeded");
        solverResults.setHasConverged(false);
        solverResults.setError(error);
        solverResults.setValue(xFinal);
        solverResults.setNumberOfIterations(currentIteration);
        return solverResults;
    }
}