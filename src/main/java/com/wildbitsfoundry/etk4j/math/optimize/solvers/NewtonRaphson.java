package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

import static com.wildbitsfoundry.etk4j.math.optimize.solvers.NewtonSolverResults.NewtonSolverStatus;


public final class NewtonRaphson {

    // TODO remove solver results?
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

    public NewtonSolverResults<Double> solve() {
        int maxNumberOfIterations = this.maxNumberOfIterations;
        UnivariateFunction func = this.func;

        double xCurrent = x0;
        double xFinal = 0.0;
        if (derivative != null) {
            while (maxNumberOfIterations-- > 0) {
                double funcValue = func.evaluateAt(xCurrent);
                if (funcValue == 0) {
                    double error = Math.abs(xFinal - xCurrent);
                    NewtonSolverResults<Double> solverResults = new NewtonSolverResults<>();
                    solverResults.setSolverStatus(NewtonSolverStatus.CONVERGED);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(maxNumberOfIterations);
                    return solverResults;
                }
                double funcDerivativeValue = derivative.evaluateAt(xCurrent);
                if (funcDerivativeValue == 0) {
                    double error = Math.abs(xFinal - xCurrent);
                    NewtonSolverResults<Double> solverResults = new NewtonSolverResults<>();
                    solverResults.setSolverStatus(NewtonSolverStatus.DERIVATIVE_WAS_ZERO);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(maxNumberOfIterations);
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
                    NewtonSolverResults<Double> solverResults = new NewtonSolverResults<>();
                    // TODO Make a boolean field converged and return the solver status as a string
                    solverResults.setSolverStatus(NewtonSolverStatus.CONVERGED);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(maxNumberOfIterations);
                    return solverResults;
                }
                xCurrent = xFinal;
            }
        } else {
            // Secant method
            if (x1 != null) {
                if (x1 == x0) {
                    double error = Math.abs(xFinal - xCurrent);
                    NewtonSolverResults<Double> solverResults = new NewtonSolverResults<>();
                    solverResults.setSolverStatus(NewtonSolverStatus
                            .SECOND_INITIAL_GUESS_EQUAL_TO_FIRST_INITIAL_GUESS);
                    solverResults.setError(error);
                    solverResults.setValue(xFinal);
                    solverResults.setNumberOfIterations(maxNumberOfIterations);
                    return solverResults;
                }
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
            while (maxNumberOfIterations-- > 0) {
                double x = 0;
                if(q1 == q0) {
                    if(xFinal != xCurrent) {
                        double error = Math.abs(xFinal - xCurrent);
                        NewtonSolverResults<Double> newtonSolverResults = new NewtonSolverResults<>();
                        newtonSolverResults.setSolverStatus(NewtonSolverStatus.TOLERANCE_WAS_REACHED);
                        newtonSolverResults.setError(error);
                        newtonSolverResults.setValue(xFinal);
                        newtonSolverResults.setNumberOfIterations(maxNumberOfIterations);
                        return newtonSolverResults;
                    }
                    x = (xFinal + xCurrent) / 2.0;
                    double error = Math.abs(xFinal - xCurrent);
                    NewtonSolverResults<Double> solverResults = new NewtonSolverResults<>();
                    solverResults.setSolverStatus(NewtonSolverStatus.CONVERGED);
                    solverResults.setError(error);
                    solverResults.setValue(x);
                    solverResults.setNumberOfIterations(maxNumberOfIterations);
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
                    NewtonSolverResults<Double> solverResults = new NewtonSolverResults<>();
                    solverResults.setSolverStatus(NewtonSolverStatus.CONVERGED);
                    solverResults.setError(error);
                    solverResults.setValue(x);
                    solverResults.setNumberOfIterations(maxNumberOfIterations);
                    return solverResults;
                }
                xCurrent = xFinal;
                q0 = q1;
                xFinal = x;
                q1 = func.evaluateAt(xFinal);
            }
        }
        double error = Math.abs(xFinal - xCurrent);
        NewtonSolverResults<Double> newtonSolverResults = new NewtonSolverResults<>();
        newtonSolverResults.setSolverStatus(NewtonSolverStatus.MAX_NUMBER_OF_ITERATIONS_EXCEEDED);
        newtonSolverResults.setError(error);
        newtonSolverResults.setValue(xFinal);
        newtonSolverResults.setNumberOfIterations(maxNumberOfIterations);
        return newtonSolverResults;
    }

    public static void main(String[] args) {
        NewtonSolverResults nr = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20)
            .derivative(x -> 3 * x * x - 2 * x)
            .solve();

        NewtonSolverResults<Double> nr2 = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20)
            .derivative(x -> 3 * x * x - 2 * x)
            .secondDerivative(x -> 6 * x)
            .solve();

        NewtonSolverResults<Double> nr3 = new NewtonRaphson(x ->  x * x * x - x * x + 2, -20)
                .solve();
        System.out.println("s");
    }
}