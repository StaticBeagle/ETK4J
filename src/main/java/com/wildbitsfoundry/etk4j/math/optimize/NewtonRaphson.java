package com.wildbitsfoundry.etk4j.math.optimize;

import com.wildbitsfoundry.etk4j.math.calculus.Derivatives;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

import static com.wildbitsfoundry.etk4j.math.optimize.SolverResults.SolverStatus;


public final class NewtonRaphson {

    private NewtonRaphson() {}

    protected int _maxIter = 100;
    protected double _absTol = 1e-9;
    protected double _relTol = 0;
    protected double _maxVal = Double.POSITIVE_INFINITY;
    private double _minVal = Double.NEGATIVE_INFINITY;
    protected double _step = 0.001;
    protected double _x0;

    protected UnivariateFunction _func;
    protected UnivariateFunction _derivative;


    public NewtonRaphson(UnivariateFunction func, double initialGuess) {
        _func = func;
        _x0 = initialGuess;
    }

    // create constructor that takes in derivative type

    public NewtonRaphson(UnivariateFunction func, UnivariateFunction derivative, double initialGuess) {
        _func = func;
        _derivative = derivative;
        _x0 = initialGuess;
    }

    public NewtonRaphson iterationLimit(int limit) {
        _maxIter = limit;
        return this;
    }

    public NewtonRaphson absTolerance(double tol) {
        _absTol = tol;
        return this;
    }

    public NewtonRaphson relTolerance(double tol) {
        _relTol = tol;
        return this;
    }

    public NewtonRaphson maxAllowedValue(double max) {
        _maxVal = max;
        return this;
    }

    public NewtonRaphson minAllowedValue(double min) {
        _minVal = min;
        return this;
    }

    public NewtonRaphson differentiationStepSize(double step) {
        _step = step;
        return this;
    }

    protected static SolverResults buildResults(double xfinal, SolverStatus status, int iterCount, double error) {
        SolverResults<Double> sr = new SolverResults<>();
        sr.setValue(xfinal);
        sr.setSolverStatus(status);
        sr.setNumberOfIterations(iterCount);
        sr.setError(error);
        return sr;
    }

    public SolverResults solve() {
        return this.solve(_derivative);
    }

    protected SolverResults solve(UnivariateFunction derivative) {
        int maxiter = _maxIter;
        double maxval = _maxVal;
        double step = _step;
        UnivariateFunction func = _func;

        double xcurrent = _x0;
        double xfinal = 0.0;
        double error = 0.0;
        while (maxiter-- > 0) {
            if (derivative != null) {
                xfinal = xcurrent - func.evaluateAt(xcurrent) / derivative.evaluateAt(xcurrent);
            } else {
                double fprime = Derivatives.centeredDifference7Points(func, xcurrent, step);
                xfinal = xcurrent - func.evaluateAt(xcurrent) / fprime;
            }

            error = Math.abs(xfinal - xcurrent);
            if (error < _absTol + _relTol * Math.min(Math.abs(xfinal), Math.abs(xcurrent))) {
                return buildResults(xfinal, SolverStatus.CONVERGED, _maxIter - maxiter, error);
            }

            if (Double.compare(xcurrent, maxval) > 0) {
                return buildResults(xfinal, SolverStatus.MAX_VALUE_EXCEEDED, _maxIter - maxiter, error);
            }
            if (Double.compare(xcurrent, _minVal) < 0) {
                return buildResults(xfinal, SolverStatus.MIN_VALUE_EXCEEDED, _maxIter - maxiter, error);
            }

            xcurrent = xfinal;
        }
        return buildResults(xfinal, SolverStatus.MAX_NUMBER_OF_ITERATIONS_EXCEEDED, Math.abs(_maxIter - maxiter), error);
    }
}