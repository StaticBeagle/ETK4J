package com.wildbitsfoundry.etk4j.math.optimize.solvers;

import com.wildbitsfoundry.etk4j.math.MathETK;
import com.wildbitsfoundry.etk4j.math.complex.Complex;
import com.wildbitsfoundry.etk4j.math.functions.ComplexUnivariateFunction;

import java.util.function.Function;

import static com.wildbitsfoundry.etk4j.math.optimize.solvers.SolverResults.SolverStatus;


public final class NewtonRaphsonComplex {

    protected int _maxIter = 100;
    protected double _absTol = 1e-9;
    protected double _relTol = 1e-6;
    protected double _maxVal = Double.POSITIVE_INFINITY;
    private double _minVal = Double.NEGATIVE_INFINITY;
    protected double _step = 0.001;
    protected Complex _x0;

    protected ComplexUnivariateFunction _func;
    protected ComplexUnivariateFunction _derivative;

    // create constructor that takes in derivative type

    public NewtonRaphsonComplex(ComplexUnivariateFunction func, ComplexUnivariateFunction derivative, Complex initialGuess) {
        _func = func;
        _derivative = derivative;
        _x0 = initialGuess;
    }

    public NewtonRaphsonComplex iterationLimit(int limit) {
        _maxIter = limit;
        return this;
    }

    public NewtonRaphsonComplex absTolerance(double tol) {
        _absTol = tol;
        return this;
    }

    public NewtonRaphsonComplex relTolerance(double tol) {
        _relTol = tol;
        return this;
    }

    public NewtonRaphsonComplex maxAllowedValue(double max) {
        _maxVal = max;
        return this;
    }

    public NewtonRaphsonComplex minAllowedValue(double min) {
        _minVal = min;
        return this;
    }

    public NewtonRaphsonComplex differentiationStepSize(double step) {
        _step = step;
        return this;
    }

    protected static SolverResults<Complex> buildResults(Complex xfinal, SolverStatus status, int iterCount, double error) {
        SolverResults<Complex> sr = new SolverResults<>();
        sr.setValue(xfinal);
        sr.setSolverStatus(status);
        sr.setNumberOfIterations(iterCount);
        sr.setError(error);
        return sr;
    }

    public SolverResults<Complex> solve() {
        return this.solve(_derivative);
    }

    protected SolverResults<Complex> solve(ComplexUnivariateFunction derivative) {
        int maxiter = _maxIter;
        double maxval = _maxVal;
        ComplexUnivariateFunction func = _func;

        Complex xcurrent = _x0;
        Complex xfinal = new Complex();
        double error = 0.0;
        while (maxiter-- > 0) {

            xfinal = xcurrent.subtract(func.evaluateAt(xcurrent).divide(derivative.evaluateAt(xcurrent)));

            error = xfinal.subtract(xcurrent).abs();
            if (error < _absTol + _relTol * xcurrent.abs()) {
                return buildResults(xfinal, SolverStatus.CONVERGED, _maxIter - maxiter, error);
            }

            if (Double.compare(xcurrent.abs(), maxval) > 0) {
                return buildResults(xfinal, SolverStatus.MAX_VALUE_EXCEEDED, _maxIter - maxiter, error);
            }
            if (Double.compare(xcurrent.abs(), _minVal) < 0) {
                return buildResults(xfinal, SolverStatus.MAX_VALUE_EXCEEDED, _maxIter - maxiter, error);
            }

            xcurrent = xfinal;
        }
        return buildResults(xfinal, SolverStatus.MAX_NUMBER_OF_ITERATIONS_EXCEEDED, Math.abs(_maxIter - maxiter), error);
    }

    public static void main(String[] args) {
        ComplexUnivariateFunction fx = z -> z.pow(2).add(1.0);
        ComplexUnivariateFunction fp = z -> z.multiply(2.0);
        SolverResults<Complex> results = new NewtonRaphsonComplex(fx, fp, new Complex(5, -2))
                .absTolerance(1e-9)
                .relTolerance(0.0)
                .iterationLimit(50)
                .solve();

        System.out.println(results);

        System.out.println(MathETK.rem(0.6, 0.2));
    }
}