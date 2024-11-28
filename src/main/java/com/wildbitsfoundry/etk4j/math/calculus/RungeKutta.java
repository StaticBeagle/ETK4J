package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.functions.BivariateFunction;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;
import com.wildbitsfoundry.etk4j.util.Tuples;

public abstract class RungeKutta extends OdeSolver {
//    import numpy as np
//    from .base import OdeSolver, DenseOutput
//    from .common import (validate_max_step, validate_tol, select_initial_step,
//    norm, warn_extraneous, validate_first_step)
//    from . import dop853_coefficients

//Multiply steps computed from asymptotic behaviour of errors by this.
    private static final double SAFETY = 0.9;

    private static double MIN_FACTOR = 0.2; // Minimum allowed decrease in a step size.
    private static double MAX_FACTOR = 10; // Maximum allowed increase in a step size.
    protected double[] C;
    protected double[][] A;
    protected double[] B;
    protected double[] E;
    protected double[][] P;
    protected double[][] K;
    protected int order;
    protected int errorEstimatorOrder;
    protected int nStages;
    protected double yOld;
    protected double maxStep;
    protected double rTol;
    protected double aTol;
    protected double f;
    protected double hAbs;
    protected double errorExponent;
    protected Double hPrevious;

    // TODO add builder

    public RungeKutta(BivariateFunction func, double t0, double y0, Double tBound,
                      double maxStep, double rTol, double aTol, Double firstStep, int errorEstimatorOrder, int nStages,
                      double[][] A, double[] B, double[] C, double[] E, double[][] P) {
        super(func, t0, y0, tBound);
        this.maxStep = validateMaxStep(maxStep);
        Tuples.Tuple2<Double, Double> tolerances = validateTol(rTol, aTol, this.n);
        this.rTol = tolerances.getItem1();
        this.aTol = tolerances.getItem2();
        this.errorEstimatorOrder = errorEstimatorOrder;
        this.nStages = nStages;
        this.A = A;
        this.B = B;
        this.C = C;
        this.E = E;
        this.P = P;
        this.tBound = tBound;

        this.f = func.evaluateAt(this.t, this.y);
        if (firstStep == null) {
            this.hAbs = selectInitialStep(this.func, this.t, this.y, tBound, maxStep, this.f,
                    this.direction, this.errorEstimatorOrder, this.rTol, this.aTol);
        } else {
            this.hAbs = validateFirstStep(firstStep, t0, tBound);
        }
        this.K = new double[this.nStages + 1][this.n];
        this.errorExponent = -1.0 / (this.errorEstimatorOrder + 1);
    }

    protected double selectInitialStep(BivariateFunction func, double t0, double y0, Double tBound, double maxStep,
                                       double f0, double direction, double order, double rTol, double aTol) {
        double internalLength = Math.abs(tBound - t0);
        if(internalLength == 0) {
            return 0;
        }

        double scale = aTol + Math.abs(y0) * rTol;
        double d0 = Math.abs(y0 / scale);
        double d1 = Math.abs(f0 / scale);
        double h0;
        if(d0 < 1e-5 || d1 < 1e-5) {
            h0 = 1e-6;
        } else {
            h0 = 0.01 * d0 / d1;
        }

        h0 = Math.min(h0, internalLength);
        double y1 = y0 + h0 * direction * f0;
        double f1 = func.evaluateAt(t0 + h0 * direction, y1);
        double d2 = Math.abs((f1 - f0) / scale) / h0;

        double h1;
        if(d1 <= 1e-15 && d2 <= 1e-15) {
            h1 = Math.max(1e-6, h0 * 1e-3);
        } else {
            h1 = Math.pow(0.01 / Math.max(d1, d2), 1 / (order + 1));
        }
        return DoubleArrays.min(100 * h0, h1, internalLength, maxStep);
    }

    protected double validateMaxStep(double maxStep) {
        if(maxStep <= 0) {
            throw new IllegalArgumentException("Max step must be positive");
        }
        return maxStep;
    }

    protected Tuples.Tuple2<Double, Double> validateTol(double rTol, double aTol, int n) {
        if(rTol < 100 * ConstantsETK.DOUBLE_EPS) {
            rTol = Math.max(rTol, 100 * ConstantsETK.DOUBLE_EPS);
        }
        if(aTol < 0) {
            throw new IllegalArgumentException("aTol must be positive");
        }

        return new Tuples.Tuple2<>(rTol, aTol);
    }

    protected double validateFirstStep(double firstStep, double t0, Double tBound) {
        if(firstStep <= 0) {
            throw new IllegalArgumentException("first step must be positive");
        }
        if(firstStep > Math.abs(tBound - t0)) {
            throw new IllegalArgumentException("first step exceeds bound.");
        }
        return firstStep;
    }

    @Override
    protected Tuples.Tuple2<Boolean, String> stepImpl() {
        double t = this.t;
        double y = this.y;

        double maxStep = this.maxStep;
        double rTol = this.rTol;
        double aTol = this.aTol;
        double h = 0;
        double tNew = 0;
        double yNew = 0;
        double fNew = 0;

        double minStep = 10 * Math.abs(Math.nextAfter(t, direction * Double.POSITIVE_INFINITY) - t);

        double hAbs;
        if(this.hAbs > maxStep) {
            hAbs = maxStep;
        } else if(this.hAbs < minStep) {
            hAbs =  minStep;
        } else {
            hAbs = this.hAbs;
        }

        boolean stepAccepted = false;
        boolean stepRejected = false;

        while(!stepAccepted) {
            if(hAbs < minStep) {
                return new Tuples.Tuple2<>(false, "Too Small Step");
            }

            h = hAbs * this.direction;
            tNew = t + h;

            if(this.direction * (tNew - this.tBound) > 0) {
                tNew = this.tBound;
            }

            h = tNew - t;
            hAbs = Math.abs(h);

            Tuples.Tuple2<Double, Double> newValues = rkStep(this.func, t, y, this.f, h, this.A,
                    this.B, this.C, this.K);

            yNew = newValues.getItem1();
            fNew = newValues.getItem2();

            double scale = aTol + Math.max(Math.abs(y), Math.abs(yNew)) * rTol;
            double errorNorm = 0;
           // double errorNorm = this.estimateErrorNorm(this.K, h, scale);

            double factor = 0;
            if(errorNorm < 1) {
                if(errorNorm == 0) {
                    factor = MAX_FACTOR;
                } else {
                    Math.min(MAX_FACTOR, SAFETY * Math.pow(errorNorm, this.errorExponent));
                }

                if(stepRejected) {
                    factor = Math.min(1, factor);
                }

                hAbs *= factor;

                stepAccepted = true;
            } else {
                hAbs *= Math.max(MIN_FACTOR, SAFETY * Math.pow(errorNorm, this.errorExponent));
                stepRejected = true;
            }
        }

        this.hPrevious = h;
        this.yOld = y;

        this.t = tNew;
        this.y = yNew;

        this.hAbs = hAbs;
        this.f = fNew;

        return new Tuples.Tuple2<>(true, null);
    }

    private Tuples.Tuple2<Double, Double> rkStep(BivariateFunction func, double t, double y, double f, double h, double[][] A,
                          double[] B, double[] C, double[][] K) {
        K[0] = new double[] { f };
        for(int i = 1; i < A.length; i++) {
            double[][] kTranspose = DoubleArrays.transpose(K);
            double[] aSub = new double[i];
            double cSub = C[i];
            double[] kSub = new double[i];

            System.arraycopy(A[i], 0, aSub, 0, i);
            System.arraycopy(kTranspose[0], 0, kSub, 0, i);

            double dy = DoubleArrays.dot(kSub, aSub) * h;
            K[i][0] = func.evaluateAt(t + cSub * h, y + dy);
        }
        double[] kSub = new double[K.length - 1];
        System.arraycopy(DoubleArrays.transpose(K)[0], 0, kSub, 0, kSub.length);

        double yNew = y + h * DoubleArrays.dot(kSub, B);
        double fNew = func.evaluateAt(t + h, yNew);

        K[K.length - 1][0] = fNew;
        return new Tuples.Tuple2<>(yNew, fNew);
    }

    @Override
    protected Object getDenseOutputImpl() {
        return null;
    }

    public static void main(String[] args) {
        BivariateFunction func = (t, x) -> -x;
        RungeKutta rungeKutta = new RungeKutta45(func, 0.0, 1.0, 10.0, 1.0, 0.001,
                Math.exp(-6), null);

        for(int i = 0; i < 100; i++) {
            rungeKutta.step();
        }
    }
//
}
