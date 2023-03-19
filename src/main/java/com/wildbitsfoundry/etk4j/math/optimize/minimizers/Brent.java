package com.wildbitsfoundry.etk4j.math.optimize.minimizers;

import com.wildbitsfoundry.etk4j.math.optimize.OptimizerStatusType;

import java.util.function.BiFunction;

/*
Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/SciPy.
 */

/**
 * The {@code Brent} class contains the Brent's algorithm to find the minimum of a real valued function.
 */
public class Brent {

    private final BiFunction<Double, Object[], Double> function;
    private final double a;
    private final double b;
    private final Object[] params;
    private double tol = 1e-5;
    private int maxNumberOfIterations = 500;

    public Brent(BiFunction<Double, Object[], Double> function, double a, double b, Object... params) {
        this.function = function;
        this.a = a;
        this.b = b;
        this.params = params;
    }

    public Brent tolerance(double tol) {
        this.tol = tol;
        return this;
    }

    public Brent iterationLimit(int limit) {
        this.maxNumberOfIterations = limit;
        return this;
    }

    //https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fminbound.html
    public MinimizerResults<Double> minimize() {
        double x0 = a;
        double x1 = b;
        double sqrtEPS = Math.sqrt(2.2e-16);
        double goldenMean = 0.5 * (3.0 - Math.sqrt(5.0));
        double fulc = x0 + goldenMean * (x1 - x0);
        double nfc = fulc;
        double xf = fulc;
        double rat = 0;
        double e = 0.0;
        double x = xf;
        double fx = function.apply(x, params);
        int num = 1;
        double fu = Double.POSITIVE_INFINITY;

        double ffulc = fx;
        double fnfc = fx;
        double xm = 0.5 * (x0 + x1);
        double tol1 = sqrtEPS * Math.abs(xf) + tol / 3.0;
        double tol2 = 2.0 * tol1;

        while (Math.abs(xf - xm) > (tol2 - 0.5 * (x1 - x0))) {
            int golden = 1;
            // Check for parabolic fit
            if (Math.abs(e) > tol1) {
                golden = 0;
                double r = (xf - nfc) * (fx - ffulc);
                double q = (xf - fulc) * (fx - fnfc);
                double p = (xf - fulc) * q - (xf - nfc) * r;
                q = 2.0 * (q - r);
                if (q > 0.0) {
                    p = -p;
                }
                q = Math.abs(q);
                r = e;
                e = rat;

                // Check for acceptability of parabola
                if ((Math.abs(p) < Math.abs(0.5 * q * r)) && (p > q * (x0 - xf)) && (p < q * (x1 - xf))) {
                    rat = (p + 0.0) / q;
                    x = xf + rat;


                    if ((x - x0) < tol2 || (x1 - x) < tol2) {
                        double si = Math.signum(xm - xf) + ((xm - xf) == 0 ? 1 : 0);
                        rat = tol1 * si;
                    }
                } else {      // #do a golden -section step
                    golden = 1;
                }
            }

            if (golden != 0) { //#do a golden -section step
                if (xf >= xm) {
                    e = x0 - xf;
                } else {
                    e = x1 - xf;
                }
                rat = goldenMean * e;
            }

            double si = Math.signum(rat) + (rat == 0 ? 1 : 0);
            x = xf + si * Math.max(Math.abs(rat), tol1);
            fu = function.apply(x, params);
            num += 1;

            if (fu <= fx) {
                if (x >= xf) {
                    x0 = xf;
                } else {
                    x1 = xf;
                }
                fulc = nfc;
                ffulc = fnfc;
                nfc = xf;
                fnfc = fx;
                xf = x;
                fx = fu;
            } else {
                if (x < xf) {
                    x0 = x;
                } else {
                    x1 = x;
                }
                if (fu <= fnfc || nfc == xf) {
                    fulc = nfc;
                    ffulc = fnfc;
                    nfc = x;
                    fnfc = fu;
                } else if (fu <= ffulc || fulc == xf || fulc == nfc) {
                    fulc = x;
                    ffulc = fu;
                }
            }

            xm = 0.5 * (x0 + x1);
            tol1 = sqrtEPS * Math.abs(xf) + tol / 3.0;
            tol2 = 2.0 * tol1;

            if (num >= maxNumberOfIterations) {
                MinimizerResults<Double> minResults = new MinimizerResults<>();
                minResults.setMinimizerStatus("Maximum number of iterations exceeded");
                minResults.setOptimizerStatusType(OptimizerStatusType.MAXIMUM_NUMBER_OF_ITERATIONS_EXCEEDED);
                minResults.setHasConverged(false);
                minResults.setValue(xf);
                minResults.setFunctionValue(fx);
                minResults.setNumberOfIterations(num);
                return minResults;
            }
        }

        if (Double.isNaN(xf) || Double.isNaN(fx) || Double.isNaN(fu)) {
            MinimizerResults<Double> minResults = new MinimizerResults<>();
            minResults.setMinimizerStatus("The minimum found was Double.NaN");
            minResults.setOptimizerStatusType(OptimizerStatusType.VALUE_FOUND_WAS_NAN);
            minResults.setHasConverged(false);
            minResults.setValue(Double.NaN);
            minResults.setValue(Double.NaN);
            minResults.setNumberOfIterations(num);
            return minResults;
        }
        MinimizerResults<Double> minResults = new MinimizerResults<>();
        minResults.setMinimizerStatus("Converged");
        minResults.setOptimizerStatusType(OptimizerStatusType.CONVERGED);
        minResults.setHasConverged(true);
        minResults.setValue(xf);
        minResults.setFunctionValue(fx);
        minResults.setNumberOfIterations(num);
        return minResults;
    }
}
