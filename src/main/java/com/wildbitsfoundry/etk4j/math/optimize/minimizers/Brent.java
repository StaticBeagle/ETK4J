package com.wildbitsfoundry.etk4j.math.optimize.minimizers;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;

import java.util.function.BiFunction;

// TODO move brent from solvers here
public class Brent {

    /*
    Copyright (c) 2001-2002 Enthought, Inc. 2003-2022, SciPy Developers.
    All rights reserved. See https://github.com/StaticBeagle/ETK4J/blob/master/COPYING.
     */
    // TODO test
    // https://github.com/scipy/scipy/blob/v1.7.1/scipy/optimize/optimize.py#L1904-L1979
    public static double brentsMinimizer(BiFunction<Double, Object[], Double> func,
                                          double a, double b, double tol, int maxIter, Object... params) {

        double sqrtEPS = Math.sqrt(2.2e-16);
        double goldenMean = 0.5 * (3.0 - Math.sqrt(5.0));
        double fulc = a + goldenMean * (b - a);
        double nfc = fulc;
        double xf = fulc;
        double rat = 0;
        double e = 0.0;
        double x = xf;
        double fx = func.apply(x, params);
        int num = 1;
        double fu = Double.POSITIVE_INFINITY;

        double ffulc = fx;
        double fnfc = fx;
        double xm = 0.5 * (a + b);
        double tol1 = sqrtEPS * Math.abs(xf) + tol / 3.0;
        double tol2 = 2.0 * tol1;


        while (Math.abs(xf - xm) > (tol2 - 0.5 * (b - a))) {
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
                if ((Math.abs(p) < Math.abs(0.5 * q * r)) && (p > q * (a - xf)) && (p < q * (b - xf))) {
                    rat = (p + 0.0) / q;
                    x = xf + rat;


                    if ((x - a) < tol2 || (b - x) < tol2) {
                        double si = Math.signum(xm - xf) + ((xm - xf) == 0 ? 1 : 0);
                        rat = tol1 * si;
                    }
                } else {      // #do a golden -section step
                    golden = 1;
                }
            }

            if (golden != 0) { //#do a golden -section step
                if (xf >= xm) {
                    e = a - xf;
                } else {
                    e = b - xf;
                }
                rat = goldenMean * e;
            }


            double si = Math.signum(rat) + (rat == 0 ? 1 : 0);
            x = xf + si * Math.max(Math.abs(rat), tol1);
            fu = func.apply(x, params);
            num += 1;

            if (fu <= fx) {
                if (x >= xf) {
                    a = xf;
                } else {
                    b = xf;
                }
                fulc = nfc;
                ffulc = fnfc;
                nfc = xf;
                fnfc = fx;
                xf = x;
                fx = fu;
            } else {
                if (x < xf) {
                    a = x;
                } else {
                    b = x;
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

            xm = 0.5 * (a + b);
            tol1 = sqrtEPS * Math.abs(xf) + tol / 3.0;
            tol2 = 2.0 * tol1;

            if (num >= maxIter)
                // max number of iterations exceeded
                break;
        }

        if(Double.isNaN(xf) || Double.isNaN(fx) || Double.isNaN(fu)) {
            // not converged
        }
        return xf;
    }
}
