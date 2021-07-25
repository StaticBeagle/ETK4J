package com.wildbitsfoundry.etk4j.math.optimize;

import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

public class Secant {

    // make this return a solver result
    public static double solve(UnivariateFunction fn, double x0, double x1, double absTol, double relTol, int maxIter) {
        if(x1 == x0) {
            // throw exception. They can't be equal
        }

        double xFinal = (x1 + x0) * 0.5;
        double fx0 = fn.evaluateAt(x0);
        double fx1 = fn.evaluateAt(x1);

        if (fx0 * fx1 > 0.0) {
            // throw root is not bracketed
            //return new RootFindingMethods.Root(false, Double.NaN, "Root is not bracketed", 0);
        }
        // swap if fx0 > fx1? let's think about it
        int iter;
        for(iter = 0; iter < maxIter; ++iter) {
            if(fx1 == fx0) {
                if(x1 != x0) {
                    // failed to converge
                }
                // converged
                return xFinal;
            } else {
                if(Math.abs(fx1) > Math.abs(fx0)) {
                    xFinal = (-fx0 / fx1 * x1 + x0) / (1 - fx0 / fx1);
                } else {
                    xFinal = (-fx1 / fx0 * x0 + x1) / (1 - fx1 / fx0);
                }
                // create an isclose function that uses rel and abstol
                double error = Math.abs(xFinal - x1);
                if(error < absTol) {
                    return xFinal;
                }
                x0 = x1;
                fx0 = fx1;
                x1 = xFinal;
                fx1 = fn.evaluateAt(x1);
            }
        }
        // max iterations exceeded
        return Double.NaN;
    }
}
