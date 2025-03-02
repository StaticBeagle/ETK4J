package com.wildbitsfoundry.etk4j.math.calculus;

import com.wildbitsfoundry.etk4j.math.function.MultivariateFunction;

public final class PartialDerivatives {

    private PartialDerivatives() {}

    /**
     * Centered difference using a 5 point stencil.
     * @param func Function to differentiate.
     * @param x Argument at which to evaluate the derivative.
     * @param h Step size of the differentiation.
     * @return The derivative evaluated at <code>x</code> with step size <code>h</code>.
     * @throws IllegalArgumentException If the step size is less than or equal to zero.
     */
    public static double[] centeredDifference5Points(MultivariateFunction func, double[] x, double h) {
        if(h <= 0) {
            throw new IllegalArgumentException("The step size h must be greater than zero");
        }
        double[] result = new double[x.length];
        for(int i = 0; i < x.length; i++) {
            final double tmp = x[i];
            x[i] = tmp + h;
            double a = func.evaluateAt(x);

            x[i] = tmp - h;
            double b = func.evaluateAt(x);

            x[i] = tmp + 2 * h;
            double c = func.evaluateAt(x);

            x[i] = tmp - 2 * h;
            double d = func.evaluateAt(x);

            x[i] = tmp;

            double num = 8 * (a - b) - c + d;
            result[i] = num / (12 * h);
        }
        return result;
    }
}
