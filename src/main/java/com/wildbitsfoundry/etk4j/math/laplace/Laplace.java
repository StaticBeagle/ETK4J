package com.wildbitsfoundry.etk4j.math.laplace;

import com.wildbitsfoundry.etk4j.math.function.UnivariateFunction;

public final class Laplace {

    private Laplace() {

    }

    /**
     * Laplace Transform.
     * @param function Function to perform the Laplace transform to.
     * @param s Frequency at which to evaluate the transform.
     * @return {@code L{y(t} = Y(s)} evaluated at s.
     */
    public static double transform(UnivariateFunction function, double s) {
        final int DefaultIntegralN = 5000;
        double du = 0.5 / (double) DefaultIntegralN;
        double y = -function.evaluateAt(0) / 2.0;
        double u = 0;
        double limit = 1.0 - 1e-10;
        while (u < limit) {
            u += du;
            y += 2.0 * Math.pow(u, s - 1) * function.evaluateAt(-Math.log(u));
            u += du;
            y += Math.pow(u, s - 1) * function.evaluateAt(-Math.log(u));
        }
        return 2.0 * y * du / 3.0;
    }
}
