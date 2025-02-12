package com.wildbitsfoundry.etk4j.math.laplace;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.functions.ComplexUnivariateFunction;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

/**
 * The {@code InverseLaplaceTransformStehfest} implements the inverse laplace transform using Stehfest's method.
 * @see <a href="https://www.codeproject.com/Articles/25189/Numerical-Laplace-Transforms-and-Inverse-Transform">Numerical-Laplace-Transforms-and-Inverse-Transform</a>
 */
public class InverseLaplaceTransformStehfest {
    private final double[] V;       //  Stehfest coefficients
    final static double ln2 = Math.log(2.0);

    /**
     * Construct an Inverse Laplace Transform using Talbot's method. The number of Stehfest's coefficients is set to 16.
     */
    InverseLaplaceTransformStehfest() {
        this(16);
    }

    /**
     * Construct an Inverse Laplace Transform using Talbot's method.
     * @param n The number of Stehfest's coefficients.
     */
    InverseLaplaceTransformStehfest(int n) {
        int N2 = n / 2;
        int NV = 2 * N2;
        V = new double[NV];
        int sign = 1;
        if ((N2 % 2) != 0)
            sign = -1;
        for (int i = 0; i < NV; i++) {
            int kmin = (i + 2) / 2;
            int kmax = i + 1;
            if (kmax > N2)
                kmax = N2;
            V[i] = 0;
            sign = -sign;
            for (int k = kmin; k <= kmax; k++) {
                V[i] = V[i] + (Math.pow(k, N2) / factorial(k)) * (factorial(2 * k)
                        / factorial(2 * k - i - 1)) / factorial(N2 - k) / factorial(k - 1)
                        / factorial(i + 1 - k);
            }
            V[i] = sign * V[i];
        }
    }

    /**
     * Perform the Inverse Laplace Transform.
     * @param function The function to apply the Inverse Laplace to.
     * @param time Argument at which to evaluate the time response.
     * @return {@code L<sup>-1</sup>{Y(s)} = y(t) evaluated at t = time.}
     */
    public double inverseTransform(UnivariateFunction function, double time) {
        if (time == 0.0) {
            time = ConstantsETK.DOUBLE_EPS;
        } else if (time == -0.0) {
            time = -ConstantsETK.DOUBLE_EPS;
        }
        double ln2t = ln2 / time;
        double x = 0;
        double y = 0;
        for (int i = 0; i < V.length; i++) {
            x += ln2t;
            y += V[i] * function.evaluateAt(x);
        }
        return ln2t * y;
    }

    /**
     * Basic factorial calculation.
     * @param n The argument used to calculate the factorial
     * @return {@code n!}
     */
    public static double factorial(int n) {
        double x = 1;
        if (n > 1) {
            for (int i = 2; i <= n; i++)
                x = i * x;
        }
        return x;
    }
}
