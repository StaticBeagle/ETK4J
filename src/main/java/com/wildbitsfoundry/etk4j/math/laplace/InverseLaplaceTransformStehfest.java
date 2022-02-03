package com.wildbitsfoundry.etk4j.math.laplace;

import com.wildbitsfoundry.etk4j.constants.ConstantsETK;
import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;

//https://www.codeproject.com/Articles/25189/Numerical-Laplace-Transforms-and-Inverse-Transform
public class InverseLaplaceTransformStehfest {
    private final double[] V;       //  Stehfest coefficients
    final static double ln2 = Math.log(2.0);

    InverseLaplaceTransformStehfest() {
        this(16);
    }

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

    public double InverseTransform(UnivariateFunction function, double time) {
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

    public static double factorial(int N) {
        double x = 1;
        if (N > 1) {
            for (int i = 2; i <= N; i++)
                x = i * x;
        }
        return x;
    }

}
