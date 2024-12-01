package com.wildbitsfoundry.etk4j.signals.stoz;

import com.wildbitsfoundry.etk4j.control.TransferFunction;
import com.wildbitsfoundry.etk4j.math.polynomials.RationalFunction;
import com.wildbitsfoundry.etk4j.signals.filters.ButterWorth;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

import static com.wildbitsfoundry.etk4j.math.MathETK.combinations;

/**
 * The {@code Bilinear Transform} class contains methods to convert a continuous time system numerator and denominator
 * to the discrete time domain
 */
public class BilinearTransform {
    /**
     * Perform bilinear transform on analog filter coefficients. No frequency pre-warping is performed in this method.
     *
     * @param num  Numerator coefficients of the analog filter
     * @param den  Denominator coefficients of the analog filter
     * @param fs Sampling frequency in Hz
     * @return A 2D array where [0] contains the discrete time numerator and [1] contains teh discrete time denominator
     * coefficients
     */
    public static double[][] transform(double[] num, double[] den, double fs) {
        // Maximum filter order
        int d = den.length - 1;
        int n = num.length - 1;
        int m = Math.max(n, d);
        int np = m;
        int dp = m;

        double[] numDigital = new double[np + 1];
        double[] denDigital = new double[dp + 1];

        for (int j = 0; j < np + 1; j++) {
            double val = 0;
            for (int i = 0; i < n + 1; i++) {
                for (int k = 0; k < i + 1; k++) {
                    for (int l = 0; l < m - i + 1; l++) {
                        if (k + l == j) {
                            val += (combinations(i, k) * combinations(m - i, l) * num[n - i]
                                    * Math.pow(2 * fs, i) * Math.pow(-1, k));
                        }
                    }
                }
            }
            numDigital[j] = val;
        }
        for (int j = 0; j < dp + 1; j++) {
            double val = 0;
            for (int i = 0; i < d + 1; i++) {
                for (int k = 0; k < i + 1; k++) {
                    for (int l = 0; l < m - i + 1; l++) {
                        if (k + l == j) {
                            val += (combinations(i, k) * combinations(m - i, l) * den[d - i]
                                    * Math.pow(2 * fs, i) * Math.pow(-1, k));
                        }
                    }
                }
            }
            denDigital[j] = val;
        }

        // Normalize by the first non-zero element in aDigital
        int i = 0;
        while(i <= denDigital.length) {
            if(i == denDigital.length) {
                throw new RuntimeException("fix this");
            }
            if(denDigital[i] != 0) {
                break;
            }
            i++;
        }
        double normalizingFactor = denDigital[i];
        DoubleArrays.divideElementWiseInPlace(numDigital, normalizingFactor);
        DoubleArrays.divideElementWiseInPlace(denDigital, normalizingFactor);
        return new double[][]{numDigital, denDigital};
    }
}