package com.wildbitsfoundry.etk4j.math.polynomials;

import com.wildbitsfoundry.etk4j.math.functions.UnivariateFunction;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

import java.util.Arrays;

import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;

public class LagrangePolynomial implements UnivariateFunction {

    private final double[] x;
    private final double[] y;
    private final double[] weights;

    public LagrangePolynomial(double[] x, double[] y) {
        checkXYDimensions(x, y);
        this.x = Arrays.copyOf(x, x.length);
        this.y = Arrays.copyOf(y, y.length);
        this.weights = computeWeights(x);
    }

    private double[] computeWeights(double[] x) {
        int n = x.length;
        double[] w = DoubleArrays.ones(n);

        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    w[i] /= (x[i] - x[j]);
                }
            }
        }
        return w;
    }

    @Override
    public double evaluateAt(double x) {
        double numerator = 0.0;
        double denominator = 0.0;

        for(int i = 0; i < this.x.length; i++) {
            if(x == this.x[i]) {
                return y[i];
            }
            double tmp = weights[i] / (x - this.x[i]);
            numerator += tmp * y[i];
            denominator += tmp;
        }
        return numerator / denominator;
    }

    public static LagrangePolynomial lagrangeFit(double[] x, double[] y) {
        return new LagrangePolynomial(x, y);
    }
}
