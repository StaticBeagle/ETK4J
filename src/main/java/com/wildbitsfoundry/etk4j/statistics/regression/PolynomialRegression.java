package com.wildbitsfoundry.etk4j.statistics.regression;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

import java.util.Arrays;

/**
 * The {@code PolynomialRegression} class implements a polynomial fit in the least square sense for a given set of (x, y) points and degree n.
 */
public class PolynomialRegression extends UnivariateRegression {

    public PolynomialRegression(double[] x, double[] y, int n) {
        // we need over-determined so more than two points
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y dimensions must match.");
        }
        int rows = x.length;
        Matrix X = Matrix.vandermonde(x, rows, n + 1);
        this.doRegression(X, x, y);
    }

    public static void main(String[] args) {
        double[] x ={1, 2, 3, 4};
        double[] y = {1, 3.5, 8.4, 17.2};
        PolynomialRegression polynomialRegression = new PolynomialRegression(x, y, 2);
        System.out.println(polynomialRegression.R2());
        System.out.println(polynomialRegression.R());
        System.out.println(Arrays.toString(polynomialRegression.beta()));
        System.out.println(Arrays.toString(polynomialRegression.residuals()));
        System.out.println(polynomialRegression.SSE());
        System.out.println(polynomialRegression.SSR());
        System.out.println(polynomialRegression.SST());
    }
}
