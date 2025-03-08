package com.wildbitsfoundry.etk4j.statistic.regression;

import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;

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
        MatrixDense X = MatrixDense.Factory.vandermonde(x, rows, n + 1);
        this.doRegression(X, x, y);
    }
}
