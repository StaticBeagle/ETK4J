package com.wildbitsfoundry.etk4j.statistics.regression;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

/**
 * The {@code LinearRegression} class implements a liner fit in the least square sense for a given set of (x, y) points.
 */
public class LinearRegression extends UnivariateRegression {

	public LinearRegression(double[] x, double[] y) {
		// we need over-determined so more than two points
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match.");
		}
		final int n = x.length;
		Matrix X = Matrix.vandermonde(x, n, 2);
		this.doRegression(X, x, y);
	}
}
