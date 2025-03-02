package com.wildbitsfoundry.etk4j.statistic.regression;

import com.wildbitsfoundry.etk4j.math.linearalgebra.MatrixDense;

/**
 * The {@code LogarithmicRegression} class implements a logarithmic fit in the least square sense for a given set of
 * (x, y) points.
 */
public class LogarithmicRegression extends UnivariateRegression {

	public LogarithmicRegression(double[] x, double[] y) {
		// we need over-determined so more than two points
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match.");
		}
		final int n = x.length;
		double[][] X = new double[n][2];
		for (int i = 0; i < n; ++i) {
			X[i][0] = 1.0;
			X[i][1] = Math.log(x[i]);
		}
		this.doRegression(new MatrixDense(X), x, y);
	}
}
