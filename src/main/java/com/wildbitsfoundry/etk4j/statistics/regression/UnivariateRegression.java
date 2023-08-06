package com.wildbitsfoundry.etk4j.statistics.regression;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.DoubleArrays;

abstract class UnivariateRegression {

	// SSE sum of squared errors
	// SSR sum of squared regression
	// SST total sum of squares = SSE + SSR
	private double _sse;
	private double _ssr;
	private double[] _beta;
	private double[] _residuals;
	private double _rnorm;	// Norm of residuals
	
	protected void doRegression(Matrix X, double[] x, double[] y) {
		final int n = x.length;

		Matrix Y = new Matrix(y, n);
		Matrix B = X.solve(Y);

		_beta = DoubleArrays.reverse(B.getArray());

		_residuals = Y.subtract(X.multiply(B)).getArray();
		_rnorm = DoubleArrays.norm2(_residuals);
		_sse = Math.pow(_rnorm, 2);

		double ybar = DoubleArrays.mean(y);

		// SST
		for (int i = 0; i < n; ++i) {
			double dev = y[i] - ybar;
			_ssr += dev * dev;
		}
	}

	/**
	 * R (a.k.a. Correlation).
	 * @return The computed {@code R}.
	 */
	public double R() {
		return Math.sqrt(this.R2());
	}

	/**
	 * R-squared (a.k.a. Coefficient of determination).
	 * @return The computed {@code R-squared}.
	 */
	public double R2() {
		return 1.0 - _sse / _ssr;
	}

	/**
	 * Sum of squared errors. The sum of the squared differences between predicted data points and observed data points.
	 * @return The computed SSE.
	 */
	public double SSE() {
		return _sse;
 	}

	/**
	 * Sum of squared regression. The sum of the squared differences between predicted data points and the mean of the
	 * response variable.
	 * @return The computed SSR.
	 */
	public double SSR() {
		return _ssr;
	}

	/**
	 * Sum of squared total.
	 * @return SSR + SSE,
	 */
	public double SST() {
		return _sse + _ssr;
	}

	/**
	 * Line of best fit coefficients.
	 * @return The coefficients of the line of best fit.
	 */
	public double[] beta() {
		return Arrays.copyOf(_beta, _beta.length);
	}

	/**
	 * Residuals. The difference between each data point and the line of best fit.
	 * @return The computed residuals.
	 */
	public double[] residuals() {
		return Arrays.copyOf(_residuals, _residuals.length);
	}

	/**
	 * Norm of residuals. The square root of the sum of the squares of the residuals.
	 * @return The {@link DoubleArrays#norm2(double[])} of the residuals}.
	 */
	public double normOfResiduals() {
		return _rnorm;
	}
}
