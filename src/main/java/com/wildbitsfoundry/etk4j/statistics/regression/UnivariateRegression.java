package com.wildbitsfoundry.etk4j.statistics.regression;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.NumArrays;

abstract class UnivariateRegression implements RegressionModel {

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

		_beta = NumArrays.reverse(B.getArray());

		_residuals = Y.subtract(X.multiply(B)).getArray();
		_rnorm = NumArrays.norm2(_residuals);
		_sse = Math.pow(_rnorm, 2);

		double ybar = NumArrays.mean(y);

		// SST
		for (int i = 0; i < n; ++i) {
			double dev = y[i] - ybar;
			_ssr += dev * dev;
		}
	}

	@Override
	public double R2() {
		return 1.0 - _sse / _ssr;
	}

	@Override
	public double SSE() {
		return _sse;
 	}
	
	@Override
	public double SSR() {
		return _ssr;
	}

	@Override
	public double SST() {
		return _sse + _ssr;
	}

	@Override
	public double[] beta() {
		return Arrays.copyOf(_beta, _beta.length);
	}

	@Override
	public double[] residuals() {
		return Arrays.copyOf(_residuals, _residuals.length);
	}
	
	@Override
	public double normOfResiduals() {
		return _rnorm;
	}

}
