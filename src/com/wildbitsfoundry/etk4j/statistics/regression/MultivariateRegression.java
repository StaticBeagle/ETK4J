package com.wildbitsfoundry.etk4j.statistics.regression;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.NumArrays;

public class MultivariateRegression{
	// SSE sum of squared errors
	// SSR sum of squared regression
	// SST total sum of squares = SSE + SSR
	private double _sse;
	private double _ssr;
	private double[] _beta;
	private double[] _residuals;
	private double _rnorm;	// Norm of residuals
	
	protected void doRegression(Matrix X, double[][] x, double[] y) {
		final int rows = y.length;

		Matrix Y = new Matrix(y, rows);
		_beta = X.solve(Y).getArray();
	}
	
	public double[] beta() {
		return Arrays.copyOf(_beta, _beta.length);
	}
}
