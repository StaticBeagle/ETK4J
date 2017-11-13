package com.wildbitsfoundry.etk4j.statistics.regression;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

public class LinearRegression extends UnivariateRegression {

	public LinearRegression(double[] x, double[] y) {
		// we need over-determined so more than two points
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}
		final int n = x.length;
		Matrix X = Matrices.Vandermonde(x, n, 2);
		this.doRegression(X, x, y);
	}
	
	public double R() {
		return Math.sqrt(this.R2());
	}
}
