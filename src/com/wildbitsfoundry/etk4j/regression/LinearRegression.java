package com.wildbitsfoundry.etk4j.regression;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;

public class LinearRegression implements RegressionModel {

	public LinearRegression(double[] x, double[] y) {
		// we need over-determined so more than two points
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}

		final int n = x.length;
		Matrix X = Matrices.Vandermonde(x, n, 1);
		Matrix Y = new Matrix(y, n);

		Matrix B = X.solve(Y);

	}

	@Override
	public double R2() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double SSE() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double SST() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double[] beta() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] residuals() {
		// TODO Auto-generated method stub
		return null;
	}

	public static void main(String[] args) {
		double[] x = new double[] { 1, 2, 3, 4, 5, 6, 7, 8 };
		double[] y = new double[] { 2, 4, 9, 3, 7, 9, 8, 7 };
		LinearRegression lr = new LinearRegression(x, y);
	}
}
