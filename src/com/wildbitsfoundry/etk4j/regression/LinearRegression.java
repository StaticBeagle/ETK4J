package com.wildbitsfoundry.etk4j.regression;

import java.util.Arrays;

import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrices;
import com.wildbitsfoundry.etk4j.math.linearalgebra.Matrix;
import com.wildbitsfoundry.etk4j.util.ArrayUtils;

public class LinearRegression implements RegressionModel {

	private double[] _beta;
	private double[] _residuals;
	private double _sst;
	private double _r2;
	private double _sse;
	
	// SSE sum of squared errors
	// SSR sum of squared regression
	// SST total sum of squares = SSE + SSR

	public LinearRegression(double[] x, double[] y) {
		// we need over-determined so more than two points
		if (x.length != y.length) {
			throw new IllegalArgumentException("x and y dimensions must match");
		}

		final int n = x.length;
		Matrix X = Matrices.Vandermonde(x, n, 1);
		Matrix Y = new Matrix(y, n);

		Matrix B = X.solve(Y);
		
		_beta = B.getArray();
		
		_residuals = B.subtract(X.multiply(B)).getArray();
		_sse = Math.pow(ArrayUtils.norm(_residuals), 2);
		
		double ymean = ArrayUtils.mean(y);
		
		// SST
		for(int i = 0; i < n; ++i) {
			double dev = y[i] - ymean;
			_sst += dev * dev;
		}

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
		return Arrays.copyOf(_beta, _beta.length);
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
