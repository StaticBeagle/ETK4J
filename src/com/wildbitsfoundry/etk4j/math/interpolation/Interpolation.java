package com.wildbitsfoundry.etk4j.math.interpolation;

public final class Interpolation {
	private Interpolation() {}
	
	public static double linear(double x0, double x1, double y0, double y1, double x) {
		double hx = x1 - x0;
		double t = (x - x0) / hx;
		return (y1 - y0) * t + y0;
	}
	
	public static double cosine(double x0, double x1, double y0, double y1, double x) {
		double hx = x1 - x0;
		double t = (x - x0) / hx;
		double t2 = (1 - Math.cos(t * Math.PI)) * 0.5;
		return (y1 - y0) * t2 + y0;
	}

	public static double neville(double[] x, double[] y, double xi) {
		if (x.length != y.length) {
			Splines.checkXYDimensions(x, y);
		}
		int length = x.length;
		double[][] N = new double[length][length];

		// Initializing first column
		for (int i = 0; i < length; ++i) {
			N[i][0] = y[i];
		}
		// Neville's method.
		for (int i = 1; i < length; ++i) {
			for (int j = 1; j <= i; ++j) {
				N[i][j] = ((xi - x[i - j]) * (N[i][j - 1]) - (xi - x[i]) * (N[i - 1][j - 1])) / (x[i] - x[i - j]);
			}
		}
		return N[length - 1][length - 1];
	}
	
	public static double spline(double[] x, double[] y, double xi) {
		if(x.length < 4) {
			return neville(x, y, xi);
		}
		Spline sp = CubicSpline.newCubicSpline(x, y);
		return sp.evaluateAt(xi);
	}
}
