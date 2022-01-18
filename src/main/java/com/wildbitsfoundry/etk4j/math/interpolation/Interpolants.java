package com.wildbitsfoundry.etk4j.math.interpolation;

import com.wildbitsfoundry.etk4j.math.curvefitting.CurveFitting;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkXYDimensions;
import static com.wildbitsfoundry.etk4j.util.validation.DimensionCheckers.checkMinXLength;

public final class Interpolants {
	private Interpolants() {
	}

	public static double linear(double x0, double x1, double y0, double y1, double x) {
		double hx = x1 - x0;
		double t = (x - x0) / hx;
		return (y1 - y0) * t + y0;
	}

	public static double neville(double[] x, double[] y, double xi) {
		checkXYDimensions(x, y);
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

	public static double quadratic(double x0, double x1, double x2, double y0, double y1, double y2, double xi) {
		double[] parabola = CurveFitting.parabola(x0, x1, x2, y0, y1, y2);
		return parabola[2] + xi * (parabola[1] + xi * parabola[0]);
	}

	public static double spline(double[] x, double[] y, double xi) {
		checkXYDimensions(x, y);
		checkMinXLength(x, 2);
		
		final int length = x.length;
		if (length == 2) {
			return linear(x[0], x[1], y[0], y[1], xi);
		}
		if(length == 3) {
			return quadratic(x[0], x[1], x[2], y[0], y[1], y[2], xi);
		}
		return CubicSpline.newCubicSpline(x, y).evaluateAt(xi);
	}
}
